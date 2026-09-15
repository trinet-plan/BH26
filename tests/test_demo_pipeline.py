import json
import unittest
import uuid
from pathlib import Path

from acmg.cli import main
from acmg.core.models import CRITERIA


ROOT = Path(__file__).resolve().parents[1]


class DemoPipelineTests(unittest.TestCase):
    def test_all_demo_records_evaluate_from_committed_cache(self):
        base = ROOT / ".work" / "demo-pipeline" / uuid.uuid4().hex
        prepared = base / "prepared"
        evaluated = base / "evaluated"
        case2_evaluated = base / "case2-pm2-evaluated"
        code = main([
            "prepare-demo-online", "--input-dir", str(ROOT / "demo-data"),
            "--cache-dir", str(ROOT / "tests" / "fixtures" / "ensembl-cache"),
            "--evidence-cache-dir", str(ROOT / "tests" / "fixtures" / "external-cache"),
            "--output-dir", str(prepared), "--ensembl-release", "116",
            "--with-gnomad", "--gnomad-release", "4.1.1",
            "--with-clinvar", "--clinvar-release", "2026-09-15",
            "--with-pm1-hotspot", "--rules", str(ROOT / "config" / "demo-rules.json"),
            "--with-dbnsfp",
            "--offline",
        ])
        self.assertEqual(code, 0)
        audit = json.loads((prepared / "audit.json").read_text(encoding="utf-8"))
        self.assertEqual(len(audit["records"]), 28)
        self.assertFalse(any(row["identity_status"] == "PENDING" for row in audit["records"]))
        manifest = json.loads((prepared / "identity-manifest.json").read_text(encoding="utf-8"))
        self.assertFalse(manifest["network_used"])
        self.assertEqual(manifest["annotation_evidence"], 17)
        # 42 uncalibrated VEP scores plus 22 version-pinned dbNSFP scores.
        self.assertEqual(manifest["computational_evidence"], 64)
        providers = {item["provider"]: item for item in manifest["external_providers"]}
        self.assertEqual(providers["gnomAD"]["queried_variants"], 17)
        self.assertEqual(providers["gnomAD"]["observed_variants"], 15)
        self.assertEqual(providers["ClinVar"]["matched_records"], 13)
        self.assertEqual(providers["ClinVar"]["pm5_residue_searches"], 10)
        self.assertEqual(providers["ClinVar"]["ps1_comparator_errors"], [])

        code = main([
            "evaluate", "--input", str(prepared / "variants.json"),
            "--evidence", str(prepared / "evidence.json"), "--offline",
            "--config", str(ROOT / "config" / "demo-rules.json"),
            "--output-dir", str(evaluated),
        ])
        self.assertEqual(code, 0)
        results = json.loads((evaluated / "results.json").read_text(encoding="utf-8"))
        self.assertEqual(len(results["records"]), 28)
        self.assertTrue(all(len(record["results"]) == len(CRITERIA)
                            for record in results["records"]))
        self.assertEqual(results["va_spec_export"], "VALIDATED")
        self.assertTrue((evaluated / "evidence-lines.json").is_file())
        pm2 = [result for record in results["records"] for result in record["results"]
               if result["criterion"] == "PM2"]
        self.assertTrue(any(result["status"] == "NOT_MET" for result in pm2))

        # Calibrated PP3/BP4 need the dbNSFP release that produced the scores.
        dbnsfp = providers["MyVariant.info dbNSFP"]
        self.assertEqual(dbnsfp["provider_version"], "4.8a")
        self.assertEqual(dbnsfp["errors"], [])
        computational = [result for record in results["records"] for result in record["results"]
                         if result["criterion"] in ("PP3", "BP4")]
        met = [result for result in computational if result["status"] == "MET"]
        self.assertEqual(len(met), 18)
        strengths = {result["evidence_outcome"] for result in met}
        self.assertEqual(strengths, {"PP3", "PP3_moderate", "PP3_strong",
                                     "BP4", "BP4_moderate"})
        for result in met:
            mechanisms = {item["mechanism"] for item in result["provenance"]["applied_calibrations"]}
            self.assertTrue(mechanisms <= {"protein", "splicing"})
            for score in [item for item in result["evidence"] if item.get("predictor") == "REVEL"]:
                self.assertEqual(score["predictor_version"], "dbNSFP-4.8a")
                self.assertTrue(score["calibration_eligible"])
        # Splicing alone carries BP4 for synonymous variants, where no protein score applies.
        splicing_only = [result for result in met
                         if [item["predictor"] for item in
                             result["provenance"]["applied_calibrations"]] == ["SpliceAI"]]
        self.assertTrue(splicing_only)
        self.assertTrue(all(result["criterion"] == "BP4" for result in splicing_only))
        # The unreported SpliceAI release travels with the result as a declared assumption.
        assertions = {item["asserted_version"] for result in met
                      for item in result["provenance"].get("version_assertions", [])}
        self.assertEqual(len(assertions), 1)

        # ClinVar missense density around each residue, under the committed PM1 policy.
        hotspot = providers["ClinVar protein hotspot density"]
        self.assertEqual(hotspot["policy_version"], "PM1-hotspot-v1")
        self.assertEqual(hotspot["region_evidence"], 10)
        self.assertEqual(hotspot["errors"], [])
        pm1 = [result for record in results["records"] for result in record["results"]
               if result["criterion"] == "PM1"]
        statuses = {status: sum(result["status"] == status for result in pm1)
                    for status in ("MET", "NOT_MET", "NOT_EVALUATED")}
        # Every demo variant that ClinVar classifies is excluded from its own density, so no
        # record is NOT_MET on the strength of its own submitted classification.
        self.assertEqual(statuses, {"MET": 2, "NOT_MET": 0, "NOT_EVALUATED": 26})
        # PM5 needs a residue-scoped search; MYH7 p.Arg719 has a pathogenic ClinVar comparator.
        pm5 = [result for record in results["records"] for result in record["results"]
               if result["criterion"] == "PM5"]
        pm5_statuses = {status: sum(result["status"] == status for result in pm5)
                        for status in ("MET", "NOT_MET", "NOT_APPLICABLE")}
        self.assertEqual(pm5_statuses, {"MET": 1, "NOT_MET": 17, "NOT_APPLICABLE": 10})
        pm5_met = next(result for result in pm5 if result["status"] == "MET")
        self.assertEqual(pm5_met["evidence_outcome"], "PM5")
        self.assertEqual(pm5_met["provenance"]["assessment_scope"], "protein_level")
        self.assertEqual(pm5_met["provenance"]["condition_assessment"], "NOT_EVALUATED")

        met = [result for result in pm1 if result["status"] == "MET"]
        for result in met:
            self.assertEqual(result["evidence_outcome"], "PM1")
            self.assertEqual(result["provenance"]["pm1_route"], "mutational_hotspot")
            self.assertEqual(result["provenance"]["assessment_method"], "automated")
            self.assertEqual(result["provenance"]["benign_count"], 0)
            self.assertGreaterEqual(result["provenance"]["pathogenic_count"], 3)
            # Disease relevance is unresolved without a condition, so it stays a review point.
            self.assertEqual(result["provenance"]["condition_assessment"], "NOT_EVALUATED")
            self.assertTrue(result["review_points"])

        # Exercise a real case2 record and the committed gnomAD response with an explicit,
        # test-only rarity threshold. This is a regression test, not a clinical policy.
        code = main([
            "evaluate", "--input", str(prepared / "variants.json"),
            "--evidence", str(prepared / "evidence.json"), "--offline",
            "--config", str(ROOT / "tests" / "fixtures" / "case2-pm2-rules.json"),
            "--criteria", "PM2", "--output-dir", str(case2_evaluated),
        ])
        self.assertEqual(code, 0)
        case2_results = json.loads(
            (case2_evaluated / "results.json").read_text(encoding="utf-8")
        )
        by_id = {record["record_id"]: record["results"][0]
                 for record in case2_results["records"]}
        self.assertEqual(by_id["case2:24:1"]["status"], "MET")
        self.assertEqual(by_id["case2:24:1"]["evidence_outcome"], "PM2_supporting")
        maximum_af = max(float(item["AF"]) for item in by_id["case2:24:1"]["evidence"])
        self.assertLessEqual(maximum_af, 0.000014)

        line = json.loads(
            (case2_evaluated / "va-spec-1.0.1" / "case2_24_1--PM2.json")
            .read_text(encoding="utf-8")
        )
        self.assertEqual(line["directionOfEvidenceProvided"], "supports")
        self.assertEqual(line["strengthOfEvidenceProvided"]["primaryCoding"]["code"],
                         "supporting")
        self.assertEqual(line["evidenceOutcome"]["primaryCoding"]["code"], "PM2_supporting")
        self.assertTrue(all(item["type"] == "CohortAlleleFrequencyStudyResult"
                            for item in line["hasEvidenceItems"]))

        # case2-var2 is not returned by gnomAD. A missing variant response has no AN or
        # callability evidence, so it must not be silently converted to AF=0.
        self.assertEqual(by_id["case2:25:1"]["status"], "NOT_EVALUATED")
        self.assertIn("population", by_id["case2:25:1"]["missing_inputs"])


if __name__ == "__main__":
    unittest.main()
