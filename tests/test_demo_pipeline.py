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
            "--offline",
        ])
        self.assertEqual(code, 0)
        audit = json.loads((prepared / "audit.json").read_text(encoding="utf-8"))
        self.assertEqual(len(audit["records"]), 28)
        self.assertFalse(any(row["identity_status"] == "PENDING" for row in audit["records"]))
        manifest = json.loads((prepared / "identity-manifest.json").read_text(encoding="utf-8"))
        self.assertFalse(manifest["network_used"])
        self.assertEqual(manifest["annotation_evidence"], 17)
        self.assertEqual(manifest["computational_evidence"], 42)
        self.assertEqual(manifest["external_providers"][0]["queried_variants"], 17)
        self.assertEqual(manifest["external_providers"][0]["observed_variants"], 15)
        self.assertEqual(manifest["external_providers"][1]["matched_records"], 13)

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

        # ClinVar missense density around each residue, under the committed PM1 policy.
        hotspot = manifest["external_providers"][2]
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
