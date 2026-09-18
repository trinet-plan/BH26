import json
import unittest
import uuid
from pathlib import Path

from acmg_pipeline.automated_cli import main
from acmg_pipeline.automated_core.models import CRITERIA


ROOT = Path(__file__).resolve().parents[1]


class DemoPipelineTests(unittest.TestCase):
    def test_all_demo_records_evaluate_from_committed_cache(self):
        base = ROOT / ".work" / "demo-pipeline" / uuid.uuid4().hex
        prepared = base / "prepared"
        evaluated = base / "evaluated"
        case2_evaluated = base / "case2-pm2-evaluated"
        mechanism_evidence = base / "gene-disease" / "evidence.json"
        mechanism_evaluated = base / "gene-disease-evaluated"
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
            "--context", str(ROOT / "config" / "curated-context.json"),
            "--output-dir", str(evaluated),
        ])
        self.assertEqual(code, 0)
        results = json.loads((evaluated / "results.json").read_text(encoding="utf-8"))
        self.assertEqual(len(results["records"]), 28)
        self.assertTrue(all(len(record["results"]) == len(CRITERIA)
                            for record in results["records"]))
        self.assertEqual(results["va_spec_export"], "VALIDATED")
        self.assertTrue((evaluated / "evidence-lines.json").is_file())
        envelope = json.loads((evaluated / "evidence-lines.json").read_text(encoding="utf-8"))
        self.assertEqual(sum(len(record["criterion_assessments"])
                             for record in envelope["records"]), 28 * len(CRITERIA))
        for record in envelope["records"]:
            self.assertEqual(
                {item["criterion"] for item in record["criterion_assessments"]},
                set(CRITERIA),
            )
            catalog = record["referenced_evidence"]
            references = {
                identifier for item in record["criterion_assessments"]
                for identifier in item["evidenceItemIds"]
            }
            self.assertTrue(references <= set(catalog))
        pm2 = [result for record in results["records"] for result in record["results"]
               if result["criterion"] == "PM2"]
        self.assertTrue(any(result["status"] == "not_met" for result in pm2))

        # Calibrated PP3/BP4 need the dbNSFP release that produced the scores.
        dbnsfp = providers["MyVariant.info dbNSFP"]
        self.assertEqual(dbnsfp["provider_version"], "4.8a")
        self.assertEqual(dbnsfp["errors"], [])
        computational = [result for record in results["records"] for result in record["results"]
                         if result["criterion"] in ("PP3", "BP4")]
        met = [result for result in computational if result["status"] == "met"]
        # PP3 6 + BP4 16. Four of the BP4 calls rest on SpliceAI 0 while REVEL (0.398, 0.398,
        # 0.398, 0.404) sits in Pejaver et al. 2022's gap between the PP3 and BP4 intervals.
        # Those four were not_met while an abstaining mechanism counted as a contradiction;
        # each now carries its abstention in provenance.uninformative_mechanisms and a review
        # point, so the benign call still says which mechanism was left open.
        self.assertEqual(len(met), 22)
        abstained = [result for result in met
                     if "uninformative_mechanisms" in result["provenance"]]
        self.assertEqual(len(abstained), 4)
        self.assertTrue(all(result["criterion"] == "BP4" and result["review_points"]
                            for result in abstained))
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
        # Unknown SpliceAI model version travels with the result as an explicit policy.
        assertions = {item["version_status"] for result in met
                      for item in result["provenance"].get("version_assertions", [])}
        self.assertEqual(assertions, {"UNKNOWN"})

        # The committed BA1 exception list resolves the high-frequency records either way.
        ba1 = [result for record in results["records"] for result in record["results"]
               if result["criterion"] == "BA1"]
        self.assertEqual(sum(result["status"] == "met" for result in ba1), 11)
        for result in ba1:
            if result["status"] == "met":
                assessment = result["provenance"]["ba1_exception_assessment"]
                self.assertFalse(assessment["is_exception"])
                self.assertEqual(assessment["list_size"], 9)
                self.assertEqual(result["evidence_outcome"], "BA1")

        # ClinVar missense density around each residue, under the committed PM1 policy.
        hotspot = providers["ClinVar protein hotspot density"]
        self.assertEqual(hotspot["policy_version"], "PM1-hotspot-v1")
        self.assertEqual(hotspot["region_evidence"], 10)
        self.assertEqual(hotspot["errors"], [])
        pm1 = [result for record in results["records"] for result in record["results"]
               if result["criterion"] == "PM1"]
        statuses = {status: sum(result["status"] == status for result in pm1)
                    for status in ("met", "not_met", "unknown")}
        # Every demo variant that ClinVar classifies is excluded from its own density, so no
        # record is NOT_MET on the strength of its own submitted classification.
        self.assertEqual(statuses, {"met": 2, "not_met": 0, "unknown": 26})
        # PM5 needs a residue-scoped search; MYH7 p.Arg719 has a pathogenic ClinVar comparator.
        pm5 = [result for record in results["records"] for result in record["results"]
               if result["criterion"] == "PM5"]
        pm5_statuses = {status: sum(result["status"] == status for result in pm5)
                        for status in ("met", "not_met", "unknown")}
        self.assertEqual(pm5_statuses, {"met": 1, "not_met": 17, "unknown": 10})
        pm5_met = next(result for result in pm5 if result["status"] == "met")
        self.assertEqual(pm5_met["evidence_outcome"], "PM5")
        self.assertEqual(pm5_met["provenance"]["assessment_scope"], "protein_level")
        # The committed record-level context now supplies HCM for this case.
        self.assertEqual(pm5_met["provenance"]["condition_assessment"], "MATCHED")

        met = [result for result in pm1 if result["status"] == "met"]
        for result in met:
            self.assertEqual(result["evidence_outcome"], "PM1")
            self.assertEqual(result["provenance"]["pm1_route"], "mutational_hotspot")
            self.assertEqual(result["provenance"]["assessment_method"], "automated")
            self.assertEqual(result["provenance"]["benign_count"], 0)
            self.assertGreaterEqual(result["provenance"]["pathogenic_count"], 3)
            # Disease relevance is unresolved without a condition, so it stays a review point.
            self.assertEqual(result["provenance"]["condition_assessment"], "NOT_EVALUATED")
            self.assertTrue(result["review_points"])

        # PP2 and BP1 need a mechanism assessment scoped to a transcript, which the online
        # providers above do not produce: ClinGen Dosage and Gene2Phenotype state whether
        # loss of function is a mechanism, which answers PVS1's gate and nothing else. The
        # transcript-scoped review that also records the missense mechanism and the variant
        # spectrum comes from build-gene-disease-evidence, and without that step PP2 and BP1
        # are unknown for every missense record in the demo set. This asserts the step is a
        # part of the pipeline rather than something a run has to know to add.
        code = main([
            "build-gene-disease-evidence",
            "--input", str(ROOT / "config" / "gene-disease-review-decisions.json"),
            "--base-evidence", str(prepared / "evidence.json"),
            "--output", str(mechanism_evidence),
        ])
        self.assertEqual(code, 0)
        code = main([
            "evaluate", "--input", str(prepared / "variants.json"),
            "--evidence", str(mechanism_evidence), "--offline",
            "--config", str(ROOT / "config" / "demo-rules.json"),
            "--context", str(ROOT / "config" / "curated-context.json"),
            "--criteria", "PP2,BP1", "--output-dir", str(mechanism_evaluated),
            "--internal-only",
        ])
        self.assertEqual(code, 0)
        mechanism_results = json.loads(
            (mechanism_evaluated / "results.json").read_text(encoding="utf-8")
        )
        mechanism = [result for record in mechanism_results["records"]
                     for result in record["results"]]
        counts = {code_: {status: sum(result["status"] == status for result in mechanism
                                      if result["criterion"] == code_)
                          for status in ("met", "not_met", "unknown")}
                  for code_ in ("PP2", "BP1")}
        # 18 of the 28 demo records are missense, but only 10 are evaluated. The other 8 are
        # missense in a gene-disease pair the review marks the criterion inapplicable to -
        # the ClinGen Cardiomyopathy VCEP's MYH7 specification (5 records), KCNJ5, whose
        # validity is for familial hyperaldosteronism rather than HCM, and MYBPC3 and TNNI3
        # against ARVC, where validity is Limited or is for HCM instead. Those join the 10
        # non-missense records as unknown, so an unknown here is not one situation.
        self.assertEqual(counts["PP2"], {"met": 1, "not_met": 9, "unknown": 18})
        self.assertEqual(counts["BP1"], {"met": 0, "not_met": 10, "unknown": 18})
        inapplicable = [result for result in mechanism if result["status"] == "unknown"
                        and result["provenance"].get("applicability_source")]
        self.assertEqual(len(inapplicable), 16)
        evaluated_records = [result for result in mechanism if result["status"] != "unknown"]
        self.assertTrue(all(result["provenance"]["assessment_scope"] == "gene_level"
                            for result in evaluated_records))
        # The step adds to the prepared evidence rather than replacing it, so the same
        # document still carries the population, annotation and ClinVar records the other
        # criteria read. Whether the gene-level records from --with-clingen-dosage can sit
        # alongside these without being counted as a competing assessment is covered by
        # tests/test_curated.py; this recipe does not request those providers.
        combined = json.loads(mechanism_evidence.read_text(encoding="utf-8"))["evidence"]
        base = json.loads((prepared / "evidence.json").read_text(encoding="utf-8"))["evidence"]
        # 7 reviewed gene-disease groups expand to 15 of the 28 demo records - a record in a
        # gene with no reviewed decision gets no assessment rather than a generalized one.
        self.assertEqual(len(combined), len(base) + 15)
        self.assertTrue(all(item.get("transcript")
                            for item in combined if item["category"] == "gene_disease"))

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
        self.assertEqual(by_id["case2:24:1"]["status"], "met")
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
        # callability evidence, so it must not be silently converted to a confirmed AF=0 -
        # but (2026-09-17 policy change) PM2 now scores that total absence as a flagged
        # supporting-strength prediction rather than withholding a call outright.
        self.assertEqual(by_id["case2:25:1"]["status"], "met")
        self.assertEqual(by_id["case2:25:1"]["evidence_outcome"], "PM2_supporting")
        self.assertTrue(by_id["case2:25:1"]["review_points"])


if __name__ == "__main__":
    unittest.main()
