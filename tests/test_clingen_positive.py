"""ClinGen-derived positive references for every supported, non-deprecated criterion.

The ERepo applied codes live only in the expected-outcome manifest.  The evaluator receives
the prepared variants and normalized independent evidence, so an expert-panel label cannot
become its own supporting evidence.
"""

import json
import unittest
import uuid
from pathlib import Path

from acmg.output import run_internal


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
SUPPORTED = (
    "PVS1", "PS1", "PM1", "PM2", "PM4", "PM5", "PP2", "PP3",
    "BA1", "BS1", "BP1", "BP3", "BP4", "BP7",
)


def load(name):
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


class ClinGenPositiveReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prepared = load("clingen-positive-prepared.json")
        cls.evidence = load("clingen-positive-evidence.json")
        cls.expected = load("clingen-positive-expected.json")
        cls.output = ROOT / ".work" / "clingen-positive" / uuid.uuid4().hex
        cls.payload = run_internal(
            FIXTURES / "clingen-positive-prepared.json",
            FIXTURES / "clingen-positive-evidence.json",
            FIXTURES / "clingen-positive-rules.json",
            cls.output,
            criteria=SUPPORTED,
            va_spec=True,
        )
        cls.by_record = {
            record["record_id"]: {item["criterion"]: item for item in record["results"]}
            for record in cls.payload["records"]
        }

    def test_manifest_covers_every_supported_non_deprecated_criterion(self):
        self.assertEqual(tuple(self.expected["criteria"]), SUPPORTED)
        self.assertNotIn("PP5", self.expected["criteria"])
        self.assertNotIn("BP6", self.expected["criteria"])

    def test_every_criterion_has_a_clingen_positive_met_result(self):
        for criterion, expectation in self.expected["criteria"].items():
            with self.subTest(criterion=criterion):
                result = self.by_record[expectation["record_id"]][criterion]
                self.assertEqual(result["status"], "MET")
                self.assertEqual(result["strength"], expectation["strength"])
                self.assertTrue(result["evidence"])

    def test_expected_labels_are_not_evaluator_inputs(self):
        evaluator_inputs = json.dumps(
            {"prepared": self.prepared, "evidence": self.evidence}, sort_keys=True
        )
        for forbidden in ("expert_panel_codes", "metCodes", "source_met_code", "expected_status"):
            self.assertNotIn(forbidden, evaluator_inputs)
        for criterion, expectation in self.expected["criteria"].items():
            self.assertTrue(expectation["source_met_code"].startswith(criterion))

    def test_every_reference_has_stable_clingen_provenance(self):
        for record in self.prepared["records"]:
            source = record["source"]
            self.assertEqual(source["dataset"], "ClinGen Evidence Repository")
            self.assertTrue(source["url"].startswith((
                "https://erepo.clinicalgenome.org/evrepo/ui/classification/",
                "https://erepo.clinicalgenome.org/evrepo/ui/interpretation/",
            )))
            self.assertTrue(source["url"].endswith(source["uuid"]))
            self.assertTrue(source["caid"].startswith("CAR:CA"))
            self.assertTrue(source["preferred_variant_title"])

    def test_va_spec_export_contains_each_positive_result(self):
        self.assertEqual(self.payload["input_errors"], [])
        self.assertEqual(self.payload["va_spec_export"], "VALIDATED")
        envelope = json.loads((self.output / "evidence-lines.json").read_text(encoding="utf-8"))
        exported = {
            (record["record_id"], wrapped["criterion"])
            for record in envelope["records"]
            for wrapped in record["evidence_lines"]
        }
        for criterion, expectation in self.expected["criteria"].items():
            self.assertIn((expectation["record_id"], criterion), exported)

        for record in envelope["records"]:
            self.assertEqual(
                {item["criterion"] for item in record["criterion_assessments"]},
                set(SUPPORTED),
            )
            catalog = record["referenced_evidence"]
            for assessment in record["criterion_assessments"]:
                self.assertTrue(set(assessment["evidenceItemIds"]) <= set(catalog))
            for wrapped in record["evidence_lines"]:
                extension = next(
                    item for item in wrapped["evidence_line"]["extensions"]
                    if item["name"] == "bh26AssessmentDetails"
                )
                self.assertEqual(extension["value"], wrapped["assessment_details"])

    def test_pvs1_reference_preserves_trace_and_rules(self):
        result = self.by_record["clingen-positive:pvs1-pah"]["PVS1"]
        self.assertEqual(result["evaluation_context"]["condition_status"], "PROVIDED")
        self.assertEqual(result["evaluation_context"]["mechanism_scope"], "CONDITION_SPECIFIC")
        self.assertTrue(result["decision_trace"])
        self.assertTrue(result["rules_used"])
        self.assertEqual(result["missing_inputs"], result["unresolved_requirements"])


if __name__ == "__main__":
    unittest.main()
