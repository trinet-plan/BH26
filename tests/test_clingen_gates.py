"""Applicability gates checked against expert-panel variants of each consequence type.

The demo VCFs contain no in-frame indel and no stop loss, so PM4 and BP3 are NOT_APPLICABLE
for all 28 records and their gates are never exercised on real data. These five variants come
from ClinGen Evidence Repository interpretations, one per consequence class, and are here to
show that each criterion opens for the variant types it covers and stays closed for the rest.
They carry no curated evidence, so an opened gate stops at NOT_EVALUATED by design.
"""

import json
import unittest
import uuid
from pathlib import Path

from acmg.output import run_internal


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
# Which criterion each record was selected to exercise, and what its consequence is.
GATES = {
    "clingen:pvs1-splice": ("PVS1", "splice_acceptor_variant"),
    "clingen:pvs1-frameshift": ("PVS1", "frameshift_variant"),
    "clingen:pm4-inframe-del": ("PM4", "inframe_deletion"),
    "clingen:bp7-synonymous-mybpc3": ("BP7", "synonymous_variant"),
    "clingen:bp7-synonymous-myh7": ("BP7", "synonymous_variant"),
}
CONSEQUENCE_GATED = ("PVS1", "PM4", "BP3", "BP7")


class ClinGenGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        output = ROOT / ".work" / "clingen-gates" / uuid.uuid4().hex
        cls.payload = run_internal(
            FIXTURES / "clingen-gate-prepared.json",
            FIXTURES / "clingen-gate-evidence.json",
            ROOT / "config" / "demo-rules.json",
            output, va_spec=True,
        )
        cls.by_record = {record["record_id"]: {item["criterion"]: item for item in record["results"]}
                         for record in cls.payload["records"]}

    def test_every_record_is_accepted(self):
        self.assertEqual(self.payload["input_errors"], [])
        self.assertEqual(len(self.payload["records"]), len(GATES))
        self.assertEqual(self.payload["va_spec_export"], "VALIDATED")

    def test_each_gate_opens_for_its_consequence(self):
        for record_id, (criterion, _) in GATES.items():
            result = self.by_record[record_id][criterion]
            self.assertNotEqual(result["status"], "NOT_APPLICABLE",
                                f"{criterion} should apply to {record_id}")
            # No curated evidence is supplied, so an opened gate stops here.
            self.assertEqual(result["status"], "NOT_EVALUATED")
            self.assertTrue(result["missing_inputs"])

    def test_an_in_frame_deletion_opens_both_length_and_repeat_criteria(self):
        results = self.by_record["clingen:pm4-inframe-del"]
        for criterion in ("PM4", "BP3"):
            self.assertEqual(results[criterion]["status"], "NOT_EVALUATED")
            self.assertEqual(results[criterion]["missing_inputs"], ["region"])
        self.assertEqual(results["PVS1"]["status"], "NOT_APPLICABLE")

    def test_gates_stay_closed_for_other_consequences(self):
        for record_id, (criterion, _) in GATES.items():
            for other in CONSEQUENCE_GATED:
                if other == criterion or (criterion == "PM4" and other == "BP3"):
                    continue
                self.assertEqual(self.by_record[record_id][other]["status"], "NOT_APPLICABLE",
                                 f"{other} should not apply to {record_id}")

    def test_loss_of_function_gate_covers_splice_and_frameshift(self):
        for record_id in ("clingen:pvs1-splice", "clingen:pvs1-frameshift"):
            result = self.by_record[record_id]["PVS1"]
            self.assertEqual(result["missing_inputs"], ["loss-of-function disease mechanism"])
            self.assertIsNone(result["strength"])

    def test_expert_panel_provenance_travels_with_each_record(self):
        records = json.loads((FIXTURES / "clingen-gate-prepared.json").read_text(encoding="utf-8"))
        for record in records["records"]:
            source = record["source"]
            self.assertTrue(source["caid"].startswith("CAR:CA"))
            self.assertTrue(source["url"].startswith("https://erepo.clinicalgenome.org/"))
            self.assertTrue(source["expert_panel_codes"])
            self.assertEqual(source["entry_method"],
                             "manual_transcription from the erepo API result")


if __name__ == "__main__":
    unittest.main()
