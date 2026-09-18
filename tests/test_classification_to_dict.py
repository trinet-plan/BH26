"""classification_to_dict() - JSON-ready ClassificationResult.

Added 2026-09-18 so a va_spec output file can carry the overall category
alongside its 28 EvidenceLines, not just the per-criterion inputs.
"""

import unittest

from acmg_pipeline.classification import (
    ClassificationCategory, CriterionEvidence, CriterionStatus, Strength,
    classification_to_dict, classify,
)


class ClassificationToDictTests(unittest.TestCase):
    def test_met_not_met_and_unevaluated_codes_all_round_trip(self):
        result = classify([
            CriterionEvidence("PVS1", CriterionStatus.MET, Strength.VERY_STRONG, "pipeline"),
            CriterionEvidence("PM2", CriterionStatus.MET, Strength.SUPPORTING, "pipeline"),
            CriterionEvidence("BA1", CriterionStatus.NOT_MET, source="pipeline"),
        ])
        as_dict = classification_to_dict(result)
        self.assertEqual(as_dict["category"], ClassificationCategory.LIKELY_PATHOGENIC.value)
        self.assertEqual(as_dict["score"], 9)
        self.assertFalse(as_dict["ba1_override"])
        self.assertEqual(
            {e["code"] for e in as_dict["met"]}, {"PVS1", "PM2"})
        self.assertEqual(
            next(e for e in as_dict["met"] if e["code"] == "PVS1")["strength"],
            "very_strong")
        self.assertEqual(
            {e["code"] for e in as_dict["not_met"]}, {"BA1"})
        self.assertIsNone(
            next(e for e in as_dict["not_met"] if e["code"] == "BA1")["strength"])
        self.assertIn("PS1", as_dict["not_evaluated_codes"])
        self.assertNotIn("PVS1", as_dict["not_evaluated_codes"])

    def test_ba1_override_is_carried_through(self):
        result = classify([
            CriterionEvidence("BA1", CriterionStatus.MET, Strength.STAND_ALONE, "pipeline"),
        ])
        as_dict = classification_to_dict(result)
        self.assertEqual(as_dict["category"], ClassificationCategory.BENIGN.value)
        self.assertTrue(as_dict["ba1_override"])


if __name__ == "__main__":
    unittest.main()
