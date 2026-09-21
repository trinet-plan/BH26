"""BP7's `synonymous_assessment` evidence, derived from VEP consequence + the SpliceAI
score already fetched for every variant (ClinGen SVI 2023, PMID:37352859) - see
synonymous_assessment.py's own module docstring."""

import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.synonymous_assessment import (
    METHOD, SPLICEAI_NO_IMPACT_THRESHOLD, get_synonymous_assessment,
)


def annotation(consequences, transcript="NM_000256.3"):
    return {"transcript": transcript, "consequences": consequences}


def spliceai(score, transcript="NM_000256.3"):
    return [{"predictor": "SpliceAI", "transcript": transcript, "score": score,
             "evidence_id": "urn:test:spliceai", "source_version": "116",
             "retrieved_at": "2026-09-21T00:00:00Z", "source": "Ensembl VEP SpliceAI"}]


class SynonymousAssessmentTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "3", 10146618, "G", "T")

    def test_low_spliceai_score_outside_splice_region_is_met_eligible(self):
        record = get_synonymous_assessment(
            self.variant, annotation(["synonymous_variant"]), spliceai(0.02), "116")[0]
        self.assertTrue(record["outside_splice_critical_region"])
        self.assertTrue(record["no_predicted_splice_impact"])
        self.assertFalse(record["contradictory_rna_evidence"])

    def test_score_at_or_above_threshold_is_not_no_impact(self):
        record = get_synonymous_assessment(
            self.variant, annotation(["synonymous_variant"]), spliceai(SPLICEAI_NO_IMPACT_THRESHOLD),
            "116")[0]
        self.assertFalse(record["no_predicted_splice_impact"])

    def test_splice_region_co_annotation_is_inside_the_critical_window(self):
        record = get_synonymous_assessment(
            self.variant, annotation(["synonymous_variant", "splice_region_variant"]),
            spliceai(0.01), "116")[0]
        self.assertFalse(record["outside_splice_critical_region"])

    def test_not_synonymous_yields_nothing(self):
        records = get_synonymous_assessment(
            self.variant, annotation(["missense_variant"]), spliceai(0.01), "116")
        self.assertEqual(records, [])

    def test_no_spliceai_score_yields_nothing(self):
        records = get_synonymous_assessment(
            self.variant, annotation(["synonymous_variant"]), [], "116")
        self.assertEqual(records, [])

    def test_spliceai_scored_for_a_different_transcript_does_not_count(self):
        records = get_synonymous_assessment(
            self.variant, annotation(["synonymous_variant"], transcript="NM_OTHER.1"),
            spliceai(0.01, transcript="NM_000256.3"), "116")
        self.assertEqual(records, [])

    def test_no_annotation_yields_nothing(self):
        self.assertEqual(get_synonymous_assessment(self.variant, None, [], "116"), [])

    def test_the_record_declares_itself_automated(self):
        record = get_synonymous_assessment(
            self.variant, annotation(["synonymous_variant"]), spliceai(0.01), "116")[0]
        self.assertEqual(record["assessment_method"], "automated")
        self.assertEqual(record["method"], METHOD)
        self.assertNotIn("curator", record)
        self.assertTrue(reviewed_or_automated(record))
        self.assertTrue(record["splice_prediction_evidence"])
        self.assertTrue(record["calibration_source"])
        self.assertTrue(record["position_rule_version"])

    def test_a_policy_version_is_required(self):
        with self.assertRaises(ValueError):
            get_synonymous_assessment(self.variant, annotation(["synonymous_variant"]),
                                      spliceai(0.01), "")


if __name__ == "__main__":
    unittest.main()
