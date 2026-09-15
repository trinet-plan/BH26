import unittest
from types import SimpleNamespace

from acmg.core.models import Status, Variant
from acmg.criteria import pp3, bp4, pp5, bp6
from acmg.services.evidence import EvidenceService


class ComputationalTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.input = {"variant": self.variant.to_dict(), "transcript": "NM_TEST.1"}
        common = {"variant_key": self.variant.key, "transcript": "NM_TEST.1",
                  "source": "synthetic", "source_version": "1", "retrieved_at": "2026-09-15",
                  "quality_status": "PASS"}
        self.annotation = {**common, "evidence_id": "test:annotation", "category": "annotation",
                           "consequences": ["missense_variant"]}
        self.prediction = {**common, "evidence_id": "test:prediction", "category": "computational",
                           "predictor": "test-predictor", "predictor_version": "1",
                           "mechanism": "protein", "score": 0.8}
        self.config = {"computational": {"selected_calibration": "test", "calibrations": {
            "test": {"predictor": "test-predictor", "predictor_version": "1", "source": "synthetic",
                     "version": "1", "mechanism": "protein", "consequences": ["missense_variant"],
                     "score_min": 0, "score_max": 1,
                     "bands": {"PP3": [{"min": 0.8, "max": 1, "strength": "supporting"},
                                       {"min": 0.9, "max": 1, "strength": "strong"}],
                               "BP4": [{"min": 0, "max": 0.2, "strength": "supporting"}]}}}}}

    def run_rule(self, module):
        services = SimpleNamespace(evidence=EvidenceService([self.annotation, self.prediction]))
        return module.evaluate(self.input, services, self.config)

    def test_calibrated_boundary(self):
        for score, expected in ((0.799, Status.NOT_MET), (0.8, Status.MET), (0.801, Status.MET)):
            self.prediction["score"] = score
            self.assertEqual(self.run_rule(pp3).status, expected)
        self.prediction["score"] = 0.95
        self.assertEqual(self.run_rule(pp3).strength, "strong")

    def test_other_predictor_cannot_use_threshold(self):
        self.prediction["predictor"] = "AlphaGenome"
        self.assertEqual(self.run_rule(pp3).status, Status.NOT_EVALUATED)

    def test_non_missense_not_applicable(self):
        self.annotation["consequences"] = ["frameshift_variant"]
        self.assertEqual(self.run_rule(pp3).status, Status.NOT_APPLICABLE)

    def test_benign_calibration(self):
        self.prediction["score"] = 0.1
        self.assertEqual(self.run_rule(bp4).status, Status.MET)
        self.assertEqual(self.run_rule(pp3).status, Status.NOT_MET)
        self.annotation["high_confidence_null_or_splice"] = True
        self.assertEqual(self.run_rule(bp4).status, Status.MANUAL_REVIEW)

    def test_no_calibration_is_not_evaluated(self):
        self.config = {}
        self.assertEqual(self.run_rule(pp3).status, Status.NOT_EVALUATED)

    def test_explicitly_uncalibrated_source_is_not_scored(self):
        self.prediction["calibration_eligible"] = False
        self.assertEqual(self.run_rule(pp3).status, Status.NOT_EVALUATED)

    def test_deprecated_independent_of_labels(self):
        self.input["CLNSIG"] = "Pathogenic"
        for module in (pp5, bp6):
            result = self.run_rule(module)
            self.assertEqual(result.status, Status.DEPRECATED)
            self.assertIsNone(result.evidence_outcome)
