from acmg_pipeline.constants import CriterionStatus
import unittest
from types import SimpleNamespace

from acmg_pipeline.automated_core.interface import inputs_from_prepared_record
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria import pp3, bp4, pp5, bp6
from acmg_pipeline.services.evidence import EvidenceService


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
        variant, clinical_note = inputs_from_prepared_record(self.input)
        return module.evaluate(variant, clinical_note, services, self.config)

    def test_calibrated_boundary(self):
        for score, expected in ((0.799, CriterionStatus.NOT_MET), (0.8, CriterionStatus.MET), (0.801, CriterionStatus.MET)):
            self.prediction["score"] = score
            self.assertEqual(self.run_rule(pp3).status, expected)
        self.prediction["score"] = 0.95
        self.assertEqual(self.run_rule(pp3).strength, "strong")

    def test_other_predictor_cannot_use_threshold(self):
        self.prediction["predictor"] = "AlphaGenome"
        self.assertEqual(self.run_rule(pp3).status, CriterionStatus.UNKNOWN)

    def test_non_missense_not_applicable(self):
        self.annotation["consequences"] = ["frameshift_variant"]
        self.assertEqual(self.run_rule(pp3).status, CriterionStatus.UNKNOWN)

    def test_benign_calibration(self):
        self.prediction["score"] = 0.1
        self.assertEqual(self.run_rule(bp4).status, CriterionStatus.MET)
        self.assertEqual(self.run_rule(pp3).status, CriterionStatus.NOT_MET)
        self.annotation["high_confidence_null_or_splice"] = True
        self.assertEqual(self.run_rule(bp4).status, CriterionStatus.UNKNOWN)

    def test_no_calibration_is_not_evaluated(self):
        self.config = {}
        self.assertEqual(self.run_rule(pp3).status, CriterionStatus.UNKNOWN)

    def test_explicitly_uncalibrated_source_is_not_scored(self):
        self.prediction["calibration_eligible"] = False
        self.assertEqual(self.run_rule(pp3).status, CriterionStatus.UNKNOWN)

    def test_deprecated_independent_of_labels(self):
        self.input["CLNSIG"] = "Pathogenic"
        for module in (pp5, bp6):
            result = self.run_rule(module)
            self.assertEqual(result.status, CriterionStatus.UNKNOWN)
            self.assertIsNone(result.evidence_outcome)


class SplicingCalibrationTests(unittest.TestCase):
    """Protein and splicing predictions are one line of evidence, not two that add up."""

    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.input = {"variant": self.variant.to_dict(), "transcript": "NM_TEST.1"}
        common = {"variant_key": self.variant.key, "transcript": "NM_TEST.1",
                  "source": "synthetic", "source_version": "1", "retrieved_at": "2026-09-15",
                  "quality_status": "PASS"}
        self.annotation = {**common, "evidence_id": "test:annotation", "category": "annotation",
                           "consequences": ["missense_variant"]}
        self.protein = {**common, "evidence_id": "test:protein", "category": "computational",
                        "predictor": "REVEL", "predictor_version": "dbNSFP-4.8a",
                        "mechanism": "protein", "score": 0.1}
        self.splicing = {**common, "evidence_id": "test:splicing", "category": "computational",
                         "source": "Ensembl VEP SpliceAI", "predictor": "SpliceAI",
                         "predictor_version": "unreported-Ensembl-116", "mechanism": "splicing",
                         "score": 0.0, "calibration_eligible": False}
        self.assertion = {
            "source": "Ensembl VEP SpliceAI", "unreported_version": "unreported-Ensembl-116",
            "asserted_version": "SpliceAI served by Ensembl VEP release 116",
            "asserted_by": "test", "justification": "test"}
        self.config = {"computational": {
            "selected_calibrations": ["protein", "splicing"],
            "calibrations": {
                "protein": {"predictor": "REVEL", "predictor_version": "dbNSFP-4.8a",
                            "source": "synthetic", "version": "1", "mechanism": "protein",
                            "consequences": ["missense_variant"], "score_min": 0, "score_max": 1,
                            "bands": {"PP3": [{"min": 0.644, "max": 1, "strength": "supporting"}],
                                      "BP4": [{"min": 0, "max": 0.29, "strength": "supporting"}]}},
                "splicing": {"predictor": "SpliceAI", "predictor_version": "SpliceAI-1.3",
                             "source": "synthetic", "version": "1", "mechanism": "splicing",
                             "consequences": ["missense_variant", "synonymous_variant"],
                             "score_min": 0, "score_max": 1,
                             "version_assertion": self.assertion,
                             "bands": {"PP3": [{"min": 0.2, "max": 1, "strength": "supporting"}],
                                       "BP4": [{"min": 0, "max": 0.1, "strength": "supporting"}]}}}}}

    def run_rule(self, module, *extra):
        services = SimpleNamespace(evidence=EvidenceService(
            [self.annotation, self.protein, self.splicing, *extra]))
        variant, clinical_note = inputs_from_prepared_record(self.input)
        return module.evaluate(variant, clinical_note, services, self.config)

    def test_both_mechanisms_are_reported(self):
        value = self.run_rule(bp4)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual({item["mechanism"] for item in value.provenance["applied_calibrations"]},
                         {"protein", "splicing"})

    def test_predicted_splice_impact_blocks_benign_evidence(self):
        """A low protein score says nothing about a disrupted splice site."""
        self.splicing["score"] = 0.9
        value = self.run_rule(bp4)
        self.assertEqual(value.status, CriterionStatus.NOT_MET)
        self.assertEqual(value.provenance["benign_blocked_by"][0]["predictor"], "SpliceAI")

    def test_splicing_alone_can_carry_pp3(self):
        self.splicing["score"] = 0.9
        value = self.run_rule(pp3)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.strength, "supporting")
        self.assertEqual(value.provenance["applied_calibrations"][1]["mechanism"], "splicing")

    def test_strengths_are_not_summed_across_mechanisms(self):
        self.protein["score"] = 0.95
        self.splicing["score"] = 0.9
        self.config["computational"]["calibrations"]["protein"]["bands"]["PP3"].append(
            {"min": 0.932, "max": 1, "strength": "strong"})
        value = self.run_rule(pp3)
        self.assertEqual(value.strength, "strong")

    def test_synonymous_variant_uses_the_splicing_calibration_only(self):
        self.annotation["consequences"] = ["synonymous_variant"]
        value = self.run_rule(bp4)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual([item["predictor"] for item in value.provenance["applied_calibrations"]],
                         ["SpliceAI"])
        # A calibration outside its consequence scope is not a missing predictor.
        self.assertNotIn("unavailable_predictors", value.provenance)

    def test_unversioned_score_needs_a_declared_assumption(self):
        del self.config["computational"]["calibrations"]["splicing"]["version_assertion"]
        value = self.run_rule(bp4)
        self.assertEqual(value.status, CriterionStatus.MET)
        # Without the assumption the splicing score is simply not used.
        self.assertEqual([item["predictor"] for item in value.provenance["applied_calibrations"]],
                         ["REVEL"])
        self.assertEqual(value.provenance["unavailable_predictors"], ["SpliceAI"])

    def test_the_assumption_is_recorded_with_the_result(self):
        value = self.run_rule(bp4)
        self.assertEqual(value.provenance["version_assertions"], [self.assertion])

    def test_an_incomplete_assumption_is_refused(self):
        self.config["computational"]["calibrations"]["splicing"]["version_assertion"] = {
            "source": "Ensembl VEP SpliceAI", "unreported_version": "unreported-Ensembl-116"}
        value = self.run_rule(bp4)
        self.assertEqual([item["predictor"] for item in value.provenance["applied_calibrations"]],
                         ["REVEL"])
'''

p = 'tests/test_computational.py'
s = io.open(p, encoding='utf-8').read()
tail = '''

if __name__ == "__main__":
    unittest.main()


if __name__ == "__main__":
    unittest.main()

