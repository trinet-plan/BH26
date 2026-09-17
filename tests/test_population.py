from acmg_pipeline.constants import CriterionStatus
import unittest
from types import SimpleNamespace

from acmg_pipeline.automated_core.interface import inputs_from_prepared_record
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria import ba1, bs1, pm2
from acmg_pipeline.providers.local import LocalPopulationProvider
from acmg_pipeline.services.population import PopulationService
from acmg_pipeline.automated_va_spec import to_evidence_line, validate_1_0_1


class PopulationTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.input = {"variant": self.variant.to_dict()}
        self.config = {code: {"minimum_an": 2000, "policy_source": "synthetic-test-policy",
                              "policy_version": "1"} for code in ("PM2", "BA1", "BS1")}
        # "Absent from controls" is stated, not defaulted: pm2.evaluate() reports an unset
        # max_af as an unconfigured policy rather than running as the strictest threshold.
        self.config["PM2"]["max_af"] = 0

    def observation(self, ac=0, an=10000, **extra):
        return {"variant_key": self.variant.key, "evidence_id": "test:frequency",
                "source": "synthetic", "source_version": "1", "population": "TEST",
                "retrieved_at": "2026-09-14T00:00:00Z", "AC": ac, "AN": an,
                "AF": ac / an, "quality_status": "PASS", "callable": True, **extra}

    def services(self, *batches):
        return SimpleNamespace(population=PopulationService([
            LocalPopulationProvider(f"test-{i}", batch) for i, batch in enumerate(batches)]))

    def evaluate(self, module, services):
        variant, clinical_note = inputs_from_prepared_record(self.input)
        return module.evaluate(variant, clinical_note, services, self.config)

    def test_absence_is_supporting(self):
        value = self.evaluate(pm2, self.services([self.observation()]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evidence_outcome, "PM2_supporting")

    def test_no_record_is_not_absence(self):
        value = self.evaluate(pm2, self.services([]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)

    def test_unset_threshold_is_not_the_strictest_threshold(self):
        """An absent max_af used to default to 0 and score PM2 on an unstated policy."""
        del self.config["PM2"]["max_af"]
        value = self.evaluate(pm2, self.services([self.observation()]))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("PM2.max_af", value.missing_inputs)

    def test_configured_zero_threshold_still_scores(self):
        """max_af = 0 is the standard "absent from controls" policy, not a missing one."""
        self.config["PM2"]["max_af"] = 0
        value = self.evaluate(pm2, self.services([self.observation()]))
        self.assertEqual(value.status, CriterionStatus.MET)

    def test_high_regional_af_is_not_hidden(self):
        value = self.evaluate(pm2, self.services(
            [self.observation()], [self.observation(ac=100, population="REGIONAL")]
        ))
        self.assertEqual(value.status, CriterionStatus.NOT_MET)

    def test_unavailable_provider_prevents_rarity(self):
        value = self.evaluate(pm2, self.services([self.observation()], []))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)

    def test_excluded_small_subgroup_does_not_invalidate_usable_population(self):
        usable = self.observation()
        too_small = self.observation(an=100, population="SMALL")
        self.assertEqual(self.evaluate(pm2, self.services([usable, too_small])).status, CriterionStatus.MET)
        self.assertEqual(self.evaluate(ba1, self.services([usable, too_small])).status, CriterionStatus.NOT_MET)

    def test_low_an_and_unknown_callability(self):
        for item in (self.observation(an=100), self.observation(callable=False),
                     self.observation(AF="NaN"), self.observation(AF=0.3)):
            with self.subTest(item=item):
                value = self.evaluate(pm2, self.services([item]))
                self.assertEqual(value.status, CriterionStatus.UNKNOWN)

    def test_ba1_boundary_and_exception(self):
        self.input["ba1_exception_assessment"] = {
            "source": "test", "source_version": "1", "reviewed_at": "2026-09-14",
            "is_exception": False,
        }
        for ac, expected in ((499, CriterionStatus.NOT_MET), (500, CriterionStatus.NOT_MET), (501, CriterionStatus.MET)):
            value = self.evaluate(ba1, self.services([self.observation(ac)]))
            self.assertEqual(value.status, expected)
        self.input["ba1_exception_assessment"]["is_exception"] = True
        self.assertEqual(
            self.evaluate(ba1, self.services([self.observation(501)])).status,
            CriterionStatus.NOT_MET,
        )

    def test_ba1_exception_list_is_policy_not_an_evidence_item(self):
        self.input["ba1_exception_assessment"] = {
            "source": "ClinGen SVI BA1 exception list", "source_version": "2018",
            "reviewed_at": "2026-09-15", "is_exception": False}
        value = self.evaluate(ba1, self.services([self.observation(501)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertTrue(all(item.get("evidence_id") for item in value.evidence))
        self.assertFalse(value.provenance["ba1_exception_assessment"]["is_exception"])
        line = to_evidence_line(value)
        validate_1_0_1(line, "BA1")
        self.assertEqual(line["directionOfEvidenceProvided"], "disputes")

    def test_bs1_requires_matching_disease(self):
        self.assertEqual(
            self.evaluate(bs1, self.services([self.observation(100)])).status,
            CriterionStatus.UNKNOWN,
        )
        self.input.update({"condition": "test:disease", "inheritance": "autosomal_dominant",
                           "disease_frequency_threshold": {
                               "condition": "test:disease", "inheritance": "autosomal_dominant",
                               "max_credible_af": 0.001, "source": "test", "source_version": "1",
                               "reviewed_at": "2026-09-14"}})
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.provenance["disease_frequency_threshold"]["max_credible_af"], 0.001)

    def test_bs1_threshold_is_policy_not_an_evidence_item(self):
        """Every exported evidence item needs an IRI; the threshold has none, so it stays out."""
        self.input.update({"condition": "test:disease", "inheritance": "autosomal_dominant",
                           "disease_frequency_threshold": {
                               "condition": "test:disease", "inheritance": "autosomal_dominant",
                               "max_credible_af": 0.001, "source": "test", "source_version": "1",
                               "reviewed_at": "2026-09-14"}})
        value = self.evaluate(bs1, self.services([self.observation(100)]))
        self.assertTrue(value.evidence)
        self.assertTrue(all(item.get("evidence_id") for item in value.evidence))
        line = to_evidence_line(value)
        validate_1_0_1(line, "BS1")
        self.assertEqual(line["evidenceOutcome"]["primaryCoding"]["code"], "BS1")
        self.assertEqual(line["directionOfEvidenceProvided"], "disputes")


if __name__ == "__main__":
    unittest.main()

