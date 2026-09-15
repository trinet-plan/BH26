import unittest
from types import SimpleNamespace

from acmg.core.models import Variant, Status
from acmg.criteria import ba1, bs1, pm2
from acmg.providers.local import LocalPopulationProvider
from acmg.services.population import PopulationService
from acmg.va_spec.mapper import to_evidence_line, validate_1_0_1


class PopulationTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.input = {"variant": self.variant.to_dict()}
        self.config = {code: {"minimum_an": 2000, "policy_source": "synthetic-test-policy",
                              "policy_version": "1"} for code in ("PM2", "BA1", "BS1")}

    def observation(self, ac=0, an=10000, **extra):
        return {"variant_key": self.variant.key, "evidence_id": "test:frequency",
                "source": "synthetic", "source_version": "1", "population": "TEST",
                "retrieved_at": "2026-09-14T00:00:00Z", "AC": ac, "AN": an,
                "AF": ac / an, "quality_status": "PASS", "callable": True, **extra}

    def services(self, *batches):
        return SimpleNamespace(population=PopulationService([
            LocalPopulationProvider(f"test-{i}", batch) for i, batch in enumerate(batches)]))

    def test_absence_is_supporting(self):
        value = pm2.evaluate(self.input, self.services([self.observation()]), self.config)
        self.assertEqual(value.status, Status.MET)
        self.assertEqual(value.evidence_outcome, "PM2_supporting")

    def test_no_record_is_not_absence(self):
        value = pm2.evaluate(self.input, self.services([]), self.config)
        self.assertEqual(value.status, Status.NOT_EVALUATED)

    def test_high_regional_af_is_not_hidden(self):
        value = pm2.evaluate(self.input, self.services([self.observation()],
                              [self.observation(ac=100, population="REGIONAL")]), self.config)
        self.assertEqual(value.status, Status.NOT_MET)

    def test_unavailable_provider_prevents_rarity(self):
        value = pm2.evaluate(self.input, self.services([self.observation()], []), self.config)
        self.assertEqual(value.status, Status.NOT_EVALUATED)

    def test_excluded_small_subgroup_does_not_invalidate_usable_population(self):
        usable = self.observation()
        too_small = self.observation(an=100, population="SMALL")
        self.assertEqual(pm2.evaluate(self.input, self.services([usable, too_small]),
                                     self.config).status, Status.MET)
        self.assertEqual(ba1.evaluate(self.input, self.services([usable, too_small]),
                                     self.config).status, Status.NOT_MET)

    def test_low_an_and_unknown_callability(self):
        for item in (self.observation(an=100), self.observation(callable=False),
                     self.observation(AF="NaN"), self.observation(AF=0.3)):
            with self.subTest(item=item):
                value = pm2.evaluate(self.input, self.services([item]), self.config)
                self.assertEqual(value.status, Status.NOT_EVALUATED)

    def test_ba1_boundary_and_exception(self):
        self.input["ba1_exception_assessment"] = {
            "source": "test", "source_version": "1", "reviewed_at": "2026-09-14",
            "is_exception": False,
        }
        for ac, expected in ((499, Status.NOT_MET), (500, Status.NOT_MET), (501, Status.MET)):
            value = ba1.evaluate(self.input, self.services([self.observation(ac)]), self.config)
            self.assertEqual(value.status, expected)
        self.input["ba1_exception_assessment"]["is_exception"] = True
        self.assertEqual(ba1.evaluate(self.input, self.services([self.observation(501)]),
                                     self.config).status, Status.NOT_MET)

    def test_ba1_exception_list_is_policy_not_an_evidence_item(self):
        self.input["ba1_exception_assessment"] = {
            "source": "ClinGen SVI BA1 exception list", "source_version": "2018",
            "reviewed_at": "2026-09-15", "is_exception": False}
        value = ba1.evaluate(self.input, self.services([self.observation(501)]), self.config)
        self.assertEqual(value.status, Status.MET)
        self.assertTrue(all(item.get("evidence_id") for item in value.evidence))
        self.assertFalse(value.provenance["ba1_exception_assessment"]["is_exception"])
        line = to_evidence_line(value)
        validate_1_0_1(line, "BA1")
        self.assertEqual(line["directionOfEvidenceProvided"], "disputes")

    def test_bs1_requires_matching_disease(self):
        self.assertEqual(bs1.evaluate(self.input, self.services([self.observation(100)]),
                                     self.config).status, Status.NOT_EVALUATED)
        self.input.update({"condition": "test:disease", "inheritance": "autosomal_dominant",
                           "disease_frequency_threshold": {
                               "condition": "test:disease", "inheritance": "autosomal_dominant",
                               "max_credible_af": 0.001, "source": "test", "source_version": "1",
                               "reviewed_at": "2026-09-14"}})
        value = bs1.evaluate(self.input, self.services([self.observation(100)]), self.config)
        self.assertEqual(value.status, Status.MET)
        self.assertEqual(value.provenance["disease_frequency_threshold"]["max_credible_af"], 0.001)

    def test_bs1_threshold_is_policy_not_an_evidence_item(self):
        """Every exported evidence item needs an IRI; the threshold has none, so it stays out."""
        self.input.update({"condition": "test:disease", "inheritance": "autosomal_dominant",
                           "disease_frequency_threshold": {
                               "condition": "test:disease", "inheritance": "autosomal_dominant",
                               "max_credible_af": 0.001, "source": "test", "source_version": "1",
                               "reviewed_at": "2026-09-14"}})
        value = bs1.evaluate(self.input, self.services([self.observation(100)]), self.config)
        self.assertTrue(value.evidence)
        self.assertTrue(all(item.get("evidence_id") for item in value.evidence))
        line = to_evidence_line(value)
        validate_1_0_1(line, "BS1")
        self.assertEqual(line["evidenceOutcome"]["primaryCoding"]["code"], "BS1")
        self.assertEqual(line["directionOfEvidenceProvided"], "disputes")


if __name__ == "__main__":
    unittest.main()
