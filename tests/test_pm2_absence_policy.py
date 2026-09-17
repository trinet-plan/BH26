"""PM2's strength, and what may count as absence.

Recorded 2026-09-17 in doc/external_data_coverage_ja.md, following ClinGen SVI's PM2
recommendation v1.0 (2020-09-04): PM2 is applied at Supporting strength, and absence from a
public database supports rarity only where coverage and callability at the locus can be
checked. TogoVar returns no locus coverage, so a variant it does not report is not absent -
it is unknown.

Both already held when the policy was written. These pin them, because the failure mode is
silent: a defaulted strength or an AF of 0 read from a missing record would both still produce
a confident-looking result.
"""

import unittest
from types import SimpleNamespace

from acmg_pipeline.automated_core.interface import inputs_from_prepared_record
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.criteria import pm2
from acmg_pipeline.criteria.common import DEFAULT_STRENGTH
from acmg_pipeline.providers.local import LocalPopulationProvider
from acmg_pipeline.services.population import PopulationService, usable_observations


class Pm2PolicyTests(unittest.TestCase):

    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.input = {"variant": self.variant.to_dict()}
        self.config = {"PM2": {"minimum_an": 2000, "policy_source": "test",
                               "policy_version": "1", "max_af": 0}}

    def observation(self, **extra):
        return {"variant_key": self.variant.key, "evidence_id": "test:frequency",
                "source": "TogoVar", "source_version": "0.9.1", "population": "global",
                "retrieved_at": "2026-09-14T00:00:00Z", "AC": 0, "AN": 10000, "AF": 0,
                "quality_status": "PASS", **extra}

    def evaluate(self, *batches):
        services = SimpleNamespace(population=PopulationService(
            [LocalPopulationProvider(f"test-{i}", batch) for i, batch in enumerate(batches)]))
        variant, note = inputs_from_prepared_record(self.input)
        return pm2.evaluate(variant, note, services, self.config)

    def test_pm2_is_applied_at_supporting_strength(self):
        value = self.evaluate([self.observation(callable=True)])
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.strength, "supporting")

    def test_the_exported_outcome_says_supporting_not_bare_pm2(self):
        """DEFAULT_STRENGTH['PM2'] is the ACMG/AMP 2015 moderate, so the SVI downgrade has to
        survive into the outcome name a consumer reads."""
        self.assertEqual(DEFAULT_STRENGTH["PM2"], "moderate")
        self.assertEqual(self.evaluate([self.observation(callable=True)]).evidence_outcome,
                         "PM2_supporting")

    def test_a_database_that_reports_nothing_is_not_absence(self):
        """No record at all cannot be read as AF 0."""
        self.assertEqual(self.evaluate([]).status, CriterionStatus.UNKNOWN)

    def test_an_ac_of_zero_without_callability_is_rejected(self):
        value = self.evaluate([self.observation()])
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.provenance["rejected_observations"],
                         [{"evidence_id": "test:frequency",
                           "reasons": ["ABSENCE_WITHOUT_CALLABILITY"]}])

    def test_callability_must_be_stated_true_not_merely_present(self):
        for value in (False, "true", 1, None):
            with self.subTest(callable=value):
                _valid, rejected = usable_observations(
                    {"observations": [self.observation(callable=value)]}, self.variant, 2000)
                self.assertEqual(rejected[0]["reasons"], ["ABSENCE_WITHOUT_CALLABILITY"])

    def test_an_observed_allele_needs_no_callability(self):
        """Callability is what makes a zero count meaningful; a seen allele speaks for itself."""
        observed = self.observation(AC=1, AF=1 / 10000)
        self.assertEqual(self.evaluate([observed]).status, CriterionStatus.NOT_MET)


if __name__ == "__main__":
    unittest.main()
