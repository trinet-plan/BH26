"""An unapproved maximum credible frequency may not carry a final judgment.

Recorded 2026-09-17 in doc/external_data_coverage_ja.md: a threshold derived from published
information, and not approved by a curator or VCEP for the disease at hand, is not a basis for
an ACMG/AMP verdict. Before this, a real API run of MYH7 c.4472C>G reached "likely benign" on
BS1 alone, and the threshold BS1 used was the configured default - labelled PLACEHOLDER in its
own source string.

The comparison still runs. Its result is handed to a curator as a draft, with the derivation
inputs it does have and the ones it does not.
"""

import unittest
from types import SimpleNamespace

from acmg_pipeline.automated_core.interface import inputs_from_prepared_record
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.criteria import bs1
from acmg_pipeline.providers.local import LocalPopulationProvider
from acmg_pipeline.services.population import PopulationService

CURATED = {"condition": "test:disease", "inheritance": "autosomal_dominant",
           "source": "ClinGen VCEP", "source_version": "1.0", "reviewed_at": "2026-01-01",
           "max_credible_af": 0.0001, "frequency_statistic": "faf95", "comparison": ">"}


class ThresholdPolicyTests(unittest.TestCase):

    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.input = {"variant": self.variant.to_dict()}
        self.config = {"BS1": {
            "minimum_an": 2000, "policy_source": "test", "policy_version": "1",
            "default_max_credible_af": 0.0001, "default_source": "test default",
            "default_source_version": "1", "default_frequency_statistic": "faf95",
            "default_comparison": ">",
            "default_derivation": {"method": "Whiffin et al. 2017", "prevalence": None,
                                   "inheritance": None, "penetrance": None,
                                   "max_genetic_contribution": None,
                                   "max_allelic_contribution": None,
                                   "target_population": None, "source_url": "https://example"},
        }}

    def observation(self, ac, an=10000):
        return {"variant_key": self.variant.key, "evidence_id": "test:frequency",
                "source": "synthetic", "source_version": "1", "population": "TEST",
                "retrieved_at": "2026-09-14T00:00:00Z", "AC": ac, "AN": an, "AF": ac / an,
                "quality_status": "PASS", "callable": True}

    def evaluate(self, ac=100):
        services = SimpleNamespace(population=PopulationService(
            [LocalPopulationProvider("test", [self.observation(ac)])]))
        variant, note = inputs_from_prepared_record(self.input)
        return bs1.evaluate(variant, note, services, self.config)

    def curated(self, **overrides):
        self.input["condition"] = "test:disease"
        self.input["inheritance"] = "autosomal_dominant"
        self.input["disease_frequency_threshold"] = {**CURATED, **overrides}

    # --- the default threshold ------------------------------------------------------

    def test_an_exceeded_default_threshold_is_not_met(self):
        value = self.evaluate()
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIsNone(value.strength)
        self.assertEqual(value.provenance["policy_status"], "DRAFT")

    def test_the_summary_says_why_it_stopped_short_of_met(self):
        self.assertIn("not approved for this disease context", self.evaluate().summary)

    def test_the_comparison_that_ran_is_handed_over_not_discarded(self):
        draft = self.evaluate().provenance["draft_threshold_candidate"]
        self.assertEqual(draft["policy_status"], "DRAFT")
        self.assertEqual(draft["max_credible_af"], "0.0001")
        self.assertEqual(draft["frequency_statistic"], "faf95")
        self.assertEqual(draft["comparison"], ">")
        self.assertTrue(draft["review_reason"])

    def test_the_draft_carries_the_observations_it_was_computed_from(self):
        observation = self.evaluate().provenance["draft_threshold_candidate"]["observations"][0]
        for field in ("evidence_id", "source", "source_version", "population", "AC", "AN", "AF"):
            self.assertIsNotNone(observation[field], field)
        self.assertIsNotNone(observation["faf95"])

    def test_missing_derivation_inputs_are_named_never_filled_in(self):
        """A blank is a question for the curator, not a value for the pipeline to choose."""
        draft = self.evaluate().provenance["draft_threshold_candidate"]
        self.assertEqual(sorted(draft["derivation_inputs_missing"]),
                         ["inheritance", "max_allelic_contribution", "max_genetic_contribution",
                          "penetrance", "prevalence", "target_population"])
        self.assertEqual(sorted(draft["derivation_inputs"]), ["method", "source_url"])

    def test_the_curator_is_asked_to_approve_before_classification(self):
        self.assertIn("Approve or replace the draft maximum credible frequency before "
                      "BS1 is used in a classification", self.evaluate().review_points)

    def test_below_a_default_threshold_is_still_not_met(self):
        """Not-met adds no evidence either way, so it does not rest on the unapproved value."""
        value = self.evaluate(ac=0)
        self.assertEqual(value.status, CriterionStatus.NOT_MET)
        self.assertNotIn("draft_threshold_candidate", value.provenance)

    # --- an approved, disease-specific threshold -------------------------------------

    def test_a_reviewed_disease_specific_threshold_still_reaches_met(self):
        self.curated()
        value = self.evaluate()
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.strength, "strong")
        self.assertEqual(value.provenance["policy_status"], "APPROVED")
        self.assertNotIn("draft_threshold_candidate", value.provenance)

    def test_a_curated_threshold_may_declare_itself_still_draft(self):
        """A curator part-way through a specification can mark it, and be honoured."""
        self.curated(policy_status="DRAFT")
        value = self.evaluate()
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.provenance["policy_status"], "DRAFT")
        self.assertIn("draft_threshold_candidate", value.provenance)


if __name__ == "__main__":
    unittest.main()
