"""A criterion with no annotation must say why it has none.

The resolver already knew: it caught the provider error and recorded it in
ResolvedEvidence.failures. But only `records` was passed on to the services, so by the time a
criterion found itself without an annotation, the reason had been discarded and all it could
report was "Transcript annotation unavailable" - which tells a curator nothing they could not
already see. A whole API run of MYH7 c.4472C>G came back with eleven criteria unknown for that
reason and no way, afterwards, to learn what had failed.
"""

import tempfile
import unittest

from acmg_pipeline.automated_engine import evaluate_prepared_record, make_services
from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.criteria.common import (
    annotation_context, curated_context, retrieval_failures,
)

VARIANT = {"assembly": "GRCh38", "chrom": "14", "pos": 23417200, "ref": "G", "alt": "C"}
INPUT = {"variant": VARIANT, "transcript": "NM_000257.4"}
ENSEMBL_FAILED = [{"provider": "Ensembl", "error": "HTTPError: 429 Too Many Requests"}]


class ReleaseLookupTests(unittest.TestCase):
    """Building the provider suite asks Ensembl which release is current, and that lookup can
    fail on its own. It sat outside resolve()'s guards, so it escaped the whole call - taking
    the population evidence already gathered with it, and reaching the API as an
    ExceptionGroup with the cause buried."""

    def resolve(self):
        from acmg_pipeline.automated_core.models import Variant
        from acmg_pipeline.services.resolve import ProviderEvidenceResolver
        resolver = ProviderEvidenceResolver(tempfile.mkdtemp(), offline=True)
        return resolver.resolve({"GENE": "MYH7", "TRANSCRIPT": "NM_000257.4",
                                 "HGVSC": "c.4472C>G"}, Variant(**VARIANT))

    def test_a_failed_release_lookup_does_not_escape(self):
        self.resolve()  # previously raised FetchError

    def test_it_is_recorded_as_a_provider_failure(self):
        failures = self.resolve().failures
        self.assertTrue(any(item["provider"] == "Ensembl" for item in failures))

    def test_the_criteria_still_run_and_say_why(self):
        from acmg_pipeline.automated_core.models import Variant
        resolved = self.resolve()
        services = make_services(resolved.records, failures=resolved.failures)
        outcome, _ = annotation_context("PM4", INPUT, services)
        self.assertIn("Ensembl", outcome.summary)


class RetrievalFailureTests(unittest.TestCase):

    def test_services_built_without_failures_still_work(self):
        """Every existing caller passes records only; absence must not become an error."""
        self.assertEqual(retrieval_failures(make_services([])), [])

    def test_a_criterion_names_the_provider_and_its_error(self):
        outcome, annotation = annotation_context(
            "PM4", INPUT, make_services([], failures=ENSEMBL_FAILED))
        self.assertIsNone(annotation)
        self.assertEqual(outcome.status, CriterionStatus.UNKNOWN)
        self.assertIn("Ensembl", outcome.summary)
        self.assertIn("429 Too Many Requests", outcome.summary)

    def test_the_failure_is_kept_in_provenance_not_only_prose(self):
        """A summary is for reading; provenance is what a later run can count."""
        outcome, _ = annotation_context("PM4", INPUT, make_services([], failures=ENSEMBL_FAILED))
        self.assertEqual(outcome.provenance["provider_failures"], ENSEMBL_FAILED)

    def test_nothing_retrieved_and_nothing_failed_is_said_to_be_different(self):
        """An empty answer from a working provider is not a provider that broke."""
        outcome, _ = annotation_context("PM4", INPUT, make_services([]))
        self.assertIn("no provider reported an error", outcome.summary)
        self.assertEqual(outcome.provenance["provider_failures"], [])

    def test_the_reason_survives_a_whole_evaluation(self):
        """The path that lost it: resolver -> make_services -> criterion -> result."""
        results = evaluate_prepared_record(
            INPUT, make_services([], failures=ENSEMBL_FAILED), {}, ["PM4", "PP2", "BP7"])
        self.assertEqual(len(results), 3)
        for outcome in results:
            self.assertEqual(outcome.status, CriterionStatus.UNKNOWN, outcome.criterion)
            self.assertIn("429 Too Many Requests", outcome.summary, outcome.criterion)

    def test_failures_can_be_narrowed_to_one_provider(self):
        services = make_services([], failures=[
            *ENSEMBL_FAILED, {"provider": "TogoVar", "error": "timeout"}])
        self.assertEqual(retrieval_failures(services, "Ensembl"), ENSEMBL_FAILED)


if __name__ == "__main__":
    unittest.main()


ANNOTATION = {
    "category": "annotation", "variant_key": "GRCh38:14:23417200:G:C",
    "evidence_id": "urn:test:annotation", "source": "Ensembl VEP HGVS", "source_version": "116",
    "retrieved_at": "2026-09-17T00:00:00+00:00", "quality_status": "PASS",
    "transcript": "NM_000257.4", "gene": "MYH7", "consequences": ["missense_variant"],
    "protein_id": "NP_000248.2", "protein_start": 1491, "protein_end": 1491,
    "ref_aa": "S", "alt_aa": "C",
}
HOTSPOT_FAILED = [{"provider": "ClinVar protein hotspot density", "error": "HTTPError: 503"}]
PS1_SEARCH_FAILED = [{"provider": "ClinVar protein comparator search:PS1",
                      "error": "HTTPError: 503"}]


class HotspotFailureTests(unittest.TestCase):
    """PM1's region assessment comes from one provider, so 'nothing was retrieved' and 'the
    provider broke' are different things to tell a curator."""

    def evaluate(self, records, failures):
        return curated_context("PM1", "region", INPUT,
                               make_services([ANNOTATION, *records], failures=failures),
                               ANNOTATION, disease_required=False)[0]

    def test_the_failed_provider_is_named(self):
        outcome = self.evaluate([], HOTSPOT_FAILED)
        self.assertIn("ClinVar protein hotspot density", outcome.summary)
        self.assertIn("503", outcome.summary)
        self.assertEqual(outcome.provenance["provider_failures"], HOTSPOT_FAILED)

    def test_nothing_retrieved_and_nothing_failed_is_unchanged(self):
        outcome = self.evaluate([], [])
        self.assertEqual(outcome.summary,
                         "No region assessment was retrieved for this variant")
        self.assertNotIn("provider_failures", outcome.provenance)

    def test_an_unrelated_failure_is_not_blamed(self):
        """Naming a provider that has nothing to do with this category would send a curator
        after the wrong thing."""
        outcome = self.evaluate([], [{"provider": "TogoVar", "error": "timeout"}])
        self.assertNotIn("TogoVar", outcome.summary)

    def test_records_that_were_retrieved_and_rejected_keep_their_own_reason(self):
        """unusable_reason already says which rejection route they took; a provider failing
        elsewhere in the run did not cause that."""
        rejected = {**ANNOTATION, "category": "region", "evidence_id": "urn:test:region",
                    "protein_id": "NP_000248.2", "start": 1, "end": 10}
        outcome = self.evaluate([rejected], HOTSPOT_FAILED)
        self.assertNotIn("ClinVar protein hotspot density", outcome.summary)


class ComparatorFailureTests(unittest.TestCase):
    """The comparator search is absorbed per criterion, so the failures must be matched
    exactly - PS1's search failing is not PM5's reason for having nothing."""

    def evaluate(self, code, failures):
        from acmg_pipeline.criteria.comparator import evaluate_comparator
        return evaluate_comparator(code, INPUT,
                                   make_services([ANNOTATION], failures=failures), {})

    def test_ps1_reports_its_own_failed_search(self):
        outcome = self.evaluate("PS1", PS1_SEARCH_FAILED)
        self.assertEqual(outcome.status, CriterionStatus.UNKNOWN)
        self.assertIn("ClinVar protein comparator search:PS1", outcome.summary)
        self.assertEqual(outcome.provenance["provider_failures"], PS1_SEARCH_FAILED)

    def test_pm5_is_not_blamed_for_ps1s_failure(self):
        outcome = self.evaluate("PM5", PS1_SEARCH_FAILED)
        self.assertEqual(outcome.summary, "Comparator search incomplete")
        self.assertNotIn("provider_failures", outcome.provenance)
