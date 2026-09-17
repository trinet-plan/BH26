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
from acmg_pipeline.criteria.common import annotation_context, retrieval_failures

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
