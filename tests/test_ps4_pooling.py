import unittest

from acmg_pipeline.common import MatchStatus, PaperContribution, VariantMatchingResult
from acmg_pipeline.criteria import ps4


def _contribution(pmid, *, affected_carriers=None, unaffected_carriers=None,
                   direction=ps4.CaseControlDirection.NOT_CLEAR, matched=True):
    judgment = ps4.PS4Judgment(
        variant_matching=VariantMatchingResult(
            match_status=MatchStatus.MATCHED if matched else MatchStatus.UNSUCCESSFUL),
        study_design=ps4.StudyDesign.CASE_SERIES,
        case_control_data=ps4.CaseControlData(
            affected_carriers=affected_carriers, unaffected_carriers=unaffected_carriers),
        overall_evidence=ps4.CaseControlEvidence(direction=direction, rationale="pooled test"),
    )
    result = ps4.finalize(judgment, pmid=pmid)
    return PaperContribution(pmid, result)


class PooledCaseSeriesTests(unittest.TestCase):
    """A real ClinGen PS4 determination is often pooled across several case-series
    papers, none of which alone is a formal case-control study - see
    ps4._pooled_case_series_aggregate()'s own docstring."""

    def test_two_papers_pooling_past_the_threshold_yields_ps4(self):
        contributions = [
            _contribution("1", affected_carriers=2),
            _contribution("2", affected_carriers=2),
        ]
        aggregated = ps4.aggregate_multi_paper_results(contributions)
        self.assertEqual(aggregated.aggregated_direction, ps4.CaseControlDirection.PS4)
        self.assertTrue(any(h.severity == "caution" for h in aggregated.aggregation_hints))

    def test_a_single_paper_alone_never_triggers_pooling(self):
        """POOLED_MIN_PAPERS=2: one case series report, even with plenty of carriers,
        must not count as pooled corroboration on its own."""
        contributions = [_contribution("1", affected_carriers=5)]
        aggregated = ps4.aggregate_multi_paper_results(contributions)
        self.assertEqual(aggregated.aggregated_direction, ps4.CaseControlDirection.NOT_CLEAR)

    def test_pooling_below_the_affected_count_threshold_stays_not_clear(self):
        contributions = [
            _contribution("1", affected_carriers=1),
            _contribution("2", affected_carriers=1),
        ]
        aggregated = ps4.aggregate_multi_paper_results(contributions)
        self.assertEqual(aggregated.aggregated_direction, ps4.CaseControlDirection.NOT_CLEAR)

    def test_a_paper_reporting_unaffected_carriers_is_excluded_from_pooling(self):
        """A paper with a real control-side count already had its own chance to
        score PS4 through the standard case-control path; pooling it here too
        would double-count evidence the primary path already rejected."""
        contributions = [
            _contribution("1", affected_carriers=2, unaffected_carriers=5),
            _contribution("2", affected_carriers=2),
        ]
        aggregated = ps4.aggregate_multi_paper_results(contributions)
        self.assertEqual(aggregated.aggregated_direction, ps4.CaseControlDirection.NOT_CLEAR)

    def test_an_unmatched_paper_is_excluded_from_pooling(self):
        contributions = [
            _contribution("1", affected_carriers=2, matched=False),
            _contribution("2", affected_carriers=2),
        ]
        aggregated = ps4.aggregate_multi_paper_results(contributions)
        self.assertEqual(aggregated.aggregated_direction, ps4.CaseControlDirection.NOT_CLEAR)

    def test_a_confident_single_paper_ps4_is_unaffected_by_pooling(self):
        """When the standard aggregation already reaches PS4 from a real
        case-control paper, pooling must never override or interfere."""
        contributions = [
            _contribution("1", affected_carriers=8, direction=ps4.CaseControlDirection.PS4),
        ]
        aggregated = ps4.aggregate_multi_paper_results(contributions)
        self.assertEqual(aggregated.aggregated_direction, ps4.CaseControlDirection.PS4)


if __name__ == "__main__":
    unittest.main()
