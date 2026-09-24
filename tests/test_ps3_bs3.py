import unittest

from acmg_pipeline.common import MatchStatus, VariantMatchingResult
from acmg_pipeline.criteria import ps3_bs3


def _judgment(rationale, directions):
    experiments = [
        ps3_bs3.ExperimentExtraction(
            assay_type=f"assay_{i}", experimental_system="cells", readout="readout",
            comparator="wild type", result_direction=d, key_findings=["a measured value of 1.0"],
        )
        for i, d in enumerate(directions)
    ]
    return ps3_bs3.PS3BS3Judgment(
        variant_matching=VariantMatchingResult(match_status=MatchStatus.MATCHED),
        experiments=experiments,
        overall_evidence=ps3_bs3.OverallEvidence(
            direction=ps3_bs3.OverallDirection.PS3, rationale=rationale,
        ),
    )


class ConflictReasoningTests(unittest.TestCase):
    """detect_experiment_conflict() still detects a paper whose extracted experiments
    disagree on direction, but finalize() no longer blanket-forces not_clear - see
    ps3_bs3._conflict_reasoned_in_rationale()'s own docstring for why the original
    "always force not_clear" behavior was itself suppressing a real MET (PTEN
    c.112C>T, ground truth PS3=MET, moderate)."""

    PTEN_C112CT_RATIONALE = (
        "Although VAMP-seq showed p.Pro38Ser has WT-like/enhanced protein "
        "abundance (not a loss-of-stability mechanism), a direct functional "
        "assay demonstrated that this variant drives increased Akt "
        "phosphorylation in the presence of endogenous wild-type PTEN, "
        "consistent with a dominant-negative pathogenic mechanism, and it "
        "is significantly enriched in melanoma. This functional (Akt "
        "phosphorylation) evidence, which most directly assays PTEN "
        "pathway activity, supports a damaging effect (PS3), even though "
        "the abundance assay alone would suggest a benign/normal result."
    )

    def test_detect_experiment_conflict_still_fires_on_mixed_directions(self):
        judgment = _judgment(
            self.PTEN_C112CT_RATIONALE,
            [ps3_bs3.ResultDirection.FUNCTIONALLY_NORMAL, ps3_bs3.ResultDirection.FUNCTIONALLY_ABNORMAL],
        )
        self.assertTrue(ps3_bs3.detect_experiment_conflict(judgment))

    def test_a_reasoned_rationale_keeps_the_direction_not_forced_to_not_clear(self):
        judgment = _judgment(
            self.PTEN_C112CT_RATIONALE,
            [ps3_bs3.ResultDirection.FUNCTIONALLY_NORMAL, ps3_bs3.ResultDirection.FUNCTIONALLY_ABNORMAL],
        )
        result = ps3_bs3.finalize(judgment, "PTEN", None, "PS3", pmid="29785012")
        self.assertEqual(result.effective_direction, ps3_bs3.OverallDirection.PS3)
        self.assertTrue(any(h.severity == "caution" for h in result.curator_hints))

    def test_an_unreasoned_rationale_still_forces_not_clear(self):
        """The original failure mode this whole check exists for: the LLM picks one
        direction and never even mentions the conflicting experiment."""
        judgment = _judgment(
            "The Akt phosphorylation assay showed a clear damaging effect, "
            "supporting PS3.",
            [ps3_bs3.ResultDirection.FUNCTIONALLY_NORMAL, ps3_bs3.ResultDirection.FUNCTIONALLY_ABNORMAL],
        )
        result = ps3_bs3.finalize(judgment, "PTEN", None, "PS3", pmid="29785012")
        self.assertEqual(result.effective_direction, ps3_bs3.OverallDirection.NOT_CLEAR)
        self.assertTrue(any(h.severity == "warning" for h in result.curator_hints))

    def test_no_conflict_at_all_is_unaffected(self):
        judgment = _judgment(
            "Both assays showed a clear damaging effect.",
            [ps3_bs3.ResultDirection.FUNCTIONALLY_ABNORMAL, ps3_bs3.ResultDirection.FUNCTIONALLY_ABNORMAL],
        )
        self.assertFalse(ps3_bs3.detect_experiment_conflict(judgment))
        result = ps3_bs3.finalize(judgment, "PTEN", None, "PS3", pmid="1")
        self.assertEqual(result.effective_direction, ps3_bs3.OverallDirection.PS3)


if __name__ == "__main__":
    unittest.main()
