import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.classification import (
    ALL_ACMG_CODES, IMPLEMENTED_CODES, PATHOGENIC_CODES, BENIGN_CODES,
    CriterionEvidence, CriterionStatus, Strength, ClassificationCategory,
    classify, from_aggregated_judgment,
)
from acmg_pipeline.criteria.registry import is_implemented, get_criterion_evidence
from acmg_pipeline.criteria import stubs
from acmg_pipeline.common import MatchStatus, VariantMatchingResult, PaperContribution
from acmg_pipeline.criteria import ps4, segregation as seg
from test_harness import Harness

h = Harness()
check = h.check


# --- 1. Code list sanity ---
print("[1] ACMG code list sanity")
check("28 total codes", len(ALL_ACMG_CODES) == 28)
check("no overlap between pathogenic/benign lists", not (set(PATHOGENIC_CODES) & set(BENIGN_CODES)))
# IMPLEMENTED_CODES history: {"PS3","BS3","PS4","PP1","BS4"} (5) through
# 2026-09-15 -> {"PS3","BS3","PS4"} (3) on 2026-09-16 when PP1/BS4/PP4 were
# handed off to another team -> +16 automated codes (19) once that team's
# evidence-cli engine was merged in -> +PP1/BS4/PP4 again (22) on
# 2026-09-17 once acmg_pipeline.criteria.pp1_bs4_pp4_engine connected that
# same team's ClinGen 2024 PP4+PP1/BS4 evaluator (segregation.py is still
# real, tested code, still exercised directly in section [5] below, but no
# longer the thing IMPLEMENTED_CODES credits for PP1/BS4 - pp1_bs4_pp4_
# engine is) -> +PS2/PM6 (24) on 2026-09-19 once acmg_pipeline.criteria.
# ps2_pm6 connected a rule-based de-novo evaluator.
check("24 implemented codes", len(IMPLEMENTED_CODES) == 24)
check("4 stub codes", len(stubs.STUB_CODES) == 4)
check("implemented + stub codes cover all 28 with no overlap",
      set(stubs.STUB_CODES) | IMPLEMENTED_CODES == set(ALL_ACMG_CODES)
      and not (set(stubs.STUB_CODES) & IMPLEMENTED_CODES))

# --- 2. Tavtigian point thresholds (Table 3-equivalent categories) ---
print("\n[2] Point-score -> classification category thresholds")
def score_of(*pairs):
    """pairs: (code, Strength) tuples, all MET."""
    return classify([CriterionEvidence(c, CriterionStatus.MET, s) for c, s in pairs]).score

check("PVS1(very_strong)+PM2(moderate) = 10 -> Pathogenic",
      classify([CriterionEvidence("PVS1", CriterionStatus.MET, Strength.VERY_STRONG),
                CriterionEvidence("PM2", CriterionStatus.MET, Strength.MODERATE)]).category
      == ClassificationCategory.PATHOGENIC)
check("PVS1(very_strong) alone = 8 -> Likely Pathogenic",
      classify([CriterionEvidence("PVS1", CriterionStatus.MET, Strength.VERY_STRONG)]).category
      == ClassificationCategory.LIKELY_PATHOGENIC)
check("no evidence at all = 0 -> Uncertain Significance",
      classify([]).category == ClassificationCategory.UNCERTAIN_SIGNIFICANCE)
check("BS3(strong)+BP4(supporting) = -5 -> Likely Benign",
      classify([CriterionEvidence("BS3", CriterionStatus.MET, Strength.STRONG),
                CriterionEvidence("BP4", CriterionStatus.MET, Strength.SUPPORTING)]).category
      == ClassificationCategory.LIKELY_BENIGN)
check("BS3(very_strong)+BS4(very_strong) = -16 -> Benign",
      classify([CriterionEvidence("BS3", CriterionStatus.MET, Strength.VERY_STRONG),
                CriterionEvidence("BS4", CriterionStatus.MET, Strength.VERY_STRONG)]).category
      == ClassificationCategory.BENIGN)

# --- 3. BA1 stand-alone override ---
print("\n[3] BA1 stand-alone override")
result = classify([
    CriterionEvidence("BA1", CriterionStatus.MET, Strength.STAND_ALONE),
    CriterionEvidence("PVS1", CriterionStatus.MET, Strength.VERY_STRONG),
])
check("BA1 met forces Benign even with strong pathogenic evidence present", result.category == ClassificationCategory.BENIGN)
check("ba1_override flag is set", result.ba1_override)
check("the overridden pathogenic evidence is still visible in `met` for curator review",
      any(e.code == "PVS1" for e in result.met))

# --- 4. Duplicate-code rejection ---
print("\n[4] Duplicate evidence for the same code is rejected")
try:
    classify([
        CriterionEvidence("PS3", CriterionStatus.MET, Strength.STRONG),
        CriterionEvidence("PS3", CriterionStatus.MET, Strength.SUPPORTING),
    ])
    check("duplicate PS3 raises ValueError", False)
except ValueError:
    check("duplicate PS3 raises ValueError", True)

# --- 5. from_aggregated_judgment(): real per-paper data -> CriterionEvidence ---
# PP1/BS4 (via segregation.py) are exercised here as a generic test of
# from_aggregated_judgment()'s direction-matching logic, independent of
# IMPLEMENTED_CODES - segregation.py is still real, working, previously-
# validated code (see acmg_pipeline.pipeline.ENGINE_BY_CRITERION), even
# though PP1/BS4 are no longer in IMPLEMENTED_CODES as of 2026-09-16.
print("\n[5] from_aggregated_judgment() (using PS4/PP1 schemas with realistic fixture data)")

j_ps4 = ps4.PS4Judgment(
    variant_matching=VariantMatchingResult(match_status=MatchStatus.MATCHED),
    study_design=ps4.StudyDesign.CASE_CONTROL,
    case_control_data=ps4.CaseControlData(odds_ratio=20.0, p_value=0.0001),
    overall_evidence=ps4.CaseControlEvidence(direction=ps4.CaseControlDirection.PS4, rationale="OR=20, p<0.0001"),
)
r_ps4 = ps4.finalize(j_ps4, pmid="31160754")
agg_ps4 = ps4.aggregate_multi_paper_results([PaperContribution("31160754", r_ps4)])
ev_ps4 = from_aggregated_judgment(agg_ps4, "PS4")
check("PS4 direction=PS4 -> status=MET", ev_ps4.status == CriterionStatus.MET)
check("1 relevant paper -> strength=supporting", ev_ps4.strength == Strength.SUPPORTING)

fam = seg.FamilySegregationData(family_id="F1", affected_with_variant=3, affected_without_variant=0,
                                  unaffected_with_variant=0, unaffected_without_variant=2)
j_pp1 = seg.SegregationJudgment(
    variant_matching=VariantMatchingResult(match_status=MatchStatus.MATCHED),
    families=[fam],
    overall_evidence=seg.SegregationEvidence(direction=seg.SegregationDirection.PP1, rationale="segregates"),
)
r_pp1 = seg.finalize(j_pp1, pmid="99999999")
agg_pp1 = seg.aggregate_multi_paper_results([PaperContribution("99999999", r_pp1)])

check("asking for BS4 when the evidence actually points PP1 -> NOT_MET (not mislabeled MET)",
      from_aggregated_judgment(agg_pp1, "BS4").status == CriterionStatus.NOT_MET)
check("asking for PP1 (the direction that was actually found) -> MET",
      from_aggregated_judgment(agg_pp1, "PP1").status == CriterionStatus.MET)

# not_clear -> UNKNOWN, not a crash and not silently MET
agg_empty = ps4.aggregate_multi_paper_results([])
check("not_clear aggregate -> UNKNOWN", from_aggregated_judgment(agg_empty, "PS4").status == CriterionStatus.UNKNOWN)

# --- 6. registry.get_criterion_evidence(): full 28-code loop ---
print("\n[6] registry.get_criterion_evidence() over all 28 codes")
fakes = {code: CriterionEvidence(code, CriterionStatus.MET, Strength.SUPPORTING) for code in IMPLEMENTED_CODES}
evidence = [get_criterion_evidence(code, fakes.get(code)) for code in ALL_ACMG_CODES]
check("produces exactly 28 CriterionEvidence entries", len(evidence) == 28)
check("implemented codes come through unchanged", all(e.status == CriterionStatus.MET for e in evidence if e.code in IMPLEMENTED_CODES))
check("stub codes come through as UNKNOWN/not-evaluated", all(e.status == CriterionStatus.UNKNOWN for e in evidence if e.code in stubs.STUB_CODES))

try:
    get_criterion_evidence("PS3", None)
    check("an implemented code with no real evidence supplied raises ValueError", False)
except ValueError:
    check("an implemented code with no real evidence supplied raises ValueError", True)

full_result = classify(evidence)
check("UNKNOWN entries remain not evaluated", set(full_result.not_evaluated_codes) == set(stubs.STUB_CODES))

h.report_and_exit()
