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

passed = 0
failed = 0


def check(label, cond):
    global passed, failed
    if cond:
        passed += 1
        print(f"  OK   {label}")
    else:
        failed += 1
        print(f"  FAIL {label}")


# --- 1. Code list sanity ---
print("[1] ACMG code list sanity")
check("28 total codes", len(ALL_ACMG_CODES) == 28)
check("no overlap between pathogenic/benign lists", not (set(PATHOGENIC_CODES) & set(BENIGN_CODES)))
check("21 implemented codes", len(IMPLEMENTED_CODES) == 21)
check("7 stub codes", len(stubs.STUB_CODES) == 7)
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

# not_clear -> NOT_MET, not a crash and not silently MET
agg_empty = ps4.aggregate_multi_paper_results([])
check("not_clear aggregate -> NOT_MET", from_aggregated_judgment(agg_empty, "PS4").status == CriterionStatus.NOT_MET)

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
check("classify() accepts the full 28-entry set without error", full_result.not_evaluated_codes == [])

print(f"\n{'='*40}\n{passed} passed, {failed} failed\n{'='*40}")
sys.exit(1 if failed else 0)
