"""
test_full_criteria_ground_truth.py

Two things, deliberately kept separate:

  [1] Integrity checks on test_data/full_criteria_ground_truth.py itself -
      this is a DATA sanity check (are the entries well-formed?), not a
      test of any judgment logic. No acmg_pipeline module beyond
      classification.py/gate.py's shared enums is exercised here.

  [2] A classify()-vs-reality comparison, run and REPORTED but never
      asserted pass/fail. Feeding classify() the criteria represented in
      the collected ground truth is expected to differ from a complete
      each variant's real, already-known classification for most of these
      variant review - that gap is
      the whole point of classification.py's `not_evaluated_codes` field,
      and printing it honestly is more useful than hiding it behind a
      passing assertion that doesn't mean what it looks like it means.
"""

import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.classification import (
    ALL_ACMG_CODES, AUTOMATED_CODES, IMPLEMENTED_CODES, LITERATURE_CODES,
    CriterionEvidence, classify,
)
from acmg_pipeline.gate import CriterionStatus
from acmg_pipeline.criteria import stubs
from test_data.full_criteria_ground_truth import GROUND_TRUTH, unique_variants, entries_for

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


# ============================================================================
# [1] Data integrity
# ============================================================================
print("[1] Dataset integrity")
check(f"at least 300 entries collected (got {len(GROUND_TRUTH)})", len(GROUND_TRUTH) >= 300)
check("every entry's criterion is a recognized ACMG code",
      all(e.criterion in ALL_ACMG_CODES for e in GROUND_TRUTH))
check("every MET entry carries a strength", all(e.strength is not None for e in GROUND_TRUTH if e.status == CriterionStatus.MET))
check("every NOT_MET entry has no strength", all(e.strength is None for e in GROUND_TRUTH if e.status == CriterionStatus.NOT_MET))
check("source is always erepo or democase", all(e.source in ("erepo", "democase") for e in GROUND_TRUTH))
check("tier is always A/B/C", all(e.tier in ("A", "B", "C") for e in GROUND_TRUTH))
check("no exact duplicate (gene, hgvsc, criterion, source) rows",
      len({(e.gene, e.hgvsc, e.criterion, e.source) for e in GROUND_TRUTH}) == len(GROUND_TRUTH))

by_criterion = Counter(e.criterion for e in GROUND_TRUTH)
covered = {c for c in ALL_ACMG_CODES if by_criterion.get(c, 0) > 0}
print(f"\n  criteria with at least 1 ground-truth entry: {len(covered)}/28")
print(f"  criteria with zero entries: {sorted(set(ALL_ACMG_CODES) - covered)}")
check("the 5 literature codes all have ground-truth coverage",
      LITERATURE_CODES.issubset(covered))
check("PP4 has real ground-truth coverage (this session's specific ask)",
      by_criterion.get("PP4", 0) > 0)
_layer1_codes = sorted(AUTOMATED_CODES)
_layer1_covered = [c for c in _layer1_codes if by_criterion.get(c, 0) > 0]
check(f"at least 14 of the 16 Layer-1 codes have coverage (got {len(_layer1_covered)})",
      len(_layer1_covered) >= 14)
if len(_layer1_covered) < len(_layer1_codes):
    missing = sorted(set(_layer1_codes) - set(_layer1_covered))
    print(f"  NOTE Layer-1 codes with zero ground-truth entries: {missing} "
          f"(PP5/BP6 = pure ClinVar/reputable-source pass-through; VCEPs that "
          f"curate their own primary evidence rarely need to cite them - not a "
          f"data-collection gap, see design doc section 15-14)")

print(f"\n{sum(1 for _ in unique_variants())} unique variants across the dataset")


# ============================================================================
# [2] classify() vs. each variant's real, already-known classification
# ============================================================================
print("\n[2] classify() (available integrated codes) vs. each variant's real classification")
print("    (mismatches are EXPECTED - see module docstring; this is a report, not a pass/fail gate)\n")

match_n = 0
compared_n = 0
for gene, hgvsc in unique_variants():
    entries = entries_for(gene, hgvsc)
    # Only variants with a known, single-answer overall classification are
    # comparable at all (democase's DSG2 entries deliberately leave
    # variant_outcome=None - see that module's docstring - because the
    # source document itself reports no single settled answer for it).
    outcomes = {e.variant_outcome for e in entries if e.variant_outcome}
    if len(outcomes) != 1:
        continue
    real_outcome = outcomes.pop()

    our_evidence = [
        CriterionEvidence(code=e.criterion, status=e.status, strength=e.strength)
        for e in entries if e.criterion in IMPLEMENTED_CODES
    ]
    if not our_evidence:
        continue  # no integrated implemented criterion exists for this variant

    result = classify(our_evidence)
    compared_n += 1
    is_match = result.category.value == real_outcome
    match_n += is_match
    codes_used = ", ".join(
        f"{e.code}({e.strength.value})" if e.strength else f"{e.code}(not_met)"
        for e in our_evidence
    )
    print(f"  {'MATCH' if is_match else 'differs'}  {gene} {hgvsc}: "
          f"our classify()={result.category.value} (score={result.score}, from {codes_used}) "
          f"vs. real={real_outcome}")

print(f"\n{compared_n} variant(s) had both a real known classification and >=1 integrated implemented "
      f"codes; {match_n}/{compared_n} matched using only those codes (the rest may be unevaluated).")


print(f"\n{'='*40}\n{passed} passed, {failed} failed\n{'='*40}")
sys.exit(1 if failed else 0)
