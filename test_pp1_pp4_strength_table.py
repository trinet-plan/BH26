"""
test_pp1_pp4_strength_table.py

Tests acmg_pipeline.criteria.pp1_pp4_strength_table.combined_pp1_pp4_strength()
against Table 4 of Biesecker et al., 2024 (AJHG 111:24-38) directly, at the
pure-function level - see test_pp1_bs4_pp4_engine.py for the same logic
exercised through the full engine.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.classification import Strength
from acmg_pipeline.criteria.pp1_pp4_strength_table import combined_pp1_pp4_strength
from test_harness import Harness

h = Harness()
check = h.check


# ============================================================================
# [1] Only one criterion has evidence -> ordinary independent floor
# ============================================================================
print("[1] Single-criterion cases")

r = combined_pp1_pp4_strength(pp1_points=0.0, pp4_points=4.0)
check("PP4 alone at 4.0 -> Strong, PP1 None", r.pp4 == Strength.STRONG and r.pp1 is None)

r = combined_pp1_pp4_strength(pp1_points=1.0, pp4_points=0.0)
check("PP1 alone at 1.0 -> Supporting, PP4 None", r.pp1 == Strength.SUPPORTING and r.pp4 is None)

r = combined_pp1_pp4_strength(pp1_points=0.0, pp4_points=0.0)
check("Neither -> both None", r.pp1 is None and r.pp4 is None)

# PP1 alone is not capped anywhere else in pp4_pp1_bs4.py's segregation
# scoring (unlike PP4, which evaluator.evaluate_pp4() caps internally) - the
# paper's Table 3 footnote states the +5.0 cap applies to "all locus
# evidence (PP1 and PP4)," so a solo PP1 value above the cap must still
# floor as if capped at 5.0, not at its own raw (higher) value.
r = combined_pp1_pp4_strength(pp1_points=10.0, pp4_points=0.0)
check("PP1 alone at 10.0 raw -> capped to 5.0 -> Strong, not VeryStrong",
      r.pp1 == Strength.STRONG and r.pp4 is None)


# ============================================================================
# [2] Both criteria contribute -> Table 4 governs the split
# ============================================================================
print("\n[2] Combined PP1+PP4 cases (Table 4)")

# Reproduces test_pp1_bs4_pp4_engine.py's [2] sibling_results fixture:
# PP4=4.0, PP1=1.0 -> total 5.0 -> bucket >=5.0 options (4,1)/(1,4).
# PP4 has the larger raw share, so PP4 gets the larger (Strong) tier.
r = combined_pp1_pp4_strength(pp1_points=1.0, pp4_points=4.0)
check("PP4=4.0,PP1=1.0 (total 5.0) -> PP4 Strong, PP1 Supporting",
      r.pp4 == Strength.STRONG and r.pp1 == Strength.SUPPORTING)

# Same total (5.0) but PP1 now has the larger raw share -> PP1 gets Strong.
r = combined_pp1_pp4_strength(pp1_points=4.0, pp4_points=1.0)
check("PP1=4.0,PP4=1.0 (total 5.0) -> PP1 Strong, PP4 Supporting",
      r.pp1 == Strength.STRONG and r.pp4 == Strength.SUPPORTING)

# PP4=3.0, PP1=3.0: raw sum 6.0, capped to 5.0 -> bucket >=5.0. Exact tie ->
# PP4 (evaluated first in the paper's own Figure 1 flow) gets the larger tier.
r = combined_pp1_pp4_strength(pp1_points=3.0, pp4_points=3.0)
check("PP4=3.0,PP1=3.0 (capped total 5.0), tie -> PP4 Strong, PP1 Supporting",
      r.pp4 == Strength.STRONG and r.pp1 == Strength.SUPPORTING)

# total 4.0-4.9 bucket with both nonzero -> the only combo option, not the
# solo-Strong option, since PP1 and PP4 both have real evidence here.
r = combined_pp1_pp4_strength(pp1_points=2.0, pp4_points=2.0)
check("PP4=2.0,PP1=2.0 (total 4.0) -> PP1 Moderate, PP4 Moderate",
      r.pp1 == Strength.MODERATE and r.pp4 == Strength.MODERATE)

# total 3.0-3.9 bucket: no solo option exists at this bucket at all (see
# module docstring) - only the two combo splits.
r = combined_pp1_pp4_strength(pp1_points=1.0, pp4_points=2.0)
check("PP4=2.0,PP1=1.0 (total 3.0) -> PP4 Moderate, PP1 Supporting",
      r.pp4 == Strength.MODERATE and r.pp1 == Strength.SUPPORTING)

r = combined_pp1_pp4_strength(pp1_points=2.0, pp4_points=1.0)
check("PP1=2.0,PP4=1.0 (total 3.0) -> PP1 Moderate, PP4 Supporting",
      r.pp1 == Strength.MODERATE and r.pp4 == Strength.SUPPORTING)

# total 2.0-2.9 bucket, both nonzero and equal -> symmetric Supporting+Supporting.
r = combined_pp1_pp4_strength(pp1_points=1.0, pp4_points=1.0)
check("PP4=1.0,PP1=1.0 (total 2.0) -> both Supporting",
      r.pp4 == Strength.SUPPORTING and r.pp1 == Strength.SUPPORTING)

# total just under 1.0 combined -> not applicable, both None, even though
# both criteria technically have some (tiny) nonzero evidence.
r = combined_pp1_pp4_strength(pp1_points=0.4, pp4_points=0.4)
check("PP4=0.4,PP1=0.4 (total 0.8, below bucket 1) -> both None",
      r.pp4 is None and r.pp1 is None)

# Bucket 1.0-1.9 has no both-nonzero option (minimum tier already fills the
# bucket) - the minor contributor (PP1's +0.4 AR-unaffected increment) rounds
# down to nothing rather than fabricating a fractional Supporting tier.
r = combined_pp1_pp4_strength(pp1_points=0.4, pp4_points=1.0)
check("PP4=1.0,PP1=0.4 (total 1.4, no both-nonzero split) -> PP4 Supporting, PP1 None",
      r.pp4 == Strength.SUPPORTING and r.pp1 is None)


# ============================================================================
# [3] Input validation
# ============================================================================
print("\n[3] Input validation")
try:
    combined_pp1_pp4_strength(pp1_points=-1.0, pp4_points=0.0)
    check("negative pp1_points raises ValueError", False)
except ValueError:
    check("negative pp1_points raises ValueError", True)


h.report_and_exit()
