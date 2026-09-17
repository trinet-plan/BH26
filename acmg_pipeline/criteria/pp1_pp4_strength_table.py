"""
pp1_pp4_strength_table.py

Table 4 from Biesecker et al., 2024, Am. J. Hum. Genet. 111:24-38
("ClinGen guidance for use of the PP1/BS4 co-segregation and PP4 phenotype
specificity criteria for sequence variant pathogenicity classification"),
DOI:10.1016/j.ajhg.2023.11.009 - "Conversion of combined phenotype
specificity and co-segregation Bayesian points to a combined PP1/PP4
strength."

Why this is a separate module from pp1_bs4_pp4_engine.py's own,
now-superseded strength floor
---------------------------------------------------------------------------
pp1_bs4_pp4_engine.py previously floored PP4's own points and PP1's own
points to a Strength INDEPENDENTLY. That is only correct when just one of
the two criteria actually has evidence for the locus (the common case,
handled the same way here). When BOTH PP4 and PP1 contribute points for the
SAME locus, flooring each independently can silently exceed the paper's own
+5.0-point shared locus-evidence cap: e.g. PP4=4.0 points (Strong) and
PP1=1.0 point (Supporting) independently floor to a REPORTED "Strong +
Supporting" pair that not only happens to already equal the paper's own
Table 4 answer for a combined total of 5.0, but a slightly different raw
split - say PP4=3.0, PP1=3.0 (identical combined total of 6.0, capped to
5.0) - would floor independently to "Moderate + Moderate" (a combined
"value" of 4, understating the actual capped total of 5) or, worse, two
raw values that independently floor above what +5.0 combined evidence can
actually support. Table 4 is ClinGen's own answer to exactly this: given
the combined, already-capped point total, which specific (PP1 tier, PP4
tier) pairs are the maximum allowed, so the two reported strengths are
consistent with the actual combined evidence rather than two independent,
uncoordinated floors.

[Table 4, as enumerated tier pairs]
  Reconstructed from the paper's Table 4 (points range -> maximum allowable
  strength of combined codes). Points follow the standard ACMG/AMP Bayesian
  point scale (Tavtigian et al., 2020): Supporting=1, Moderate=2, Strong=4,
  VeryStrong=8 (VeryStrong never appears here - the shared +5.0 locus cap
  makes it unreachable for PP1/PP4). Every listed pair sums to exactly that
  bucket's own floor value, consistent with Table 2/Table 3's own "round
  down between rows" convention:
    total >= 5.0 : (PP1_Strong, PP4_Supporting) or (PP1_Supporting, PP4_Strong)
    4.0-4.9      : (PP1_Strong, -) or (-, PP4_Strong) or (PP1_Moderate, PP4_Moderate)
    3.0-3.9      : (PP1_Supporting, PP4_Moderate) or (PP1_Moderate, PP4_Supporting)
    2.0-2.9      : (PP1_Moderate, -) or (-, PP4_Moderate) or (PP1_Supporting, PP4_Supporting)
    1.0-1.9      : (PP1_Supporting, -) or (-, PP4_Supporting)
    < 1.0        : not applicable (neither criterion met)
  Note bucket 3.0-3.9 has no single-criterion-alone option: no single
  Tavtigian tier lands exactly there (the next tier above Moderate=2 is
  Strong=4, which is bucket 4.0-4.9), so a solo-criterion total of 3.x
  never occurs; it can only arise as a genuine PP1+PP4 combination.

[Which option to pick when Table 4 lists more than one]
  The paper explicitly leaves this kind of allocation to "professional
  judgment" (its own words, discussing the analogous PP4/PS4 split). For an
  automated pipeline we need a deterministic rule instead: when a bucket
  offers an option that credits BOTH criteria (both tiers nonzero), this
  module prefers that over an "alone" option, since both genuinely
  contributed evidence here (bucket 1.0-1.9 has no both-nonzero option at
  all - its minimum tier already fills the whole bucket - so there the
  minor contributor's share rounds down to nothing, the same as this
  project's other floor/round-down conventions). Among the remaining
  candidates, the module gives the larger tier to whichever criterion (PP4
  or PP1) actually earned the larger RAW point share for this locus, so the
  reported split reflects where the real evidence came from rather than an
  arbitrary convention. On an exact tie with both sides non-zero, PP4 gets
  the larger tier (PP4 phenotype-specificity is evaluated first in the
  paper's own Figure 1 flow diagram). This tie-break, like
  pp1_bs4_pp4_engine's own former per-criterion floor, is this project's
  interpretation, not an authoritative ClinGen reading - flagged here the
  same way.

[The +5.0 cap also applies to a single criterion alone]
  Table 3's own footnote states the +5.0 cap applies to "all locus evidence
  (PP1 and PP4)," not only to a PP1+PP4 combination. A single criterion's
  raw points (in particular PP1, which pp4_pp1_bs4.py's segregation scoring
  does not itself cap - unlike PP4, which evaluator.evaluate_pp4() already
  caps internally at LOCUS_EVIDENCE_CAP) are therefore also capped here
  before flooring, even when the other criterion has no evidence at all.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from acmg_pipeline.classification import Strength
from acmg_pipeline.criteria.pp4_pp1_bs4 import LOCUS_EVIDENCE_CAP

# (pp1_tier, pp4_tier) pairs, keyed by each bucket's own floor value.
_TABLE_4_OPTIONS: tuple[tuple[float, tuple[tuple[int, int], ...]], ...] = (
    (5.0, ((4, 1), (1, 4))),
    (4.0, ((4, 0), (0, 4), (2, 2))),
    (3.0, ((1, 2), (2, 1))),
    (2.0, ((2, 0), (0, 2), (1, 1))),
    (1.0, ((1, 0), (0, 1))),
)

_TIER_TO_STRENGTH: dict[int, Optional[Strength]] = {
    0: None,
    1: Strength.SUPPORTING,
    2: Strength.MODERATE,
    4: Strength.STRONG,
    8: Strength.VERY_STRONG,
}


@dataclass
class CombinedStrength:
    pp1: Optional[Strength]
    pp4: Optional[Strength]


def _floor_strength(points: float) -> Optional[Strength]:
    if points >= 8:
        return Strength.VERY_STRONG
    if points >= 4:
        return Strength.STRONG
    if points >= 2:
        return Strength.MODERATE
    if points >= 1:
        return Strength.SUPPORTING
    return None


def _table_4_bucket(total: float) -> tuple[tuple[int, int], ...]:
    for floor, options in _TABLE_4_OPTIONS:
        if total >= floor:
            return options
    return ()


def combined_pp1_pp4_strength(*, pp1_points: float, pp4_points: float) -> CombinedStrength:
    """Split a locus's PP1+PP4 evidence into per-criterion Strengths.

    `pp1_points`/`pp4_points` are each criterion's own RAW point
    contribution for this locus - e.g.
    LocusEvidenceResult.segregation.pp1_points_used and
    LocusEvidenceResult.pp4.points - not pre-floored to a Strength, and not
    necessarily pre-capped. Pass 0 for a criterion with no evidence at all
    (not evaluable, or not applicable).
    """
    if pp1_points < 0 or pp4_points < 0:
        raise ValueError("pp1_points and pp4_points must be >= 0")

    if pp1_points <= 0 or pp4_points <= 0:
        # Only one criterion has evidence for this locus - Table 4 governs
        # splitting evidence BETWEEN PP1 and PP4, not this ordinary case.
        return CombinedStrength(
            pp1=_floor_strength(min(pp1_points, LOCUS_EVIDENCE_CAP)) if pp1_points > 0 else None,
            pp4=_floor_strength(min(pp4_points, LOCUS_EVIDENCE_CAP)) if pp4_points > 0 else None,
        )

    total = min(pp1_points + pp4_points, LOCUS_EVIDENCE_CAP)
    options = _table_4_bucket(total)
    if not options:
        return CombinedStrength(pp1=None, pp4=None)

    # Prefer a pair that credits BOTH criteria (both tiers nonzero) since
    # both genuinely contributed evidence here; a bucket like 1.0-1.9 offers
    # no such pair (its only options are (1,0)/(0,1) - the minimum tier
    # already consumes the whole bucket), so the minor contributor's share
    # rounds down to nothing, same as this project's other floor/round-down
    # conventions.
    both_nonzero = [pair for pair in options if pair[0] > 0 and pair[1] > 0]
    candidates = both_nonzero or options

    major_is_pp4 = pp4_points >= pp1_points
    pp1_tier, pp4_tier = max(candidates, key=lambda pair: pair[1] if major_is_pp4 else pair[0])
    return CombinedStrength(pp1=_TIER_TO_STRENGTH[pp1_tier], pp4=_TIER_TO_STRENGTH[pp4_tier])
