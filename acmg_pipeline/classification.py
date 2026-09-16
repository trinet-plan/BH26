"""
classification.py
Combines per-criterion evidence (met / not_met / not_evaluated, plus a
strength for met calls) into a single final ACMG/AMP variant classification.

[Method: Tavtigian et al. 2018 point-based system]
  Uses the published point-based extension of the original Richards et al.
  2015 ACMG/AMP categorical combining rules (Tavtigian SV et al.,
  "Modeling the ACMG/AMP variant classification guidelines as a Bayesian
  classification framework", Genetics in Medicine 2018) - already cited in
  this project's own design doc and in the curator UI requirements doc
  (doc/recs for expert board.docx, its Table 2/Table 3 images). Each met
  criterion contributes points by its strength (Supporting=1, Moderate=2,
  Strong=4, Very Strong=8 for pathogenic evidence; the same magnitudes
  negative for benign evidence), and the summed score maps to one of the
  five classification categories via published thresholds. This numeric
  method is a well-documented, standard extension of the categorical rules
  (it reproduces Richards et al. Table 5's combinations, not an alternative
  to them) - implementing it here is applying a public scientific method as
  code, not reproducing any particular publication's text/figures.

[Why this exists now]
  Per the user's decision (2026-09-15): since the eventual external API
  interface with other teams isn't settled yet, build the classification
  engine first - it can be built and tested today using only this
  project's own implemented criteria, and needs no changes when other
  teams' criteria (PVS1, PM1, PM2, ...) are added later, since it only
  cares about (code, status, strength) tuples, not how they were produced.

[Where the other ACMG codes come from]
  This module owns the canonical list of all 28 ACMG/AMP 2015 codes
  (ALL_ACMG_CODES) and which ones this project actually implements
  (IMPLEMENTED_CODES = {"PS3", "BS3", "PS4"} as of 2026-09-16 - see the
  comment on IMPLEMENTED_CODES below for the PP1/BS4/PP4 handoff history).
  acmg_pipeline/criteria/stubs.py uses these to generate a NOT_EVALUATED
  CriterionEvidence for every other code, so classify() can be called with
  a genuinely complete evidence set (missing evidence for an unimplemented
  code is explicit and visible, not silently absent).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from acmg_pipeline.common import AggregatedJudgment, is_not_clear, strength_tier_from_paper_count
from acmg_pipeline.gate import CriterionStatus, MatchStatus

# ============================================================================
# 1. The full ACMG/AMP 2015 code list
# ============================================================================

PATHOGENIC_CODES = [
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
]
BENIGN_CODES = [
    "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
]
ALL_ACMG_CODES = PATHOGENIC_CODES + BENIGN_CODES

# The codes this project has a real judgment engine for (see
# acmg_pipeline/criteria/{ps3_bs3,ps4}.py and design doc section 15-5 for
# why these and not the other 25).
#
# [2026-09-16] Was {"PS3", "BS3", "PS4", "PP1", "BS4"} through the first
# several months of this project - PP1/BS4 (acmg_pipeline/criteria/
# segregation.py) were fully implemented and validated here (real LLM
# runs, ERepo-sourced ground truth, VA-Spec export). Per the user's
# explicit decision, PP1/BS4's judgment logic (and PP4's design research,
# see stubs.py) have since been handed off to another team, with this
# project's history transferred at both the doc and code level. This
# project's own implemented scope is now PS3/BS3/PS4 only.
# segregation.py itself is NOT deleted (still real, working, tested code -
# see acmg_pipeline.pipeline.ENGINE_BY_CRITERION, which still maps PP1/BS4
# to it) - only this constant, and everything derived from it (stubs.py's
# STUB_CODES, the default criteria this project's own pipeline entry
# points evaluate), changed to reflect that PP1/BS4 aren't this project's
# own scope to claim credit for going forward.
IMPLEMENTED_CODES = {"PS3", "BS3", "PS4"}


class Strength(str, Enum):
    SUPPORTING = "supporting"
    MODERATE = "moderate"
    STRONG = "strong"
    VERY_STRONG = "very_strong"
    STAND_ALONE = "stand_alone"


# Tavtigian et al. 2018 point values. STAND_ALONE has no formal point value
# in the original table (it is a categorical override, not a point tier -
# see the BA1 special case in classify() below) but is included here at the
# same magnitude as VERY_STRONG so a caller that passes it by mistake still
# gets a sane (if slightly wrong) score rather than a KeyError.
_POINTS_PATHOGENIC = {
    Strength.SUPPORTING: 1,
    Strength.MODERATE: 2,
    Strength.STRONG: 4,
    Strength.VERY_STRONG: 8,
    Strength.STAND_ALONE: 8,
}
_POINTS_BENIGN = {
    Strength.SUPPORTING: -1,
    Strength.MODERATE: -2,
    Strength.STRONG: -4,
    Strength.VERY_STRONG: -8,
    Strength.STAND_ALONE: -8,
}


class ClassificationCategory(str, Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    UNCERTAIN_SIGNIFICANCE = "Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


# ============================================================================
# 2. Per-criterion evidence + the classification result
# ============================================================================

@dataclass
class CriterionEvidence:
    """
    One ACMG code's status for a single variant, from whatever source
    produced it (this project's LLM pipeline, a stub, ERepo, or eventually
    another team's rule-based module - classify() does not care).
    """
    code: str
    status: CriterionStatus
    strength: Optional[Strength] = None  # required iff status == MET
    source: str = ""  # free text provenance, e.g. "llm_pipeline:PMID27532257,PMID...", "stub", "erepo_gate"

    def __post_init__(self):
        if self.code not in ALL_ACMG_CODES:
            raise ValueError(f"{self.code!r} is not a recognized ACMG/AMP 2015 code")
        if self.status == CriterionStatus.MET and self.strength is None:
            raise ValueError(f"{self.code}: status=MET requires a strength")


@dataclass
class ClassificationResult:
    category: ClassificationCategory
    score: int
    met: list[CriterionEvidence] = field(default_factory=list)
    not_met: list[CriterionEvidence] = field(default_factory=list)
    not_evaluated_codes: list[str] = field(default_factory=list)
    ba1_override: bool = False  # True if classified Benign via the BA1 stand-alone rule, not the point score


# ============================================================================
# 3. The classification function itself
# ============================================================================

def classify(evidence: list[CriterionEvidence]) -> ClassificationResult:
    """
    Combines a list of CriterionEvidence into one ClassificationResult.

    Duplicate codes are not allowed (a caller with e.g. two different
    strength estimates for PS3 must resolve that itself first - classify()
    has no basis to prefer one over the other).

    Special case: BA1 met is a stand-alone Benign call per Richards et al.
    2015 (this specific override is a categorical rule from the original
    guidelines, not part of Tavtigian's point extension) - it short-circuits
    the point score entirely, matching the standard "BA1 that isn't
    contradicted by other strong evidence is a definitive Benign" reading.
    This implementation does NOT check for a PVS1/PS-level contradiction
    before applying it (Richards et al. note BA1 can be overridden by
    conflicting strong evidence in specific circumstances) - that
    adjudication is left to the human curator, flagged via `met` still
    listing every other met criterion for their review.
    """
    by_code: dict[str, CriterionEvidence] = {}
    for e in evidence:
        if e.code in by_code:
            raise ValueError(f"duplicate evidence for {e.code} - resolve to one CriterionEvidence before calling classify()")
        by_code[e.code] = e

    met = [e for e in by_code.values() if e.status == CriterionStatus.MET]
    not_met = [e for e in by_code.values() if e.status == CriterionStatus.NOT_MET]
    not_evaluated_codes = sorted(set(ALL_ACMG_CODES) - set(by_code))

    ba1 = by_code.get("BA1")
    if ba1 is not None and ba1.status == CriterionStatus.MET:
        return ClassificationResult(
            category=ClassificationCategory.BENIGN,
            score=sum(_score_one(e) for e in met),
            met=met, not_met=not_met, not_evaluated_codes=not_evaluated_codes,
            ba1_override=True,
        )

    score = sum(_score_one(e) for e in met)
    category = _category_for_score(score)
    return ClassificationResult(category, score, met, not_met, not_evaluated_codes)


def _score_one(e: CriterionEvidence) -> int:
    table = _POINTS_PATHOGENIC if e.code in PATHOGENIC_CODES else _POINTS_BENIGN
    return table[e.strength]


def _category_for_score(score: int) -> ClassificationCategory:
    if score >= 10:
        return ClassificationCategory.PATHOGENIC
    if 6 <= score <= 9:
        return ClassificationCategory.LIKELY_PATHOGENIC
    if 0 <= score <= 5:
        return ClassificationCategory.UNCERTAIN_SIGNIFICANCE
    if -6 <= score <= -1:
        return ClassificationCategory.LIKELY_BENIGN
    return ClassificationCategory.BENIGN


# ============================================================================
# 4. Building CriterionEvidence from this project's own pipeline output
# ============================================================================

def from_aggregated_judgment(aggregated: AggregatedJudgment, code: str) -> CriterionEvidence:
    """
    Converts this pipeline's AggregatedJudgment (see acmg_pipeline/pipeline.
    py's judge_variant()) into a CriterionEvidence for `code`.

    `code` is the ACMG code this run was actually evaluating (e.g. "PS3",
    "BS3", "PP1", "BS4") - NOT necessarily the same as
    `aggregated.aggregated_direction.value`, which is whichever direction
    the evidence actually pointed (they can legitimately differ: evaluating
    BS3 but finding PS3-direction evidence means `code`="BS3" but
    aggregated_direction.value=="PS3" - a NOT_MET result for BS3, not a
    MET result mislabeled as BS3). This mirrors the same distinction
    export.py's _strength_blocks() already makes for VA-Spec output.

    not_clear -> NOT_MET (this project's design already treats not_clear as
    "the safe non-committal answer", i.e. it does not contribute pathogenic
    OR benign points - encoding it as NOT_MET rather than a separate
    "inconclusive" CriterionStatus reuses the existing gate.CriterionStatus
    enum instead of introducing a third status value that classify() would
    then have to special-case identically to NOT_MET anyway).
    """
    direction = aggregated.aggregated_direction
    relevant_pmids = [
        c.pmid for c in aggregated.contributions
        if c.result.judgment.variant_matching.match_status != MatchStatus.UNSUCCESSFUL
        and not is_not_clear(c.result.effective_direction)
    ]
    source = f"llm_pipeline:{','.join(f'PMID{p}' for p in relevant_pmids)}" if relevant_pmids else "llm_pipeline:no_usable_paper"

    if is_not_clear(direction) or direction.value != code:
        return CriterionEvidence(code=code, status=CriterionStatus.NOT_MET, source=source)

    tier = strength_tier_from_paper_count(len(relevant_pmids))
    strength = Strength(tier) if tier else Strength.SUPPORTING
    return CriterionEvidence(code=code, status=CriterionStatus.MET, strength=strength, source=source)
