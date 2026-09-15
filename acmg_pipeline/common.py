"""
evidence_common.py
Infrastructure shared across every criterion-specific judgment module in this
project (ps3_bs3_judgment.py for PS3/BS3, ps4_judgment.py for PS4,
segregation_judgment.py for PP1/BS4, and any future ones).

[Why this module exists]
  ps3_bs3_judgment.py originally defined all of this itself, since PS3/BS3
  was the only criterion pair implemented. Extracted here (2026-09-15) when
  PS4 and PP1/BS4 were added, so those criteria don't have to duplicate the
  same per-paper result bookkeeping, human-curator-hint plumbing, and
  multi-paper aggregation logic under a different name. None of what's here
  is specific to functional-assay evidence - it only assumes:
    - a paper can either be matched to the target variant or not
      (MatchStatus / VariantMatchingResult)
    - a per-paper judgment produces some "effective direction" value with a
      literal string member "not_clear" as its neutral/inconclusive case
      (see aggregate_multi_paper_results' `not_clear` parameter - this is
      why it takes the sentinel as an argument instead of importing a
      specific enum: OverallDirection for PS3/BS3, a PS4-specific direction
      enum, and a PP1/BS4 direction enum are all different Python Enum
      classes, but every one of them follows this "not_clear" convention)
    - that judgment can carry a list of free-text CuratorHint objects

[Not here]
  Anything about what the LLM is actually asked to extract from a paper
  (the structured JSON schema, the prompt template, the criterion-specific
  safety-net heuristics like PS3/BS3's detect_experiment_conflict) stays in
  each criterion-specific module, since that part genuinely differs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional


class MatchStatus(str, Enum):
    MATCHED = "matched"
    HEURISTIC = "heuristic"
    SINGLE_VARIANT_STUDY = "single_variant_study"
    UNSUCCESSFUL = "unsuccessful"


@dataclass
class VariantMatchingResult:
    match_status: MatchStatus
    match_type: Optional[str] = None  # protein_notation / hgvsc / rsid / abbreviation, etc.
    confidence: str = "low"  # high/medium/low
    notes: str = ""


@dataclass
class CuratorHint:
    severity: str  # "info" / "caution" / "warning"
    message: str


@dataclass
class FinalResult:
    """
    `judgment` and `effective_direction` are intentionally left as `Any`:
    which structured-output dataclass and which direction Enum they hold
    depends on which criterion-specific module produced this FinalResult
    (PS3BS3Judgment/OverallDirection, CaseControlJudgment/CaseControlDirection,
    SegregationJudgment/SegregationDirection, ...). Every caller that reads
    `judgment` needs to already know which module it came from (it received
    the FinalResult from that module's own finalize() in the first place);
    every caller that reads `effective_direction` only needs `.value` and
    equality against that module's own NOT_CLEAR member, both of which any
    of these Enums provide.
    """
    judgment: Any
    curator_hints: list[CuratorHint] = field(default_factory=list)
    # Separate from the LLM's raw judgment: this is "the judgment that should
    # actually be used", reflecting any forced override from a criterion-
    # specific safety net. Downstream processing and scoring against ground
    # truth should use this field, not the raw one.
    effective_direction: Any = None


@dataclass
class PaperContribution:
    pmid: str
    result: FinalResult


@dataclass
class AggregatedJudgment:
    contributions: list[PaperContribution]
    aggregated_direction: Any
    aggregation_hints: list[CuratorHint]


def aggregate_multi_paper_results(contributions: list[PaperContribution], not_clear: Any) -> AggregatedJudgment:
    """
    Combines the per-paper FinalResult for each cited/searched PMID into a
    single variant-level direction, for whichever criterion this is being
    run for. `not_clear` must be that criterion module's own "no confident
    direction" Enum member (e.g. OverallDirection.NOT_CLEAR) - passed in
    explicitly, rather than imported, so this one implementation works for
    every criterion-specific direction Enum instead of each module needing
    its own copy of this logic.

    A paper is counted as "relevant" only if its own variant-matching
    succeeded (match_status != unsuccessful) AND its own effective_direction
    (after that module's own per-paper safety nets) is not itself not_clear.
    Papers that turned out to be about something else, or whose own
    judgment was already forced to not_clear, do not count toward either
    the agreement or the conflict below.

      - 0 relevant papers -> not_clear (no usable evidence found in any of
        the papers examined)
      - all relevant papers agree on one direction -> that direction is
        adopted; the number of agreeing papers is reported so a human
        curator can judge strength (this replaces any LLM self-report of
        "is this a single study" with an actual count)
      - relevant papers disagree -> not_clear, flagged as a cross-study
        conflict requiring human adjudication
    """
    hints: list[CuratorHint] = []
    relevant = [
        c for c in contributions
        if c.result.judgment.variant_matching.match_status != MatchStatus.UNSUCCESSFUL
        and c.result.effective_direction != not_clear
    ]

    if not relevant:
        hints.append(CuratorHint(
            "info",
            f"None of the {len(contributions)} paper(s) examined yielded a "
            "usable evidence-based judgment for this variant (either the "
            "variant could not be identified in the text, or the per-paper "
            "judgment was itself inconclusive).",
        ))
        return AggregatedJudgment(contributions, not_clear, hints)

    directions = {c.result.effective_direction for c in relevant}
    if len(directions) > 1:
        conflicting_pmids = ", ".join(f"{c.pmid}={c.result.effective_direction.value}" for c in relevant)
        hints.append(CuratorHint(
            "warning",
            f"Across {len(relevant)} independent papers with usable "
            f"evidence, the direction conflicts ({conflicting_pmids}). "
            "Forced to not_clear pending human adjudication of which "
            "paper's evidence should take precedence.",
        ))
        return AggregatedJudgment(contributions, not_clear, hints)

    agreed_direction = next(iter(directions))
    if len(relevant) == 1:
        hints.append(CuratorHint(
            "caution",
            f"Only 1 of {len(contributions)} paper(s) examined yielded "
            "usable evidence (this count is computed directly from the "
            "papers actually fetched and judged, not from any self-report "
            "on a single paper).",
        ))
    else:
        hints.append(CuratorHint(
            "info",
            f"{len(relevant)} independent papers agree on "
            f"{agreed_direction.value}, which is stronger support than a "
            "single-paper judgment alone.",
        ))
    return AggregatedJudgment(contributions, agreed_direction, hints)


def is_not_clear(direction: Any) -> bool:
    """
    True iff `direction` is that criterion module's own neutral/inconclusive
    member. Compares by value (not identity/type) so this works for any of
    the direction Enums in this project (every one of them uses the literal
    string "not_clear" for its neutral member - see this module's own
    docstring). Shared by export.py and classification.py (2026-09-15) so
    the same check isn't duplicated in both places.
    """
    return direction.value == "not_clear"


def strength_tier_from_paper_count(n_agreeing_papers: int) -> Optional[str]:
    """
    1 agreeing paper -> "supporting", 2 -> "moderate", 3+ -> "strong" (the
    real ga4gh.va_spec StrengthOfEvidenceProvided vocabulary, minus
    "very strong"/"standalone" which this heuristic never reaches). Shared
    by export.py (VA-Spec strengthOfEvidenceProvided) and classification.py
    (Tavtigian point-scale strength) so both use the exact same heuristic -
    see export.py's own docstring for why this is a heuristic and not a
    real VCEP-calibrated strength determination.
    """
    if n_agreeing_papers <= 0:
        return None
    if n_agreeing_papers == 1:
        return "supporting"
    if n_agreeing_papers == 2:
        return "moderate"
    return "strong"
