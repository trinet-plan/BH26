"""
va_spec_export.py
Converts this project's judgment output - PS3/BS3 (ps3_bs3_judgment.py),
PS4 (ps4_judgment.py), or PP1/BS4 (segregation_judgment.py) - into a GA4GH
VA-Spec EvidenceLine (as a plain dict, ready for json.dump()).

[v2, 2026-09-15: now built on the real ga4gh.va_spec Pydantic library]
  The first version of this module (see git history / design doc section
  15-2) hand-built plain dicts shaped to match the VA-Spec docs and the
  sample fixture (va_spec_example.yaml), since neither `ga4gh.va_spec` nor
  any other VA-Spec schema library was installed or imported. That worked,
  but was never actually validated against the real schema - it could have
  had undetected drift (wrong field shapes, wrong enum values) anywhere.
  This version imports and constructs the real Pydantic models from
  `ga4gh.va_spec` (`pip install ga4gh.va_spec`) instead, so every
  EvidenceLine this module produces is guaranteed to satisfy the actual
  schema (Pydantic validation raises immediately on construction if not).
  Two real discrepancies were caught by switching to the real models:
    - `Contribution.activityType` is a plain `str`, not a MappableConcept
      (the old code nested it as `{"name": "..."}`, which was wrong).
    - `Contribution.contributor` is an `Agent` model, not an arbitrary dict.
  Everything else (Direction enum values "supports"/"neutral"/"disputes",
  the "ACMG Guidelines, 2015" System string, the Document/Method/
  MappableConcept/Coding/Extension field names) matched what had already
  been hand-built, confirming the original reading of the docs was
  otherwise accurate.

[Criterion-agnostic by design]
  This module only imports the generic shapes from evidence_common.py
  (MatchStatus, CuratorHint, PaperContribution, AggregatedJudgment) - never
  a criterion-specific direction Enum (OverallDirection, or whatever PS4 /
  PP1&BS4 define). "Is this the neutral/inconclusive case" is checked as
  `direction.value == "not_clear"` rather than `direction == SomeEnum.
  NOT_CLEAR`, since every direction Enum in this project uses that literal
  string for its neutral member (see evidence_common.py's own docstring).
  This is what lets build_evidence_line() work unmodified for any of the
  three judgment modules above, and any future one that follows the same
  convention.

[Scope]
  Builds one EvidenceLine per (variant, criterion) - aggregated across every
  cited paper examined for that criterion - not a full VariantPathogenicity
  Statement/Proposition for the whole variant. That wrapping is a separate,
  later concern; the immediate goal agreed with the user is the EvidenceLine
  itself.

[Design decisions, finalized in conversation with the user (2026-09-15),
 superseding the informal v3-era mapping in design doc section 5]
  - Granularity: one top-level EvidenceLine per (variant, criterion),
    corresponding to this pipeline's AggregatedJudgment. Each paper actually
    examined (PaperContribution) becomes a nested EvidenceLine under
    hasEvidenceItems, so the aggregated conclusion and the individual
    per-paper evidence both survive in the same document (the real schema's
    `hasEvidenceItems` explicitly allows nesting another EvidenceLine, not
    just StudyResult/Statement - confirmed via EvidenceLine.model_fields).
  - not_clear -> still emit an EvidenceLine (do not drop it / route it to a
    separate "Evidence Gap" queue as design doc v3 proposed). Use
    directionOfEvidenceProvided = Direction.NEUTRAL, one of exactly three
    real enum members (supports/neutral/disputes - confirmed via
    ga4gh.va_spec.base.enums.Direction) for exactly this case. reportedIn
    (the papers actually examined) and description (why no confident
    direction could be reached) are still populated; strength fields are
    omitted since no strength claim is being made.
  - Curator hints (structured severity+message, no dedicated VA-Spec field)
    go into `extensions`. A short free-text summary goes into the inherited
    `description` field (InformationEntity.description).
  - Paper references (PMID/URL) go into `reportedIn` (Document).
  - Strength: the LLM never actually judges strength (see
    ps3_bs3_judgment.OverallEvidence.strength_hint, which the prompt fixes to
    the literal string "not_clear" - see PROMPT_TEMPLATE - so it carries no
    real signal and is deliberately NOT used here). Instead,
    strengthOfEvidenceProvided/evidenceOutcome are populated from a rough,
    explicitly-flagged heuristic: the count of independently agreeing papers
    that aggregate_multi_paper_results() already computes (1->Supporting,
    2->Moderate, 3+->Strong - all three are real values from
    ga4gh.va_spec.base.enums.STRENGTH_OF_EVIDENCE_PROVIDED_VALUES). This is
    NOT a VCEP-calibrated ACMG strength determination - every EvidenceLine
    that carries a strength value also carries a `strengthEstimationMethod`
    extension saying so, so a human curator does not mistake it for an
    authoritative call.
  - structuredEvidenceItems: a checklist-style extension (see
    _structured_items_extension below) added 2026-09-15 in response to a
    real curator UI requirements doc (doc/recs for expert board.docx).
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

from ga4gh.core.models import Coding, Extension, MappableConcept
from ga4gh.va_spec.base.core import Agent, Contribution, Direction, Document, EvidenceLine, Method
from ga4gh.va_spec.base.enums import System

from acmg_pipeline.common import (
    AggregatedJudgment, CuratorHint, MatchStatus, PaperContribution,
    is_not_clear, strength_tier_from_paper_count,
)
from acmg_pipeline.inputs import variant_identity
from acmg_pipeline.vcf_record import VariantRecord

PIPELINE_AGENT = Agent(
    id="acmg-literature-llm-pipeline",
    name=(
        "ACMG literature-evidence LLM pipeline (gemma-4 via vLLM) - "
        "AI-generated draft evidence, requires human curator review before use"
    ),
)


def _pubmed_url(pmid: str) -> str:
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def _document(pmid: str) -> Document:
    return Document(pmid=pmid, urls=[_pubmed_url(pmid)])


def _direction_of_evidence(direction, criterion: str) -> Direction:
    if is_not_clear(direction):
        return Direction.NEUTRAL
    return Direction.SUPPORTS if direction.value == criterion else Direction.DISPUTES


def _hints_extension(hints: list[CuratorHint]) -> Optional[Extension]:
    if not hints:
        return None
    return Extension(
        name="curatorHints",
        value=[{"severity": h.severity, "message": h.message} for h in hints],
    )


def _structured_items_extension(judgment) -> Optional[Extension]:
    """
    Pulls the criterion-specific checklist items off `judgment` via duck
    typing - `structured_evidence_items()` is an OPTIONAL method each
    criterion module's Judgment class may implement (PS3BS3Judgment,
    PS4Judgment, SegregationJudgment all do); this module never imports any
    of those classes, so a getattr probe is how it stays criterion-agnostic
    the same way common.is_not_clear() does for direction.

    Each item is {"label": str, "checked": bool, "detail": str} - see
    design doc section 15-4-1 for the full convention.
    """
    getter = getattr(judgment, "structured_evidence_items", None)
    if not callable(getter):
        return None
    items = getter()
    if not items:
        return None
    return Extension(name="structuredEvidenceItems", value=items)


def _strength_blocks(
    direction, n_agreeing_papers: int,
) -> tuple[Optional[MappableConcept], Optional[MappableConcept], Optional[Extension]]:
    """
    Returns (strengthOfEvidenceProvided, evidenceOutcome, heuristic-disclosure
    extension), any of which may be None.

    The ACMG code is built from the direction actually established
    (`direction.value`, e.g. "PS3" or "BS3"), NOT from whichever criterion
    the caller happened to be testing for. Using the tested criterion
    unconditionally would mislabel a disputing result - e.g. testing PS3 but
    finding BS3-direction evidence must be coded "BS3_moderate: Criterion
    Met", not "PS3_moderate: Criterion Met" (which would claim the opposite
    of what was found).
    """
    if is_not_clear(direction):
        return None, None, None
    tier = strength_tier_from_paper_count(n_agreeing_papers)
    if tier is None:
        return None, None, None
    code_base = direction.value
    strength = MappableConcept(primaryCoding=Coding(code=tier, system=System.ACMG.value))
    outcome = MappableConcept(
        primaryCoding=Coding(code=f"{code_base}_{tier}", system=System.ACMG.value),
        name=f"ACMG 2015 {code_base} {tier} Criterion Met (heuristic estimate, unconfirmed)",
    )
    disclosure = Extension(
        name="strengthEstimationMethod",
        value=(
            f"heuristic: {n_agreeing_papers} independent paper(s) with agreeing "
            "direction (1=supporting, 2=moderate, 3+=strong). This is a rough "
            "proxy, NOT a VCEP-calibrated ACMG strength determination - the "
            "underlying LLM judgment only assesses direction, never strength "
            "- and requires human curator confirmation."
        ),
    )
    return strength, outcome, disclosure


def build_paper_evidence_line(contribution: PaperContribution, criterion: str, index: int) -> EvidenceLine:
    """Builds the nested, per-paper EvidenceLine for one PaperContribution."""
    result = contribution.result
    hints_ext = _hints_extension(result.curator_hints)
    items_ext = _structured_items_extension(result.judgment)
    extensions = [e for e in (hints_ext, items_ext) if e] or None

    return EvidenceLine(
        id=f"evline:paper-{index}-{contribution.pmid}",
        directionOfEvidenceProvided=_direction_of_evidence(result.effective_direction, criterion),
        reportedIn=[_document(contribution.pmid)],
        description=result.judgment.overall_evidence.rationale or None,
        extensions=extensions,
    )


def build_evidence_line(
    aggregated: AggregatedJudgment,
    variant: VariantRecord,
    criterion: str,
    vcep_name: Optional[str] = None,
) -> dict:
    """
    Builds the top-level, aggregated VA-Spec EvidenceLine for one
    (variant, criterion) pair from this pipeline's AggregatedJudgment, and
    returns it as a JSON-ready dict (`.model_dump(mode="json")` on a real,
    schema-validated ga4gh.va_spec.base.core.EvidenceLine instance - not a
    hand-built dict).
    """
    gene, hgvsc, _, _ = variant_identity(variant)
    direction = aggregated.aggregated_direction
    all_pmids = [c.pmid for c in aggregated.contributions]
    relevant_pmids = [
        c.pmid for c in aggregated.contributions
        if c.result.judgment.variant_matching.match_status != MatchStatus.UNSUCCESSFUL
        and not is_not_clear(c.result.effective_direction)
    ]

    strength, outcome, strength_disclosure = _strength_blocks(direction, len(relevant_pmids))
    hints_ext = _hints_extension(aggregated.aggregation_hints)
    extensions = [e for e in (hints_ext, strength_disclosure) if e] or None

    description = " ".join(h.message for h in aggregated.aggregation_hints) or None

    safe_hgvsc = (
        hgvsc.replace(">", "_").replace(".", "_").replace("+", "p").replace("-", "m")
    )

    evidence_line = EvidenceLine(
        id=f"evline:{gene}_{safe_hgvsc}_{criterion}",
        directionOfEvidenceProvided=_direction_of_evidence(direction, criterion),
        reportedIn=[_document(pmid) for pmid in all_pmids] or None,
        description=description,
        hasEvidenceItems=[
            build_paper_evidence_line(c, criterion, i)
            for i, c in enumerate(aggregated.contributions, start=1)
        ] or None,
        strengthOfEvidenceProvided=strength,
        evidenceOutcome=outcome,
        specifiedBy=Method(
            methodType=criterion,
            name=(
                f"{vcep_name or 'Unspecified VCEP'} ACMG/AMP {criterion} "
                "evidence assessment (LLM-assisted draft)"
            ),
        ),
        contributions=[Contribution(
            contributor=PIPELINE_AGENT,
            activityType="automated evidence evaluation",
            date=datetime.now(timezone.utc),
        )],
        extensions=extensions,
    )
    return evidence_line.model_dump(mode="json", exclude_none=True)
