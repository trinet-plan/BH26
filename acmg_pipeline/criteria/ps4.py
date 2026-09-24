"""
ps4_judgment.py
PS4 (case-control / cohort prevalence) judgment logic (prompt design +
structured output + human-curator hints), modeled directly on
ps3_bs3_judgment.py but for a completely different extraction target:
case-control/cohort counts instead of functional-assay results.

[Status - carried over from design doc section 4]
  No existing tool (InterVar, AutoPVS1, AutoGVP, BIAS-2015, AcmGENTIC,
  Diablo_ACMG - see doc/acmg_tool_comparison_v1_en.md) has independently
  validated PS4-specific automation; AcmGENTIC only covers PS3/BS3. This
  module's prompt and structured schema are a first-pass, UNVALIDATED
  design that carries over the PS3/BS3 pattern (variant-matching gate +
  structured extraction + multi-paper aggregation) to case-control data,
  per design doc section 4's own draft. Treat accordingly - this has not
  been through the multiple rounds of real-data-driven refinement that
  ps3_bs3_judgment.py's safety nets (detect_experiment_conflict etc.) went
  through.

[Key design-doc lesson this module is built around, section 4/3-3]
  A real ERepo PS4 determination is very often the aggregate of MANY papers
  (a real example: 9+ papers bundled into "20+ affected reported"), not one
  paper read in isolation - reading a single cited PMID and concluding
  "insufficient" from it alone previously produced a wrong rejection. This
  is exactly what the multi-paper aggregation in evidence_common.py already
  handles (built for PS3/BS3, reused here unchanged) - judge_variant() in
  the pipeline should be called with ALL of ERepo's evidence_pmids for a
  PS4 lookup, never just one.

[No BS4-like counterpart in this module]
  Unlike PS3/BS3 or PP1/BS4, ACMG 2015 has no dedicated "opposite" code for
  case-control evidence pointing the other way (BS2, "observed in healthy
  individuals", is a related but distinct criterion with its own different
  data source - population/health databases, not case-control literature -
  and is explicitly out of scope per the Layer 3 "manual" classification in
  doc/BH26_participant_briefing_v3_en.md). So CaseControlDirection only has
  PS4 / not_clear, not a two-sided pair like OverallDirection or
  SegregationDirection.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.inputs import variant_identity
from acmg_pipeline.vcf_record import VariantRecord

from acmg_pipeline.common import (
    MatchStatus, VariantMatchingResult, CuratorHint, FinalResult,
    PaperContribution, AggregatedJudgment,
    aggregate_multi_paper_results as _generic_aggregate_multi_paper_results,
)


def _clean_enum_token(raw: str) -> str:
    """Same defensive normalization as ps3_bs3_judgment._clean_enum_token."""
    if not isinstance(raw, str):
        return raw
    return re.split(r"[\(（]", raw.strip())[0].strip()


# ============================================================================
# 1. Structured output schema
# ============================================================================

class StudyDesign(str, Enum):
    CASE_CONTROL = "case_control"
    FAMILY_COHORT = "family_cohort"
    CASE_SERIES = "case_series"
    NOT_CASE_CONTROL = "not_case_control"


class CaseControlDirection(str, Enum):
    PS4 = "PS4"
    NOT_CLEAR = "not_clear"


@dataclass
class CaseControlData:
    """
    All fields are nullable: per design doc section 4's own prompt draft, a
    paper reporting case-control evidence in prose ("detected in 8 of 8
    affected individuals across 6 unrelated families") does not always give
    a clean denominator for the control/unaffected side, and the LLM should
    leave a field null rather than guess.

    odds_ratio/p_value were added after a real gemma-4 run (2026-09-15,
    GJB2 c.109G>A / PMID:31160754) reported strong case-control evidence
    entirely as inferential statistics ("OR of 20, 95%CI 17-24, Z=31,
    p<0.0001") with NO raw affected_carriers/total_affected counts stated
    anywhere in the text. detect_definitive_without_numbers() originally
    only checked the raw-count fields and forced this genuinely strong,
    quantified result to not_clear - a real case-control paper reporting
    significance via OR/p-value instead of raw 2x2 counts is normal
    practice, not missing data, so the raw counts alone were the wrong
    signal for "was this actually quantified".
    """
    affected_carriers: Optional[int] = None
    unaffected_carriers: Optional[int] = None
    total_affected: Optional[int] = None
    total_unaffected: Optional[int] = None
    odds_ratio: Optional[float] = None
    p_value: Optional[float] = None


@dataclass
class CaseControlEvidence:
    direction: CaseControlDirection
    strength_hint: str = "not_clear"  # fixed; this design does not target strength (same as PS3/BS3)
    rationale: str = ""


@dataclass
class PS4Judgment:
    """Data class mirroring the LLM's structured output directly."""
    variant_matching: VariantMatchingResult
    study_design: StudyDesign
    case_control_data: CaseControlData
    overall_evidence: CaseControlEvidence
    single_study_only: bool = True

    @staticmethod
    def from_json(data: dict) -> "PS4Judgment":
        vm = data["variant_matching"]
        if isinstance(vm, str):
            vm = {"match_status": vm}
        ccd = data.get("case_control_data") or {}
        return PS4Judgment(
            variant_matching=VariantMatchingResult(
                match_status=MatchStatus(_clean_enum_token(vm["match_status"])),
                match_type=vm.get("match_type"),
                confidence=vm.get("confidence", "low"),
                notes=vm.get("notes", ""),
            ),
            study_design=StudyDesign(_clean_enum_token(data.get("study_design", "not_case_control"))),
            case_control_data=CaseControlData(
                affected_carriers=ccd.get("affected_carriers"),
                unaffected_carriers=ccd.get("unaffected_carriers"),
                total_affected=ccd.get("total_affected"),
                total_unaffected=ccd.get("total_unaffected"),
                odds_ratio=ccd.get("odds_ratio"),
                p_value=ccd.get("p_value"),
            ),
            # `data.get("overall_evidence") or {}`, not `data["overall_evidence"]`:
            # same null-vs-missing quirk as ps3_bs3.py's PS3BS3Judgment.from_json()
            # (see its own comment) - Claude can emit "overall_evidence": null
            # outright, which `data["overall_evidence"]["direction"]` then crashes
            # on with "'NoneType' object is not subscriptable".
            overall_evidence=CaseControlEvidence(
                direction=CaseControlDirection(_clean_enum_token(
                    (data.get("overall_evidence") or {}).get("direction", "not_clear"))),
                strength_hint=(data.get("overall_evidence") or {}).get("strength_hint", "not_clear"),
                rationale=(data.get("overall_evidence") or {}).get("rationale", ""),
            ),
            single_study_only=data.get("single_study_only", True),
        )

    def structured_evidence_items(self) -> list[dict]:
        """
        Checklist-style items for a curator-facing UI, mirroring
        ps3_bs3_judgment.PS3BS3Judgment.structured_evidence_items() - see
        that docstring for why this exists (doc/recs for expert board.docx).
        One item per distinct piece of quantitative backing actually
        extracted (raw counts, control-side counts, summary statistic),
        plus one for the study design itself, so a curator can see exactly
        which of PS4's requirements were and weren't satisfied by this
        paper without re-reading the rationale text.
        """
        items = []
        ccd = self.case_control_data
        if ccd.affected_carriers is not None or ccd.total_affected is not None:
            has_both = ccd.affected_carriers is not None and ccd.total_affected is not None
            items.append({
                "label": "Affected-side counts",
                "checked": has_both,
                "detail": (
                    f"{ccd.affected_carriers}/{ccd.total_affected} affected individuals carry the variant"
                    if has_both else
                    f"affected_carriers={ccd.affected_carriers}, total_affected={ccd.total_affected}"
                ),
            })
        if ccd.unaffected_carriers is not None or ccd.total_unaffected is not None:
            has_both = ccd.unaffected_carriers is not None and ccd.total_unaffected is not None
            items.append({
                "label": "Control-side counts",
                "checked": has_both,
                "detail": (
                    f"{ccd.unaffected_carriers}/{ccd.total_unaffected} unaffected/control individuals carry the variant"
                    if has_both else
                    f"unaffected_carriers={ccd.unaffected_carriers}, total_unaffected={ccd.total_unaffected}"
                ),
            })
        if ccd.odds_ratio is not None or ccd.p_value is not None:
            items.append({
                "label": "Reported summary statistic",
                "checked": True,
                "detail": f"odds_ratio={ccd.odds_ratio}, p_value={ccd.p_value}",
            })
        items.append({
            "label": f"Study design: {self.study_design.value}",
            "checked": self.study_design == StudyDesign.CASE_CONTROL,
            "detail": "",
        })
        return items


# ============================================================================
# 2. Prompt template (from design doc section 4, "3-2. 抽出プロンプト(PS4専用)")
# ============================================================================

PROMPT_TEMPLATE = """\
You are a clinical genetics curator supporting ACMG/AMP PS4 (case-control \
study) evaluation.

Target variant: {gene} {hgvsc} ({hgvsp})
Also accept these equivalent notations as the same variant: {equivalents}

Full text of the paper:
{full_text}

Perform the following steps and output ONLY the JSON described at the end \
(no other explanatory text):

1. Variant identification (variant_matching):
   match_status: matched / heuristic / single_variant_study / unsuccessful
   - matched: the exact target variant (or a notational equivalent from the \
list above) is explicitly tested/reported in this paper
   - unsuccessful: the paper does not test or report the target variant at \
all (if so, stop here; the rest of the fields should reflect "no evidence")
   match_type: protein_notation / hgvsc / rsid / abbreviation, etc.
   confidence: high / medium / low
   notes: brief note on how the match was made

2. Study design (study_design):
   case_control: a defined affected group vs a defined control/unaffected \
group, both with known denominators
   family_cohort: multiple unrelated families/pedigrees pooled together \
(with or without a formal control group)
   case_series: a series of affected individuals reported, but with no \
denominator or comparison group at all
   not_case_control: the paper is not a case-control/cohort study of this \
variant (e.g. it is a pure functional-assay paper) - if so, this paper \
cannot contribute PS4 evidence; set overall_evidence.direction to \
"not_clear" and stop there

3. Case-control data extraction (case_control_data) - leave any field null \
if the paper does not give a number for it; do NOT guess or infer a number \
that is not explicitly stated. Many case-control papers report only \
inferential statistics (odds ratio, p-value) rather than the raw 2x2 \
counts - fill in whichever of these the paper actually states; it is fine \
and expected for only some of these six fields to be non-null:
   affected_carriers: number of affected individuals who carry the target \
variant
   unaffected_carriers: number of unaffected/control individuals who carry \
the target variant
   total_affected: total number of affected individuals examined
   total_unaffected: total number of unaffected/control individuals examined
   odds_ratio: the odds ratio (or risk ratio) reported for this variant, if \
stated
   p_value: the p-value (or a bound on it, e.g. 0.0001 for "p<0.0001") \
reported for the enrichment/association test, if stated

4. Self-report (single_study_only): based on citations/mentions within this \
paper alone (not an independent literature search), does this paper appear \
to be the only source of case-control evidence for this variant, or does it \
reference/build on other studies?

5. Overall judgment (overall_evidence):
   direction: "PS4" if the variant is significantly enriched among affected \
individuals compared to controls/background (the paper's own text supports \
this conclusion, not just raw counts you calculate yourself); "not_clear" \
otherwise (including: not_case_control design, or a case-control design but \
the paper's own text does not draw a significance conclusion, or you cannot \
determine this from the text)
   strength_hint: always the fixed string "not_clear" (strength is not \
targeted by this design)
   rationale: 1-3 sentences explaining the basis for the direction above, \
citing the actual reported numbers/statistics where available

Output exactly this JSON structure:
{{
  "variant_matching": {{"match_status": "matched", "match_type": "protein_notation", "confidence": "high", "notes": "..."}},
  "study_design": "case_control",
  "case_control_data": {{"affected_carriers": 8, "unaffected_carriers": 0, "total_affected": 8, "total_unaffected": 500, "odds_ratio": null, "p_value": null}},
  "single_study_only": true,
  "overall_evidence": {{"direction": "PS4", "strength_hint": "not_clear", "rationale": "..."}}
}}
"""


def build_prompt(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    full_text: str,
) -> str:
    """Build the literature prompt from the shared criterion input types."""
    del clinical_note
    gene, hgvsc, hgvsp, equivalents = variant_identity(variant)
    return PROMPT_TEMPLATE.format(
        gene=gene, hgvsc=hgvsc, hgvsp=hgvsp,
        equivalents=", ".join(equivalents), full_text=full_text,
    )


# ============================================================================
# 3. Post-processing: human-curator hint generation + safety nets
# ============================================================================

def generate_curator_hints(judgment: PS4Judgment, pmid: Optional[str] = None) -> list[CuratorHint]:
    """
    Analogous to ps3_bs3_judgment.generate_curator_hints, but for the PS4
    schema. `pmid` is optional but should be passed whenever available, for
    the same reason as in ps3_bs3_judgment.py: an unattributed "this paper"/
    "this judgment" hint becomes ambiguous the moment it is read outside the
    one FinalResult it was generated for.
    """
    hints: list[CuratorHint] = []
    paper_ref = f" (PMID:{pmid})" if pmid else ""

    if judgment.variant_matching.match_status == MatchStatus.UNSUCCESSFUL:
        return hints

    ccd = judgment.case_control_data
    if judgment.overall_evidence.direction == CaseControlDirection.PS4:
        if ccd.affected_carriers is None and ccd.total_affected is None and ccd.odds_ratio is None and ccd.p_value is None:
            hints.append(CuratorHint(
                "warning",
                f"This paper{paper_ref} was judged PS4-supporting, but "
                "neither raw counts (affected_carriers/total_affected) nor "
                "a summary statistic (odds_ratio/p_value) could be "
                "extracted. A significance claim with no extractable "
                "quantitative backing at all is hard for a curator to "
                "independently verify from this record alone.",
            ))
        if ccd.unaffected_carriers is None and ccd.total_unaffected is None and ccd.odds_ratio is None:
            hints.append(CuratorHint(
                "caution",
                f"This paper{paper_ref} does not report any control/"
                "unaffected-side numbers, nor an odds/risk ratio that would "
                "itself imply a comparison group. PS4 requires a comparison "
                "against a control/background rate, not just a count of "
                "affected carriers; confirm the comparison basis manually.",
            ))

    if judgment.study_design != StudyDesign.CASE_CONTROL and judgment.overall_evidence.direction == CaseControlDirection.PS4:
        hints.append(CuratorHint(
            "caution",
            f"This paper{paper_ref} was judged PS4-supporting despite a "
            f"study design of \"{judgment.study_design.value}\" rather than "
            "a formal case-control design. Case series / family-cohort "
            "evidence without a defined control group is generally weaker "
            "support for PS4 and may need a lower strength than a true "
            "case-control study would.",
        ))

    if judgment.single_study_only and judgment.overall_evidence.direction != CaseControlDirection.NOT_CLEAR:
        hints.append(CuratorHint(
            "caution",
            f"This judgment{paper_ref} is based on a single study only "
            "(self-reported). Design doc section 4 found a real case "
            "(MYH7 c.2155C>T, p.Arg719Trp, an HCM variant) where the true "
            "ERepo PS4 basis was 9+ pooled papers, not the single PMID "
            "initially checked - confirm via ERepo's full evidenceLinks "
            "list, not this paper alone.",
        ))

    return hints


def detect_definitive_without_numbers(judgment: PS4Judgment) -> bool:
    """
    Mirrors ps3_bs3_judgment.detect_empty_experiments_with_definitive_
    direction: a confident PS4 call with NO extractable quantitative
    backing at all - neither raw counts nor a summary statistic - is
    internally inconsistent and should be forced to not_clear, not just
    flagged.

    Deliberately checks all four of affected_carriers/total_affected/
    odds_ratio/p_value (any one being non-null is enough to pass) rather
    than just the two raw-count fields - see the CaseControlData docstring
    for the real gemma-4 run (GJB2 c.109G>A / PMID:31160754) that reported
    strong evidence purely as "OR of 20, p<0.0001" with zero raw counts
    stated, which the raw-count-only version of this check wrongly forced
    to not_clear.
    """
    if judgment.variant_matching.match_status == MatchStatus.UNSUCCESSFUL:
        return False
    if judgment.overall_evidence.direction != CaseControlDirection.PS4:
        return False
    ccd = judgment.case_control_data
    return (
        ccd.affected_carriers is None and ccd.total_affected is None
        and ccd.odds_ratio is None and ccd.p_value is None
    )


def detect_not_case_control_with_ps4(judgment: PS4Judgment) -> bool:
    """A "not_case_control" study design cannot itself support PS4 - if the LLM said PS4 anyway, that is a direct contradiction of its own study_design field."""
    return (
        judgment.study_design == StudyDesign.NOT_CASE_CONTROL
        and judgment.overall_evidence.direction == CaseControlDirection.PS4
    )


def finalize(judgment: PS4Judgment, pmid: Optional[str] = None) -> FinalResult:
    hints = generate_curator_hints(judgment, pmid=pmid)
    effective_direction = judgment.overall_evidence.direction
    paper_ref = f" (PMID:{pmid})" if pmid else ""

    if detect_not_case_control_with_ps4(judgment):
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, study_design was reported as "
            "\"not_case_control\" yet the LLM still concluded PS4. These "
            "two fields directly contradict each other, so the judgment "
            "has been forced to not_clear.",
        ))
        effective_direction = CaseControlDirection.NOT_CLEAR
    elif detect_definitive_without_numbers(judgment):
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, the LLM concluded PS4 but extracted "
            "no quantitative backing at all - no raw counts "
            "(affected_carriers/total_affected) and no summary statistic "
            "(odds_ratio/p_value). A definitive enrichment claim with "
            "nothing concrete behind it is internally inconsistent, so the "
            "judgment has been forced to not_clear.",
        ))
        effective_direction = CaseControlDirection.NOT_CLEAR

    return FinalResult(judgment=judgment, curator_hints=hints, effective_direction=effective_direction)


# ============================================================================
# 4. Multi-paper aggregation - thin wrapper, see evidence_common.py, plus a
#    pooled-case-series fallback specific to PS4 (see
#    _pooled_case_series_aggregate() below)
# ============================================================================

# A single small case series alone is exactly what the per-paper prompt
# already refuses to call PS4 on its own (no defined control/background
# comparison) - see PROMPT_TEMPLATE step 5. Requiring at least this many
# INDEPENDENT papers, each with its own stated affected-carrier count,
# before pooling is a conservative floor against one paper's case series
# masquerading as corroborated evidence; POOLED_MIN_PAPERS papers times a
# typical single-digit case series each is what actually produces a
# ClinGen-style "9+ affected reported across N papers" determination (see
# this module's own docstring, design doc section 4/3-3).
POOLED_MIN_PAPERS = 2
# ClinGen VCEPs commonly treat a handful of independently reported affected
# carriers, pooled across sources, as PS4_Supporting-level evidence (see
# ps3_bs3_ps4_gate's own note on RUNX1 c.601C>T: a real ERepo PS4 basis
# pooled 9+ papers into "20+ affected reported"). 3 is a deliberately
# conservative floor for this project's own pooling (not a VCEP-published
# number), disclosed as such in the caution hint below.
POOLED_MIN_AFFECTED = 3


def _pooled_case_series_aggregate(contributions: list[PaperContribution]) -> Optional[AggregatedJudgment]:
    """
    Real PS4 determinations are very often built by pooling several case
    reports/series across papers, none of which alone is a formal
    case-control study with a defined denominator (see this module's own
    docstring, design doc section 4/3-3: a real ERepo PS4 basis pooled 9+
    papers into "20+ affected reported"). The per-paper prompt correctly
    refuses to call a single case series PS4 on its own - it has no
    control/background comparison - so the standard aggregation above
    (which only counts papers whose OWN direction was already PS4) finds
    nothing to aggregate and returns not_clear, even when several papers
    each independently reported real affected carriers for this variant.

    This pools `case_control_data.affected_carriers` across every paper
    whose variant match succeeded, regardless of that paper's own
    study_design or direction, and returns a PS4 call only when: at least
    POOLED_MIN_PAPERS distinct papers each report a positive
    affected_carriers count, none of the pooled papers reports any
    unaffected/control carriers (a paper with both affected AND control
    counts already had every chance to be scored as a real case-control
    study through the normal path above - pooling it here too would double
    count evidence the primary path already had and rejected), and the sum
    reaches POOLED_MIN_AFFECTED. Returns None (caller keeps the original
    not_clear) when these conditions are not met, so this only ever adds a
    signal the standard aggregation missed - it never overrides one the
    standard path already reached.
    """
    pooled = [
        c for c in contributions
        if c.result.judgment.variant_matching.match_status != MatchStatus.UNSUCCESSFUL
        and (c.result.judgment.case_control_data.affected_carriers or 0) > 0
        and not c.result.judgment.case_control_data.unaffected_carriers
    ]
    if len(pooled) < POOLED_MIN_PAPERS:
        return None
    total_affected = sum(c.result.judgment.case_control_data.affected_carriers for c in pooled)
    if total_affected < POOLED_MIN_AFFECTED:
        return None
    per_paper = ", ".join(
        f"PMID:{c.pmid}={c.result.judgment.case_control_data.affected_carriers}" for c in pooled
    )
    hints = [CuratorHint(
        "caution",
        f"No single paper reached a formal case-control PS4 conclusion on "
        f"its own, but {len(pooled)} independent papers each reported real "
        f"affected carriers with no contradicting control-side count "
        f"({per_paper}), pooling to {total_affected} affected individuals "
        "total. This is this project's OWN pooling heuristic (a threshold "
        f"of {POOLED_MIN_AFFECTED}, not a VCEP-published number), not a "
        "verified ClinGen aggregate determination - confirm this pooling "
        "is appropriate (e.g. the papers are not reporting overlapping "
        "cohorts/families) before relying on it.",
    )]
    return AggregatedJudgment(contributions, CaseControlDirection.PS4, hints)


def aggregate_multi_paper_results(contributions: list[PaperContribution]) -> AggregatedJudgment:
    """
    Same generic aggregation as ps3_bs3_judgment.py's wrapper, but this is
    the criterion this design-doc lesson was actually about (section 4/3-3:
    a real PS4 determination pooled 9+ papers) - see evidence_common.
    aggregate_multi_paper_results for the logic itself. When that generic
    pass finds no single paper confidently reaching PS4 on its own, tries
    _pooled_case_series_aggregate() as a fallback before giving up to
    not_clear - see that function's own docstring for why and how.
    """
    result = _generic_aggregate_multi_paper_results(contributions, not_clear=CaseControlDirection.NOT_CLEAR)
    if result.aggregated_direction != CaseControlDirection.NOT_CLEAR:
        return result
    pooled = _pooled_case_series_aggregate(contributions)
    return pooled if pooled is not None else result
