"""
segregation_judgment.py
PP1/BS4 (family segregation) judgment logic (prompt design + structured
output + human-curator hints), modeled on ps3_bs3_judgment.py but for a
third extraction target: per-family pedigree genotype/phenotype counts,
instead of functional-assay results (PS3/BS3) or case-control counts (PS4).

[Why PP1 and BS4 share one module]
  ACMG 2015 treats these as a direction pair over the same underlying data:
  a family showing the variant present in every affected member and absent
  in every unaffected member is PP1 (cosegregation, pathogenic-supporting);
  a family showing the variant in an unaffected member, or missing from an
  affected member, is BS4 (lack of segregation, benign-supporting). This is
  structurally the same PS3/BS3 pattern (one paper's pedigree data can point
  either way), so - exactly like OverallDirection for PS3/BS3 -
  SegregationDirection is PP1 / BS4 / not_clear, not two separate modules.

[Status - new, no design-doc prior art]
  Unlike PS4 (which had a draft prompt in design doc section 4) or PS3/BS3
  (extensively refined against real ERepo/CGBench data), this schema and
  prompt are a first pass with no prior validation against real segregation
  papers. The only real data point available in this project's local CGBench
  pool (clingen_vci_pubmed_fulltext_dedup_pmid_CORRECTED.csv) was 4 BS4 rows
  (USH2A x2, GJB2, RUNX1), all met_status=not_met - not enough to validate a
  design against. Treat this module's safety nets as reasonable first
  guesses to be revised once run against real data, the same way PS3/BS3's
  detect_experiment_conflict etc. were only added after real gemma-4 runs
  surfaced actual failure modes.

[Scope note]
  Per the user's decision (2026-09-15): only the criteria whose supporting
  evidence is actually literature (papers on PubMed) are in scope for this
  LLM pipeline - PS3, BS3, PS4, PP1, BS4. The other "manual"/Layer-3 ACMG
  codes (PS2, PM3, PM6, BS2, BP2, BP5 - see doc/BH26_participant_briefing_
  v3_en.md line 66) depend on the patient's own clinical/genetic-testing
  records, not on papers, and are out of scope here.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

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

class SegregationDirection(str, Enum):
    PP1 = "PP1"
    BS4 = "BS4"
    NOT_CLEAR = "not_clear"


@dataclass
class FamilySegregationData:
    """
    One family/pedigree's genotype-phenotype counts. All four counts are
    nullable (a paper may describe a family narratively - "the variant was
    found in all three affected siblings but neither unaffected parent" -
    without ever stating denominators as bare numbers); leave null rather
    than infer.
    """
    family_id: str = ""
    affected_with_variant: Optional[int] = None
    affected_without_variant: Optional[int] = None  # key BS4 signal
    unaffected_with_variant: Optional[int] = None    # key BS4 signal
    unaffected_without_variant: Optional[int] = None
    informative_meioses: Optional[int] = None  # if the paper reports a LOD score/meiosis count directly
    notes: str = ""


@dataclass
class SegregationEvidence:
    direction: SegregationDirection
    strength_hint: str = "not_clear"  # fixed; this design does not target strength (same as PS3/BS3)
    rationale: str = ""


@dataclass
class SegregationJudgment:
    """Data class mirroring the LLM's structured output directly."""
    variant_matching: VariantMatchingResult
    families: list[FamilySegregationData]
    overall_evidence: SegregationEvidence
    single_study_only: bool = True

    @staticmethod
    def from_json(data: dict) -> "SegregationJudgment":
        vm = data["variant_matching"]
        if isinstance(vm, str):
            vm = {"match_status": vm}
        return SegregationJudgment(
            variant_matching=VariantMatchingResult(
                match_status=MatchStatus(_clean_enum_token(vm["match_status"])),
                match_type=vm.get("match_type"),
                confidence=vm.get("confidence", "low"),
                notes=vm.get("notes", ""),
            ),
            families=[
                FamilySegregationData(
                    family_id=fam.get("family_id", ""),
                    affected_with_variant=fam.get("affected_with_variant"),
                    affected_without_variant=fam.get("affected_without_variant"),
                    unaffected_with_variant=fam.get("unaffected_with_variant"),
                    unaffected_without_variant=fam.get("unaffected_without_variant"),
                    informative_meioses=fam.get("informative_meioses"),
                    notes=fam.get("notes", ""),
                )
                for fam in data.get("families", [])
            ],
            overall_evidence=SegregationEvidence(
                direction=SegregationDirection(_clean_enum_token(data["overall_evidence"]["direction"])),
                strength_hint=data["overall_evidence"].get("strength_hint", "not_clear"),
                rationale=data["overall_evidence"].get("rationale", ""),
            ),
            single_study_only=data.get("single_study_only", True),
        )

    def structured_evidence_items(self) -> list[dict]:
        """
        One checklist-style item per family, for a curator-facing UI -
        mirrors ps3_bs3_judgment.PS3BS3Judgment.structured_evidence_items()
        (see that docstring for why this exists: doc/recs for expert
        board.docx). `checked` means "this family's own extracted counts
        are consistent with full segregation (no affected-without-variant
        or unaffected-with-variant individuals)" - a fact about that one
        family, independent of the paper's overall adopted direction (a
        single non-segregating family is exactly the kind of per-item
        detail a curator needs to see directly, not buried in prose).
        """
        items = []
        for fam in self.families:
            non_segregating = (fam.affected_without_variant or 0) > 0 or (fam.unaffected_with_variant or 0) > 0
            detail = (
                f"affected_with_variant={fam.affected_with_variant}, "
                f"affected_without_variant={fam.affected_without_variant}, "
                f"unaffected_with_variant={fam.unaffected_with_variant}, "
                f"unaffected_without_variant={fam.unaffected_without_variant}"
            )
            if fam.informative_meioses is not None:
                detail += f", informative_meioses={fam.informative_meioses}"
            if fam.notes:
                detail += f" ({fam.notes})"
            items.append({
                "label": f"Family {fam.family_id}" if fam.family_id else "Family (unnamed)",
                "checked": not non_segregating,
                "detail": detail,
            })
        return items


# ============================================================================
# 2. Prompt template
# ============================================================================

PROMPT_TEMPLATE = """\
You are a clinical genetics curator supporting ACMG/AMP PP1 (cosegregation) \
and BS4 (lack of segregation) evaluation. Both criteria are assessed from \
the SAME underlying family/pedigree data - PP1 if the variant segregates \
perfectly with disease, BS4 if it clearly does not.

Target variant: {gene} {hgvsc} ({hgvsp})
Also accept these equivalent notations as the same variant: {equivalents}

Full text of the paper:
{full_text}

Perform the following steps and output ONLY the JSON described at the end \
(no other explanatory text):

1. Variant identification (variant_matching):
   match_status: matched / heuristic / single_variant_study / unsuccessful
   - matched: the exact target variant (or a notational equivalent from the \
list above) is explicitly reported as segregating (or not) within at least \
one family/pedigree in this paper
   - unsuccessful: the paper does not report any family/pedigree \
segregation data for the target variant at all (if so, stop here; the rest \
of the fields should reflect "no evidence")
   match_type: protein_notation / hgvsc / rsid / abbreviation, etc.
   confidence: high / medium / low
   notes: brief note on how the match was made

2. Per-family extraction (families) - one entry per distinct family/\
pedigree reported, each with counts left null if the paper does not state \
that specific number (do NOT guess or infer a number that is not \
explicitly stated):
   family_id: however the paper identifies this family (e.g. "Family 1", \
"kindred A", or a short description if unnamed)
   affected_with_variant: affected individuals in this family who carry \
the target variant
   affected_without_variant: affected individuals in this family who do \
NOT carry the target variant (this is the single most important BS4 \
signal - a nonzero value here directly indicates non-segregation)
   unaffected_with_variant: unaffected individuals in this family who DO \
carry the target variant (also a key BS4 signal, especially past the \
age of onset for the condition)
   unaffected_without_variant: unaffected individuals in this family who \
do not carry the target variant
   informative_meioses: only if the paper explicitly states a meiosis \
count or LOD score for this family; otherwise null
   notes: anything relevant not captured above (e.g. reduced penetrance \
caveats, phenocopies, incomplete testing of some family members)

3. Self-report (single_study_only): based on citations/mentions within this \
paper alone (not an independent literature search), does this paper appear \
to be the only source of segregation evidence for this variant, or does it \
reference/build on other published families?

4. Overall judgment (overall_evidence):
   direction:
     "PP1" if every informative family member examined is consistent with \
full cosegregation (every affected individual carries the variant, every \
unaffected individual does not), across all families reported
     "BS4" if at least one family clearly shows non-segregation (an \
affected individual without the variant, or an unaffected individual with \
it, not explained away by an established phenocopy/reduced-penetrance \
mechanism the paper itself describes)
     "not_clear" otherwise (including: no families reported, data too \
sparse/ambiguous to call either way, or you cannot determine this from the \
text)
   strength_hint: always the fixed string "not_clear" (strength is not \
targeted by this design)
   rationale: 1-3 sentences explaining the basis for the direction above, \
citing the actual per-family counts where available

Output exactly this JSON structure:
{{
  "variant_matching": {{"match_status": "matched", "match_type": "protein_notation", "confidence": "high", "notes": "..."}},
  "families": [
    {{"family_id": "Family 1", "affected_with_variant": 3, "affected_without_variant": 0, "unaffected_with_variant": 0, "unaffected_without_variant": 2, "informative_meioses": null, "notes": "..."}}
  ],
  "single_study_only": true,
  "overall_evidence": {{"direction": "PP1", "strength_hint": "not_clear", "rationale": "..."}}
}}
"""


def build_prompt(gene: str, hgvsc: str, hgvsp: str, equivalents: list[str], full_text: str) -> str:
    return PROMPT_TEMPLATE.format(
        gene=gene, hgvsc=hgvsc, hgvsp=hgvsp,
        equivalents=", ".join(equivalents), full_text=full_text,
    )


# ============================================================================
# 3. Post-processing: human-curator hint generation + safety nets
# ============================================================================

def _has_non_segregation_signal(families: list[FamilySegregationData]) -> bool:
    return any(
        (fam.affected_without_variant or 0) > 0 or (fam.unaffected_with_variant or 0) > 0
        for fam in families
    )


def _has_full_segregation_signal(families: list[FamilySegregationData]) -> bool:
    """
    At least one family has some informative count AND no family shows a
    non-segregation signal.
    """
    has_any_count = any(
        (fam.affected_with_variant or 0) > 0 or (fam.unaffected_without_variant or 0) > 0
        for fam in families
    )
    return has_any_count and not _has_non_segregation_signal(families)


def generate_curator_hints(judgment: SegregationJudgment, pmid: Optional[str] = None) -> list[CuratorHint]:
    """Analogous to ps3_bs3_judgment.generate_curator_hints, for the PP1/BS4 schema."""
    hints: list[CuratorHint] = []
    paper_ref = f" (PMID:{pmid})" if pmid else ""

    if judgment.variant_matching.match_status == MatchStatus.UNSUCCESSFUL:
        return hints

    if len(judgment.families) > 1:
        hints.append(CuratorHint(
            "info",
            f"This paper{paper_ref} reports {len(judgment.families)} "
            "separate families. If they disagree on segregation, that "
            "should already be reflected in the overall direction below - "
            "check the per-family notes if the direction is not_clear.",
        ))

    for fam in judgment.families:
        if fam.notes and re.search(r"phenocop|reduced penetrance|incomplete", fam.notes, re.IGNORECASE):
            hints.append(CuratorHint(
                "caution",
                f"Family \"{fam.family_id}\"{paper_ref} has a noted caveat "
                f"({fam.notes}) that may explain an apparent non-"
                "segregation observation without it actually counting "
                "against the variant - confirm the paper's own reasoning "
                "before treating this family's counts as BS4 evidence.",
            ))

    if judgment.single_study_only and judgment.overall_evidence.direction != SegregationDirection.NOT_CLEAR:
        hints.append(CuratorHint(
            "caution",
            f"This judgment{paper_ref} is based on a single study only "
            "(self-reported). As with PS4 (design doc section 4), a "
            "variant's true segregation evidence can be pooled across "
            "several published families - confirm via ERepo's full "
            "evidenceLinks list, not this paper alone.",
        ))

    return hints


def detect_segregation_direction_conflict(judgment: SegregationJudgment) -> bool:
    """
    Mirrors ps3_bs3_judgment.detect_experiment_conflict: the LLM concluded
    PP1 (full cosegregation) despite at least one family's own extracted
    counts showing a non-segregation signal (or vice versa: concluded BS4
    despite every family's counts being consistent with full segregation).
    This is the single most important safety net for this schema, since PP1
    vs BS4 is decided by the sign of exactly these numbers.
    """
    direction = judgment.overall_evidence.direction
    if direction == SegregationDirection.PP1 and _has_non_segregation_signal(judgment.families):
        return True
    if direction == SegregationDirection.BS4 and _has_full_segregation_signal(judgment.families) and not _has_non_segregation_signal(judgment.families):
        return True
    return False


def detect_no_families_with_definitive_direction(judgment: SegregationJudgment) -> bool:
    """Mirrors ps3_bs3_judgment.detect_empty_experiments_with_definitive_direction."""
    if judgment.variant_matching.match_status == MatchStatus.UNSUCCESSFUL:
        return False
    return (
        len(judgment.families) == 0
        and judgment.overall_evidence.direction != SegregationDirection.NOT_CLEAR
    )


def finalize(judgment: SegregationJudgment, pmid: Optional[str] = None) -> FinalResult:
    hints = generate_curator_hints(judgment, pmid=pmid)
    effective_direction = judgment.overall_evidence.direction
    paper_ref = f" (PMID:{pmid})" if pmid else ""

    if detect_segregation_direction_conflict(judgment):
        counts = "; ".join(
            f"{fam.family_id}: affected_without_variant={fam.affected_without_variant}, "
            f"unaffected_with_variant={fam.unaffected_with_variant}"
            for fam in judgment.families
        )
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, the per-family counts ({counts}) "
            f"are inconsistent with the concluded direction "
            f"({effective_direction.value}). This conflict was detected "
            "automatically and the judgment has been forced to not_clear.",
        ))
        effective_direction = SegregationDirection.NOT_CLEAR
    elif detect_no_families_with_definitive_direction(judgment):
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, match_status was reported as "
            f"{judgment.variant_matching.match_status.value}, but zero "
            f"families were extracted, and the LLM still concluded "
            f"{effective_direction.value}. A definitive segregation "
            "direction with no concrete supporting family is internally "
            "inconsistent, so the judgment has been forced to not_clear.",
        ))
        effective_direction = SegregationDirection.NOT_CLEAR

    return FinalResult(judgment=judgment, curator_hints=hints, effective_direction=effective_direction)


# ============================================================================
# 4. Multi-paper aggregation - thin wrapper, see evidence_common.py
# ============================================================================

def aggregate_multi_paper_results(contributions: list[PaperContribution]) -> AggregatedJudgment:
    return _generic_aggregate_multi_paper_results(contributions, not_clear=SegregationDirection.NOT_CLEAR)
