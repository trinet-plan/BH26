"""
bp5.py
BP5 ("variant found in a case with an alternate molecular basis for
disease") judgment logic - modeled directly on ps4.py, but with a
completely different extraction target: does this paper report that the
patient carrying the target variant was ALSO found to carry a
pathogenic/likely pathogenic variant in a DIFFERENT gene that the paper's
own authors consider to (better) explain the patient's phenotype.

[Why this is feasible where PM3/BP2 are not]
  PM3 ("in trans with a pathogenic variant") and BP2 ("in trans/cis with a
  pathogenic variant") both need case-level PHASING data this project has no
  source for at all - ERepo's own API only ever returns the VCEP's final
  MET/NOT_MET verdict for those codes (confirmed live, 2026-09-22: querying
  ERepo for a real PM3-met variant returns just
  {"status": "Met", "label": "PM3_Strong"}, no trans-variant detail),
  and the pipeline's own request contract is a single variant per case
  (ApiCaseInput / VariantRecord), so there is no second variant to check
  phase against even in principle.

  BP5 is different: "an alternate molecular basis was found" is exactly the
  kind of narrative signal a case-report paper states in prose (e.g. "the
  patient's phenotype was ultimately attributed to a de novo pathogenic
  variant in GENE2; the GENE1 variant reported here is considered
  incidental/not causative"), which is the same shape of signal PS3/BS3/PS4
  already extract from full-text via judge_single_paper() - so this reuses
  that proven infrastructure unchanged, just with a different prompt and
  schema.

[Unlike PS4, no BS-pair]
  Same situation as ps4.py's own note: ACMG 2015 has no dedicated "opposite"
  of BP5 in the direction sense - a paper either does or doesn't report an
  alternate molecular basis for THIS proband. AlternateBasisDirection only
  has BP5 / not_clear.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
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
    """Same defensive normalization as ps3_bs3.py/ps4.py's own _clean_enum_token."""
    if not isinstance(raw, str):
        return raw
    return re.split(r"[\(（]", raw.strip())[0].strip()


# ============================================================================
# 1. Structured output schema
# ============================================================================

class AlternateBasisDirection(str, Enum):
    BP5 = "BP5"
    NOT_CLEAR = "not_clear"


@dataclass
class AlternateDiagnosisData:
    """
    All fields nullable, same convention as ps4.py's CaseControlData: a
    paper reporting an alternate diagnosis in prose does not always name
    the gene precisely or state an exact variant, and the LLM should leave
    a field null rather than guess.
    """
    alternate_gene: Optional[str] = None
    alternate_variant: Optional[str] = None
    # Does the paper's own text say the alternate finding explains/accounts
    # for the patient's phenotype (not just an incidental secondary
    # finding of uncertain relevance)?
    phenotype_explained_by_alternate: Optional[bool] = None


@dataclass
class AlternateBasisEvidence:
    direction: AlternateBasisDirection
    strength_hint: str = "not_clear"  # fixed; this design does not target strength (same as PS4)
    rationale: str = ""


@dataclass
class BP5Judgment:
    """Data class mirroring the LLM's structured output directly."""
    variant_matching: VariantMatchingResult
    alternate_diagnosis_data: AlternateDiagnosisData
    overall_evidence: AlternateBasisEvidence
    single_study_only: bool = True

    @staticmethod
    def from_json(data: dict) -> "BP5Judgment":
        vm = data["variant_matching"]
        if isinstance(vm, str):
            vm = {"match_status": vm}
        add = data.get("alternate_diagnosis_data") or {}
        return BP5Judgment(
            variant_matching=VariantMatchingResult(
                match_status=MatchStatus(_clean_enum_token(vm["match_status"])),
                match_type=vm.get("match_type"),
                confidence=vm.get("confidence", "low"),
                notes=vm.get("notes", ""),
            ),
            alternate_diagnosis_data=AlternateDiagnosisData(
                alternate_gene=add.get("alternate_gene"),
                alternate_variant=add.get("alternate_variant"),
                phenotype_explained_by_alternate=add.get("phenotype_explained_by_alternate"),
            ),
            # `data.get("overall_evidence") or {}`, not `data["overall_evidence"]`:
            # same null-vs-missing quirk as ps3_bs3.py/ps4.py's own from_json() -
            # the LLM can emit "overall_evidence": null outright.
            overall_evidence=AlternateBasisEvidence(
                direction=AlternateBasisDirection(_clean_enum_token(
                    (data.get("overall_evidence") or {}).get("direction", "not_clear"))),
                strength_hint=(data.get("overall_evidence") or {}).get("strength_hint", "not_clear"),
                rationale=(data.get("overall_evidence") or {}).get("rationale", ""),
            ),
            single_study_only=data.get("single_study_only", True),
        )

    def structured_evidence_items(self) -> list[dict]:
        """Checklist-style items for a curator-facing UI, same pattern as ps4.PS4Judgment."""
        items = []
        add = self.alternate_diagnosis_data
        if add.alternate_gene is not None:
            items.append({
                "label": "Alternate gene named",
                "checked": True,
                "detail": add.alternate_gene,
            })
        if add.alternate_variant is not None:
            items.append({
                "label": "Alternate variant named",
                "checked": True,
                "detail": add.alternate_variant,
            })
        items.append({
            "label": "Alternate finding explains the phenotype (per the paper's own text)",
            "checked": bool(add.phenotype_explained_by_alternate),
            "detail": "" if add.phenotype_explained_by_alternate is None
            else str(add.phenotype_explained_by_alternate),
        })
        return items


# ============================================================================
# 2. Prompt template
# ============================================================================

PROMPT_TEMPLATE = """\
You are a clinical genetics curator supporting ACMG/AMP BP5 (alternate \
molecular basis) evaluation.

Target variant: {gene} {hgvsc} ({hgvsp})
Also accept these equivalent notations as the same variant: {equivalents}

Full text of the paper:
{full_text}

Perform the following steps and output ONLY the JSON described at the end \
(no other explanatory text):

1. Variant identification (variant_matching):
   match_status: matched / heuristic / single_variant_study / unsuccessful
   - matched: the exact target variant (or a notational equivalent from the \
above list) is explicitly reported in this paper, in a proband whose case is \
described
   - unsuccessful: the paper does not report the target variant in a \
described proband at all (if so, stop here; the rest of the fields should \
reflect "no evidence")
   match_type: protein_notation / hgvsc / rsid / abbreviation, etc.
   confidence: high / medium / low
   notes: brief note on how the match was made

2. Alternate diagnosis extraction (alternate_diagnosis_data) - leave any \
field null if the paper does not state it; do NOT guess or infer something \
not explicitly written:
   alternate_gene: the name of a DIFFERENT gene (not {gene}) in which a \
pathogenic or likely pathogenic variant was ALSO identified in the SAME \
proband who carries the target variant, if the paper reports one
   alternate_variant: the HGVS notation or other description of that \
alternate variant, if stated
   phenotype_explained_by_alternate: true if the paper's own text states or \
concludes that this alternate finding explains/accounts for the patient's \
phenotype (e.g. "the patient's phenotype was attributed to the GENE2 \
variant" / "the {gene} variant is considered incidental" / "re-classified \
as the true molecular diagnosis"); false if the paper reports an alternate \
finding but explicitly treats it as NOT explaining the phenotype (e.g. an \
incidental secondary finding of uncertain significance mentioned in \
passing); null if this cannot be determined from the text

3. Self-report (single_study_only): based on citations/mentions within this \
paper alone (not an independent literature search), does this paper appear \
to be the only source reporting this alternate-diagnosis finding for this \
proband, or does it reference/build on other reports?

4. Overall judgment (overall_evidence):
   direction: "BP5" if the paper's own text reports that the proband \
carrying the target variant was found to have a pathogenic/likely \
pathogenic variant in a DIFFERENT gene that the paper's authors consider to \
(better) explain the patient's phenotype - i.e. the target variant's own \
causal role for this patient's phenotype is thereby doubted; "not_clear" \
otherwise (including: the target variant was not matched, no alternate \
finding is reported at all, an alternate finding is reported but not stated \
to explain the phenotype, or you cannot determine this from the text). Do \
NOT infer "BP5" merely from the target variant's own classification being \
uncertain/benign in the paper - BP5 requires a POSITIVELY IDENTIFIED \
alternate causal variant in the same patient, not just doubt about the \
target variant on its own.
   strength_hint: always the fixed string "not_clear" (strength is not \
targeted by this design)
   rationale: 1-3 sentences explaining the basis for the direction above, \
citing what the paper actually states about the alternate finding

Output exactly this JSON structure:
{{
  "variant_matching": {{"match_status": "matched", "match_type": "protein_notation", "confidence": "high", "notes": "..."}},
  "alternate_diagnosis_data": {{"alternate_gene": "GENE2", "alternate_variant": "c.123A>G", "phenotype_explained_by_alternate": true}},
  "single_study_only": true,
  "overall_evidence": {{"direction": "BP5", "strength_hint": "not_clear", "rationale": "..."}}
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

def generate_curator_hints(judgment: BP5Judgment, pmid: Optional[str] = None) -> list[CuratorHint]:
    """Analogous to ps4.generate_curator_hints, but for the BP5 schema."""
    hints: list[CuratorHint] = []
    paper_ref = f" (PMID:{pmid})" if pmid else ""

    if judgment.variant_matching.match_status == MatchStatus.UNSUCCESSFUL:
        return hints

    add = judgment.alternate_diagnosis_data
    if judgment.overall_evidence.direction == AlternateBasisDirection.BP5:
        if add.alternate_gene is None:
            hints.append(CuratorHint(
                "warning",
                f"This paper{paper_ref} was judged BP5-supporting, but no "
                "alternate gene was extracted at all. A BP5 call with no "
                "named alternate finding is hard for a curator to "
                "independently verify from this record alone.",
            ))
        if add.phenotype_explained_by_alternate is not True:
            hints.append(CuratorHint(
                "caution",
                f"This paper{paper_ref} was judged BP5-supporting, but "
                "phenotype_explained_by_alternate was not extracted as "
                "true. BP5 requires the alternate finding to actually "
                "explain the phenotype, not just be present alongside the "
                "target variant; confirm this manually.",
            ))

    if judgment.single_study_only and judgment.overall_evidence.direction != AlternateBasisDirection.NOT_CLEAR:
        hints.append(CuratorHint(
            "caution",
            f"This judgment{paper_ref} is based on a single study only "
            "(self-reported). Confirm via ERepo's full evidenceLinks list, "
            "not this paper alone, the same caution PS4's own judgment "
            "carries.",
        ))

    return hints


def detect_definitive_without_alternate_gene(judgment: BP5Judgment) -> bool:
    """
    Mirrors ps4.detect_definitive_without_numbers: a confident BP5 call with
    no named alternate gene at all is internally inconsistent - there is
    nothing concrete for a curator (or a downstream re-read) to check - and
    should be forced to not_clear, not just flagged.
    """
    if judgment.variant_matching.match_status == MatchStatus.UNSUCCESSFUL:
        return False
    if judgment.overall_evidence.direction != AlternateBasisDirection.BP5:
        return False
    return judgment.alternate_diagnosis_data.alternate_gene is None


def detect_bp5_without_phenotype_explained(judgment: BP5Judgment) -> bool:
    """
    A BP5 call where the paper's own text does NOT say the alternate
    finding explains the phenotype (phenotype_explained_by_alternate is
    False or null) contradicts BP5's actual definition - an alternate
    finding that is merely present, e.g. an incidental secondary finding,
    does not itself doubt the target variant's causal role.
    """
    if judgment.overall_evidence.direction != AlternateBasisDirection.BP5:
        return False
    return judgment.alternate_diagnosis_data.phenotype_explained_by_alternate is not True


def finalize(judgment: BP5Judgment, pmid: Optional[str] = None) -> FinalResult:
    hints = generate_curator_hints(judgment, pmid=pmid)
    effective_direction = judgment.overall_evidence.direction
    paper_ref = f" (PMID:{pmid})" if pmid else ""

    if detect_definitive_without_alternate_gene(judgment):
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, the LLM concluded BP5 but named no "
            "alternate gene at all. A definitive alternate-basis claim "
            "with nothing concrete behind it is internally inconsistent, "
            "so the judgment has been forced to not_clear.",
        ))
        effective_direction = AlternateBasisDirection.NOT_CLEAR
    elif detect_bp5_without_phenotype_explained(judgment):
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, the LLM concluded BP5 but did not "
            "state that the alternate finding explains the patient's "
            "phenotype (phenotype_explained_by_alternate was not true). "
            "BP5 requires the alternate finding to actually account for "
            "the phenotype, not merely be present, so the judgment has "
            "been forced to not_clear.",
        ))
        effective_direction = AlternateBasisDirection.NOT_CLEAR

    return FinalResult(judgment=judgment, curator_hints=hints, effective_direction=effective_direction)


# ============================================================================
# 4. Multi-paper aggregation - thin wrapper, see common.py
# ============================================================================

def aggregate_multi_paper_results(contributions: list[PaperContribution]) -> AggregatedJudgment:
    """Same generic aggregation as ps3_bs3.py/ps4.py's own wrapper - see common.aggregate_multi_paper_results."""
    return _generic_aggregate_multi_paper_results(contributions, not_clear=AlternateBasisDirection.NOT_CLEAR)
