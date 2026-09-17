"""Integrated PP4 + PP1/BS4 locus-evidence evaluator for BH26.

This module implements the ClinGen 2024 idea that phenotype specificity (PP4)
and co-segregation / non-segregation (PP1 / BS4) are coupled locus-level
evidence rather than three independent signals.

Expected upstream flow
----------------------
clinical_note
    -> clinical_extraction.py
    -> ClinicalNoteExtraction
    -> hpo_extraction.py
    -> ClinicalNoteExtraction with proband ClinicalFeature.hpo_id populated
    -> evaluate_locus_evidence(...)

Important boundaries
--------------------
* HPO IDs are used only to decide whether the patient's phenotype constellation
  matches the curated phenotype definition attached to the diagnostic-yield
  reference. The HPO match itself does not create PP4 points.
* Diagnostic yield must come from curated/citable gene-phenotype data. It is
  never inferred from HPO terms or by an LLM.
* Family data come from ClinicalNoteExtraction. Missing family data do not block
  PP4, but they make patient-level PP1/BS4 unevaluable.
* The PP1/BS4 logic here is intentionally conservative. Reduced penetrance,
  phenocopies, large/complex pedigrees, and incomplete genotype data should be
  escalated to manual/formal segregation review.
* PP4 + PP1 positive locus evidence is capped at +5 points per variant.

The current ClinicalFeature model stores the source-grounded ``label`` and a
verified ``hpo_id``. It does not currently store the official HPO label. This
module therefore reports both the source label and HPO ID for auditability.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from enum import Enum
from typing import Optional

from acmg_pipeline.clinical_note import ClinicalNoteExtraction, Relative
from acmg_pipeline.criteria.evaluator import (
    LOCUS_EVIDENCE_CAP,
    PP4EvaluationResult,
    evaluate_pp4,
)

@dataclass(frozen=True)
class PP4ReferenceRecord:
    reference_id: str
    gene: str
    phenotype_label: str
    phenotype_hpo: tuple[str, ...]
    locus_model: str
    diagnostic_yield: float
    testing_method: str
    source_citation: str
    source_note: str = ""

class InheritanceMode(str, Enum):
    AUTOSOMAL_DOMINANT = "autosomal_dominant"
    AUTOSOMAL_RECESSIVE = "autosomal_recessive"
    X_LINKED_RECESSIVE = "x_linked_recessive"
    UNKNOWN = "unknown"


_PARENT_RELATIONSHIPS = {"father", "mother", "parent"}
_MALE_RELATIONSHIPS = {
    "father",
    "son",
    "brother",
    "uncle",
    "grandfather",
    "grandson",
    "nephew",
}


@dataclass
class PatientHpoTerm:
    hpo_id: str
    source_label: str


@dataclass
class PhenotypeMatchResult:
    matched: Optional[bool]
    patient_terms: list[PatientHpoTerm] = field(default_factory=list)
    required_hpo_ids: list[str] = field(default_factory=list)
    matched_hpo_ids: list[str] = field(default_factory=list)
    missing_hpo_ids: list[str] = field(default_factory=list)
    reason: str = ""


@dataclass
class RelativeSegregationObservation:
    relationship: str
    affected_status: Optional[bool]
    variant_status: Optional[bool]
    zygosity: Optional[str]
    role: str
    points: float = 0.0
    note: str = ""


@dataclass
class SegregationAssessment:
    inheritance_mode: InheritanceMode
    pp1_evaluable: bool = False
    pp1_points_raw: float = 0.0
    pp1_points_used: float = 0.0
    bs4_evaluable: bool = False
    bs4_met: bool = False
    bs4_points: float = 0.0  # benign points are represented as negative here
    observations: list[RelativeSegregationObservation] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass
class LocusEvidenceResult:
    gene: str
    reference_id: str
    phenotype_match: PhenotypeMatchResult
    pp4: Optional[PP4EvaluationResult]
    segregation: SegregationAssessment
    positive_locus_points_before_cap: float
    combined_positive_locus_points: float
    bs4_points: float
    net_locus_points: float
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        data = asdict(self)
        data["segregation"]["inheritance_mode"] = self.segregation.inheritance_mode.value
        return data


def _normalize_inheritance_mode(value: Optional[str]) -> InheritanceMode:
    if not value:
        return InheritanceMode.UNKNOWN

    token = value.strip().lower().replace("-", " ").replace("_", " ")
    token = " ".join(token.split())

    mapping = {
        "autosomal dominant": InheritanceMode.AUTOSOMAL_DOMINANT,
        "ad": InheritanceMode.AUTOSOMAL_DOMINANT,
        "autosomal recessive": InheritanceMode.AUTOSOMAL_RECESSIVE,
        "ar": InheritanceMode.AUTOSOMAL_RECESSIVE,
        "x linked recessive": InheritanceMode.X_LINKED_RECESSIVE,
        "x-linked recessive": InheritanceMode.X_LINKED_RECESSIVE,
        "xlr": InheritanceMode.X_LINKED_RECESSIVE,
    }
    return mapping.get(token, InheritanceMode.UNKNOWN)


def patient_hpo_terms(extraction: ClinicalNoteExtraction) -> list[PatientHpoTerm]:
    """Return deduplicated proband HPO IDs with their source labels."""
    seen: set[str] = set()
    terms: list[PatientHpoTerm] = []
    for feature in extraction.proband.phenotype.clinical_features:
        if not feature.hpo_id:
            continue
        hpo_id = feature.hpo_id.strip()
        if not hpo_id or hpo_id in seen:
            continue
        seen.add(hpo_id)
        terms.append(PatientHpoTerm(hpo_id=hpo_id, source_label=feature.label.strip()))
    return terms


def match_phenotype_constellation(
    extraction: ClinicalNoteExtraction,
    reference: PP4ReferenceRecord,
) -> PhenotypeMatchResult:
    """Conservative HPO-set match for a curated diagnostic-yield phenotype.

    For the MVP, every HPO ID listed in ``reference.phenotype_hpo`` is treated
    as required. This deliberately avoids similarity thresholds or LLM-based
    guessing. If the curated reference has no HPO definition, the result is
    ``matched=None`` rather than silently treating it as a match.

    A future version can replace this exact/subset rule with a curated
    phenotype-definition object (required/optional/forbidden HPO terms,
    ontology ancestry, formal clinical criteria, etc.) without changing the
    PP4/PP1 integration below.
    """
    terms = patient_hpo_terms(extraction)
    patient_ids = {term.hpo_id for term in terms}
    required = list(dict.fromkeys(reference.phenotype_hpo))

    if not terms:
        return PhenotypeMatchResult(
            matched=None,
            patient_terms=terms,
            required_hpo_ids=required,
            reason="no_normalized_proband_hpo",
        )

    if not required:
        return PhenotypeMatchResult(
            matched=None,
            patient_terms=terms,
            required_hpo_ids=[],
            reason="reference_has_no_curated_hpo_definition",
        )

    matched_ids = [hpo_id for hpo_id in required if hpo_id in patient_ids]
    missing_ids = [hpo_id for hpo_id in required if hpo_id not in patient_ids]
    return PhenotypeMatchResult(
        matched=not missing_ids,
        patient_terms=terms,
        required_hpo_ids=required,
        matched_hpo_ids=matched_ids,
        missing_hpo_ids=missing_ids,
        reason=(
            "all_required_reference_hpo_present"
            if not missing_ids
            else "reference_phenotype_constellation_not_fully_matched"
        ),
    )


def _is_parent(relative: Relative) -> bool:
    return relative.relationship.strip().lower() in _PARENT_RELATIONSHIPS


def _is_known_male_relationship(relative: Relative) -> bool:
    return relative.relationship.strip().lower() in _MALE_RELATIONSHIPS


def _score_family_segregation(
    extraction: ClinicalNoteExtraction,
    *,
    inheritance_mode: InheritanceMode,
    fully_penetrant: Optional[bool],
    low_phenocopy: Optional[bool],
    ar_case_mode: Optional[str],
) -> SegregationAssessment:
    """Score patient-family PP1/BS4 evidence from ClinicalNoteExtraction.

    PP1 points follow ClinGen 2024 Table 3 where the current schema can support
    the calculation:
      * autosomal dominant: +1.0 per informative co-segregating relative
      * autosomal recessive affected: +2.0 per affected co-segregating relative
      * autosomal recessive unaffected: +0.4 per informative unaffected relative
      * X-linked recessive male: +1.0 per informative male relative

    Conservative restrictions:
      * unaffected relatives are counted for PP1 only when fully_penetrant=True
      * unaffected parents are not counted for PP1 (they may be used for phase)
      * AR scoring requires ar_case_mode='homozygous' or 'compound_heterozygous'
      * BS4 from an unaffected carrier requires fully_penetrant=True
      * BS4 from non-segregation requires low_phenocopy=True
      * AR compound-heterozygous relative non-segregation is not auto-BS4
    """
    result = SegregationAssessment(inheritance_mode=inheritance_mode)

    relatives = extraction.family.relatives
    if not relatives:
        result.warnings.append("no_family_relatives_available")
        return result

    if inheritance_mode == InheritanceMode.UNKNOWN:
        result.warnings.append("inheritance_mode_unknown_pp1_bs4_not_scored")
        for relative in relatives:
            result.observations.append(
                RelativeSegregationObservation(
                    relationship=relative.relationship,
                    affected_status=relative.affected_status,
                    variant_status=relative.variant_status,
                    zygosity=relative.zygosity,
                    role="not_scored",
                    note="inheritance mode unknown",
                )
            )
        return result

    ar_mode = (ar_case_mode or "").strip().lower()
    if inheritance_mode == InheritanceMode.AUTOSOMAL_RECESSIVE and ar_mode not in {
        "homozygous",
        "compound_heterozygous",
    }:
        result.warnings.append(
            "autosomal_recessive_requires_ar_case_mode_homozygous_or_compound_heterozygous"
        )

    result.pp1_evaluable = True
    result.bs4_evaluable = True

    for relative in relatives:
        affected = relative.affected_status
        variant = relative.variant_status
        obs = RelativeSegregationObservation(
            relationship=relative.relationship,
            affected_status=affected,
            variant_status=variant,
            zygosity=relative.zygosity,
            role="uninformative",
        )

        if affected is None or variant is None:
            obs.note = "affected or target-variant status is unknown"
            result.observations.append(obs)
            continue

        # ------------------------------
        # Non-segregation / BS4 signals
        # ------------------------------
        affected_without_variant = affected is True and variant is False
        unaffected_with_variant = affected is False and variant is True

        if affected_without_variant or unaffected_with_variant:
            # ClinGen notes that AR compound-heterozygous non-segregation in
            # relatives may provide little/no benign evidence for one variant.
            if (
                inheritance_mode == InheritanceMode.AUTOSOMAL_RECESSIVE
                and ar_mode == "compound_heterozygous"
            ):
                obs.role = "possible_nonsegregation_not_auto_bs4"
                obs.note = "AR compound-heterozygous non-segregation is not auto-scored as BS4"
                result.observations.append(obs)
                continue

            if unaffected_with_variant and fully_penetrant is not True:
                obs.role = "possible_nonsegregation_not_auto_bs4"
                obs.note = "unaffected carrier not counted as BS4 without confirmed full penetrance"
                result.observations.append(obs)
                continue

            if low_phenocopy is not True:
                obs.role = "possible_nonsegregation_not_auto_bs4"
                obs.note = "non-segregation not counted as BS4 without low-phenocopy assumption"
                result.observations.append(obs)
                continue

            if (
                inheritance_mode == InheritanceMode.X_LINKED_RECESSIVE
                and not _is_known_male_relationship(relative)
            ):
                obs.role = "possible_nonsegregation_not_auto_bs4"
                obs.note = "X-linked recessive scoring needs an informative male relationship"
                result.observations.append(obs)
                continue

            result.bs4_met = True
            result.bs4_points = -4.0
            obs.role = "BS4"
            obs.points = -4.0
            obs.note = "clear non-segregation under supplied penetrance/phenocopy assumptions"
            result.observations.append(obs)
            continue

        # ------------------------------
        # Co-segregation / PP1 signals
        # ------------------------------
        if affected is True and variant is True:
            if inheritance_mode == InheritanceMode.AUTOSOMAL_DOMINANT:
                obs.role = "PP1"
                obs.points = 1.0
            elif inheritance_mode == InheritanceMode.AUTOSOMAL_RECESSIVE:
                if ar_mode in {"homozygous", "compound_heterozygous"}:
                    obs.role = "PP1"
                    obs.points = 2.0
                else:
                    obs.note = "AR case mode missing; affected co-segregation not scored"
            elif inheritance_mode == InheritanceMode.X_LINKED_RECESSIVE:
                if _is_known_male_relationship(relative):
                    obs.role = "PP1"
                    obs.points = 1.0
                else:
                    obs.note = "X-linked recessive scoring needs an informative male relationship"

        elif affected is False and variant is False:
            if fully_penetrant is not True:
                obs.note = "unaffected non-carrier not counted for PP1 without confirmed full penetrance"
            elif _is_parent(relative):
                obs.note = "unaffected parent not counted for PP1; may be used to establish phase"
            elif inheritance_mode == InheritanceMode.AUTOSOMAL_DOMINANT:
                obs.role = "PP1"
                obs.points = 1.0
            elif inheritance_mode == InheritanceMode.AUTOSOMAL_RECESSIVE:
                if ar_mode in {"homozygous", "compound_heterozygous"}:
                    obs.role = "PP1"
                    obs.points = 0.4
                else:
                    obs.note = "AR case mode missing; unaffected co-segregation not scored"
            elif inheritance_mode == InheritanceMode.X_LINKED_RECESSIVE:
                if _is_known_male_relationship(relative):
                    obs.role = "PP1"
                    obs.points = 1.0
                else:
                    obs.note = "X-linked recessive scoring needs an informative male relationship"

        result.pp1_points_raw += max(obs.points, 0.0)
        result.observations.append(obs)

    return result


def evaluate_locus_evidence(
    extraction: ClinicalNoteExtraction,
    reference: PP4ReferenceRecord,
    *,
    method_comparable: bool,
    inheritance_mode: Optional[str] = None,
    fully_penetrant: Optional[bool] = None,
    low_phenocopy: Optional[bool] = None,
    ar_case_mode: Optional[str] = None,
    candidate_variants_on_allele: int = 1,
    adjusted_diagnostic_yield: Optional[float] = None,
    phenotype_match_override: Optional[PhenotypeMatchResult] = None,
) -> LocusEvidenceResult:
    """Evaluate PP4 and patient-family PP1/BS4 in one coordinated pass.

    Parameters
    ----------
    extraction:
        ClinicalNoteExtraction after HPO normalization.
    reference:
        Curated gene-phenotype/testing-method diagnostic-yield record.
    method_comparable:
        True only when the current test methodology is sufficiently comparable
        to the methodology underlying the reference diagnostic yield.
    inheritance_mode:
        Explicit override such as ``autosomal_dominant``. If omitted, the
        function uses ``extraction.family.inheritance_pattern``. No inheritance
        mode is inferred from the pedigree.
    fully_penetrant / low_phenocopy:
        Conservative gates for using unaffected relatives / BS4 automatically.
        Leave as None when not established.
    ar_case_mode:
        ``homozygous`` or ``compound_heterozygous`` for AR cases.
    adjusted_diagnostic_yield:
        Optional externally curated/recalculated yield after other loci have
        been excluded. The function does not invent this value.
    phenotype_match_override:
        When supplied, used instead of calling ``match_phenotype_constellation``
        against ``reference.phenotype_hpo``. This is the seam
        ``acmg_pipeline.pubcasefinder`` hooks into: PubCaseFinder's HPO-based
        gene ranking (Layer 2 in the BH26 ACMG criteria definition, v2 2026-09-14)
        answers "is this gene the best phenotype match" directly from the
        patient's own HPO profile, so callers with that signal available no
        longer need a hand-curated ``phenotype_hpo`` required-term list per
        reference record.
    """
    if adjusted_diagnostic_yield is not None and not 0.0 <= adjusted_diagnostic_yield <= 1.0:
        raise ValueError("adjusted_diagnostic_yield must be between 0 and 1")

    phenotype = (
        phenotype_match_override
        if phenotype_match_override is not None
        else match_phenotype_constellation(extraction, reference)
    )

    pp4_reference = (
        replace(reference, diagnostic_yield=adjusted_diagnostic_yield)
        if adjusted_diagnostic_yield is not None
        else reference
    )

    pp4_result: Optional[PP4EvaluationResult]
    warnings: list[str] = []

    if phenotype.matched is None:
        pp4_result = None
        warnings.append(f"pp4_not_evaluable:{phenotype.reason}")
    else:
        pp4_result = evaluate_pp4(
            pp4_reference,
            phenotype_match=phenotype.matched,
            method_comparable=method_comparable,
            candidate_variants_on_allele=candidate_variants_on_allele,
        )

    mode = _normalize_inheritance_mode(
        inheritance_mode or extraction.family.inheritance_pattern
    )
    segregation = _score_family_segregation(
        extraction,
        inheritance_mode=mode,
        fully_penetrant=fully_penetrant,
        low_phenocopy=low_phenocopy,
        ar_case_mode=ar_case_mode,
    )

    pp4_points = pp4_result.points if pp4_result and pp4_result.applicable else 0.0
    raw_pp1 = segregation.pp1_points_raw

    # High-yield locus-homogeneous phenotype: PP1 pathogenic support overlaps
    # with PP4 and should not be added again.
    if pp4_result is not None and pp4_result.pp1_support_allowed is False:
        segregation.pp1_points_used = 0.0
        if raw_pp1 > 0:
            warnings.append("pp1_not_added_high_yield_locus_homogeneity")
    else:
        segregation.pp1_points_used = raw_pp1

    positive_before_cap = pp4_points + segregation.pp1_points_used
    combined_positive = min(positive_before_cap, LOCUS_EVIDENCE_CAP)

    # A clear BS4 observation at the target locus is contradictory to using
    # positive PP4/PP1 locus evidence automatically. Preserve the raw values
    # for audit, but suppress automatic positive locus points pending review.
    if segregation.bs4_met:
        if combined_positive > 0:
            warnings.append("bs4_conflicts_with_positive_locus_evidence_manual_review_required")
        combined_for_net = 0.0
    else:
        combined_for_net = combined_positive

    net = combined_for_net + segregation.bs4_points

    return LocusEvidenceResult(
        gene=reference.gene,
        reference_id=reference.reference_id,
        phenotype_match=phenotype,
        pp4=pp4_result,
        segregation=segregation,
        positive_locus_points_before_cap=positive_before_cap,
        combined_positive_locus_points=combined_for_net,
        bs4_points=segregation.bs4_points,
        net_locus_points=net,
        warnings=warnings + segregation.warnings,
    )
