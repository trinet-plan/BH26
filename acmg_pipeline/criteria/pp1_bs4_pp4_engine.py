"""
pp1_bs4_pp4_engine.py

Wires acmg_pipeline.criteria.pp4_pp1_bs4.evaluate_locus_evidence() (pulled
into main 2026-09-17 from r-kobayashi's pp4_pp1_bs4 branch - see that
module's own docstring) into this project's classify()/VA-Spec pipeline,
per the user's explicit direction (2026-09-17): "pipelineと繋げて。テスト
も流れるようにして" (connect it to the pipeline, make sure tests run too).

[Why this is a separate module from pp4_pp1_bs4.py]
  pp4_pp1_bs4.py is the reviewed judgment logic itself (ClinGen 2024
  Bayesian points, phenotype matching, family segregation scoring) and is
  left exactly as pulled in, unmodified. Everything specific to THIS
  project's own integration - the points-to-Strength mapping, the curated
  reference-record lookup, and the CriterionEvidence/EvidenceLine
  conversion - lives here instead, the same separation segregation.py
  (judgment logic) already keeps from export.py (VA-Spec conversion).

[The points-to-Strength mapping is a provisional design choice, NOT
 reviewed by the team that wrote pp4_pp1_bs4.py/evaluator.py]
  evaluator.py's DIAGNOSTIC_YIELD_POINT_TABLE and pp4_pp1_bs4.py's
  segregation scoring both produce points on the same Tavtigian-compatible
  scale acmg_pipeline.classification already uses (1/2/4/8 = Supporting/
  Moderate/Strong/Very_Strong), but neither module says how a continuous
  point value (e.g. 3.5, or a value above 8) should map onto this
  project's discrete Strength enum for classify(). _strength_for_points()
  below floors to the nearest tier at or below the point value (matching
  evaluator.py's own "for values between rows, use the lower point value"
  convention for its own table) - this is this project's own
  interpretation, flagged here for the pp4_pp1_bs4 branch's own author to
  confirm or correct, not an authoritative ClinGen reading.

[Why PP1/BS4/PP4 need a curated reference record, and what happens
 without one]
  evaluate_locus_evidence() requires a PP4ReferenceRecord (curated gene-
  phenotype diagnostic-yield data) and several explicit conservative gates
  (method_comparable, fully_penetrant, low_phenocopy, ar_case_mode) that
  the function deliberately never infers. This project has no such
  curated database yet - config/pp4_reference_records.json is a real,
  empty DRAFT registry (same "registry_status: DRAFT, only APPROVED
  entries load" convention as config/bs1_thresholds_draft.json), not a
  fabricated one. Until entries are curated and approved there, every
  real gene correctly evaluates to UNKNOWN for PP1/BS4/PP4 - the same
  "honest gap, not a guess" behavior every other stub in this project
  already uses. See test_pp1_bs4_pp4_engine.py for both the "no curated
  reference" path and (using a synthetic in-test reference record, the
  same pattern test_pp4_pp1_bs4.py's own fixtures already use) the "real
  evidence produced" path.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from acmg_pipeline.classification import CriterionEvidence, Strength
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.criteria.pp4_pp1_bs4 import (
    LocusEvidenceResult,
    PP4ReferenceRecord,
    evaluate_locus_evidence,
)
from acmg_pipeline.vcf_record import VariantRecord

PP1_BS4_PP4_CODES = frozenset({"PP1", "BS4", "PP4"})

_DEFAULT_REGISTRY_PATH = Path(__file__).resolve().parents[2] / "config" / "pp4_reference_records.json"


def load_reference_records(path: Path = _DEFAULT_REGISTRY_PATH) -> dict[str, dict]:
    """gene -> {"reference": PP4ReferenceRecord, "gates": {...}} for APPROVED entries only.

    An entry's `gates` dict carries the conservative parameters evaluate_
    locus_evidence() requires and never infers (method_comparable,
    inheritance_mode, fully_penetrant, low_phenocopy, ar_case_mode) - these
    are as much curated, per-gene/per-VCEP-specification facts as the
    diagnostic yield itself, so they are curated in the same registry
    entry rather than guessed at evaluation time.
    """
    if not path.is_file():
        return {}
    registry = json.loads(path.read_text(encoding="utf-8"))
    result: dict[str, dict] = {}
    for entry in registry.get("entries", []):
        if entry.get("status") != "APPROVED":
            continue
        reference = PP4ReferenceRecord(
            reference_id=entry["id"],
            gene=entry["gene"],
            phenotype_label=entry["phenotype_label"],
            phenotype_hpo=tuple(entry.get("phenotype_hpo", [])),
            locus_model=entry["locus_model"],
            diagnostic_yield=entry["diagnostic_yield"],
            testing_method=entry["testing_method"],
            source_citation=entry["source_citation"],
            source_note=entry.get("source_note", ""),
        )
        result[entry["gene"]] = {"reference": reference, "gates": entry.get("gates", {})}
    return result


def _strength_for_points(points: float) -> Optional[Strength]:
    """Floor `points` to the highest Tavtigian tier it reaches, or None below Supporting(1)."""
    if points >= 8:
        return Strength.VERY_STRONG
    if points >= 4:
        return Strength.STRONG
    if points >= 2:
        return Strength.MODERATE
    if points >= 1:
        return Strength.SUPPORTING
    return None


def _unknown_all(reason: str) -> dict[str, CriterionEvidence]:
    return {
        code: CriterionEvidence(code=code, status=CriterionStatus.UNKNOWN, source=reason)
        for code in PP1_BS4_PP4_CODES
    }


def to_criterion_evidence(result: LocusEvidenceResult) -> dict[str, CriterionEvidence]:
    """Map one LocusEvidenceResult onto separate PP1/BS4/PP4 CriterionEvidence.

    PP4 is scored from result.pp4 alone (independent of segregation). PP1
    is scored from result.segregation.pp1_points_used - already zeroed by
    evaluate_locus_evidence() itself when high-yield locus homogeneity
    makes it redundant with PP4, so reading it directly (rather than
    pp1_points_raw) correctly reports PP1 as NOT_MET in that case, not as
    a fabricated MET that double-counts the same evidence. BS4 is scored
    from segregation.bs4_met/bs4_points (always exactly -4.0 in the
    current segregation scoring logic - see pp4_pp1_bs4.py - so no
    points-to-strength floor is needed for it).
    """
    evidence: dict[str, CriterionEvidence] = {}

    if result.pp4 is not None and result.pp4.applicable:
        strength = _strength_for_points(result.pp4.points)
        evidence["PP4"] = (
            CriterionEvidence(code="PP4", status=CriterionStatus.MET, strength=strength,
                              source=f"pp1_bs4_pp4_engine: {result.pp4.reference_id}")
            if strength is not None
            else CriterionEvidence(code="PP4", status=CriterionStatus.NOT_MET,
                                    source=f"pp1_bs4_pp4_engine: {result.pp4.reference_id} (below Supporting threshold)")
        )
    elif result.pp4 is not None:
        evidence["PP4"] = CriterionEvidence(
            code="PP4", status=CriterionStatus.NOT_MET,
            source=f"pp1_bs4_pp4_engine: {result.pp4.reason}",
        )
    else:
        evidence["PP4"] = CriterionEvidence(
            code="PP4", status=CriterionStatus.UNKNOWN,
            source="pp1_bs4_pp4_engine: phenotype not evaluable against curated reference",
        )

    if not result.segregation.pp1_evaluable:
        evidence["PP1"] = CriterionEvidence(
            code="PP1", status=CriterionStatus.UNKNOWN,
            source="pp1_bs4_pp4_engine: segregation not evaluable (no family data or unknown inheritance mode)",
        )
    else:
        pp1_strength = _strength_for_points(result.segregation.pp1_points_used)
        evidence["PP1"] = (
            CriterionEvidence(code="PP1", status=CriterionStatus.MET, strength=pp1_strength,
                              source="pp1_bs4_pp4_engine: family segregation")
            if pp1_strength is not None
            else CriterionEvidence(code="PP1", status=CriterionStatus.NOT_MET,
                                    source="pp1_bs4_pp4_engine: family segregation below Supporting threshold")
        )

    if not result.segregation.bs4_evaluable:
        evidence["BS4"] = CriterionEvidence(
            code="BS4", status=CriterionStatus.UNKNOWN,
            source="pp1_bs4_pp4_engine: segregation not evaluable (no family data or unknown inheritance mode)",
        )
    elif result.segregation.bs4_met:
        evidence["BS4"] = CriterionEvidence(
            code="BS4", status=CriterionStatus.MET, strength=Strength.STRONG,
            source="pp1_bs4_pp4_engine: clear non-segregation",
        )
    else:
        evidence["BS4"] = CriterionEvidence(
            code="BS4", status=CriterionStatus.NOT_MET,
            source="pp1_bs4_pp4_engine: no non-segregation observed",
        )

    return evidence


def build_evidence_line(code: str, evidence: CriterionEvidence, variant: VariantRecord) -> dict:
    """Convert one PP1/BS4/PP4 CriterionEvidence into a real VA-Spec EvidenceLine.

    Mirrors acmg_pipeline.export.build_evidence_line() (the literature
    engine's own builder) closely enough to slot into the same 28-line
    document, but works from a plain CriterionEvidence instead of an
    AggregatedJudgment - there is no per-paper contribution list here, so
    `hasEvidenceItems`/`reportedIn` are omitted rather than fabricated.
    """
    from acmg_pipeline.export import (
        ACMG_2015_METHOD_DOCUMENT, _variant_identity, build_reference_extensions,
        evidence_line_id, validate_integrated_line,
    )
    from acmg_pipeline.classification import PATHOGENIC_CODES
    from ga4gh.core.models import Extension
    from ga4gh.va_spec.base.core import Direction, EvidenceLine, Method

    gene, _hgvsc, safe_hgvsc = _variant_identity(variant)
    status = evidence.status

    if status == CriterionStatus.MET:
        direction = Direction.SUPPORTS if code in PATHOGENIC_CODES else Direction.DISPUTES
        outcome_code = code if evidence.strength and evidence.strength.value == _default_strength(code) else f"{code}_{evidence.strength.value}"
        evidence_outcome = {
            "name": f"ACMG 2015 {code} criterion met",
            "primaryCoding": {"system": "ACMG Guidelines, 2015", "code": outcome_code},
        }
        strength_block = {"primaryCoding": {"system": "ACMG Guidelines, 2015", "code": evidence.strength.value}}
    elif status == CriterionStatus.NOT_MET:
        direction = Direction.NEUTRAL
        evidence_outcome = {
            "name": f"ACMG 2015 {code} criterion not met",
            "primaryCoding": {"system": "ACMG Guidelines, 2015", "code": f"{code}_not_met"},
        }
        strength_block = None
    else:
        direction = Direction.NEUTRAL
        evidence_outcome = None
        strength_block = None

    assessment = {"criterion": code, "status": status.value, "summary": evidence.source or ""}
    extensions = [Extension(name="bh26AssessmentDetails", value=assessment)]
    extensions.extend(build_reference_extensions(code, variant))

    line = EvidenceLine(
        id=evidence_line_id(code, gene, safe_hgvsc),
        directionOfEvidenceProvided=direction,
        description=evidence.source or None,
        specifiedBy=Method(
            methodType=code,
            name=f"ACMG/AMP {code} assessment (family segregation / phenotype specificity)",
            reportedIn=ACMG_2015_METHOD_DOCUMENT,
        ),
        extensions=extensions,
    )
    payload = line.model_dump(mode="json", exclude_none=True)
    if evidence_outcome is not None:
        payload["evidenceOutcome"] = evidence_outcome
    if strength_block is not None:
        payload["strengthOfEvidenceProvided"] = strength_block
    return validate_integrated_line(payload, code)


def _default_strength(code: str) -> str:
    from acmg_pipeline.export import default_strength
    return default_strength(code)


async def evaluate(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    config: dict,
) -> dict[str, CriterionEvidence]:
    """Return {code: CriterionEvidence} for PP1, BS4, PP4 - the entry point
    acmg_pipeline.pipeline.evaluate_variant_evidence_lines() calls.

    `config` follows the same convention as the automated engine's
    `automated_config` (acmg_pipeline.pipeline_interface.load_automated_config()) -
    an optional "pp4_reference_records_path" key overrides the default
    registry location, mainly for tests.

    HPO normalization (acmg_pipeline.hpo_extraction.normalize_hpo(), a
    TogoMCP + LLM round trip) only runs once a curated reference record is
    actually found for this gene - with no APPROVED entry the result is
    UNKNOWN regardless of phenotype, so normalizing first would just spend
    real network/LLM cost on an answer that's already decided. The import
    is local (not top-level) because hpo_extraction requires VLLM_BASE_URL/
    VLLM_API_KEY at import time (same reason clinical_note.py lazily
    imports clinical_extraction.py instead of importing it at module load).
    """
    gene = str(variant.info.get("GENE", ""))
    registry_path = Path(config.get("pp4_reference_records_path", _DEFAULT_REGISTRY_PATH))
    registry = load_reference_records(registry_path)
    entry = registry.get(gene)
    if entry is None:
        return _unknown_all(f"pp1_bs4_pp4_engine: no APPROVED curated PP4 reference record for gene {gene!r}")

    from acmg_pipeline import hpo_extraction
    clinical_note = await hpo_extraction.normalize_hpo(clinical_note)

    gates = entry["gates"]
    result = evaluate_locus_evidence(
        clinical_note,
        entry["reference"],
        method_comparable=bool(gates.get("method_comparable", False)),
        inheritance_mode=gates.get("inheritance_mode"),
        fully_penetrant=gates.get("fully_penetrant"),
        low_phenocopy=gates.get("low_phenocopy"),
        ar_case_mode=gates.get("ar_case_mode"),
    )
    return to_criterion_evidence(result)
