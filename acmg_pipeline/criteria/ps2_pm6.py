"""
ps2_pm6.py

Real, rule-based evaluator for PS2 ("de novo, both parentage confirmed, no
family history") and PM6 ("assumed de novo, parentage not confirmed") - the
"PS2/PM6's de-novo rule" call site pipeline.py's own comments (search_
candidate_pmids()'s docstring, evaluate_variant_evidence_lines()'s call
site) already anticipated as "not-yet-wired", per acmg_pipeline.criteria.
stubs.py's own instructions for replacing a stub: "build a module shaped
like ps4.py or segregation.py, then register it in registry.py in place of
the stub."

Unlike ps3_bs3/ps4/segregation, this needs no literature search or LLM
call at all: acmg_pipeline.clinical_note.ClinicalNoteExtraction.de_novo
(father/mother variant_status, paternity/maternity confirmed) is already
extracted from the clinical note by clinical_extraction.py - see that
dataclass's own docstring, "exactly the data PS2 ... and PM6 ... need".
This module only has to apply the ACMG/AMP 2015 rule to fields that are
already there, the same no-literature-needed shape segregation.
from_clinical_note_family() already uses for the patient's own family data.

PS2 and PM6 are mutually exclusive here (same de novo signal, differing
only in whether parentage was confirmed) - the same one-pair-at-a-time
shape ps3_bs3.py already uses for PS3/BS3.
"""

from __future__ import annotations

from acmg_pipeline.classification import CriterionEvidence, Strength
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.vcf_record import VariantRecord


def evaluate(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    config: dict,
) -> dict[str, CriterionEvidence]:
    """Return {"PS2": ..., "PM6": ...} from clinical_note.de_novo alone.

    `variant`/`config` are accepted for call-site parity with
    pp1_bs4_pp4_engine.evaluate() (acmg_pipeline.pipeline.evaluate_variant_
    evidence_lines() dispatches every phenotype/family-derived engine the
    same way) but are not read - this rule needs only the clinical note.
    """
    dn = clinical_note.de_novo
    affected = clinical_note.proband.phenotype.affected_status
    family_history = clinical_note.family.family_history_status

    if dn.father_variant_status is None or dn.mother_variant_status is None:
        reason = (
            "father/mother variant_status not reported in the clinical note - "
            "de novo occurrence cannot be determined"
        )
        return {
            "PS2": CriterionEvidence(code="PS2", status=CriterionStatus.UNKNOWN, source=reason),
            "PM6": CriterionEvidence(code="PM6", status=CriterionStatus.UNKNOWN, source=reason),
        }

    if dn.father_variant_status or dn.mother_variant_status:
        reason = (
            "variant is present in at least one biological parent - "
            "not consistent with a de novo occurrence"
        )
        return {
            "PS2": CriterionEvidence(code="PS2", status=CriterionStatus.NOT_MET, source=reason),
            "PM6": CriterionEvidence(code="PM6", status=CriterionStatus.NOT_MET, source=reason),
        }

    # Both biological parents confirmed variant-negative from here on - a
    # genuine de novo occurrence, pending the patient/family checks below.
    if affected is not True:
        reason = "proband affected_status is not confirmed - PS2/PM6 require a patient with the disease"
        return {
            "PS2": CriterionEvidence(code="PS2", status=CriterionStatus.UNKNOWN, source=reason),
            "PM6": CriterionEvidence(code="PM6", status=CriterionStatus.UNKNOWN, source=reason),
        }

    if family_history:
        # PS2's own ACMG/AMP 2015 definition requires "no family history" -
        # a reported positive family history directly contradicts an
        # isolated de novo occurrence for this patient (regardless of
        # parentage confirmation), so neither code applies.
        reason = "family_history_status is positive, which contradicts an isolated de novo occurrence"
        return {
            "PS2": CriterionEvidence(code="PS2", status=CriterionStatus.NOT_MET, source=reason),
            "PM6": CriterionEvidence(code="PM6", status=CriterionStatus.NOT_MET, source=reason),
        }

    if dn.paternity_confirmed and dn.maternity_confirmed:
        return {
            "PS2": CriterionEvidence(
                code="PS2", status=CriterionStatus.MET, strength=Strength.STRONG,
                source="clinical_note.de_novo: variant absent in both biological parents, "
                       "paternity and maternity confirmed, patient affected, no family history",
            ),
            "PM6": CriterionEvidence(
                code="PM6", status=CriterionStatus.NOT_MET,
                source="paternity/maternity both confirmed - evaluated as PS2 instead",
            ),
        }

    return {
        "PS2": CriterionEvidence(
            code="PS2", status=CriterionStatus.NOT_MET,
            source="paternity/maternity not both confirmed - evaluated as PM6 instead",
        ),
        "PM6": CriterionEvidence(
            code="PM6", status=CriterionStatus.MET, strength=Strength.MODERATE,
            source="clinical_note.de_novo: variant absent in both biological parents "
                   "(assumed de novo), paternity/maternity not both confirmed",
        ),
    }


def build_evidence_line(code: str, evidence: CriterionEvidence, variant: VariantRecord) -> dict:
    """Convert one PS2/PM6 CriterionEvidence into a real VA-Spec EvidenceLine.

    Mirrors pp1_bs4_pp4_engine.build_evidence_line() closely (same 28-line
    document, same UNKNOWN-reported-as-not_met-with-a-curatorHint
    convention - see that function's own docstring for why), differing
    only in the method name/description, which describes de novo status
    rather than family segregation / phenotype specificity.
    """
    from acmg_pipeline.export import (
        ACMG_2015_METHOD_DOCUMENT, _variant_identity, build_reference_extensions,
        evidence_line_id, validate_integrated_line,
    )
    from acmg_pipeline.automated_va_spec import details_as_extensions
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

    curator_hints = []
    reported_status = status
    if status == CriterionStatus.UNKNOWN:
        reported_status = CriterionStatus.NOT_MET
        curator_hints.append({
            "severity": "caution", "category": "unevaluated",
            "message": f"{code} could not actually be evaluated from the "
                       "available data (reported as not_met rather than left "
                       "unknown).",
        })
    extensions = [Extension(**d) for d in details_as_extensions({"status": reported_status.value})]
    if curator_hints:
        extensions.append(Extension(name="curatorHints", value=curator_hints))
    ref_extensions, reference_evidence_items = build_reference_extensions(code, variant)
    extensions.extend(ref_extensions)

    line = EvidenceLine(
        id=evidence_line_id(code, gene, safe_hgvsc),
        directionOfEvidenceProvided=direction,
        description=evidence.source or None,
        hasEvidenceItems=reference_evidence_items or None,
        specifiedBy=Method(
            methodType=code,
            name=f"ACMG/AMP {code} assessment (de novo occurrence)",
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
