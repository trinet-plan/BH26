from acmg_pipeline.constants import CriterionStatus
"""Gene-disease mechanism and variant spectrum, not constraint alone."""

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import (
    NOT_APPLICABLE, annotation_context, curated_context, require_boolean_fields, result,
)


def _gene_disease_draft_record(input_data, services, gene):
    """The gene_disease_draft record for this gene, if resolve.py's optional
    ClinGen Gene-Disease Validity + gnomAD constraint step produced one.

    Bypasses the standard PASS-only evidence filter
    (EvidenceService.get_candidates()) deliberately:
    acmg_pipeline.providers.gene_disease_draft.GeneDiseaseDraftProvider marks
    its own output quality_status="DRAFT" specifically so the normal
    evidence path never picks it up by accident. Reading it here is this
    module's own explicit, flagged fallback (see _mechanism_from_draft()),
    not a bypass of that design.
    """
    variant = Variant(**input_data["variant"])
    records = getattr(services.evidence, "records", [])
    return next(
        (item for item in records
         if item.get("category") == "gene_disease_draft"
         and item.get("variant_key") == variant.key and item.get("gene") == gene),
        None,
    )


def _mechanism_from_draft(code, input_data, services, annotation):
    """A flagged MET/NOT_MET from the statistical gene-disease suggestion.

    Only for CANDIDATE/NOT_SUGGESTED - the suggestion is direct enough there
    to report a prediction, always with the caveat that it is a statistical
    suggestion (ClinGen Gene-Disease Validity + gnomAD constraint), not a
    curator's own reviewed gene_disease record. INSUFFICIENT/
    CURATED_NEGATIVE returns None so the caller falls through to the
    original UNKNOWN - no call is fabricated from an inconclusive or
    disputed signal.
    """
    draft = _gene_disease_draft_record(input_data, services, annotation.get("gene"))
    if draft is None:
        return None
    suggestion = (draft.get("suggestions") or {}).get(code)
    if not suggestion or suggestion.get("status") not in ("CANDIDATE", "NOT_SUGGESTED"):
        return None
    reasons = "; ".join(suggestion.get("reasons") or []) or "no reasons recorded"
    evidence = [annotation, draft]
    review = ([f"Statistical suggestion, not a curator-reviewed gene_disease record: {reasons}"]
              + list(suggestion.get("requires_review_of") or []))
    provenance = {"assessment_method": "gene_disease_draft", "suggestion": suggestion,
                  "validity_state": draft.get("validity_state")}
    status = (CriterionStatus.MET if suggestion["status"] == "CANDIDATE"
              else CriterionStatus.NOT_MET)
    return result(code, input_data, status,
                  f"Statistical gene-disease mechanism suggestion (not a curator review): {reasons}",
                  strength="supporting" if status == CriterionStatus.MET else None,
                  evidence=evidence, review=review, provenance=provenance)


def evaluate_mechanism(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    if "missense_variant" not in annotation["consequences"]:
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      f"{code} is not applicable: it evaluates missense variants, and the "
                      f"annotation reports {', '.join(sorted(annotation['consequences']))}.",
                      evidence=[annotation], provenance=NOT_APPLICABLE)
    # The mechanism statement is about the gene, so a condition-agnostic assessment is
    # usable and the disease relevance is reported instead of being required up front.
    early, mechanism = curated_context(code, "gene_disease", input_data, services, annotation,
                                       disease_required=False)
    if early:
        draft_result = _mechanism_from_draft(code, input_data, services, annotation)
        return draft_result if draft_result is not None else early
    evidence = [annotation, mechanism]
    if not annotation.get("gene") or annotation["gene"] != mechanism.get("gene"):
        return result(code, input_data, CriterionStatus.UNKNOWN, "Gene identity mismatch",
                      evidence=evidence, review=["Resolve gene-disease mapping"])
    applicability_field = f"{code.lower()}_applicable"
    applicable = mechanism.get(applicability_field)
    if applicable is False:
        return result(
            code,
            input_data,
            CriterionStatus.UNKNOWN,
            f"{code} is not applicable: the reviewed gene-disease specification for "
            f"{mechanism.get('gene')} marks it as not applicable",
            evidence=evidence,
            provenance={
                **NOT_APPLICABLE,
                "assessment_scope": "condition_specific" if mechanism.get("condition") else "gene_level",
                "condition_assessment": (
                    "MATCHED" if mechanism.get("condition") == input_data.get("condition")
                    else "NOT_EVALUATED"
                ),
                "applicability_source": mechanism.get("applicability_source"),
            },
        )
    if applicable is not None and applicable is not True:
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      "Criterion applicability assessment is invalid", evidence=evidence,
                      missing=[applicability_field])
    fields = ["missense_mechanism_established", "spectrum_review_complete"]
    fields += ["low_benign_missense_variation"] if code == "PP2" else ["predominantly_truncating"]
    early = require_boolean_fields(code, input_data, evidence, mechanism, fields)
    if early:
        return early
    if not mechanism["spectrum_review_complete"]:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Variant spectrum review incomplete",
                      evidence=evidence, review=["Review pathogenic and benign variant spectrum"])
    if code == "PP2":
        met = mechanism["missense_mechanism_established"] and mechanism["low_benign_missense_variation"]
    else:
        met = mechanism["predominantly_truncating"] and not mechanism["missense_mechanism_established"]
    condition = input_data.get("condition")
    matched = bool(condition) and mechanism.get("condition") == condition
    review = [] if matched or not met else [
        "Confirm the gene-disease mechanism for the disease context before final classification"]
    return result(code, input_data, CriterionStatus.MET if met else CriterionStatus.NOT_MET,
                  "Reviewed disease mechanism and variant spectrum evaluated",
                  strength="supporting" if met else None, evidence=evidence, review=review,
                  provenance={"assessment_scope": "gene_level",
                              "condition_assessment": "MATCHED" if matched else "NOT_EVALUATED"})
