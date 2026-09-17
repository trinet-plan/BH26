from acmg_pipeline.constants import CriterionStatus
"""Gene-disease mechanism and variant spectrum, not constraint alone."""

from acmg_pipeline.criteria.common import annotation_context, curated_context, require_boolean_fields, result


def evaluate_mechanism(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    if "missense_variant" not in annotation["consequences"]:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Requires a missense variant",
                      evidence=[annotation])
    # The mechanism statement is about the gene, so a condition-agnostic assessment is
    # usable and the disease relevance is reported instead of being required up front.
    early, mechanism = curated_context(code, "gene_disease", input_data, services, annotation,
                                       disease_required=False)
    if early:
        return early
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
            "Criterion is not applicable under the reviewed gene-disease specification",
            evidence=evidence,
            provenance={
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
