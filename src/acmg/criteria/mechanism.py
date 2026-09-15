"""Gene-disease mechanism and variant spectrum, not constraint alone."""

from acmg.core.models import Status
from acmg.criteria.common import annotation_context, curated_context, require_boolean_fields, result


def evaluate_mechanism(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    if "missense_variant" not in annotation["consequences"]:
        return result(code, input_data, Status.NOT_APPLICABLE, "Requires a missense variant",
                      evidence=[annotation])
    early, mechanism = curated_context(code, "gene_disease", input_data, services, annotation)
    if early:
        return early
    evidence = [annotation, mechanism]
    if not annotation.get("gene") or annotation["gene"] != mechanism.get("gene"):
        return result(code, input_data, Status.MANUAL_REVIEW, "Gene identity mismatch",
                      evidence=evidence, review=["Resolve gene-disease mapping"])
    fields = ["missense_mechanism_established", "spectrum_review_complete"]
    fields += ["low_benign_missense_variation"] if code == "PP2" else ["predominantly_truncating"]
    early = require_boolean_fields(code, input_data, evidence, mechanism, fields)
    if early:
        return early
    if not mechanism["spectrum_review_complete"]:
        return result(code, input_data, Status.MANUAL_REVIEW, "Variant spectrum review incomplete",
                      evidence=evidence, review=["Review pathogenic and benign variant spectrum"])
    if code == "PP2":
        met = mechanism["missense_mechanism_established"] and mechanism["low_benign_missense_variation"]
    else:
        met = mechanism["predominantly_truncating"] and not mechanism["missense_mechanism_established"]
    return result(code, input_data, Status.MET if met else Status.NOT_MET,
                  "Reviewed disease mechanism and variant spectrum evaluated",
                  strength="supporting" if met else None, evidence=evidence)
