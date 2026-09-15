"""Protein length / repeat / critical-region assessments with explicit missingness."""

from acmg.core.models import Status
from acmg.criteria.common import (
    annotation_context, curated_context, require_boolean_fields, result,
)


def evaluate_region(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    consequences = set(annotation["consequences"])
    indel = bool(consequences & {"inframe_insertion", "inframe_deletion"})
    if code == "PM4" and not (indel or "stop_lost" in consequences):
        return result(code, input_data, Status.NOT_APPLICABLE, "Requires in-frame indel or stop loss",
                      evidence=[annotation])
    if code == "BP3" and not indel:
        return result(code, input_data, Status.NOT_APPLICABLE, "Requires in-frame indel",
                      evidence=[annotation])
    # PM1 is a protein-level statement about the region itself, so a condition-agnostic
    # reviewed assessment is usable; the disease relevance is reported separately.
    early, region = curated_context(code, "region", input_data, services, annotation,
                                    disease_required=False)
    if early:
        return early
    evidence = [annotation, region]
    protein = annotation.get("protein_id")
    position = annotation.get("protein_start")
    end = annotation.get("protein_end", position)
    if not protein or region.get("protein_id") != protein:
        return result(code, input_data, Status.NOT_EVALUATED, "Protein reference mapping unavailable",
                      evidence=evidence, missing=["protein_id"])
    start, stop = region.get("start"), region.get("end")
    if not all(type(value) is int and value > 0 for value in (position, end, start, stop)) or start > stop or position > end:
        return result(code, input_data, Status.NOT_EVALUATED, "Protein coordinates unavailable",
                      evidence=evidence, missing=["protein_interval"])
    # A partial intersection is insufficient to classify the full altered region as non-functional.
    contained = start <= position <= end <= stop
    if not contained:
        return result(code, input_data, Status.MANUAL_REVIEW, "Assessment does not cover altered protein interval",
                      evidence=evidence, review=["Review complete affected interval"])
    fields = {
        "PM1": ["critical_region", "pathogenic_enrichment", "benign_depletion"],
        "PM4": ["nonfunctional_repeat", "functional_review_complete"],
        "BP3": ["repetitive", "functional_importance", "functional_review_complete"],
    }[code]
    early = require_boolean_fields(code, input_data, evidence, region, fields)
    if early:
        return early
    if code == "PM1":
        met = all(region[field] for field in fields)
    else:
        if not region["functional_review_complete"]:
            return result(code, input_data, Status.MANUAL_REVIEW, "Functional region review incomplete",
                          evidence=evidence, review=["Confirm functional relevance of repeat/region"])
        if code == "PM4":
            delta = annotation.get("protein_length_change")
            if type(delta) is not int:
                return result(code, input_data, Status.NOT_EVALUATED, "Protein length change unavailable",
                              evidence=evidence, missing=["protein_length_change"])
            met = delta != 0 and not region["nonfunctional_repeat"]
        else:
            met = region["repetitive"] and not region["functional_importance"]
    strength = "supporting" if code == "BP3" else "moderate"
    extra = {}
    review = []
    if code == "PM1":
        condition = input_data.get("condition")
        matched = bool(condition) and region.get("condition") == condition
        extra = {"assessment_scope": "protein_level",
                 "condition_assessment": "MATCHED" if matched else "NOT_EVALUATED"}
        if met and not matched:
            review = ["Confirm region criticality for the disease context before final classification"]
    return result(code, input_data, Status.MET if met else Status.NOT_MET,
                  "Altered protein interval and reviewed region evidence evaluated",
                  strength=strength if met else None, evidence=evidence,
                  review=review, provenance=extra)
