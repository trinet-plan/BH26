"""Provisional PVS1 decision path; never exports an applied pathogenic strength."""

from acmg.core.models import Status
from acmg.criteria.common import annotation_context, curated_context, result


ELIGIBLE = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost"}


def evaluate(input_data, services, config):
    early, annotation = annotation_context("PVS1", input_data, services)
    if early:
        return early
    consequences = set(annotation["consequences"])
    if not consequences & ELIGIBLE:
        return result("PVS1", input_data, Status.NOT_APPLICABLE, "Not an eligible loss-of-function consequence",
                      evidence=[annotation])
    early, mechanism = curated_context("PVS1", "gene_disease", input_data, services, annotation)
    if early:
        return early
    evidence = [annotation, mechanism]
    if not annotation.get("gene") or mechanism.get("gene") != annotation["gene"]:
        return result("PVS1", input_data, Status.MANUAL_REVIEW, "Gene-disease mapping unresolved",
                      evidence=evidence, review=["Confirm gene-disease mechanism"])
    if mechanism.get("lof_mechanism_established") is False:
        return result("PVS1", input_data, Status.NOT_MET, "Reviewed disease mechanism does not support LoF",
                      evidence=evidence)
    if mechanism.get("lof_mechanism_established") is not True:
        return result("PVS1", input_data, Status.NOT_EVALUATED, "LoF mechanism unknown",
                      evidence=evidence, missing=["lof_mechanism_established"])
    early, assessment = curated_context("PVS1", "lof_assessment", input_data, services, annotation)
    if early:
        early.evidence = evidence + early.evidence[1:]
        return early
    evidence.append(assessment)
    trace = ["Eligible LoF consequence", "Gene-disease LoF mechanism established"]
    recommendation = None
    points = []
    if assessment.get("biologically_relevant_transcript") is not True or assessment.get("exon_relevant") is not True:
        points.append("Confirm biologically relevant transcript and exon")
    elif "start_lost" in consequences:
        trace.append("Initiation loss: alternative initiation requires review")
        points.append("Assess alternative initiation and pathogenic variants upstream of next start")
    elif consequences & {"splice_donor_variant", "splice_acceptor_variant"} and assessment.get("splice_null_effect_confirmed") is not True:
        trace.append("Splice site: null outcome not established")
        points.append("Confirm exon skipping/cryptic splice effect and reading-frame outcome")
    else:
        nmd = assessment.get("nmd_predicted")
        if type(nmd) is not bool or not assessment.get("nmd_rule_source"):
            points.append("Determine NMD using documented transcript-specific rule")
        elif nmd:
            trace.append("NMD predicted in relevant transcript/exon")
            recommendation = "very_strong"
        else:
            trace.append("Predicted NMD escape")
            if assessment.get("critical_region_disrupted") is True:
                recommendation = "strong"
                trace.append("Reviewed critical region disrupted")
            else:
                points.append("Assess extent and importance of lost protein region")
    # Remaining decision-tree branches are not approximated; a reviewed branch may suggest a strength.
    if recommendation is None and assessment.get("decision_tree_branch") and assessment.get("decision_tree_source"):
        reviewed_strength = assessment.get("reviewed_recommended_strength")
        if reviewed_strength in {"very_strong", "strong", "moderate", "supporting"}:
            recommendation = reviewed_strength
            trace.append(f"Curated branch: {assessment['decision_tree_branch']}")
    return result("PVS1", input_data, Status.MANUAL_REVIEW, "Provisional PVS1 candidate; curator confirmation required",
                  evidence=evidence, review=[*points, "Confirm PVS1 decision path and proposed strength"],
                  provenance={"provisional": True, "decision_path": trace, "recommended_strength": recommendation})
