from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import annotation_context, curated_context, require_boolean_fields, result as _base_result
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def result(code, input_data, status, summary, **kwargs):
    assert code == "BP7"
    outcome = ("the available evidence satisfies BP7" if status == CriterionStatus.MET else
               "BP7 was evaluated but its requirements were not satisfied" if status == CriterionStatus.NOT_MET else
               "no MET or NOT_MET judgment is made from the available information")
    message = (f"{summary.rstrip('.')}. BP7 evaluates synonymous variants outside splice-critical regions "
               f"when splice prediction and conservation evidence support a benign interpretation. "
               f"Outcome: {status.value.upper()}; {outcome}.")
    return _base_result(code, input_data, status, message, **kwargs)


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    early, annotation = annotation_context("BP7", input_data, services)
    if early:
        return early
    if "synonymous_variant" not in annotation["consequences"]:
        noncoding = any("intron" in c or "non_coding" in c or "UTR" in c for c in annotation["consequences"])
        return result("BP7", input_data, CriterionStatus.UNKNOWN if noncoding else CriterionStatus.UNKNOWN,
                      "BP7 noncoding extension requires review: this noncoding consequence needs "
                      "a dedicated splice/position assessment."
                      if noncoding else
                      "BP7 is not applicable: the annotation does not indicate a synonymous variant; "
                      "BP7 evaluates synonymous changes outside splice-critical regions.",
                      evidence=[annotation], review=["Review noncoding BP7 applicability"] if noncoding else [])
    early, assessment = curated_context("BP7", "synonymous_assessment", input_data, services,
                                       annotation, disease_required=False)
    if early:
        return early
    evidence = [annotation, assessment]
    fields = ("outside_splice_critical_region", "no_predicted_splice_impact", "not_conserved",
              "contradictory_rna_evidence")
    early = require_boolean_fields("BP7", input_data, evidence, assessment, fields)
    if early:
        return early
    if not all(assessment.get(key) for key in ("splice_prediction_evidence", "calibration_source", "conservation_evidence", "position_rule_version")):
        return result("BP7", input_data, CriterionStatus.UNKNOWN, "Prediction/position/conservation provenance missing",
                      evidence=evidence, missing=["BP7_prediction_and_position_policy"])
    if assessment["contradictory_rna_evidence"]:
        return result("BP7", input_data, CriterionStatus.UNKNOWN, "RNA evidence contradicts prediction",
                      evidence=evidence, review=["Resolve RNA contradiction"])
    met = all(assessment[field] for field in fields[:3])
    return result("BP7", input_data, CriterionStatus.MET if met else CriterionStatus.NOT_MET,
                  "Reviewed synonymous position, splice prediction and conservation",
                  strength="supporting" if met else None, evidence=evidence)
