from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import citable, population_context, result as _base_result
from acmg_pipeline.services.population import number
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def result(code, input_data, status, summary, **kwargs):
    assert code == "BA1"
    outcome = "satisfies BA1" if status == CriterionStatus.MET else (
        "was evaluated but does not satisfy BA1" if status == CriterionStatus.NOT_MET
        else "cannot be assigned MET or NOT_MET from the available information")
    message = (f"{summary.rstrip('.')}. BA1 evaluates whether population allele frequency exceeds "
               f"the stand-alone benign threshold after exception review. Outcome: {status.value.upper()}; {outcome}.")
    return _base_result(code, input_data, status, message, **kwargs)


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    early, context = population_context("BA1", input_data, services, config)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    high = [item for item in observations if number(item["AF"]) > number("0.05")]
    if not high:
        status = CriterionStatus.UNKNOWN if failures else CriterionStatus.NOT_MET
        return result("BA1", input_data, status, "No reliable AF above 5%",
                      evidence=observations, provenance=provenance)
    exception = input_data.get("ba1_exception_assessment", {})
    if not all(exception.get(key) for key in ("source", "source_version", "reviewed_at")):
        return result("BA1", input_data, CriterionStatus.UNKNOWN, "BA1 exception check unavailable",
                      evidence=observations, missing=["ba1_exception_assessment"], provenance=provenance)
    # The exception list is curated policy: it is always recorded, and cited as an evidence
    # item only when it carries a retrievable identifier.
    provenance = {**provenance, "ba1_exception_assessment": exception}
    cited = citable(exception)
    if exception.get("is_exception") is True:
        return result("BA1", input_data, CriterionStatus.NOT_MET, "Known BA1 exception",
                      evidence=observations + cited, provenance=provenance)
    if exception.get("is_exception") is not False:
        return result("BA1", input_data, CriterionStatus.UNKNOWN, "BA1 exception status unresolved",
                      evidence=observations + cited, review=["Confirm exception status"],
                      provenance=provenance)
    return result("BA1", input_data, CriterionStatus.MET, "Reliable AF above 5% and exception check negative",
                  strength="stand_alone", evidence=observations + cited, provenance=provenance)
