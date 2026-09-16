from acmg.core.models import Status
from acmg.core.interface import criterion_input
from acmg.criteria.common import citable, population_context, result
from acmg.services.population import number
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    early, context = population_context("BA1", input_data, services, config)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    high = [item for item in observations if number(item["AF"]) > number("0.05")]
    if not high:
        status = Status.NOT_EVALUATED if failures else Status.NOT_MET
        return result("BA1", input_data, status, "No reliable AF above 5%",
                      evidence=observations, provenance=provenance)
    exception = input_data.get("ba1_exception_assessment", {})
    if not all(exception.get(key) for key in ("source", "source_version", "reviewed_at")):
        return result("BA1", input_data, Status.NOT_EVALUATED, "BA1 exception check unavailable",
                      evidence=observations, missing=["ba1_exception_assessment"], provenance=provenance)
    # The exception list is curated policy: it is always recorded, and cited as an evidence
    # item only when it carries a retrievable identifier.
    provenance = {**provenance, "ba1_exception_assessment": exception}
    cited = citable(exception)
    if exception.get("is_exception") is True:
        return result("BA1", input_data, Status.NOT_MET, "Known BA1 exception",
                      evidence=observations + cited, provenance=provenance)
    if exception.get("is_exception") is not False:
        return result("BA1", input_data, Status.MANUAL_REVIEW, "BA1 exception status unresolved",
                      evidence=observations + cited, review=["Confirm exception status"],
                      provenance=provenance)
    return result("BA1", input_data, Status.MET, "Reliable AF above 5% and exception check negative",
                  strength="stand_alone", evidence=observations + cited, provenance=provenance)
