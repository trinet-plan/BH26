from acmg.core.models import Status
from acmg.criteria.common import population_context, result
from acmg.services.population import number


def evaluate(input_data, services, config):
    early, context = population_context("BA1", input_data, services, config)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    high = [item for item in observations if number(item["AF"]) > number("0.05")]
    if not high:
        status = Status.NOT_EVALUATED if rejected or failures else Status.NOT_MET
        return result("BA1", input_data, status, "No reliable AF above 5%",
                      evidence=observations, provenance=provenance)
    exception = input_data.get("ba1_exception_assessment", {})
    if not all(exception.get(key) for key in ("source", "source_version", "reviewed_at")):
        return result("BA1", input_data, Status.NOT_EVALUATED, "BA1 exception check unavailable",
                      evidence=observations, missing=["ba1_exception_assessment"], provenance=provenance)
    if exception.get("is_exception") is True:
        return result("BA1", input_data, Status.NOT_MET, "Known BA1 exception",
                      evidence=observations + [exception], provenance=provenance)
    if exception.get("is_exception") is not False:
        return result("BA1", input_data, Status.MANUAL_REVIEW, "BA1 exception status unresolved",
                      evidence=observations + [exception], review=["Confirm exception status"])
    return result("BA1", input_data, Status.MET, "Reliable AF above 5% and exception check negative",
                  strength="stand_alone", evidence=observations + [exception], provenance=provenance)
