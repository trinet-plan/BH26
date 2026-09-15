from acmg.core.models import Status
from acmg.criteria.common import population_context, result
from acmg.services.population import number


def evaluate(input_data, services, config):
    assessment = input_data.get("disease_frequency_threshold", {})
    threshold = number(assessment.get("max_credible_af"))
    condition = input_data.get("condition")
    if not condition or assessment.get("condition") != condition or threshold is None or not 0 < threshold < 1:
        return result("BS1", input_data, Status.NOT_EVALUATED, "Disease-specific frequency threshold missing",
                      missing=["condition", "disease_frequency_threshold"])
    if not all(assessment.get(key) for key in ("source", "source_version", "reviewed_at", "inheritance")):
        return result("BS1", input_data, Status.NOT_EVALUATED, "Disease threshold provenance missing",
                      missing=["disease_frequency_threshold.provenance"])
    if input_data.get("inheritance") != assessment["inheritance"]:
        return result("BS1", input_data, Status.MANUAL_REVIEW, "Threshold inheritance context mismatch",
                      review=["Confirm disease inheritance"])
    early, context = population_context("BS1", input_data, services, config)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    above = any(number(item["AF"]) > threshold for item in observations)
    status = Status.MET if above else Status.NOT_EVALUATED if rejected or failures else Status.NOT_MET
    return result("BS1", input_data, status, "Observed AF compared with curated disease threshold",
                  strength="strong" if above else None, evidence=observations + [assessment],
                  provenance=provenance)
