from acmg.core.models import Status
from acmg.criteria.common import population_context, result
from acmg.services.population import number


def evaluate(input_data, services, config):
    early, context = population_context("PM2", input_data, services, config)
    if early:
        return early
    rule, observations, rejected, failures, provenance = context
    threshold = number(rule.get("max_af", 0))
    if threshold is None or not 0 <= threshold < 1:
        return result("PM2", input_data, Status.NOT_EVALUATED, "Invalid rarity threshold",
                      missing=["PM2.max_af"])
    maximum = max(number(item["AF"]) for item in observations)
    # A reliable counterexample defeats rarity even if another source is unavailable.
    if maximum > threshold:
        return result("PM2", input_data, Status.NOT_MET, "Observed AF exceeds configured rarity cutoff",
                      evidence=observations, provenance=provenance)
    if failures:
        return result("PM2", input_data, Status.NOT_EVALUATED, "Population search incomplete",
                      evidence=observations, missing=["complete_population_evidence"],
                      provenance=provenance)
    return result("PM2", input_data, Status.MET, "All resolved observations satisfy rarity policy",
                  strength="supporting", evidence=observations, provenance=provenance)
