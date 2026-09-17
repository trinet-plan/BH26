from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import population_context, result
from acmg_pipeline.services.population import number
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    early, context = population_context("PM2", input_data, services, config)
    if early:
        return early
    rule, observations, rejected, failures, provenance = context
    threshold = number(rule.get("max_af", 0))
    if threshold is None or not 0 <= threshold < 1:
        return result("PM2", input_data, CriterionStatus.UNKNOWN, "Invalid rarity threshold",
                      missing=["PM2.max_af"])
    maximum = max(number(item["AF"]) for item in observations)
    # A reliable counterexample defeats rarity even if another source is unavailable.
    if maximum > threshold:
        return result("PM2", input_data, CriterionStatus.NOT_MET, "Observed AF exceeds configured rarity cutoff",
                      evidence=observations, provenance=provenance)
    if failures:
        return result("PM2", input_data, CriterionStatus.UNKNOWN, "Population search incomplete",
                      evidence=observations, missing=["complete_population_evidence"],
                      provenance=provenance)
    return result("PM2", input_data, CriterionStatus.MET, "All resolved observations satisfy rarity policy",
                  strength="supporting", evidence=observations, provenance=provenance)
