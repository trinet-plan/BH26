from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import population_context, result as _base_result
from acmg_pipeline.services.population import number
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def result(code, input_data, status, summary, **kwargs):
    assert code == "PM2"
    outcome = "satisfies PM2" if status == CriterionStatus.MET else (
        "was evaluated but does not satisfy PM2" if status == CriterionStatus.NOT_MET
        else "cannot be assigned MET or NOT_MET from the available information")
    message = (f"{summary.rstrip('.')}. PM2 evaluates whether reliable population observations satisfy "
               f"the configured rarity threshold. Outcome: {status.value.upper()}; {outcome}.")
    return _base_result(code, input_data, status, message, **kwargs)


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
        return result("PM2", input_data, CriterionStatus.NOT_MET,
                      f"The highest reliable observed AF ({maximum}) exceeds the configured "
                      f"rarity cutoff ({threshold})",
                      evidence=observations, provenance=provenance)
    if failures:
        return result("PM2", input_data, CriterionStatus.UNKNOWN,
                      f"Every resolved observation is at or below the rarity cutoff ({threshold}), "
                      f"but {len(failures)} population source(s) could not be queried, so the "
                      f"search is incomplete",
                      evidence=observations, missing=["complete_population_evidence"],
                      provenance=provenance)
    return result("PM2", input_data, CriterionStatus.MET,
                  f"Every population source resolved and the highest observed AF ({maximum}) is at "
                  f"or below the configured rarity cutoff ({threshold})",
                  strength="supporting", evidence=observations, provenance=provenance)
