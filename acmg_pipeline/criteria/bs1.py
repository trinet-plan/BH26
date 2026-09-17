from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import citable, population_context, result as _base_result
from acmg_pipeline.services.population import number
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def result(code, input_data, status, summary, **kwargs):
    assert code == "BS1"
    outcome = "satisfies BS1" if status == CriterionStatus.MET else (
        "was evaluated but does not satisfy BS1" if status == CriterionStatus.NOT_MET
        else "cannot be assigned MET or NOT_MET from the available information")
    message = (f"{summary.rstrip('.')}. BS1 evaluates whether population allele frequency exceeds the "
               f"reviewed maximum credible frequency for the disease context. Outcome: {status.value.upper()}; {outcome}.")
    return _base_result(code, input_data, status, message, **kwargs)


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    assessment = input_data.get("disease_frequency_threshold", {})
    threshold = number(assessment.get("max_credible_af"))
    condition = input_data.get("condition")
    # Four different situations reach this gate; a curator supplies a disease context,
    # retargets an existing threshold, or fixes a bad value, so each is reported on its own.
    if not condition:
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      "No disease context was supplied, and BS1's maximum credible frequency is "
                      "defined per disease", missing=["condition"])
    if not assessment:
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      f"No curated maximum credible allele frequency was supplied for {condition}",
                      missing=["disease_frequency_threshold"])
    if assessment.get("condition") != condition:
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      f"The curated frequency threshold is defined for "
                      f"{assessment.get('condition')!r}, not the requested {condition!r}",
                      review=["Resolve the disease context of the frequency threshold"])
    if threshold is None or not 0 < threshold < 1:
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      f"The curated maximum credible allele frequency "
                      f"({assessment.get('max_credible_af')!r}) is not a fraction between 0 and 1",
                      missing=["disease_frequency_threshold.max_credible_af"])
    absent = [key for key in ("source", "source_version", "reviewed_at", "inheritance")
              if not assessment.get(key)]
    if absent:
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      f"The curated frequency threshold does not record {', '.join(absent)}",
                      missing=[f"disease_frequency_threshold.{key}" for key in absent])
    if input_data.get("inheritance") != assessment["inheritance"]:
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      f"The threshold was derived for {assessment['inheritance']!r} inheritance, "
                      f"but this variant is assessed as {input_data.get('inheritance')!r}",
                      review=["Confirm disease inheritance"])
    early, context = population_context("BS1", input_data, services, config)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    exceeding = [item for item in observations if number(item["AF"]) > threshold]
    # One summary for all three outcomes would state the verdict without its reason, so each
    # says what was actually compared - and an incomplete search is not a negative result.
    if exceeding:
        status, strength = CriterionStatus.MET, "strong"
        summary = (f"{len(exceeding)} of {len(observations)} resolved observation(s) report an AF "
                   f"above the curated maximum credible frequency for {condition} ({threshold})")
    elif failures:
        status, strength = CriterionStatus.UNKNOWN, None
        summary = (f"No resolved observation exceeds the curated threshold for {condition} "
                   f"({threshold}), but {len(failures)} population source(s) could not be "
                   f"queried, so the search is incomplete")
    else:
        status, strength = CriterionStatus.NOT_MET, None
        summary = (f"Every population source resolved and none reports an AF above the curated "
                   f"maximum credible frequency for {condition} ({threshold})")
    # The threshold is the policy the observations are judged against: always recorded, and
    # cited as an evidence item only when it carries a retrievable identifier.
    return result("BS1", input_data, status, summary,
                  strength=strength, evidence=observations + citable(assessment),
                  missing=["complete_population_evidence"] if status == CriterionStatus.UNKNOWN else [],
                  provenance={**provenance, "disease_frequency_threshold": assessment})
