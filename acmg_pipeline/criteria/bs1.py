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


def disease_specific_threshold(input_data):
    """Why the curated threshold is or is not usable for this record's disease context.

    Returns (threshold, assessment, reason) - `reason` is None when the threshold applies,
    and otherwise says, for a curator, what would make it apply. A threshold that is absent,
    aimed at another disease, missing its provenance or derived for another inheritance mode
    is not a disease-specific threshold for THIS record, so all four route to the configured
    default rather than each stopping the criterion. A value outside (0, 1) is different: the
    threshold is present and wrong, so it is reported instead of being papered over.
    """
    condition = input_data.get("condition")
    assessment = input_data.get("disease_frequency_threshold", {})
    if not condition:
        return None, assessment, "no disease context was supplied"
    if not assessment:
        return None, assessment, f"no curated maximum credible allele frequency was supplied for {condition}"
    if assessment.get("condition") != condition:
        return None, assessment, (f"the curated threshold is defined for "
                                  f"{assessment.get('condition')!r}, not the requested {condition!r}")
    absent = [key for key in ("source", "source_version", "reviewed_at", "inheritance")
              if not assessment.get(key)]
    if absent:
        return None, assessment, f"the curated threshold does not record {', '.join(absent)}"
    if input_data.get("inheritance") != assessment["inheritance"]:
        return None, assessment, (f"the curated threshold was derived for "
                                  f"{assessment['inheritance']!r} inheritance, but this variant is "
                                  f"assessed as {input_data.get('inheritance')!r}")
    return number(assessment["max_credible_af"]), assessment, None


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    condition = input_data.get("condition")
    threshold, assessment, reason = disease_specific_threshold(input_data)
    rule = config.get("BS1", {})
    if reason is None:
        scope, applied, fallback_note = "disease_specific", assessment, ""
        label = f"the curated maximum credible frequency for {condition}"
    else:
        # A default threshold keeps BS1 evaluable without a per-disease review, but it is a
        # weaker claim than the criterion describes, so the substitution travels with the
        # result instead of disappearing into a number.
        applied = {key: rule.get(f"default_{key}") for key in
                   ("max_credible_af", "source", "source_version")}
        threshold = number(applied["max_credible_af"])
        if not all(applied.values()):
            return result("BS1", input_data, CriterionStatus.UNKNOWN,
                          f"No disease-specific threshold applies ({reason}) and no default "
                          f"threshold is configured for BS1",
                          missing=[f"BS1.default_{key}" for key, value in applied.items()
                                   if not value])
        scope = "default"
        fallback_note = (f" (applied because {reason}; this is not a disease-specific "
                         f"threshold and has to be confirmed before final classification)")
        label = (f"the configured default maximum credible frequency "
                 f"[{applied['source']} {applied['source_version']}]")
    if threshold is None or not 0 < threshold < 1:
        source = "curated" if scope == "disease_specific" else "configured default"
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      f"The {source} maximum credible allele frequency "
                      f"({applied.get('max_credible_af')!r}) is not a fraction between 0 and 1",
                      missing=["disease_frequency_threshold.max_credible_af"
                               if scope == "disease_specific" else "BS1.default_max_credible_af"])
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
                   f"above {label} ({threshold})")
    elif failures:
        status, strength = CriterionStatus.UNKNOWN, None
        summary = (f"No resolved observation exceeds {label} ({threshold}), but "
                   f"{len(failures)} population source(s) could not be queried, so the search "
                   f"is incomplete")
    else:
        status, strength = CriterionStatus.NOT_MET, None
        summary = (f"Every population source resolved and none reports an AF above {label} "
                   f"({threshold})")
    review = [] if scope == "disease_specific" else [
        "Confirm BS1 against a disease-specific maximum credible frequency"]
    # The threshold is the policy the observations are judged against: always recorded, and
    # cited as an evidence item only when it carries a retrievable identifier.
    return result("BS1", input_data, status, summary + fallback_note,
                  strength=strength, evidence=observations + citable(assessment), review=review,
                  missing=["complete_population_evidence"] if status == CriterionStatus.UNKNOWN else [],
                  provenance={**provenance, "threshold_scope": scope,
                              "disease_frequency_threshold": applied})
