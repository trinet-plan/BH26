from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import citable, population_context, result as _base_result
from acmg_pipeline.services.population import (
    FAF_METHOD, STATISTICS, article, number, observed_frequencies,
)
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


def frequency_statistic(input_data, applied, scope):
    """Which quantity this threshold was calibrated against - AF or FAF - never a free choice.

    A threshold and the statistic it is compared with are one unit. The ClinGen Cardiomyopathy
    VCEP's 0.02% for MYH7 is a point-estimate MAF cutoff; comparing a filtering allele
    frequency with it silently makes the criterion stricter than the specification says,
    because FAF is always the lower number. gnomAD's FAF and ClinGen SVI's later advice to
    prefer it do not retroactively recalibrate specifications written before them.

    So the policy names its own statistic. A curated assessment that does not is read as the
    point-estimate cutoff its wording implies, and the assumption is reported rather than
    made quietly. The configured default is ours to state, so it has to state it.

    Returns (early_result, statistic, assumed).
    """
    declared = applied.get("frequency_statistic")
    if declared is None and scope == "disease_specific":
        return None, "af", True
    if declared not in STATISTICS:
        key = ("disease_frequency_threshold.frequency_statistic" if scope == "disease_specific"
               else "BS1.default_frequency_statistic")
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      f"The threshold names {declared!r} as the frequency statistic to compare "
                      f"against; BS1 supports {' and '.join(sorted(STATISTICS))}",
                      missing=[key]), None, False
    return None, declared, False


COMPARISONS = {">": lambda value, threshold: value > threshold,
               ">=": lambda value, threshold: value >= threshold}


def comparison(input_data, applied, scope):
    """Whether the specification says "above" the threshold or "at or above" it.

    VCEP specifications differ on this and the difference is not cosmetic: ClinGen's PAH and
    Lysosomal Storage Disorders panels write BS1 as AF > their threshold, while the Hearing
    Loss, PTEN, Glaucoma and Familial Hypercholesterolemia panels write >=. Evaluating a >=
    specification with > drops exactly the variants that sit on the boundary the panel chose,
    which are the ones it was drawing a line at.

    Absent, ">" is used: it is the weaker claim of the two, so assuming it can only withhold
    BS1, never assert it on a variant the specification would have left out. Like the
    statistic, the assumption is reported rather than made quietly.

    Returns (early_result, operator, symbol, assumed).
    """
    declared = applied.get("comparison")
    if declared is None:
        return None, COMPARISONS[">"], ">", True
    if declared not in COMPARISONS:
        key = ("disease_frequency_threshold.comparison" if scope == "disease_specific"
               else "BS1.default_comparison")
        return result("BS1", input_data, CriterionStatus.UNKNOWN,
                      f"The threshold names {declared!r} as its comparison; BS1 supports "
                      f"{' and '.join(sorted(COMPARISONS))}", missing=[key]), None, None, False
    return None, COMPARISONS[declared], declared, False


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
                   ("max_credible_af", "source", "source_version", "frequency_statistic")}
        applied["comparison"] = rule.get("default_comparison")
        threshold = number(applied["max_credible_af"])
        if not all(applied[key] for key in
                   ("max_credible_af", "source", "source_version", "frequency_statistic")):
            return result("BS1", input_data, CriterionStatus.UNKNOWN,
                          f"No disease-specific threshold applies ({reason}) and no default "
                          f"threshold is configured for BS1",
                          missing=[f"BS1.default_{key}" for key in
                                   ("max_credible_af", "source", "source_version",
                                    "frequency_statistic") if not applied[key]])
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
    early, statistic, assumed = frequency_statistic(input_data, applied, scope)
    if early:
        return early
    early, exceeds, symbol, assumed_comparison = comparison(input_data, applied, scope)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    measure = STATISTICS[statistic]
    scored = observed_frequencies(observations, statistic)
    if statistic == "faf95":
        provenance = {**provenance, "faf_method": FAF_METHOD}
    exceeding = [item for item, value in scored if value is not None and exceeds(value, threshold)]
    highest = max((value for _, value in scored if value is not None), default=None)
    provenance = {**provenance, "frequency_statistic": statistic, "comparison": symbol,
                  "highest_observed": str(highest) if highest is not None else None}
    # One summary for all three outcomes would state the verdict without its reason, so each
    # says what was actually compared - and an incomplete search is not a negative result.
    if exceeding:
        status, strength = CriterionStatus.MET, "strong"
        summary = (f"{len(exceeding)} of {len(observations)} resolved observation(s) have {article(measure)} "
                   f"{symbol} {label} ({threshold}); highest {measure} {highest}")
    elif failures:
        status, strength = CriterionStatus.UNKNOWN, None
        summary = (f"No resolved observation has {article(measure)} {symbol} {label} "
                   f"({threshold}; highest {measure} {highest}), but "
                   f"{len(failures)} population source(s) could not be queried, so the search "
                   f"is incomplete")
    else:
        status, strength = CriterionStatus.NOT_MET, None
        summary = (f"Every population source resolved and none has {article(measure)} {symbol} {label} "
                   f"({threshold}; highest {measure} {highest})")
    review = [] if scope == "disease_specific" else [
        "Confirm BS1 against a disease-specific maximum credible frequency"]
    if assumed:
        review = review + ["Confirm the curated threshold is a point-estimate AF cutoff, "
                           "not a filtering allele frequency"]
    if assumed_comparison:
        review = review + [f"Confirm the specification says {measure} > the threshold, "
                           f"not >="]
    # The threshold is the policy the observations are judged against: always recorded, and
    # cited as an evidence item only when it carries a retrievable identifier.
    return result("BS1", input_data, status, summary + fallback_note,
                  strength=strength, evidence=observations + citable(assessment), review=review,
                  missing=["complete_population_evidence"] if status == CriterionStatus.UNKNOWN else [],
                  provenance={**provenance, "threshold_scope": scope,
                              "disease_frequency_threshold": applied})
