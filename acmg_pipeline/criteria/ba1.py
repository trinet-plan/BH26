from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import citable, population_context, result as _base_result
from acmg_pipeline.services.population import (
    FAF_METHOD, STATISTICS, article, number, observed_frequencies,
)
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


def threshold_policy(input_data, config):
    """BA1's stand-alone cutoff and the statistic it is compared with, both from policy.

    5% is the ACMG/AMP 2015 general figure, but a VCEP specification replaces it - the
    ClinGen Cardiomyopathy panel uses 0.1% - so hardcoding it made this criterion silently
    disagree with whichever specification the run was meant to follow. There is no default:
    an unconfigured threshold is reported, the way PM2's max_af stopped defaulting to zero.

    Returns (early_result, threshold, statistic).
    """
    rule = config.get("BA1", {})
    absent = [key for key in ("max_af", "frequency_statistic") if rule.get(key) is None]
    if absent:
        return result("BA1", input_data, CriterionStatus.UNKNOWN,
                      f"BA1's stand-alone threshold is not configured "
                      f"({', '.join(f'BA1.{key}' for key in absent)})",
                      missing=[f"BA1.{key}" for key in absent]), None, None
    threshold = number(rule["max_af"])
    if threshold is None or not 0 < threshold < 1:
        return result("BA1", input_data, CriterionStatus.UNKNOWN,
                      f"The configured BA1 threshold ({rule['max_af']!r}) is not a fraction "
                      f"between 0 and 1", missing=["BA1.max_af"]), None, None
    if rule["frequency_statistic"] not in STATISTICS:
        return result("BA1", input_data, CriterionStatus.UNKNOWN,
                      f"The BA1 policy names {rule['frequency_statistic']!r} as the frequency "
                      f"statistic to compare against; BA1 supports "
                      f"{' and '.join(sorted(STATISTICS))}",
                      missing=["BA1.frequency_statistic"]), None, None
    return None, threshold, rule["frequency_statistic"]


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    early, threshold, statistic = threshold_policy(input_data, config)
    if early:
        return early
    early, context = population_context("BA1", input_data, services, config)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    measure = STATISTICS[statistic]
    scored = observed_frequencies(observations, statistic)
    high = [item for item, value in scored if value is not None and value > threshold]
    highest = max((value for _, value in scored if value is not None), default=None)
    provenance = {**provenance, "frequency_statistic": statistic,
                  "stand_alone_threshold": str(threshold),
                  "highest_observed": str(highest) if highest is not None else None}
    if statistic == "faf95":
        provenance = {**provenance, "faf_method": FAF_METHOD}
    if not high:
        # "Nothing above the threshold was found" and "nothing above it exists" are different
        # claims: an unreachable provider could still hold the observation that carries BA1.
        if failures:
            return result("BA1", input_data, CriterionStatus.UNKNOWN,
                          f"No resolved observation has {article(measure)} above the configured "
                          f"stand-alone threshold ({threshold}; highest {measure} {highest}), "
                          f"but {len(failures)} population source(s) could not be queried, so "
                          f"the search is incomplete",
                          evidence=observations, missing=["complete_population_evidence"],
                          provenance=provenance)
        return result("BA1", input_data, CriterionStatus.NOT_MET,
                      f"Every population source resolved and none has {article(measure)} above the "
                      f"configured stand-alone threshold ({threshold}; highest {measure} "
                      f"{highest})",
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
    return result("BA1", input_data, CriterionStatus.MET,
                  f"{len(high)} of {len(observations)} resolved observation(s) have {article(measure)} "
                  f"above the configured stand-alone threshold ({threshold}), and the BA1 "
                  f"exception check is negative",
                  strength="stand_alone", evidence=observations + cited, provenance=provenance)
