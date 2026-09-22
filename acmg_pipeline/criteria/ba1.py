from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import citable, population_context, result as _base_result
from acmg_pipeline.services.population import (
    COMPARISONS, FAF_METHOD, STATISTICS, article, number, observed_frequencies,
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
    """BA1's stand-alone cutoff, the statistic and comparison it is checked with, and where
    they came from - all from policy, never hardcoded.

    5% is the ACMG/AMP 2015 general figure, but a VCEP specification replaces it - e.g. the
    ClinGen Cardiomyopathy panel's real MYH7 specification uses 0.1% FAF popmax, ClinGen
    RASopathy's RAF1 specification 0.05% - so hardcoding either the number or the statistic
    made this criterion silently disagree with whichever specification the run was meant to
    follow. input_data["ba1_threshold_override"] (see automated_core/context.py's
    gene_frequency_thresholds, populated from real ClinGen CSpec data - not derived here)
    is checked first and, when present and complete, replaces config["BA1"]'s global
    default entirely for this variant's gene; there is no partial override.

    Neither threshold has a default: an unconfigured one is reported, the way PM2's max_af
    stopped defaulting to zero.

    Returns (early_result, threshold, statistic, comparator, symbol, scope, rule).
    """
    override = input_data.get("ba1_threshold_override")
    if override and all(
        override.get(key) for key in
        ("max_af", "frequency_statistic", "source", "source_version", "reviewed_at")
    ):
        rule, scope = override, "gene_specific"
    else:
        rule, scope = config.get("BA1", {}), "default"
    label = "The gene-specific BA1 threshold" if scope == "gene_specific" else "BA1's stand-alone threshold"
    absent = [key for key in ("max_af", "frequency_statistic") if rule.get(key) is None]
    if absent:
        missing_key = "ba1_threshold_override" if scope == "gene_specific" else "BA1"
        return result("BA1", input_data, CriterionStatus.UNKNOWN,
                      f"{label} is not configured "
                      f"({', '.join(f'{missing_key}.{key}' for key in absent)})",
                      missing=[f"{missing_key}.{key}" for key in absent]), None, None, None, None, None, None
    threshold = number(rule["max_af"])
    if threshold is None or not 0 < threshold < 1:
        return result("BA1", input_data, CriterionStatus.UNKNOWN,
                      f"{label} ({rule['max_af']!r}) is not a fraction between 0 and 1",
                      missing=["ba1_threshold_override.max_af" if scope == "gene_specific"
                              else "BA1.max_af"]), None, None, None, None, None, None
    if rule["frequency_statistic"] not in STATISTICS:
        return result("BA1", input_data, CriterionStatus.UNKNOWN,
                      f"{label} names {rule['frequency_statistic']!r} as the frequency "
                      f"statistic to compare against; BA1 supports "
                      f"{' and '.join(sorted(STATISTICS))}",
                      missing=["ba1_threshold_override.frequency_statistic" if scope == "gene_specific"
                              else "BA1.frequency_statistic"]), None, None, None, None, None, None
    symbol = rule.get("comparison") or ">"
    if symbol not in COMPARISONS:
        return result("BA1", input_data, CriterionStatus.UNKNOWN,
                      f"{label} names {symbol!r} as its comparison; BA1 supports "
                      f"{' and '.join(sorted(COMPARISONS))}",
                      missing=["ba1_threshold_override.comparison" if scope == "gene_specific"
                              else "BA1.comparison"]), None, None, None, None, None, None
    return None, threshold, rule["frequency_statistic"], COMPARISONS[symbol], symbol, scope, rule


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    early, threshold, statistic, exceeds, symbol, scope, rule = threshold_policy(input_data, config)
    if early:
        return early
    early, context = population_context("BA1", input_data, services, config)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    measure = STATISTICS[statistic]
    scored = observed_frequencies(observations, statistic)
    high = [item for item, value in scored if value is not None and exceeds(value, threshold)]
    highest = max((value for _, value in scored if value is not None), default=None)
    provenance = {**provenance, "frequency_statistic": statistic, "comparison": symbol,
                  "threshold_scope": scope,
                  "stand_alone_threshold": str(threshold),
                  "highest_observed": str(highest) if highest is not None else None}
    if scope == "gene_specific":
        provenance = {**provenance, "gene_threshold": rule}
    if statistic == "faf95":
        provenance = {**provenance, "faf_method": FAF_METHOD}
    cited_threshold = citable(rule) if scope == "gene_specific" else []
    if not high:
        # "Nothing above the threshold was found" and "nothing above it exists" are different
        # claims: an unreachable provider could still hold the observation that carries BA1.
        if failures:
            return result("BA1", input_data, CriterionStatus.UNKNOWN,
                          f"No resolved observation has {article(measure)} {symbol} the configured "
                          f"stand-alone threshold ({threshold}; highest {measure} {highest}), "
                          f"but {len(failures)} population source(s) could not be queried, so "
                          f"the search is incomplete",
                          evidence=observations, missing=["complete_population_evidence"],
                          provenance=provenance)
        return result("BA1", input_data, CriterionStatus.NOT_MET,
                      f"Every population source resolved and none has {article(measure)} {symbol} the "
                      f"configured stand-alone threshold ({threshold}; highest {measure} "
                      f"{highest})",
                      evidence=observations + cited_threshold, provenance=provenance)
    exception = input_data.get("ba1_exception_assessment", {})
    if not all(exception.get(key) for key in ("source", "source_version", "reviewed_at")):
        return result("BA1", input_data, CriterionStatus.UNKNOWN, "BA1 exception check unavailable",
                      evidence=observations + cited_threshold, missing=["ba1_exception_assessment"],
                      provenance=provenance)
    # The exception list is curated policy: it is always recorded, and cited as an evidence
    # item only when it carries a retrievable identifier.
    provenance = {**provenance, "ba1_exception_assessment": exception}
    cited = citable(exception)
    if exception.get("is_exception") is True:
        return result("BA1", input_data, CriterionStatus.NOT_MET, "Known BA1 exception",
                      evidence=observations + cited_threshold + cited, provenance=provenance)
    if exception.get("is_exception") is not False:
        return result("BA1", input_data, CriterionStatus.UNKNOWN, "BA1 exception status unresolved",
                      evidence=observations + cited_threshold + cited, review=["Confirm exception status"],
                      provenance=provenance)
    return result("BA1", input_data, CriterionStatus.MET,
                  f"{len(high)} of {len(observations)} resolved observation(s) have {article(measure)} "
                  f"{symbol} the configured stand-alone threshold ({threshold}), and the BA1 "
                  f"exception check is negative",
                  strength="stand_alone", evidence=observations + cited_threshold + cited,
                  provenance=provenance)
