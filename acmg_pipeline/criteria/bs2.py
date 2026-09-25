"""BS2: observed in a healthy adult for a fully penetrant, early-onset disorder.

[Why this is a reinterpretation of the project's own earlier BS2 scope decision]
  acmg_pipeline.criteria.stubs's own docstring (2026-09-xx) classified BS2 alongside PM3/
  BP2/BS4-style "patient-record" criteria - evidence that lives in the PATIENT's own
  clinical/genetic-testing records (a specific documented healthy carrier), out of scope for
  this project's literature/population-provider pipeline by design. This module instead
  automates the narrower, common real-world proxy many ClinGen VCEP specifications and
  clinical labs actually use: a gnomAD population database's own genotype counts stand in
  for "observed in [a] healthy adult individual" (gnomAD's own cohorts are curated to
  exclude early-onset severe pediatric disease, but are not confirmed disease-free for any
  other condition - the same caveat every VCEP threshold based on gnomAD counts carries).
  This is a distinct, explicit interpretive choice - like bs1.py's frequency_statistic()
  choosing AF over FAF when a specification does not say - not a discovery that the earlier
  stub classification was wrong.

[Why only autosomal-recessive/semidominant homozygotes and X-linked hemizygotes are
 automated here - autosomal-dominant/X-linked-dominant/mitochondrial are left UNKNOWN, per
 the user's explicit direction (2026-09-25)]
  For a fully penetrant, early-onset AUTOSOMAL DOMINANT disorder, a single HETEROZYGOUS
  observation in a presumed-general-population database already means the variant has a
  non-zero population frequency - which is exactly BA1/BS1's own allele-frequency signal,
  not a distinct one. Automating an AD path here would double-count the same underlying
  observation under a second criterion code. BS2's own distinct signal - a genotype dose
  (homozygous for a recessive disorder, hemizygous for an X-linked one) that a fully
  penetrant, early-onset disorder should make vanishingly rare even at real-world population
  allele frequencies - only exists for the recessive/hemizygous case. X-linked DOMINANT is
  excluded for the same double-counting reason as autosomal dominant (heterozygous females
  and hemizygous males are both directly assessed by allele frequency); mitochondrial
  inheritance has no autosomal/X homozygote-hemizygote analogue this module's data supports.

  SEMIDOMINANT (added 2026-09-25, found via live ERepo verification: LDLR/familial
  hypercholesterolemia is curated by ClinGen Gene-Disease Validity as "SD", not AD or AR) is
  read the same way as autosomal_recessive, not excluded like autosomal dominant: a
  semidominant disorder's defining feature is that heterozygotes have a milder/later-onset
  phenotype while HOMOZYGOTES have the severe, fully penetrant, early-onset one FH itself is
  the textbook example of. That is exactly BS2's own "fully penetrant... early age" question,
  answered by the homozygote count the same way autosomal_recessive's is - not a case where
  the heterozygous allele-frequency signal BA1/BS1 already covers would be double-counted,
  since a heterozygous observation says nothing about the severe homozygous presentation a
  semidominant disorder's own BS2 evidence has to be about.

[Threshold pattern mirrors ba1.py's threshold_policy(), not bs1.py's
 disease_specific_threshold()]
  BS1's per-disease threshold needs a matching disease/condition because a maximum credible
  allele frequency is derived FROM disease prevalence and penetrance. BS2's question - "is
  this genotype count too high for ANY fully penetrant, early-onset presentation of this
  gene's disease to be plausible" - does not need a condition match the way BS1's number
  does, so this mirrors BA1's simpler gene-specific-override-or-configured-default shape
  instead: input_data["bs2_threshold_override"] (automated_core/context.py's
  gene_frequency_thresholds, mirroring ba1_threshold_override) first, config["BS2"]'s
  global default otherwise - never a silent default when neither is configured.
"""

from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import (
    citable, normalize_inheritance, population_context, reviewed_or_automated,
    result as _base_result,
)
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def result(code, input_data, status, summary, **kwargs):
    assert code == "BS2"
    outcome = "satisfies BS2" if status == CriterionStatus.MET else (
        "was evaluated but does not satisfy BS2" if status == CriterionStatus.NOT_MET
        else "cannot be assigned MET or NOT_MET from the available information")
    message = (f"{summary.rstrip('.')}. BS2 evaluates whether a population database reports this "
               f"variant's genotype count (homozygotes for a recessive disorder, hemizygotes for "
               f"an X-linked one) too high for a fully penetrant, early-onset presentation to be "
               f"plausible. Outcome: {status.value.upper()}; {outcome}.")
    return _base_result(code, input_data, status, message, **kwargs)


# See this module's own docstring for why AD/XLD/mitochondrial are excluded rather than
# scored UNKNOWN-by-omission: it is a deliberate scope boundary, not a missing case.
# semidominant reads the same field as autosomal_recessive - see module docstring.
RECESSIVE_MODES = {"autosomal_recessive", "semidominant"}
HEMIZYGOUS_MODES = {"x_linked", "x_linked_recessive"}


def _resolve_inheritance(input_data, services):
    """A best-effort inheritance mode from curated gene-disease sources, tried only when
    the caller supplied none directly - added 2026-09-25 per the user's explicit direction,
    after confirming that most of this project's ERepo-sourced records have no clinical
    note to extract input_data["inheritance"] from at all, which otherwise leaves BS2
    UNKNOWN for nearly every real case regardless of how good its genotype-count logic is.

    Unanimous agreement only, across Gene2Phenotype ("gene_disease") and ClinGen Gene-
    Disease Validity ("gene_disease_validity") records for this gene - the same "offered,
    never guessed" caution pvs1.py's own _candidate_conditions() documents for the same
    real risk: some genes (that module names ABCC9, ACTB, ATP1A2, CACNA1D) are curated for
    more than one disease with DIFFERENT inheritance modes (confirmed directly in this
    project's own cached data, 2026-09-25: BRCA1 carries both an autosomal dominant cancer-
    predisposition record and an autosomal recessive Fanconi-anemia-like one). Auto-picking
    one mode for a multi-mode gene could silently pick the wrong disease's inheritance
    rather than none.

    [Why this does NOT also require the record's own condition to match
     input_data["condition"], despite that being the obvious next precision improvement]
      Tried first, and reverted after live verification (2026-09-25) against real ERepo
      ground truth: exact condition matching caused PDHA1 (unanimous "x_linked" across
      6 real records spanning 4 different MONDO IDs for what is, in every curated source,
      the same underlying disease at different granularity) to resolve to NOTHING, because
      none of those 4 MONDO IDs happened to equal input_data["condition"]'s own value - the
      same real MONDO-granularity mismatch bs1.py's own disease_specific_threshold()
      docstring already documents for BRCA2/RAF1/CDH1. A gene where every source that states
      a mode agrees is exactly the safe case this function exists to resolve, regardless of
      which of several equally-valid MONDO IDs each source happened to use; a gene where
      sources genuinely disagree (BRCA1 above) stays unresolved either way, condition
      matching or not, since nothing here ever picks between disagreeing modes.
    """
    # get_candidates() already filters to this variant's own variant_key (services/
    # evidence.py), which is gene-scoped by construction (one variant belongs to one gene) -
    # no separate gene match is needed or available here (unlike pvs1.py's annotation-
    # derived gene, criterion_input() carries no gene field of its own for this pipeline).
    variant = Variant(**input_data["variant"])
    modes = set()

    def votes(category, field, *, require_reviewed_or_automated):
        for item in services.evidence.get_candidates(category, variant):
            if require_reviewed_or_automated and not reviewed_or_automated(item):
                continue
            mode = normalize_inheritance(item.get(field))
            if mode is not None:
                modes.add(mode)

    # gene_disease (Gene2Phenotype) mechanism assessments are held to the same
    # reviewed_or_automated() bar pvs1.py's own _candidates() applies. gene_disease_validity
    # (ClinGen Gene-Disease Validity) is not - see pvs1.py's _candidate_conditions()
    # docstring for why: it is a weaker "gene and disease are related" statement, not an
    # individually reviewed assessment, so the same bar would only drop real context.
    votes("gene_disease", "inheritance", require_reviewed_or_automated=True)
    votes("gene_disease_validity", "moi", require_reviewed_or_automated=False)

    return modes.pop() if len(modes) == 1 else None


def zygosity_field(inheritance):
    """Which gnomAD genotype-count field BS2 checks for this inheritance mode.

    Returns the observation field name to read, or None when this inheritance mode is out
    of this module's automated scope (see module docstring).
    """
    if inheritance in RECESSIVE_MODES:
        return "homozygote_count"
    if inheritance in HEMIZYGOUS_MODES:
        return "hemizygote_count"
    return None


APPROVED, DRAFT = "APPROVED", "DRAFT"


def policy_status(scope, rule):
    """Whether this threshold may carry a final judgment, or only a draft for review.

    Same convention as bs1.py's policy_status(): a gene-specific threshold is approved
    unless it says of itself that it is not (a curator part-way through transcribing a VCEP
    specification can mark it DRAFT and have that honoured); the configured default is not
    reviewed for any particular gene, so it is never approved. Today config["BS2"]'s default
    is the ONLY threshold this project has (no curated bs2_threshold_override entries exist
    yet - see config/demo-rules.json's own BS2.source note), so every BS2 MET this module
    produces is currently a flagged prediction, not a final call, until real per-gene VCEP
    numbers are curated.
    """
    if scope != "gene_specific":
        return DRAFT
    return DRAFT if str(rule.get("policy_status", "")).upper() == DRAFT else APPROVED


def threshold_policy(input_data, config):
    """BS2's max permissible genotype count, and where it came from.

    Mirrors ba1.py's threshold_policy(): input_data["bs2_threshold_override"] (a real
    gene-specific VCEP number, via automated_core/context.py's gene_frequency_thresholds)
    is checked first and, when present and complete, replaces config["BS2"]'s global default
    entirely for this variant's gene. Neither has an implicit default - an unconfigured
    value is reported rather than assumed to be 0.

    Returns (early_result, max_count, scope, rule).
    """
    override = input_data.get("bs2_threshold_override")
    if override and all(
        override.get(key) for key in ("source", "source_version", "reviewed_at")
    ) and isinstance(override.get("max_count"), int) and not isinstance(override.get("max_count"), bool):
        rule, scope = override, "gene_specific"
    else:
        rule, scope = config.get("BS2", {}), "default"
    label = "The gene-specific BS2 threshold" if scope == "gene_specific" else "BS2's configured default threshold"
    max_count = rule.get("max_count")
    if not isinstance(max_count, int) or isinstance(max_count, bool) or max_count < 0:
        missing_key = "bs2_threshold_override" if scope == "gene_specific" else "BS2"
        return result("BS2", input_data, CriterionStatus.UNKNOWN,
                      f"{label} is not configured with a non-negative integer max_count",
                      missing=[f"{missing_key}.max_count"]), None, None, None
    return None, max_count, scope, rule


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    inheritance = normalize_inheritance(input_data.get("inheritance"))
    inheritance_source = "supplied" if inheritance else None
    if inheritance is None:
        inheritance = _resolve_inheritance(input_data, services)
        if inheritance is not None:
            inheritance_source = "resolved_from_gene_disease_evidence"
    field = zygosity_field(inheritance)
    if field is None:
        return result("BS2", input_data, CriterionStatus.UNKNOWN,
                      "BS2 is automated for autosomal-recessive (homozygote count) and X-linked "
                      "(hemizygote count) inheritance only" +
                      (f"; this variant is assessed as {inheritance!r} inheritance"
                       if inheritance else "; no inheritance mode was supplied or unanimously "
                       "resolved from curated gene-disease evidence"),
                      missing=["inheritance"])
    if field == "hemizygote_count" and input_data["variant"]["chrom"] != "X":
        return result("BS2", input_data, CriterionStatus.UNKNOWN,
                      f"This variant is assessed as {inheritance!r} inheritance but its chromosome "
                      f"({input_data['variant']['chrom']}) is not X, so a hemizygote count would "
                      f"not be biologically meaningful",
                      missing=["inheritance"])
    early, max_count, scope, rule = threshold_policy(input_data, config)
    if early:
        return early
    early, context = population_context("BS2", input_data, services, config)
    if early:
        return early
    _, observations, rejected, failures, provenance = context
    counted = [(item, item.get(field)) for item in observations]
    invalid = [item for item, value in counted
              if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 0)]
    if invalid:
        return result("BS2", input_data, CriterionStatus.UNKNOWN,
                      f"{len(invalid)} of {len(observations)} resolved observation(s) report an "
                      f"invalid {field}",
                      evidence=observations, missing=[field], provenance=provenance)
    exceeding = [item for item, value in counted if value is not None and value > max_count]
    countable = [value for _, value in counted if value is not None]
    highest = max(countable, default=None)
    status_of_policy = policy_status(scope, rule)
    provenance = {**provenance, "zygosity_field": field, "inheritance": inheritance,
                  "inheritance_source": inheritance_source,
                  "threshold_scope": scope, "max_count": max_count, "highest_observed": highest,
                  "policy_status": status_of_policy}
    if scope == "gene_specific":
        provenance["gene_threshold"] = rule
    cited_threshold = citable(rule) if scope == "gene_specific" else []
    mode_label = inheritance.replace("_", " ")
    if exceeding:
        summary = (f"{len(exceeding)} of {len(observations)} resolved observation(s) report "
                   f"{field} > {max_count} (highest {highest}), more than expected in a "
                   f"presumed-general population for a fully penetrant, early-onset "
                   f"{mode_label} disorder")
        review = []
        if status_of_policy == DRAFT:
            # An unapproved threshold still reflects a real comparison against the best
            # available number, so it is reported as a (flagged) prediction rather than
            # withheld outright - same "report, don't withhold" choice as bs1.py's own MET-
            # while-DRAFT branch.
            summary += (", but this threshold is not approved for this gene, so this MET is "
                       "a prediction pending curator sign-off, not a final call")
            review = ["Approve or replace the draft BS2 genotype-count threshold before "
                      "BS2 is used in a classification"]
        return result("BS2", input_data, CriterionStatus.MET, summary,
                      strength="strong", evidence=observations + cited_threshold,
                      review=review, provenance=provenance)
    if not countable:
        return result("BS2", input_data, CriterionStatus.UNKNOWN,
                      f"No resolved observation reports a {field}, so BS2 cannot be evaluated "
                      f"from the available population data",
                      evidence=observations, missing=[field], provenance=provenance)
    if failures:
        return result("BS2", input_data, CriterionStatus.UNKNOWN,
                      f"No resolved observation reports {field} > {max_count} (highest {highest}), "
                      f"but {len(failures)} population source(s) could not be queried, so the "
                      f"search is incomplete",
                      evidence=observations, missing=["complete_population_evidence"],
                      provenance=provenance)
    return result("BS2", input_data, CriterionStatus.NOT_MET,
                  f"Every population source resolved and none reports {field} > {max_count} "
                  f"(highest {highest})",
                  evidence=observations + cited_threshold, provenance=provenance)
