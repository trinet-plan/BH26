from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.models import CriterionResult, Variant
from acmg_pipeline.services.population import number, usable_observations


DEFAULT_STRENGTH = {
    "PVS1": "very_strong", "PS1": "strong", "PM1": "moderate", "PM2": "moderate",
    "PM4": "moderate", "PM5": "moderate", "PP2": "supporting", "PP3": "supporting",
    "BA1": "stand_alone", "BS1": "strong", "BP1": "supporting", "BP3": "supporting",
    "BP4": "supporting", "BP7": "supporting",
}


def result(code, input_data, status, summary, *, strength=None, evidence=None,
           missing=None, review=None, provenance=None, evaluation_context=None,
           decision_trace=None, rules_used=None, warnings=None):
    direction, outcome = None, None
    if status == CriterionStatus.MET:
        direction = "disputes" if code.startswith("B") else "supports"
        outcome = code if strength == DEFAULT_STRENGTH.get(code) else f"{code}_{strength}"
    elif status == CriterionStatus.NOT_MET:
        direction, outcome = "none", f"{code}_not_met"
    return CriterionResult(
        code, status, input_data["variant"], summary, strength, direction, outcome,
        evidence=evidence or [], missing_inputs=missing or [], review_points=review or [],
        provenance={"rule_version": f"{code}-v1", **(provenance or {})},
        evaluation_context=evaluation_context, decision_trace=decision_trace or [],
        rules_used=rules_used or [], warnings=warnings or [],
    )


# "This criterion does not apply to this variant" and "this criterion could not be decided"
# are both UNKNOWN - the status enum has no fourth value, and adding one would change the
# result schema every downstream consumer reads.  Marking the first kind inside the existing
# provenance dict keeps the two separable when results are counted, without that change:
# an inapplicable criterion needs no further evidence, an indeterminate one does.
NOT_APPLICABLE = {"assessment_outcome": "not_applicable"}


def citable(assessment):
    """Curated policy is exported as an evidence item only when it has a retrievable IRI.

    An assessment typed inline in the prepared input has no identifier, so it is recorded as
    provenance instead; every exported evidence item must be resolvable.
    """
    return [assessment] if assessment.get("evidence_id") else []


def population_context(code, input_data, services, config, *, treat_total_absence_as_evidence=False):
    """`treat_total_absence_as_evidence`: only PM2 passes this. For a rarity-seeking
    criterion, a variant that no queried population source returned ANYTHING for (not
    even a rejected/low-quality record) is itself a weak signal of rarity, not merely an
    unanswerable gap - unlike BA1/BS1, where the same silence cannot argue a variant is
    common. When true and every provider resolved cleanly with nothing to report (no real
    fetch error - see the NO_OBSERVATION check below), the early UNKNOWN below is skipped
    and an empty-but-valid context is returned instead, so the caller can score the absence
    itself (with its own caveat) rather than being forced into UNKNOWN here. A provider that
    could not be reached at all (a real error, not a confirmed empty result) still blocks
    this path - that is a search gap, not an observation.
    """
    rule = config.get(code, {})
    minimum_an = number(rule.get("minimum_an"))
    if minimum_an is None or minimum_an < 1 or minimum_an != minimum_an.to_integral_value():
        return result(code, input_data, CriterionStatus.UNKNOWN, "Population quality policy missing",
                      missing=[f"{code}.minimum_an"]), None
    if not all(rule.get(key) for key in ("policy_source", "policy_version")):
        return result(code, input_data, CriterionStatus.UNKNOWN, "Population policy provenance missing",
                      missing=[f"{code}.policy_source", f"{code}.policy_version"]), None
    variant = Variant(**input_data["variant"])
    resolved = services.population.get_resolved_evidence(variant, input_data)
    valid, rejected = usable_observations(resolved, variant, minimum_an)
    provenance = {"policy_source": rule["policy_source"], "policy_version": rule["policy_version"],
                  "rejected_observations": rejected, "provider_failures": resolved["failures"]}
    if not valid:
        total_absence = (
            treat_total_absence_as_evidence
            and not resolved["observations"]
            and all(f.get("reason") == "NO_OBSERVATION" for f in resolved["failures"])
        )
        if total_absence:
            return None, (rule, [], rejected, resolved["failures"], provenance)
        return result(code, input_data, CriterionStatus.UNKNOWN, "No reliable population observation",
                      evidence=resolved["observations"], missing=["population"],
                      provenance=provenance), None
    return None, (rule, valid, rejected, resolved["failures"], provenance)


def get_evidence(category, input_data, services):
    return services.evidence.get(category, Variant(**input_data["variant"]), input_data)


# Which provider supplies each evidence category, so a criterion left with nothing can name
# the one that failed rather than every failure of the run. The comparator search is absorbed
# per criterion ("ClinVar protein comparator search:PS1"), so these match on the name before
# the colon as well as on the whole.
CATEGORY_PROVIDERS = {
    "region": ("ClinVar protein hotspot density",),
    "comparator_search": ("ClinVar protein comparator search",),
}


def retrieval_failures(services, provider=None, *, exact=False):
    """What the resolver could not retrieve, so a criterion can say why it has nothing.

    Without this, a provider that errored and a provider that legitimately returned nothing
    are indistinguishable downstream, and the curator is told only that evidence is absent -
    which is the one thing they could already see.

    `provider` is a name or names. By default a failure also matches that provider's
    per-criterion entries ("<name>:PS1"); with `exact`, only the names given match - PS1's
    search failing is not PM5's reason for having nothing.
    """
    failures = getattr(services, "failures", None) or []
    if provider is None:
        return list(failures)
    names = (provider,) if isinstance(provider, str) else tuple(provider)
    return [item for item in failures
            if any(str(item.get("provider")) == name
                   or (not exact and str(item.get("provider")).startswith(f"{name}:"))
                   for name in names)]


def failure_detail(failures):
    """The failures as one readable clause, or "" when nothing failed."""
    return "; ".join(f"{item.get('provider')}: {item.get('error')}" for item in failures)


def annotation_context(code, input_data, services):
    annotations = get_evidence("annotation", input_data, services)
    if not annotations:
        failures = retrieval_failures(services)
        detail = failure_detail(failures)
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      f"Transcript annotation unavailable - {detail}" if detail
                      else "Transcript annotation unavailable; no provider reported an error, "
                           "so the annotation source returned nothing for this variant",
                      missing=["annotation", "transcript"],
                      provenance={"provider_failures": failures}), None
    if len(annotations) != 1:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Multiple transcript annotations require resolution",
                      evidence=annotations, review=["Select disease-relevant transcript"]), None
    annotation = annotations[0]
    # Only the fields actually absent are reported: naming a field the annotation does carry
    # sends a curator looking for evidence that is already there.
    absent = [field for field in ("consequences", "transcript") if not annotation.get(field)]
    if absent:
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      f"Incomplete transcript annotation: {', '.join(absent)} missing",
                      evidence=annotations, missing=absent), None
    return None, annotation


def reviewed_or_automated(record):
    """Human curation, or an automated assessment that names its versioned policy."""
    if record.get("assessment_method") == "automated":
        return bool(record.get("method") and record.get("policy_version"))
    return bool(record.get("curator") and record.get("reviewed_at"))


def transcript_compatible(record, annotation):
    """Whether this record's scope covers the annotated transcript.

    A record that names no transcript is not a record about a different one: a gene-level
    mechanism statement (ClinGen Dosage, Gene2Phenotype) has no transcript to name, and
    reading its absence as a mismatch told curators to "resolve the transcript" of a record
    that never had one. This is the same rule PVS1's own _candidates() already applies -
    the two matchers disagreeing was accidental, not a policy.

    What actually keeps a gene-level record out of a criterion that cannot use it is
    `required_fields`, which asks whether the record answers this criterion's question.
    """
    return not record.get("transcript") or record.get("transcript") == annotation["transcript"]


def answers_criterion(record, required_fields, answered_by=()):
    """Whether this record answers the question the criterion is about to ask.

    One evidence category can have several producers with different scopes. `gene_disease`
    carries both a gene-level LoF statement (which answers PVS1's mechanism gate and nothing
    else) and a transcript-scoped mechanism and spectrum review (which answers PP2/BP1 too).
    Selecting on the fields a criterion actually reads says which records can answer it,
    instead of leaning on a proxy - and lets both producers be supplied at once without the
    ones that cannot answer being counted as competing assessments.

    `answered_by` names the fields that answer the criterion on their own. An expert panel
    that marks a criterion not applicable to its gene-disease pair has answered it outright,
    and deliberately records no mechanism or spectrum finding to go with that - requiring
    the fields of a judgment the panel said not to make would discard the panel's decision.

    Presence is the test here, not type: require_boolean_fields() still rules on whether a
    present value is an explicit true/false, so a half-filled record keeps its own, more
    specific message rather than being silently dropped here.
    """
    if any(record.get(field) is not None for field in answered_by):
        return True
    return all(record.get(field) is not None for field in required_fields)


def unusable_reason(category, records, annotation, condition, disease_required,
                    required_fields=(), answered_by=()):
    """Why every retrieved `category` record was rejected - one cause per rejection route.

    "No usable assessment" covers five situations a curator has to act on differently:
    nothing was retrieved at all (produce the assessment), records exist but name no
    reviewer or policy version (record the provenance), records were reviewed against a
    different transcript (resolve the transcript), against a different disease (resolve the
    disease context), or they are a different kind of assessment that does not carry what
    this criterion reads (produce the assessment this criterion needs - the transcript and
    the disease are not the problem).  Reporting them as one message would tell a curator to
    create evidence that already exists, or to fix a field that is not what stopped the run.
    Returns (summary, missing, review).
    """
    if not records:
        return (f"No {category} assessment was retrieved for this variant", [category], [])
    unreviewed = [r for r in records if not reviewed_or_automated(r)]
    reviewed = [r for r in records if reviewed_or_automated(r)]
    transcript_mismatch = [r for r in reviewed if not transcript_compatible(r, annotation)]
    in_scope = [r for r in reviewed if transcript_compatible(r, annotation)]
    condition_mismatch = [r for r in in_scope
                          if disease_required and r.get("condition") != condition]
    # Reported before the two mismatches: these records are in the right transcript and
    # disease context and were still rejected, so they got furthest and theirs is the
    # specific cause. A curator sent after a transcript here would find nothing wrong.
    wrong_kind = [r for r in in_scope
                  if (not disease_required or r.get("condition") == condition)
                  and not answers_criterion(r, required_fields, answered_by)]
    if wrong_kind:
        absent = sorted({field for r in wrong_kind for field in required_fields
                         if r.get(field) is None})
        sources = sorted({str(r.get("source")) for r in wrong_kind})
        return (f"{len(wrong_kind)} {category} assessment(s) apply to this variant "
                f"({', '.join(sources)}) but are a different kind of assessment: they do not "
                f"record {', '.join(absent)}, which this criterion reads",
                list(absent),
                [f"Supply a {category} assessment that records {', '.join(absent)}"])
    if condition_mismatch:
        seen = sorted({str(r.get("condition")) for r in condition_mismatch})
        return (f"{len(condition_mismatch)} reviewed {category} assessment(s) cover the requested "
                f"transcript but a different disease context ({', '.join(seen)}, not {condition})",
                [], [f"Resolve the disease context of the {category} assessment"])
    if transcript_mismatch:
        seen = sorted({str(r.get("transcript")) for r in transcript_mismatch})
        return (f"{len(transcript_mismatch)} reviewed {category} assessment(s) were made against "
                f"a different transcript ({', '.join(seen)}, not {annotation['transcript']})",
                [], [f"Resolve the transcript of the {category} assessment"])
    return (f"{len(unreviewed)} {category} assessment(s) were retrieved but none names a curator "
            f"and review date, or an assessment method and policy version",
            [f"{category}.review_provenance"],
            [f"Record review provenance for the {category} assessment"])


def curated_context(code, category, input_data, services, annotation, *, disease_required=True,
                   required_fields=(), answered_by=()):
    """One reviewed assessment in this disease/transcript scope that answers this criterion.

    `required_fields` are the fields the caller's judgment reads. They select, not just
    validate: a category can hold assessments of different kinds, and one that does not
    record what the criterion reads is not a competing assessment of the same question - it
    is an assessment of a different one. Passing them here keeps such a record from being
    counted as a conflicting second opinion, and makes the rejection say what is actually
    missing (see answers_criterion() and unusable_reason()).
    """
    if disease_required and not input_data.get("condition"):
        return result(code, input_data, CriterionStatus.UNKNOWN, "Disease context required",
                      evidence=[annotation], missing=["condition"]), None
    condition = input_data.get("condition")
    retrieved = get_evidence(category, input_data, services)
    records = [r for r in retrieved if reviewed_or_automated(r)
               and transcript_compatible(r, annotation)
               and (not disease_required or r.get("condition") == condition)
               and answers_criterion(r, required_fields, answered_by)]
    if not records:
        summary, missing, review = unusable_reason(category, retrieved, annotation, condition,
                                                   disease_required, required_fields,
                                                   answered_by)
        # Only when nothing came back at all: if records were retrieved and rejected,
        # unusable_reason already says which rejection route they took, and a provider that
        # failed elsewhere in the run did not cause that.
        failures = (retrieval_failures(services, CATEGORY_PROVIDERS.get(category, ()))
                    if not retrieved else [])
        detail = failure_detail(failures)
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      f"{summary} - {detail}" if detail else summary,
                      evidence=[annotation, *retrieved], missing=missing, review=review,
                      provenance={"provider_failures": failures} if failures else None), None
    if len(records) != 1:
        return result(code, input_data, CriterionStatus.UNKNOWN, f"Multiple {category} assessments",
                      evidence=[annotation, *records], review=[f"Resolve {category} assessments"]), None
    return None, records[0]


def require_boolean_fields(code, input_data, evidence, assessment, fields):
    """Every field must be an explicit true/false - absent and null are not "false"."""
    missing = [field for field in fields if type(assessment.get(field)) is not bool]
    if missing:
        # Six criteria share this gate, so the summary names the criterion and the fields:
        # "Assessment incomplete" alone left a curator to guess which of them stalled.
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      f"The reviewed assessment records no explicit true/false value for "
                      f"{', '.join(missing)}, which {code} requires",
                      evidence=evidence, missing=missing)
    return None


# The inheritance-mode vocabularies this pipeline has to reconcile. Curated specifications
# write "AD"/"AR", clinical notes write "autosomal dominant", and prepared records have been
# seen carrying "autosomal_recessive". Comparing those as raw strings silently fails to
# match, and a silent non-match here is worse than a loud one: it drops a disease-specific
# assessment and falls back to a weaker default without saying so.
INHERITANCE_MODES = {
    "autosomal_dominant": ("ad", "autosomal dominant", "autosomal dominant inheritance"),
    "autosomal_recessive": ("ar", "autosomal recessive", "autosomal recessive inheritance"),
    "x_linked_dominant": ("xld", "x linked dominant", "x linked dominant inheritance"),
    "x_linked_recessive": ("xlr", "x linked recessive", "x linked recessive inheritance"),
    "x_linked": ("xl", "x linked", "x linked inheritance"),
    "mitochondrial": ("mt", "mitochondrial", "mitochondrial inheritance"),
}

_INHERITANCE_ALIASES = {alias: canonical
                       for canonical, aliases in INHERITANCE_MODES.items()
                       for alias in (canonical.replace("_", " "), *aliases)}


def normalize_inheritance(value):
    """Canonical inheritance-mode token, or None when absent or unrecognized.

    Absent and unrecognized deliberately collapse to None at this level: both mean "this
    string does not name a mode we can compare". Callers that must tell them apart check the
    raw value first, because an unrecognized mode is a curation defect a curator should see,
    while an absent one just means the record is not scoped to a mode.
    """
    if not isinstance(value, str):
        return None
    token = " ".join(value.strip().lower().replace("-", " ").replace("_", " ").split())
    return _INHERITANCE_ALIASES.get(token)


def resolved_condition(holder):
    """The disease a case or record is scoped to, and how that identifier was arrived at.

    A case and a curation can name the same disease in different vocabularies, so comparing
    the identifiers as written answers only when both happen to use the same one. When a
    `condition_mapping` is attached, the identifier it resolved to is what the two are
    compared on, and the mapping type it came from says whether that was an identifier match
    or an equivalence. Without one, the condition stands for itself.

    Returns (identifier, how) where `how` is "identity", the mapping's own type, or None when
    no disease is named at all.
    """
    mapping = holder.get("condition_mapping")
    if isinstance(mapping, dict) and mapping.get("normalized_condition"):
        return mapping["normalized_condition"], mapping.get("mapping_type") or "equivalent"
    condition = holder.get("condition")
    return condition, "identity" if condition else None


def condition_ancestors(holder):
    """The MONDO terms this holder's disease sits under, as far as they were resolved."""
    value = holder.get("condition_ancestors")
    if isinstance(value, dict):
        value = value.get("ancestors")
    return set(value) if isinstance(value, (list, set, tuple)) else set()


# MONDO condition pairs known, by manual review, to name the same real-world disease
# concept despite carrying no formal ontology link (no is_a relationship in either
# direction, no shared database_cross_reference) - a MONDO curation gap that no
# ancestor-based or xref-based check can ever find on its own. Verified directly against
# OLS4 (2026-09-22), found on the 679-variant ground truth run:
#   MONDO:0007648 "hereditary diffuse gastric adenocarcinoma" is the term ClinGen's own
#     CDH1 gene-disease validity curation is filed under (its "curated content resource"
#     annotation points at ClinGen's own condition page for this exact ID).
#   MONDO:0100488 "CDH1-related diffuse gastric and lobular breast cancer syndrome" is a
#     newer, broader MONDO term - cross-referencing OMIM:137215, which 0007648 does not -
#     covering the same CDH1-driven cancer predisposition, now naming the lobular breast
#     cancer component too. Its full ancestor chain (fetched live via OLS4) does not
#     include 0007648, and neither term's database_cross_reference list overlaps the
#     other's, so no existing automated check can connect them.
# Several CDH1 variants in this project's ground truth are exactly this: the gene's
# curated LoF mechanism is filed under 0007648, the case is filed under 0100488, and
# PVS1's mechanism gate (the only caller of ontology_related()) had no way to see they
# are the same disease. Add a pair here only after the same manual verification (both an
# OLS4 ancestor check and a cross-reference check showing no formal link) - this is a
# documented, reviewed exception list, not a general substitute for the live ontology
# data ontology_related() otherwise relies on.
_MANUALLY_MAPPED_EQUIVALENT_CONDITIONS = frozenset({
    frozenset({"MONDO:0007648", "MONDO:0100488"}),
})


def _manually_mapped_equivalent(case_condition, record_condition):
    return frozenset({case_condition, record_condition}) in _MANUALLY_MAPPED_EQUIVALENT_CONDITIONS


def ontology_related(case, case_condition, record, record_condition):
    """Whether two diseases are parent and child in MONDO, in either direction - or a
    manually reviewed equivalent pair, see _MANUALLY_MAPPED_EQUIVALENT_CONDITIONS above.

    Deliberately not equivalence. The same gene can lose function in one subtype and gain it
    in another, and subtypes can differ in inheritance mode, so this says only that the two
    terms are on one path - which is a reason to ask a curator, never a reason to decide.
    """
    if not case_condition or not record_condition or case_condition == record_condition:
        return False
    return (record_condition in condition_ancestors(case)
            or case_condition in condition_ancestors(record)
            or _manually_mapped_equivalent(case_condition, record_condition))


def condition_scope(record, key):
    """The phenotypes a curated disease was recorded as covering, or as keeping out.

    An absent scope is an empty set, which is not the same as a scope that says the set is
    empty: the provider that supplies this never reports an unreadable curation as an empty
    one, so nothing here has to tell the two apart.
    """
    scope = record.get("condition_scope")
    value = scope.get(key) if isinstance(scope, dict) else None
    return set(value) if isinstance(value, (list, set, tuple)) else set()
