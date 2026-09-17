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
           decision_trace=None, rules_used=None, warnings=None, unresolved_requirements=None):
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
        unresolved_requirements=unresolved_requirements or [],
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


def population_context(code, input_data, services, config):
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
        return result(code, input_data, CriterionStatus.UNKNOWN, "No reliable population observation",
                      evidence=resolved["observations"], missing=["population"],
                      provenance=provenance), None
    return None, (rule, valid, rejected, resolved["failures"], provenance)


def get_evidence(category, input_data, services):
    return services.evidence.get(category, Variant(**input_data["variant"]), input_data)


def retrieval_failures(services, provider=None):
    """What the resolver could not retrieve, so a criterion can say why it has nothing.

    Without this, a provider that errored and a provider that legitimately returned nothing
    are indistinguishable downstream, and the curator is told only that evidence is absent -
    which is the one thing they could already see.
    """
    failures = getattr(services, "failures", None) or []
    if provider is None:
        return list(failures)
    return [item for item in failures if item.get("provider") == provider]


def annotation_context(code, input_data, services):
    annotations = get_evidence("annotation", input_data, services)
    if not annotations:
        failures = retrieval_failures(services)
        detail = "; ".join(f"{item.get('provider')}: {item.get('error')}" for item in failures)
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


def unusable_reason(category, records, annotation, condition, disease_required):
    """Why every retrieved `category` record was rejected - one cause per rejection route.

    "No usable assessment" covers four situations a curator has to act on differently:
    nothing was retrieved at all (produce the assessment), records exist but name no
    reviewer or policy version (record the provenance), records were reviewed against a
    different transcript (resolve the transcript), or against a different disease (resolve
    the disease context).  Reporting them as one message would tell a curator to create
    evidence that already exists.  Returns (summary, missing, review).
    """
    if not records:
        return (f"No {category} assessment was retrieved for this variant", [category], [])
    unreviewed = [r for r in records if not reviewed_or_automated(r)]
    reviewed = [r for r in records if reviewed_or_automated(r)]
    transcript_mismatch = [r for r in reviewed if r.get("transcript") != annotation["transcript"]]
    condition_mismatch = [r for r in reviewed if r.get("transcript") == annotation["transcript"]
                          and disease_required and r.get("condition") != condition]
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


def curated_context(code, category, input_data, services, annotation, *, disease_required=True):
    """One reviewed assessment in the exact disease/transcript context, never a DB label."""
    if disease_required and not input_data.get("condition"):
        return result(code, input_data, CriterionStatus.UNKNOWN, "Disease context required",
                      evidence=[annotation], missing=["condition"]), None
    condition = input_data.get("condition")
    retrieved = get_evidence(category, input_data, services)
    records = [r for r in retrieved if reviewed_or_automated(r)
               and r.get("transcript") == annotation["transcript"]
               and (not disease_required or r.get("condition") == condition)]
    if not records:
        summary, missing, review = unusable_reason(category, retrieved, annotation, condition,
                                                   disease_required)
        return result(code, input_data, CriterionStatus.UNKNOWN, summary,
                      evidence=[annotation, *retrieved], missing=missing, review=review), None
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
