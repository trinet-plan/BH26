from acmg.core.models import CriterionResult, Status, Variant
from acmg.services.population import number, usable_observations


DEFAULT_STRENGTH = {
    "PVS1": "very_strong", "PS1": "strong", "PM1": "moderate", "PM2": "moderate",
    "PM4": "moderate", "PM5": "moderate", "PP2": "supporting", "PP3": "supporting",
    "BA1": "stand_alone", "BS1": "strong", "BP1": "supporting", "BP3": "supporting",
    "BP4": "supporting", "BP7": "supporting",
}


def result(code, input_data, status, summary, *, strength=None, evidence=None,
           missing=None, review=None, provenance=None):
    direction, outcome = None, None
    if status == Status.MET:
        direction = "disputes" if code.startswith("B") else "supports"
        outcome = code if strength == DEFAULT_STRENGTH.get(code) else f"{code}_{strength}"
    elif status == Status.NOT_MET:
        direction, outcome = "none", f"{code}_not_met"
    return CriterionResult(
        code, status, input_data["variant"], summary, strength, direction, outcome,
        evidence=evidence or [], missing_inputs=missing or [], review_points=review or [],
        provenance={"rule_version": f"{code}-v1", **(provenance or {})},
    )


def population_context(code, input_data, services, config):
    rule = config.get(code, {})
    minimum_an = number(rule.get("minimum_an"))
    if minimum_an is None or minimum_an < 1 or minimum_an != minimum_an.to_integral_value():
        return result(code, input_data, Status.NOT_EVALUATED, "Population quality policy missing",
                      missing=[f"{code}.minimum_an"]), None
    if not all(rule.get(key) for key in ("policy_source", "policy_version")):
        return result(code, input_data, Status.NOT_EVALUATED, "Population policy provenance missing",
                      missing=[f"{code}.policy_source", f"{code}.policy_version"]), None
    variant = Variant(**input_data["variant"])
    resolved = services.population.get_resolved_evidence(variant, input_data)
    valid, rejected = usable_observations(resolved, variant, minimum_an)
    provenance = {"policy_source": rule["policy_source"], "policy_version": rule["policy_version"],
                  "rejected_observations": rejected, "provider_failures": resolved["failures"]}
    if not valid:
        return result(code, input_data, Status.NOT_EVALUATED, "No reliable population observation",
                      evidence=resolved["observations"], missing=["population"],
                      provenance=provenance), None
    return None, (rule, valid, rejected, resolved["failures"], provenance)


def get_evidence(category, input_data, services):
    return services.evidence.get(category, Variant(**input_data["variant"]), input_data)


def annotation_context(code, input_data, services):
    annotations = get_evidence("annotation", input_data, services)
    if not annotations:
        return result(code, input_data, Status.NOT_EVALUATED, "Transcript annotation unavailable",
                      missing=["annotation", "transcript"]), None
    if len(annotations) != 1:
        return result(code, input_data, Status.MANUAL_REVIEW, "Multiple transcript annotations require resolution",
                      evidence=annotations, review=["Select disease-relevant transcript"]), None
    annotation = annotations[0]
    if not annotation.get("consequences") or not annotation.get("transcript"):
        return result(code, input_data, Status.NOT_EVALUATED, "Incomplete annotation",
                      evidence=annotations, missing=["consequences", "transcript"]), None
    return None, annotation


def curated_context(code, category, input_data, services, annotation, *, disease_required=True):
    """One reviewed assessment in the exact disease/transcript context, never a DB label."""
    if disease_required and not input_data.get("condition"):
        return result(code, input_data, Status.NOT_EVALUATED, "Disease context required",
                      evidence=[annotation], missing=["condition"]), None
    records = get_evidence(category, input_data, services)
    records = [r for r in records if r.get("curator") and r.get("reviewed_at")
               and r.get("transcript") == annotation["transcript"]
               and (not disease_required or r.get("condition") == input_data["condition"])]
    if not records:
        return result(code, input_data, Status.NOT_EVALUATED, f"Reviewed {category} evidence unavailable",
                      evidence=[annotation], missing=[category]), None
    if len(records) != 1:
        return result(code, input_data, Status.MANUAL_REVIEW, f"Multiple {category} assessments",
                      evidence=[annotation, *records], review=[f"Resolve {category} assessments"]), None
    return None, records[0]


def require_boolean_fields(code, input_data, evidence, assessment, fields):
    missing = [field for field in fields if type(assessment.get(field)) is not bool]
    if missing:
        return result(code, input_data, Status.NOT_EVALUATED, "Assessment incomplete",
                      evidence=evidence, missing=missing)
    return None
