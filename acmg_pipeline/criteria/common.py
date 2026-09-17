from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.models import CriterionResult, Variant
from acmg_pipeline.services.population import number, usable_observations


DEFAULT_STRENGTH = {
    "PVS1": "very_strong", "PS1": "strong", "PM1": "moderate", "PM2": "moderate",
    "PM4": "moderate", "PM5": "moderate", "PP2": "supporting", "PP3": "supporting",
    "BA1": "stand_alone", "BS1": "strong", "BP1": "supporting", "BP3": "supporting",
    "BP4": "supporting", "BP7": "supporting",
}


# Human-readable scope is included in every automated-criterion summary.  The
# structured fields (status, missing_inputs, review_points and decision_trace)
# remain the machine-readable source of truth; this text makes an exported
# EvidenceLine understandable without requiring a curator to reconstruct the
# ACMG definition from the criterion code alone.
CRITERION_SCOPE = {
    "PVS1": "PVS1 evaluates predicted loss-of-function variants only when loss of function is an established disease mechanism for the gene.",
    "PS1": "PS1 compares a missense variant with an independently established pathogenic variant producing the same amino-acid substitution.",
    "PM1": "PM1 evaluates a variant in a reviewed mutational hotspot or critical functional domain with no conflicting benign variation.",
    "PM2": "PM2 evaluates whether reliable population observations satisfy the configured rarity threshold.",
    "PM4": "PM4 evaluates protein-length changes caused by an in-frame insertion, in-frame deletion, or stop-loss variant.",
    "PM5": "PM5 compares a missense variant with an independently established pathogenic variant causing a different amino-acid substitution at the same residue.",
    "PP2": "PP2 evaluates missense variants in genes where missense is an established disease mechanism and benign missense variation is constrained.",
    "PP3": "PP3 evaluates calibrated computational evidence supporting a damaging effect.",
    "PP5": "PP5 is not scored from external assertions alone; primary evidence must be reviewed instead.",
    "BA1": "BA1 evaluates whether population allele frequency exceeds the stand-alone benign threshold after exception review.",
    "BS1": "BS1 evaluates whether population allele frequency exceeds the reviewed maximum credible frequency for the disease context.",
    "BP1": "BP1 evaluates missense variants in genes where disease is predominantly caused by truncating variants rather than missense variation.",
    "BP3": "BP3 evaluates in-frame insertions or deletions in a repetitive region without known function.",
    "BP4": "BP4 evaluates calibrated computational evidence supporting a benign effect.",
    "BP6": "BP6 is not scored from an external assertion alone; primary evidence must be reviewed instead.",
    "BP7": "BP7 evaluates synonymous variants outside splice-critical regions when splice prediction and conservation evidence support a benign interpretation.",
}


def _summary_with_context(code, status, summary):
    """Add a stable criterion definition and an explicit three-state outcome."""
    text = summary.rstrip(".") + "."
    scope = CRITERION_SCOPE.get(code, "")
    outcome = {
        CriterionStatus.MET: "Outcome: MET; the available evidence satisfies this criterion.",
        CriterionStatus.NOT_MET: "Outcome: NOT_MET; the criterion was evaluated but its requirements were not satisfied.",
        CriterionStatus.UNKNOWN: "Outcome: UNKNOWN; no MET or NOT_MET judgment is made from the available information.",
    }[status]
    return " ".join(part for part in (text, scope, outcome) if part)


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
        code, status, input_data["variant"], _summary_with_context(code, status, summary), strength, direction, outcome,
        evidence=evidence or [], missing_inputs=missing or [], review_points=review or [],
        provenance={"rule_version": f"{code}-v1", **(provenance or {})},
        evaluation_context=evaluation_context, decision_trace=decision_trace or [],
        rules_used=rules_used or [], warnings=warnings or [],
        unresolved_requirements=unresolved_requirements or [],
    )


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


def annotation_context(code, input_data, services):
    annotations = get_evidence("annotation", input_data, services)
    if not annotations:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Transcript annotation unavailable",
                      missing=["annotation", "transcript"]), None
    if len(annotations) != 1:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Multiple transcript annotations require resolution",
                      evidence=annotations, review=["Select disease-relevant transcript"]), None
    annotation = annotations[0]
    if not annotation.get("consequences") or not annotation.get("transcript"):
        return result(code, input_data, CriterionStatus.UNKNOWN, "Incomplete annotation",
                      evidence=annotations, missing=["consequences", "transcript"]), None
    return None, annotation


def reviewed_or_automated(record):
    """Human curation, or an automated assessment that names its versioned policy."""
    if record.get("assessment_method") == "automated":
        return bool(record.get("method") and record.get("policy_version"))
    return bool(record.get("curator") and record.get("reviewed_at"))


def curated_context(code, category, input_data, services, annotation, *, disease_required=True):
    """One reviewed assessment in the exact disease/transcript context, never a DB label."""
    if disease_required and not input_data.get("condition"):
        return result(code, input_data, CriterionStatus.UNKNOWN, "Disease context required",
                      evidence=[annotation], missing=["condition"]), None
    records = get_evidence(category, input_data, services)
    records = [r for r in records if reviewed_or_automated(r)
               and r.get("transcript") == annotation["transcript"]
               and (not disease_required or r.get("condition") == input_data["condition"])]
    if not records:
        return result(code, input_data, CriterionStatus.UNKNOWN, f"Reviewed {category} evidence unavailable",
                      evidence=[annotation], missing=[category]), None
    if len(records) != 1:
        return result(code, input_data, CriterionStatus.UNKNOWN, f"Multiple {category} assessments",
                      evidence=[annotation, *records], review=[f"Resolve {category} assessments"]), None
    return None, records[0]


def require_boolean_fields(code, input_data, evidence, assessment, fields):
    missing = [field for field in fields if type(assessment.get(field)) is not bool]
    if missing:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Assessment incomplete",
                      evidence=evidence, missing=missing)
    return None
