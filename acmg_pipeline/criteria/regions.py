from dataclasses import replace

from acmg_pipeline.constants import CriterionStatus
"""Protein length / repeat / critical-region assessments with explicit missingness."""

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import (
    NOT_APPLICABLE, annotation_context, curated_context, require_boolean_fields, result,
)

# ACMG PM1 reads "mutational hot spot and/or critical and well-established functional
# domain, without benign variation", so the two routes are alternatives and only the
# absence of benign variation is shared. A hotspot is derived from observed variant
# density, a critical domain from independently established function; requiring both would
# demand the same evidence twice and reject established active sites with few reports.
PM1_ROUTES = ("mutational_hotspot", "critical_functional_domain")
HOTSPOT_POLICY_FIELDS = ("window_aa", "min_pathogenic", "max_benign", "method",
                         "policy_source", "policy_version")


def evaluate_region(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    consequences = set(annotation["consequences"])
    indel = bool(consequences & {"inframe_insertion", "inframe_deletion"})
    if code == "PM4" and not (indel or "stop_lost" in consequences):
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      "PM4 is not applicable: the annotation does not indicate an in-frame insertion, "
                      f"in-frame deletion, or stop-loss consequence (annotated: "
                      f"{', '.join(sorted(consequences))}).",
                      evidence=[annotation], provenance=NOT_APPLICABLE)
    if code == "BP3" and not indel:
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      "BP3 is not applicable: the annotation does not indicate an in-frame insertion "
                      f"or deletion, which BP3 requires before repeat-region assessment "
                      f"(annotated: {', '.join(sorted(consequences))}).",
                      evidence=[annotation], provenance=NOT_APPLICABLE)
    if code == "PM1":
        curated = _pm1_curated_critical_domain(input_data, annotation, config)
        if curated is not None:
            return curated
    # PM1 is a protein-level statement about the region itself, so a condition-agnostic
    # reviewed assessment is usable; the disease relevance is reported separately.
    early, region = curated_context(code, "region", input_data, services, annotation,
                                    disease_required=False)
    if early:
        return early
    evidence = [annotation, region]
    protein = annotation.get("protein_id")
    position = annotation.get("protein_start")
    end = annotation.get("protein_end", position)
    if not protein or region.get("protein_id") != protein:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Protein reference mapping unavailable",
                      evidence=evidence, missing=["protein_id"])
    start, stop = region.get("start"), region.get("end")
    # An absent coordinate and a reversed interval need different fixes - one is evidence to
    # retrieve, the other is evidence to correct - so they are not reported as one cause.
    labelled = {"protein_start": position, "protein_end": end,
                "region.start": start, "region.end": stop}
    unusable = [name for name, value in labelled.items()
                if type(value) is not int or value <= 0]
    if unusable:
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      f"Protein coordinates unavailable: {', '.join(unusable)} is not a positive "
                      f"residue position",
                      evidence=evidence, missing=["protein_interval"])
    if start > stop or position > end:
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      f"Protein coordinates are inconsistent: the assessed region is "
                      f"{start}-{stop} and the altered interval is {position}-{end}; each must "
                      f"start at or before it ends",
                      evidence=evidence, review=["Correct the protein interval bounds"])
    # A partial intersection is insufficient to classify the full altered region as non-functional.
    contained = start <= position <= end <= stop
    if not contained:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Assessment does not cover altered protein interval",
                      evidence=evidence, review=["Review complete affected interval"])
    if code == "PM1":
        return evaluate_pm1(input_data, region, evidence, config)
    fields = {
        "PM4": ["nonfunctional_repeat", "functional_review_complete"],
        "BP3": ["repetitive", "functional_importance", "functional_review_complete"],
    }[code]
    early = require_boolean_fields(code, input_data, evidence, region, fields)
    if early:
        return early
    if not region["functional_review_complete"]:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Functional region review incomplete",
                      evidence=evidence, review=["Confirm functional relevance of repeat/region"])
    if code == "PM4":
        delta = annotation.get("protein_length_change")
        if type(delta) is not int:
            return result(code, input_data, CriterionStatus.UNKNOWN, "Protein length change unavailable",
                          evidence=evidence, missing=["protein_length_change"])
        met = delta != 0 and not region["nonfunctional_repeat"]
    else:
        met = region["repetitive"] and not region["functional_importance"]
    strength = "supporting" if code == "BP3" else "moderate"
    return result(code, input_data, CriterionStatus.MET if met else CriterionStatus.NOT_MET,
                  "Altered protein interval and reviewed region evidence evaluated",
                  strength=strength if met else None, evidence=evidence)


def _pm1_curated_critical_domain(input_data, annotation, config):
    """A gene's real, VCEP-published critical-domain codon range, hand-transcribed into
    curated-context.json's gene_critical_domains (see automated_core/context.py) - None
    when this gene has no such table, or the altered interval falls outside every listed
    range, so evaluate_region() falls through to the density-based hotspot route unchanged.

    [Why this bypasses curated_context()'s "region" category entirely, 2026-09-22]
      pm1_critical_domain() below has always refused assessment_method="automated" region
      evidence for this exact route ("Automated evidence cannot establish functional
      criticality") - a deliberate protection against inferring domain importance from
      variant density, which this data is not: it is transcribed directly from a gene's own
      published ClinGen VCEP specification (e.g. MECP2's real spec names "Methyl-DNA binding
      (MBD): aa 90-162" verbatim), the same manual_transcription-with-citation convention
      already used for curated-context.json's ba1_exceptions and gene_frequency_thresholds.
      Routing it through the shared "region" evidence category instead would collide with
      ClinVarHotspotProvider's own hotspot search for the same missense variant (both would
      answer curated_context()'s single-record lookup, tripping "Multiple region
      assessments") - a synthesized region dict, evaluated directly, avoids that collision
      and the automated-only gate is simply not reached for genes covered here.

    [benign_depletion defaulted True - a documented assumption, not a citation]
      Real VCEP specifications name the critical residues/ranges but do not always restate
      "and no benign variation exists there" in the same sentence - that is the premise a
      VCEP publishing the rule as Moderate/Strong is understood to have already weighed, not
      a separate number this project's own tooling can re-derive. Flagged in review_points,
      same as every other flagged default in this project (splice_default.py's
      alternative_rescue, region_repeat.py's repetitive=False).
    """
    domains = input_data.get("gene_critical_domains")
    if not domains or not domains.get("ranges"):
        return None
    protein = annotation.get("protein_id")
    position = annotation.get("protein_start")
    end = annotation.get("protein_end", position)
    if not protein or type(position) is not int or position <= 0 or type(end) is not int or end < position:
        return None
    matched = next(
        (item for item in domains["ranges"] if item["start"] <= position <= end <= item["end"]),
        None,
    )
    if matched is None:
        return None
    region = {
        "category": "region", "variant_key": Variant(**input_data["variant"]).key,
        "evidence_id": (f"urn:bh26:pm1-critical-domain:{domains.get('source_version', 'v1')}:"
                        f"{protein}:{matched['start']}-{matched['end']}"),
        "source": domains["source"], "source_version": domains["source_version"],
        "retrieved_at": domains["reviewed_at"], "quality_status": "PASS",
        "curator": "BH26 project (manual transcription from ClinGen CSpec)",
        "reviewed_at": domains["reviewed_at"],
        "protein_id": protein, "start": matched["start"], "end": matched["end"],
        "region_type": "critical_functional_domain",
        "critical_functional_region": True,
        "benign_depletion": True,
        "label": matched.get("label"),
        "policy_note": (
            "Codon range transcribed from a real ClinGen VCEP CSpec specification "
            f"({domains['source']}), not inferred from variant density - see this "
            "function's own docstring for why benign_depletion is a flagged default "
            "rather than a value the specification itself states."
        ),
    }
    evidence = [annotation, region]
    outcome = evaluate_pm1(input_data, region, evidence, config)
    if outcome.status == CriterionStatus.MET:
        outcome = replace(outcome, review_points=[
            *outcome.review_points,
            "benign_depletion is defaulted True for this gene's curated critical domain, "
            "not itself stated by the VCEP specification - confirm no benign variation is "
            "reported in this interval before relying on it",
        ])
    return outcome


def evaluate_pm1(input_data, region, evidence, config):
    route = region.get("region_type")
    if route not in PM1_ROUTES:
        return result("PM1", input_data, CriterionStatus.UNKNOWN, "Region evidence declares no PM1 route",
                      evidence=evidence, missing=["region_type"])
    if route == "mutational_hotspot":
        met, early, extra = pm1_hotspot(input_data, region, evidence, config)
    else:
        met, early, extra = pm1_critical_domain(input_data, region, evidence)
    if early:
        return early
    condition = input_data.get("condition")
    matched = bool(condition) and region.get("condition") == condition
    provenance = {"assessment_scope": "protein_level",
                  "condition_assessment": "MATCHED" if matched else "NOT_EVALUATED",
                  "assessment_method": region.get("assessment_method", "curated"), **extra}
    review = [] if matched or not met else [
        "Confirm region criticality for the disease context before final classification"]
    return result("PM1", input_data, CriterionStatus.MET if met else CriterionStatus.NOT_MET,
                  f"Reviewed {route} evidence covering the altered protein interval",
                  strength="moderate" if met else None, evidence=evidence,
                  review=review, provenance=provenance)


def pm1_critical_domain(input_data, region, evidence):
    """Criticality must be established independently, never inferred from variant density."""
    if region.get("assessment_method") == "automated":
        return False, result("PM1", input_data, CriterionStatus.UNKNOWN,
                             "Automated evidence cannot establish functional criticality",
                             evidence=evidence, missing=["curated_criticality"]), None
    fields = ("critical_functional_region", "benign_depletion")
    early = require_boolean_fields("PM1", input_data, evidence, region, fields)
    if early:
        return False, early, None
    return all(region[field] for field in fields), None, {"pm1_route": "critical_functional_domain"}


def pm1_hotspot(input_data, region, evidence, config):
    """Density-derived hotspot: thresholds come from a versioned policy, never from code."""
    policy = config.get("PM1", {}).get("hotspot", {})
    missing_policy = [field for field in HOTSPOT_POLICY_FIELDS if policy.get(field) is None or policy[field] == ""]
    if missing_policy:
        return False, result("PM1", input_data, CriterionStatus.UNKNOWN, "Hotspot density policy unavailable",
                             evidence=evidence,
                             missing=[f"PM1.hotspot.{field}" for field in missing_policy]), None
    if region.get("policy_version") != policy["policy_version"] or region.get("method") != policy["method"]:
        return False, result("PM1", input_data, CriterionStatus.UNKNOWN,
                             "Hotspot evidence policy does not match configuration",
                             evidence=evidence, missing=["policy_version", "method"]), None
    counts = {field: region.get(field) for field in ("pathogenic_count", "benign_count")}
    if any(type(value) is not int or value < 0 for value in counts.values()):
        return False, result("PM1", input_data, CriterionStatus.UNKNOWN, "Hotspot variant counts unavailable",
                             evidence=evidence, missing=sorted(counts)), None
    provenance = {"pm1_route": "mutational_hotspot", "hotspot_policy": policy, **counts}
    depleted = counts["benign_count"] <= policy["max_benign"]
    enriched = counts["pathogenic_count"] >= policy["min_pathogenic"]
    # Benign variation inside the window argues against PM1, so that is a real negative.
    # Too few pathogenic reports is missing evidence, not proof that no hotspot exists:
    # well-studied genes reach the threshold far more easily than rarely reported ones.
    if depleted and not enriched:
        return False, result("PM1", input_data, CriterionStatus.UNKNOWN,
                             "Reported pathogenic density is below the hotspot policy threshold",
                             evidence=evidence, missing=["pathogenic_density"],
                             provenance=provenance), None
    met = enriched and depleted
    derived = {"mutational_hotspot": met, "pathogenic_enrichment": enriched,
               "benign_depletion": depleted}
    conflicting = sorted(field for field, value in derived.items()
                         if type(region.get(field)) is bool and region[field] is not value)
    if conflicting:
        return False, result("PM1", input_data, CriterionStatus.UNKNOWN,
                             "Hotspot assertion conflicts with the counts it is derived from",
                             evidence=evidence,
                             review=[f"Reconcile {field}" for field in conflicting],
                             provenance=provenance), None
    return met, None, provenance
