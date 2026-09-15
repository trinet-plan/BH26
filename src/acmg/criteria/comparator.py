"""Same-residue comparison requires independent, reviewed pathogenic evidence."""

from acmg.core.models import Status, Variant
from acmg.criteria.common import annotation_context, get_evidence, result


def evaluate_comparator(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    if "missense_variant" not in annotation["consequences"]:
        splice = any("splice" in item for item in annotation["consequences"])
        return result(code, input_data, Status.MANUAL_REVIEW if code == "PS1" and splice else Status.NOT_APPLICABLE,
                      "Splice-equivalence requires review" if splice else "Requires missense substitution",
                      evidence=[annotation], review=["Assess splice-effect equivalence"] if splice else [])
    fields = ("protein_id", "protein_start", "ref_aa", "alt_aa")
    disease_required = code == "PM5"
    if not all(annotation.get(field) for field in fields) or (
            disease_required and not input_data.get("condition")):
        missing = [field for field in fields if not annotation.get(field)]
        if disease_required and not input_data.get("condition"):
            missing.append("condition")
        return result(code, input_data, Status.NOT_EVALUATED, "Protein/disease context missing",
                      evidence=[annotation], missing=missing)
    candidates = get_evidence("comparator", input_data, services)
    matching, qualified = [], []
    for candidate in candidates:
        if any(candidate.get(field) != annotation[field] for field in ("protein_id", "protein_start", "ref_aa")):
            continue
        if not candidate.get("alt_aa") or candidate["alt_aa"] == annotation["ref_aa"]:
            continue
        same_aa = candidate["alt_aa"] == annotation["alt_aa"]
        if (code == "PS1") != same_aa:
            continue
        matching.append(candidate)
        try:
            other = Variant(**candidate["comparator_variant"])
        except (ValueError, TypeError, KeyError):
            continue
        if other == Variant(**input_data["variant"]):
            continue
        if candidate.get("transcript") != annotation["transcript"]:
            continue
        if candidate.get("classification") != "Pathogenic":
            continue
        automated = (
            code == "PS1"
            and candidate.get("exact_protein_match") is True
            and candidate.get("different_nucleotide_variant") is True
            and candidate.get("review_status_eligible") is True
            and candidate.get("splice_effect_checked") is True
            and candidate.get("splice_conflict") is False
        )
        reviewed = (
            candidate.get("condition") == input_data.get("condition")
            and all(candidate.get(field) is True for field in
                    ("pathogenic_evidence_reviewed", "independent_evidence",
                     "mechanism_matches", "splice_effect_checked"))
            and candidate.get("different_splice_mechanism") is False
            and all(candidate.get(field) for field in
                    ("curator", "reviewed_at", "primary_evidence"))
        )
        if not (automated or reviewed):
            continue
        qualified.append(candidate)
    if qualified:
        condition = input_data.get("condition")
        if condition:
            disease_matches = [condition in candidate.get("conditions", [])
                               or candidate.get("condition") == condition for candidate in qualified]
            condition_status = "MATCHED" if any(disease_matches) else "MISMATCHED"
            if not any(disease_matches):
                return result(
                    code, input_data, Status.MANUAL_REVIEW,
                    "Protein-level comparator found but disease context differs",
                    evidence=[annotation, *qualified],
                    review=["Review comparator disease relevance"],
                    provenance={"assessment_scope": "protein_level",
                                "condition_assessment": condition_status},
                )
        else:
            condition_status = "NOT_EVALUATED"
        review = [] if condition_status == "MATCHED" else [
            "Confirm disease relevance before final classification",
            "Review comparator classification for circular PS1 use",
        ]
        return result(
            code, input_data, Status.MET,
            "Pathogenic ClinVar comparator produces the same amino acid substitution",
            strength="strong" if code == "PS1" else "moderate",
            evidence=[annotation, *qualified], review=review,
            provenance={"assessment_scope": "protein_level",
                        "condition_assessment": condition_status},
        )
    if matching:
        return result(code, input_data, Status.MANUAL_REVIEW, "Comparator evidence requires confirmation",
                      evidence=[annotation, *matching], review=["Verify comparator pathogenicity, independence and splice mechanism"])
    searches = get_evidence("comparator_search", input_data, services)
    complete = any(
        search.get("complete") is True
        and search.get("protein_id") == annotation["protein_id"]
        and search.get("protein_start") == annotation["protein_start"]
        and (code == "PS1" or search.get("condition") == input_data.get("condition"))
        for search in searches
    )
    return result(code, input_data, Status.NOT_MET if complete else Status.NOT_EVALUATED,
                  "No eligible comparator found" if complete else "Comparator search incomplete",
                  evidence=[annotation, *searches], missing=[] if complete else ["complete_comparator_search"])
