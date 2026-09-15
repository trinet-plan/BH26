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
    if not all(annotation.get(field) for field in fields) or not input_data.get("condition"):
        missing = [field for field in fields if not annotation.get(field)]
        if not input_data.get("condition"):
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
        if candidate.get("condition") != input_data["condition"]:
            continue
        if candidate.get("transcript") != annotation["transcript"]:
            continue
        if candidate.get("classification") != "Pathogenic":
            continue
        if not all(candidate.get(field) is True for field in
                   ("pathogenic_evidence_reviewed", "independent_evidence", "mechanism_matches", "splice_effect_checked")):
            continue
        if candidate.get("different_splice_mechanism") is not False:
            continue
        if not all(candidate.get(field) for field in ("curator", "reviewed_at", "primary_evidence")):
            continue
        qualified.append(candidate)
    if qualified:
        return result(code, input_data, Status.MET, "Reviewed independent pathogenic comparator at same residue",
                      strength="strong" if code == "PS1" else "moderate", evidence=[annotation, *qualified])
    if matching:
        return result(code, input_data, Status.MANUAL_REVIEW, "Comparator evidence requires confirmation",
                      evidence=[annotation, *matching], review=["Verify comparator pathogenicity, independence and splice mechanism"])
    searches = get_evidence("comparator_search", input_data, services)
    complete = any(search.get("complete") is True and search.get("protein_id") == annotation["protein_id"]
                   and search.get("protein_start") == annotation["protein_start"]
                   and search.get("condition") == input_data["condition"] for search in searches)
    return result(code, input_data, Status.NOT_MET if complete else Status.NOT_EVALUATED,
                  "No eligible comparator found" if complete else "Comparator search incomplete",
                  evidence=[annotation, *searches], missing=[] if complete else ["complete_comparator_search"])
