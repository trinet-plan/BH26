from acmg_pipeline.constants import CriterionStatus
"""Same-residue comparison requires independent, reviewed pathogenic evidence."""

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import (
    CATEGORY_PROVIDERS, NOT_APPLICABLE, annotation_context, failure_detail, get_evidence,
    result, retrieval_failures,
)


def evaluate_comparator(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    if "missense_variant" not in annotation["consequences"]:
        # A splice consequence is not outside the criterion the way a synonymous one is: PS1's
        # same-effect argument can extend to an equivalent splice effect, which a curator has
        # to assess.  So only the non-splice case is reported as inapplicable.
        splice = any("splice" in item for item in annotation["consequences"])
        listed = ", ".join(sorted(annotation["consequences"]))
        if splice:
            return result(code, input_data, CriterionStatus.UNKNOWN,
                          f"{code} compares amino-acid substitutions, and the annotation reports a "
                          f"splice consequence ({listed}); equivalence of the splice effect has to "
                          f"be assessed before the comparison can be made",
                          evidence=[annotation], review=["Assess splice-effect equivalence"])
        return result(code, input_data, CriterionStatus.UNKNOWN,
                      f"{code} is not applicable: it compares missense substitutions, and the "
                      f"annotation reports {listed}.",
                      evidence=[annotation], provenance=NOT_APPLICABLE)
    # Both codes compare amino acid changes at one residue, so they are protein-level
    # statements; the disease relevance is reported separately instead of being required.
    fields = ("protein_id", "protein_start", "ref_aa", "alt_aa")
    if not all(annotation.get(field) for field in fields):
        return result(code, input_data, CriterionStatus.UNKNOWN, "Protein context missing",
                      evidence=[annotation],
                      missing=[field for field in fields if not annotation.get(field)])
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
        # PS1 needs the identical protein change; PM5 needs a different change at the same
        # residue, which the provider confirms against the same protein reference.
        scoped = candidate.get("exact_protein_match" if code == "PS1" else "residue_match") is True
        automated = (
            scoped
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
                    code, input_data, CriterionStatus.UNKNOWN,
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
            code, input_data, CriterionStatus.MET,
            "Pathogenic ClinVar comparator produces the same amino acid substitution",
            strength="strong" if code == "PS1" else "moderate",
            evidence=[annotation, *qualified], review=review,
            provenance={"assessment_scope": "protein_level",
                        "condition_assessment": condition_status},
        )
    # A residue search surfaces every reported change at the residue. Only a pathogenic
    # comparator can carry the criterion, so a VUS or benign one is recorded as evidence of
    # what was searched, not raised as something a curator must resolve.
    pathogenic_matches = [candidate for candidate in matching
                          if candidate.get("classification") in ("Pathogenic", "Likely pathogenic")]
    if pathogenic_matches:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Comparator evidence requires confirmation",
                      evidence=[annotation, *pathogenic_matches],
                      review=["Verify comparator pathogenicity, independence and splice mechanism"])
    searches = get_evidence("comparator_search", input_data, services)
    # An exact protein-change search cannot show that no *other* change at the residue is
    # pathogenic, so PM5 only accepts a residue-scoped search.
    complete = any(
        search.get("complete") is True
        and search.get("protein_id") == annotation["protein_id"]
        and search.get("protein_start") == annotation["protein_start"]
        and (search.get("search_scope", "exact_protein_change") == "residue"
             or (code == "PS1" and search.get("search_scope", "exact_protein_change")
                 == "exact_protein_change"))
        for search in searches
    )
    if complete:
        return result(code, input_data, CriterionStatus.NOT_MET, "No eligible comparator found",
                      evidence=[annotation, *searches, *matching], missing=[])
    # An incomplete search has two causes a curator acts on differently: the search ran and
    # covered the wrong scope, or it never ran because the provider failed.
    # The search is absorbed per criterion, so PS1 reports PS1's failure and PM5 reports PM5's.
    failures = retrieval_failures(
        services,
        [name for base in CATEGORY_PROVIDERS["comparator_search"] for name in (base, f"{base}:{code}")],
        exact=True)
    detail = failure_detail(failures)
    return result(code, input_data, CriterionStatus.UNKNOWN,
                  f"Comparator search incomplete - {detail}" if detail
                  else "Comparator search incomplete",
                  evidence=[annotation, *searches, *matching],
                  missing=["complete_comparator_search"],
                  provenance={"provider_failures": failures} if failures else None)
