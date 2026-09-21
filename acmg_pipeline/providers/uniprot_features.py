"""Shared UniProt feature-fetch and interval-overlap helpers.

Two providers answer "does UniProt record something over this protein interval?"
against the same UniProt entry: protein_region.py (PVS1's NF04/NF06, over the
LOST interval [protein_start, total_protein_length]) and region_repeat.py
(PM4/BP3's repeat-region question, over the ALTERED interval
[protein_start, protein_end]). The fetch and overlap logic is identical between
them - only the interval and the feature types each one reads differ - so it is
factored once here rather than duplicated per criterion.
"""

from __future__ import annotations

_UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
# "Chain" is UniProt's whole-protein-span feature - callers exclude it on purpose,
# since it always overlaps any interval and would make an overlap check meaningless.
CRITICAL_FEATURE_TYPES = {"Domain", "Region", "Binding site", "Active site", "Motif", "Coiled coil"}


def is_disordered(feature: dict) -> bool:
    return "disordered" in str(feature.get("description", "")).lower()


def feature_interval(feature: dict) -> tuple[int, int] | None:
    """(start, end) as a real integer pair, or None - UniProt records an unresolved or
    fuzzy boundary (e.g. "?" or an unknown terminus) with a null `value`, which parses
    without a KeyError/TypeError and must be rejected explicitly rather than compared."""
    try:
        start = feature["location"]["start"]["value"]
        end = feature["location"]["end"]["value"]
    except (KeyError, TypeError):
        return None
    if type(start) is not int or type(end) is not int:
        return None
    return start, end


def overlaps(feature_start: int, feature_end: int, query_start: int, query_end: int) -> bool:
    return feature_start <= query_end and query_start <= feature_end


def fetch_features(client, accession: str):
    """(features, response), or (None, None) if UniProt could not be resolved."""
    try:
        response = client.fetch(_UNIPROT_ENTRY_URL.format(accession=accession), response_format="json")
    except ValueError:
        return None, None
    features = (response["body"] or {}).get("features")
    if not isinstance(features, list):
        return None, None
    return features, response


__all__ = [
    "CRITICAL_FEATURE_TYPES", "is_disordered", "feature_interval", "overlaps", "fetch_features",
]
