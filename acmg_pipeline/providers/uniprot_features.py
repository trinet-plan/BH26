"""Shared UniProt feature-fetch, sequence-complexity, and isoform helpers.

Two providers answer "does UniProt record something over this protein interval?"
against the same UniProt entry: protein_region.py (PVS1's NF04/NF06, over the
LOST interval [protein_start, total_protein_length]) and region_repeat.py
(PM4/BP3's repeat-region question, over the ALTERED interval
[protein_start, protein_end]). The fetch and overlap logic is identical between
them - only the interval and the feature types each one reads differ - so it is
factored once here rather than duplicated per criterion.

sequence_complexity() and resolve_covering_entry() exist for region_repeat.py
alone (protein_region.py never needs a sequence, only features), but live here
because they operate on the same UniProt entry shape fetch_entry() returns.
"""

from __future__ import annotations

import math
from collections import Counter

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


def fetch_entry(client, accession: str):
    """(body, response), or (None, None) if UniProt could not be resolved.

    Keeps the whole parsed entry - sequence and comments as well as features -
    for callers (resolve_covering_entry(), sequence_complexity()'s callers) that
    need more than the feature list fetch_features() returns.
    """
    try:
        response = client.fetch(_UNIPROT_ENTRY_URL.format(accession=accession), response_format="json")
    except ValueError:
        return None, None
    body = response["body"]
    if not isinstance(body, dict):
        return None, None
    return body, response


def fetch_features(client, accession: str):
    """(features, response), or (None, None) if UniProt could not be resolved."""
    body, response = fetch_entry(client, accession)
    if body is None:
        return None, None
    features = body.get("features")
    if not isinstance(features, list):
        return None, None
    return features, response


def isoform_accessions(canonical_body: dict) -> list[str]:
    """This entry's non-canonical isoform accessions, in UniProt's own listed order.

    Read from the ALTERNATIVE PRODUCTS comment - the canonical ("Displayed") isoform
    is excluded since the caller already has it (it is the entry just fetched).
    """
    for comment in canonical_body.get("comments", []) or []:
        if comment.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        return [
            iso["isoformIds"][0]
            for iso in comment.get("isoforms", []) or []
            if iso.get("isoformSequenceStatus") != "Displayed" and iso.get("isoformIds")
        ]
    return []


def resolve_covering_entry(client, accession: str, min_length: int):
    """(body, response, accession) for the first UniProt entry whose own sequence
    is at least `min_length` residues long - the canonical entry, or (only when
    that is too short) the first of its isoforms that is - or (None, None, None)
    when nothing covers it.

    Why this is needed: a gene's canonical UniProt entry is one specific isoform,
    and a transcript's own protein numbering can run past its length when the
    evaluated transcript corresponds to a longer alternatively-spliced isoform
    instead (found on RPGR: NM_001034853.2/RPGR-ORF15 reaches residue ~1040, but
    the canonical entry Q92834 is only 1020 residues; UniProt's own isoform 6,
    accession Q92834-6, is named "ORF15" and is 1152 residues - its sequence is
    identical to the canonical entry's up to residue 584, which is what makes
    reading the ORF15 transcript's numbering against it - rather than at a guess
    - reasonable). Untested past that one confirmed case, so this is a first pass
    like the rest of this module: it picks the first isoform long enough, not
    necessarily the biologically correct one, and a caller should record which
    accession was actually used.
    """
    body, response = fetch_entry(client, accession)
    if body is None:
        return None, None, None
    length = (body.get("sequence") or {}).get("length")
    # An unreadable length is not evidence the entry is too short - only a length we
    # can actually compare tells us that, so escalation is attempted only then.
    if not isinstance(length, int) or length >= min_length:
        return body, response, accession
    for candidate in isoform_accessions(body):
        iso_body, iso_response = fetch_entry(client, candidate)
        if iso_body is None:
            continue
        iso_length = (iso_body.get("sequence") or {}).get("length")
        if isinstance(iso_length, int) and iso_length >= min_length:
            return iso_body, iso_response, candidate
    return None, None, None


def sequence_complexity(window: str) -> float | None:
    """Wootton-Federhen (SEG) local complexity of a residue window, in bits.

    K = (1/L) * log2( L! / prod(n_i!) ) over the window's own amino-acid counts
    n_i - the measure SEG (Wootton & Federhen, Methods Enzymol. 1993;266:554-71;
    parameters corroborated in Wootton, Comput Chem 1994;18(3):269-85) scans a
    sequence with to find low-complexity regions. Low K means a small number of
    residue types dominate the window (a homopolymeric or near-homopolymeric
    run); high K means the window looks like a typical, diverse stretch of
    protein. None for an empty window - there is nothing to measure.

    This computes K directly over one caller-chosen window rather than running
    SEG's own trigger/extend/merge scan across a whole sequence: region_repeat.py
    only ever asks about one specific altered interval, not "where are all the
    low-complexity segments in this protein", so the scanning half of the
    original algorithm has nothing to add here - only the complexity measure and
    its published default trigger threshold (K1=2.2 bits, see region_repeat.py)
    are reused.
    """
    if not window:
        return None
    length = len(window)
    counts = Counter(window)
    log_length_factorial = math.lgamma(length + 1)
    log_denominator = sum(math.lgamma(n + 1) for n in counts.values())
    return (log_length_factorial - log_denominator) / math.log(2) / length


__all__ = [
    "CRITICAL_FEATURE_TYPES", "is_disordered", "feature_interval", "overlaps",
    "fetch_entry", "fetch_features", "isoform_accessions", "resolve_covering_entry",
    "sequence_complexity",
]
