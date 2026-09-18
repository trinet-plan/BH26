"""MONDO ancestry, read only to tell a related disease apart from the same disease.

[Why this exists]
  Curations and cases name diseases at different depths. G2P curates
  "MYBPC3-related hypertrophic cardiomyopathy" (MONDO:0007268) while a case is as likely to
  be recorded as "hypertrophic cardiomyopathy" (MONDO:0005045), and the two identifiers are
  not equal. Without ancestry the gate reports that no mechanism exists for this disease,
  which sends a curator looking for evidence that is on file under the neighbouring term.

[Why an ancestry relation is not a match]
  It is deliberately not equivalence. The same gene can lose function in one subtype and gain
  it in another, and subtypes can differ in inheritance mode, so a mechanism curated for a
  parent or a child is not automatically a mechanism for this disease. ACMG routes this to a
  person, and so does PVS1: the relation turns "no mechanism" into "a related curation
  exists, decide whether it applies", which is a question a curator can answer and the
  criterion cannot.

[Why ancestors and not descendants]
  Asking for one term's descendants is unbounded - a broad term has hundreds - while its
  ancestors are a short path to the root. Both directions are covered by comparing ancestor
  sets: the case is a descendant of the record when the record's term is among the case's
  ancestors, and the record is a descendant of the case when the case's term is among the
  record's. Neither lookup grows with how broad the other term is.
"""

from __future__ import annotations

OLS_TERMS = "https://www.ebi.ac.uk/ols4/api/ontologies/mondo/terms"
METHOD = "mondo_hierarchical_ancestors_ols4"
# One term's ancestors are a path to the root, tens of entries at most. A term that somehow
# exceeds this is not silently truncated into a shorter ancestry - it resolves to nothing.
MAX_ANCESTORS = 500


def _quoted(mondo_id):
    """OLS identifies a term by its IRI, twice URL-encoded in the path."""
    return (f"http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F"
            f"{mondo_id.replace(':', '_')}")


class MondoHierarchyProvider:
    """Answers which MONDO terms a term sits under, and nothing else about it."""

    name = "MONDO hierarchy (OLS4)"

    def __init__(self, client, *, terms_url=OLS_TERMS):
        self.client = client
        self.terms_url = terms_url

    def ancestors(self, mondo_id):
        """Every MONDO term `mondo_id` is a descendant of, or None if that cannot be settled.

        Non-MONDO ancestors are dropped: BFO's upper-level classes sit above every disease and
        relating two diseases through them would make every pair of diseases related.
        """
        if not isinstance(mondo_id, str) or not mondo_id.startswith("MONDO:"):
            return None
        url = (f"{self.terms_url}/{_quoted(mondo_id)}/hierarchicalAncestors"
               f"?size={MAX_ANCESTORS}")
        response = self.client.fetch(url)
        body = response["body"]
        total = (body.get("page") or {}).get("totalElements")
        terms = ((body.get("_embedded") or {}).get("terms")) or []
        if not isinstance(total, int) or total > MAX_ANCESTORS or len(terms) < total:
            return None
        return {
            "term": mondo_id,
            "ancestors": sorted({item.get("obo_id") for item in terms
                                 if str(item.get("obo_id", "")).startswith("MONDO:")}),
            "source": self.name,
            "source_version": None,
            "retrieved_at": response["retrieved_at"],
            "assessment_method": "automated",
            "method": METHOD,
        }
