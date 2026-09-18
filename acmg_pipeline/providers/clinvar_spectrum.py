"""Gene-wide ClinVar pathogenic/benign missense and truncating counts, for
GeneDiseaseDraftProvider's BP1 suggestion (a gene where disease is
predominantly caused by truncating rather than missense variation).

[Why this exists separately from providers/clinvar.py's own hotspot search]
  ClinVarHotspotProvider's gene-wide missense search (PM1) is bounded to a
  small window around one variant's position. BP1's question is unrelated
  to position: across the WHOLE gene, how many pathogenic missense
  variants are there compared with pathogenic truncating ones, and how
  many missense variants are outright benign? Two separate searches
  (missense; nonsense+frameshift), each counted for its own sake rather
  than filtered to a window.

[Bounded search or nothing, same convention as gene_missense_search()]
  Both searches are capped (retmax) and treated as usable only when
  ClinVar returned every hit it counted - a truncated list cannot support
  a real ratio, so an incomplete search on either side yields no record at
  all rather than a count that silently excludes an unknown number of
  variants.
"""

from __future__ import annotations

from acmg_pipeline.providers.clinvar import gene_consequence_search, summary_documents

METHOD = "clinvar_gene_classification_counts"
PATHOGENIC = ("pathogenic", "likely pathogenic")
BENIGN = ("benign", "likely benign")
MISSENSE_TERM = '"missense variant"[molecular consequence]'
TRUNCATING_TERM = '("nonsense"[molecular consequence] OR "frameshift variant"[molecular consequence])'
RETMAX = 5000
SUMMARY_BATCH = 200


def _matches_gene(document, gene):
    symbol = document.get("gene_sort") or ""
    symbols = {value.get("symbol") for value in document.get("genes", [])
              if isinstance(value, dict)}
    return symbol == gene or gene in symbols


def _bucket(document):
    classification = (document.get("germline_classification")
                      or document.get("clinical_significance") or {})
    description = (classification.get("description") or "").strip().lower()
    if "conflicting" in description:
        return "conflicting"
    if description in PATHOGENIC:
        return "pathogenic"
    if description in BENIGN:
        return "benign"
    return None


class ClinvarSpectrumProvider:
    name = "ClinVar gene classification counts"

    def __init__(self, client, release, *, retmax=RETMAX, summary_batch=SUMMARY_BATCH):
        if not release:
            raise ValueError("ClinVar release must be recorded")
        self.client = client
        self.release = release
        self.retmax = retmax
        self.summary_batch = summary_batch

    def _buckets(self, gene, consequence_term):
        """(buckets, retrieved_at), or (None, None) if the search was incomplete."""
        found = gene_consequence_search(self.client, self.release, gene,
                                        consequence_term, self.retmax)
        if not found["complete"]:
            return None, None
        if not found["ids"]:
            return [], found["retrieved_at"]
        documents, summarised_at = summary_documents(
            self.client, self.release, found["ids"], self.summary_batch)
        buckets = [_bucket(document) for document in documents.values()
                  if isinstance(document, dict) and _matches_gene(document, gene)]
        return buckets, max(summarised_at or found["retrieved_at"], found["retrieved_at"])

    def get_spectrum(self, gene):
        """The clinvar_spectrum record GeneDiseaseDraftProvider.build() takes, or None."""
        if not gene:
            return None
        missense, missense_at = self._buckets(gene, MISSENSE_TERM)
        if missense is None:
            return None
        truncating, truncating_at = self._buckets(gene, TRUNCATING_TERM)
        if truncating is None:
            return None
        retrieved_at = max(missense_at, truncating_at)
        return {
            "gene": gene,
            "source": self.name,
            "source_version": self.release,
            "retrieved_at": retrieved_at,
            "quality_status": "PASS",
            "evidence_id": f"urn:bh26:clinvar-spectrum:{gene}:{self.release}",
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": self.release,
            "complete": True,
            "pathogenic_missense_count": sum(1 for bucket in missense if bucket == "pathogenic"),
            "benign_missense_count": sum(1 for bucket in missense if bucket == "benign"),
            "pathogenic_truncating_count": sum(1 for bucket in truncating if bucket == "pathogenic"),
            "policy_note": (
                "Gene-wide ClinVar classification counts, split by missense vs "
                "nonsense/frameshift molecular consequence - not position-bounded, unlike "
                "PM1's hotspot density. Feeds GeneDiseaseDraftProvider's BP1 suggestion only; "
                "see acmg_pipeline.criteria.mechanism for how that suggestion is flagged."
            ),
        }


__all__ = ["METHOD", "ClinvarSpectrumProvider"]
