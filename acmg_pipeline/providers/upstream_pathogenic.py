"""Reported pathogenic variants upstream of a restart codon, for PVS1's IC03 gate.

[What it decides]
  When an initiation codon is lost and translation could restart at a downstream in-frame
  ATG, PVS1 asks whether anything pathogenic has been reported in the stretch that would be
  skipped. A reported pathogenic variant there says the lost N-terminus matters, and the
  criterion is applied at moderate rather than supporting.

[Why a negative is harder than a positive, and is only sometimes given]
  ClinVar's summary records position a protein change only when it is written as a plain
  single-residue substitution - R143W and the like. A nonsense or frameshift change carries
  no such position, and those are exactly the variants likely to sit near the N-terminus.

  So the two answers are not symmetric and are not treated as such. One placed pathogenic
  substitution upstream of the restart codon settles the question as yes. Saying no requires
  that every pathogenic variant the gene search returned could be placed; while any remain
  unplaced, one of them could be the upstream evidence, and no record is emitted rather than
  a no that the data does not support.

[Bounded search or nothing]
  The gene search is emitted as evidence only when ClinVar returned every hit it counted. A
  truncated list cannot support either answer, so it produces no record.
"""

from __future__ import annotations

import hashlib
from urllib.parse import urlencode

from acmg_pipeline.providers.clinvar import (
    SUMMARY_BATCH, canonical_json, protein_change_positions, summary_documents,
)

METHOD = "clinvar_upstream_pathogenic_report"
PATHOGENIC_TERM = ('({gene}[gene]) AND ("pathogenic"[Clinical significance] '
                   'OR "likely pathogenic"[Clinical significance])')
RETMAX = 5000


def gene_pathogenic_search(client, release, gene, retmax=RETMAX):
    """Every pathogenic and likely-pathogenic record ClinVar holds for `gene`.

    Not restricted to missense, unlike the PM1/PM5 searches: the question is whether
    anything pathogenic is reported upstream, and a nonsense variant answers it as well as a
    substitution does.
    """
    term = PATHOGENIC_TERM.format(gene=gene)
    query = urlencode({"db": "clinvar", "term": term, "retmode": "json", "retmax": retmax})
    response = client.fetch(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{query}",
        dataset_version=release,
    )
    body = response["body"]
    found = body.get("esearchresult") if isinstance(body, dict) else None
    ids = found.get("idlist") if isinstance(found, dict) else None
    try:
        count = int(found["count"])
        limit = int(found["retmax"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Unexpected ClinVar gene search response") from exc
    if not isinstance(ids, list) or any(not str(value).isdigit() for value in ids):
        raise ValueError("Unexpected ClinVar gene search ID list")
    return {
        "term": term, "ids": [str(value) for value in ids], "count": count,
        "complete": count <= limit and len(ids) == count,
        "retrieved_at": response["retrieved_at"],
        "digest": hashlib.sha256(canonical_json(body).encode()).hexdigest(),
    }


class UpstreamPathogenicProvider:
    """Answers IC03, as fields merged into the initiation record rather than a second one.

    PVS1 selects one `initiation_assessment` per transcript and reports two of them as a
    conflict, so this returns the fields for the caller to fold into the record
    InitiationProvider already built - the same reason the MANE provider carries NF03's
    answer instead of a provider of its own.
    """

    name = "ClinVar upstream pathogenic report"

    def __init__(self, client, release, *, summary_batch=SUMMARY_BATCH):
        if not release:
            raise ValueError("ClinVar release must be recorded")
        self.client = client
        self.release = release
        self.summary_batch = summary_batch

    def _placed(self, documents, gene):
        """(positions of placeable records, count that could not be placed)."""
        placed, unplaced = [], 0
        for document in documents.values():
            if not isinstance(document, dict):
                continue
            symbol = document.get("gene_sort") or ""
            symbols = {value.get("symbol") for value in document.get("genes", [])
                       if isinstance(value, dict)}
            if symbol != gene and gene not in symbols:
                continue
            positions = protein_change_positions(document.get("protein_change"))
            if positions:
                placed.append((document, positions))
            else:
                unplaced += 1
        return placed, unplaced

    def get_upstream_evidence(self, variant, gene, transcript, downstream_start_codon):
        """Fields to merge into the initiation record when the reports settle IC03, else {}.

        Returns a dict rather than a record: see the class docstring.
        """
        if not (gene and transcript) or not isinstance(downstream_start_codon, int):
            return {}
        if downstream_start_codon < 2:
            return {}
        found = gene_pathogenic_search(self.client, self.release, gene)
        if not found["complete"]:
            return {}
        documents, retrieved_at = summary_documents(
            self.client, self.release, found["ids"], self.summary_batch)
        placed, unplaced = self._placed(documents, gene)
        upstream = [
            {"variation_id": document.get("uid"), "accession": document.get("accession"),
             "protein_change": document.get("protein_change"),
             "positions": [position for position in positions
                           if position < downstream_start_codon],
             "classification": ((document.get("germline_classification")
                                 or document.get("clinical_significance") or {})
                                .get("description"))}
            for document, positions in placed
            if any(position < downstream_start_codon for position in positions)
        ]
        if not upstream and unplaced:
            # Any unplaced record could be the upstream evidence, so "no" is not supported.
            return {}
        digest = hashlib.sha256(
            f"{found['digest']}:{downstream_start_codon}".encode("utf-8")).hexdigest()
        return {
            "upstream_pathogenic_evidence": bool(upstream),
            "upstream_method": METHOD,
            "upstream_source": self.name,
            "upstream_source_version": self.release,
            "upstream_retrieved_at": retrieved_at or found["retrieved_at"],
            "upstream_query": found["term"],
            "upstream_returned_count": found["count"],
            "upstream_placed_count": len(placed),
            "upstream_unplaced_count": unplaced,
            "upstream_variants": upstream[:20],
            "upstream_policy_note": (
                "Positions come from ClinVar's reported protein change, which carries one "
                "only for plain single-residue substitutions. A yes needs one placed "
                "pathogenic substitution upstream of the restart codon; a no additionally "
                "needs every returned record to have been placeable, because an unplaced "
                "nonsense or frameshift change could be the upstream evidence."
            ),
            "upstream_search_sha256": found["digest"],
            "upstream_digest": digest,
        }
