"""ClinGen Gene-Disease Validity bulk CSV adapter.

ClinGen publishes every GCEP's completed Gene-Disease Validity curation as one
CSV snapshot (no API key, no per-gene query) at the URL below. Each row is one
(gene, disease) curation: a classification (Definitive/Strong/Moderate/
Limited/Disputed/Refuted/No Known Disease Relationship), the mode of
inheritance the GCEP curated it under, and a citable report URL - see
https://search.clinicalgenome.org/kb/downloads.
"""

from __future__ import annotations

import csv
import hashlib
import io
import re

DOWNLOAD_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
_FILE_CREATED = re.compile(r"FILE CREATED:\s*([0-9]{4}-[0-9]{2}-[0-9]{2})")


class ClinGenGeneValidityProvider:
    """Fetches and caches the whole snapshot once, then filters it per gene."""

    name = "ClinGen Gene-Disease Validity"

    def __init__(self, client):
        self.client = client
        self._rows = None
        self._source_version = None

    def _load(self):
        if self._rows is not None:
            return
        response = self.client.fetch(DOWNLOAD_URL, response_format="text")
        body = response["body"]
        if not isinstance(body, str) or not body.strip():
            raise ValueError("Unexpected ClinGen Gene-Disease Validity response")
        lines = body.splitlines()
        created = next((_FILE_CREATED.search(line) for line in lines[:5]
                        if _FILE_CREATED.search(line)), None)
        if created is None:
            raise ValueError("ClinGen Gene-Disease Validity CSV has no FILE CREATED header")
        header_index = next(
            (i for i, line in enumerate(lines) if line.startswith('"GENE SYMBOL"')), None,
        )
        if header_index is None:
            raise ValueError("ClinGen Gene-Disease Validity CSV has no GENE SYMBOL header row")
        reader = csv.DictReader(io.StringIO("\n".join(lines[header_index:])))
        self._rows = list(reader)
        self._source_version = created.group(1)
        self._retrieved_at = response["retrieved_at"]
        self._digest = hashlib.sha256(body.encode()).hexdigest()

    def get_validity(self, gene: str) -> list[dict]:
        """All curated (gene, disease) validity rows for `gene`, or [] if none."""
        self._load()
        return [
            {
                "category": "gene_disease_validity",
                "gene": row.get("GENE SYMBOL"),
                "hgnc_id": row.get("GENE ID (HGNC)"),
                "condition": row.get("DISEASE ID (MONDO)"),
                "condition_label": row.get("DISEASE LABEL"),
                "moi": row.get("MOI"),
                "classification": row.get("CLASSIFICATION"),
                "classification_date": row.get("CLASSIFICATION DATE"),
                "gcep": row.get("GCEP"),
                "source": self.name,
                "source_version": self._source_version,
                "retrieved_at": self._retrieved_at,
                "evidence_id": row.get("ONLINE REPORT"),
                "quality_status": "PASS",
            }
            for row in self._rows
            if row.get("GENE SYMBOL") == gene
        ]
