"""The downstream start codon PVS1's IC02 gate asks about, from the CDS sequence.

[What it supplies]
  When the initiation codon is lost, PVS1 asks whether translation could restart at an
  in-frame ATG further along. That is a property of the coding sequence and nothing else:
  read the CDS, step three bases at a time, look for ATG. Ensembl publishes the sequence, so
  the answer is derivable and exact.

[What it deliberately does not supply]
  The same `initiation_assessment` record carries two more fields, and this fills neither.

  IC01 asks whether an intact, biologically relevant alternative transcript exists. "Biologically
  relevant" is the same judgment NF06 turns on, and PVS1 reads it first - so the path still
  stops there, with the downstream-start question already answered and one left for a curator.

  IC03 asks whether a pathogenic variant has been reported upstream of that downstream start,
  which decides moderate against supporting. It needs a ClinVar region search bounded by the
  codon found here, so it is a separate piece of work rather than something this can infer.

[A caveat worth stating]
  ClinGen's general framework caps initiation-codon loss well below very_strong, but the
  expert panels that curated the two real cases here recorded PVS1 for GJB2 c.2T>C and
  PVS1_Strong for HNF4A c.3G>A. Those are VCEP specifications overriding the general rule,
  not a defect in it: a run of this path will reach moderate or supporting and disagree with
  those labels for a reason that is about which specification applies.
"""

from __future__ import annotations

import hashlib
from urllib.parse import quote

METHOD = "ensembl_cds_downstream_in_frame_start"
START_CODON = "ATG"
START_LOST = "start_lost"


def downstream_in_frame_start(cds):
    """The 1-based codon number of the first in-frame ATG after the initiation codon.

    None when the sequence does not support the question - it is not a multiple of three,
    it is too short to hold a second codon, or it simply has no further ATG. The last of
    those is a real answer, which is why the caller distinguishes it from the others.
    """
    if not isinstance(cds, str):
        return None, "unusable"
    sequence = cds.strip().upper()
    if len(sequence) < 6 or len(sequence) % 3:
        return None, "unusable"
    for offset in range(3, len(sequence) - 2, 3):
        if sequence[offset:offset + 3] == START_CODON:
            return offset // 3 + 1, "found"
    return None, "absent"


class InitiationProvider:
    """Builds `initiation_assessment` evidence carrying the downstream-start answer."""

    name = "Ensembl CDS sequence"

    def __init__(self, client, release):
        if not release:
            raise ValueError("Ensembl release must be recorded")
        self.client = client
        self.release = release

    def _consequence(self, hgvsc, gene, transcript):
        url = ("https://rest.ensembl.org/vep/human/hgvs/"
               f"{quote(hgvsc, safe='')}?hgvs=1&numbers=1&mane=1")
        response = self.client.fetch(url, response_format="json",
                                     dataset_version=self.release)
        body = response["body"]
        if not isinstance(body, list) or not body:
            raise ValueError("Unexpected VEP response")
        candidates = [item for item in body[0].get("transcript_consequences", [])
                      if item.get("gene_symbol") == gene
                      and START_LOST in (item.get("consequence_terms") or [])
                      and item.get("transcript_id")]
        preferred = [item for item in candidates
                     if (item.get("mane_select") or "").split(".")[0] == transcript.split(".")[0]]
        chosen = preferred or candidates
        identifiers = {item["transcript_id"] for item in chosen}
        return identifiers.pop() if len(identifiers) == 1 else None

    def get_initiation_assessment(self, variant, gene, transcript, hgvsc):
        """One automated record when the CDS settles the downstream start, else none."""
        if not (gene and transcript and hgvsc):
            return []
        ensembl_transcript = self._consequence(hgvsc, gene, transcript)
        if not ensembl_transcript:
            return []
        response = self.client.fetch(
            f"https://rest.ensembl.org/sequence/id/{ensembl_transcript}?type=cds",
            response_format="json", dataset_version=self.release)
        codon, outcome = downstream_in_frame_start((response["body"] or {}).get("seq"))
        if outcome == "unusable":
            return []
        digest = hashlib.sha256(response["body_sha256"].encode("utf-8")).hexdigest()
        return [{
            "category": "initiation_assessment",
            "variant_key": variant.key,
            "evidence_id": f"urn:sha256:{digest}:initiation:{transcript}:{variant.key}",
            "source": self.name,
            "source_version": self.release,
            "retrieved_at": response["retrieved_at"],
            "quality_status": "PASS",
            "transcript": transcript,
            "gene": gene,
            "downstream_in_frame_start": outcome == "found",
            # Left unset on purpose - see the module docstring. IC01 and IC03 still ask.
            "downstream_start_codon": codon,
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": self.release,
            "ensembl_transcript": ensembl_transcript,
            "policy_note": (
                "Read from the coding sequence: the first in-frame ATG after the initiation "
                "codon, or its absence. Whether an intact alternative transcript exists "
                "(IC01) and whether a pathogenic variant is reported upstream of this codon "
                "(IC03) are not derived here."
            ),
            "response_sha256": response["body_sha256"],
        }]
