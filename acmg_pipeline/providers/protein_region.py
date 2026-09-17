"""The protein-loss measurement PVS1's NF07 gate weighs, from Ensembl.

[What it supplies]
  When a premature termination codon escapes NMD, PVS1 asks how much of the protein is lost
  and splits strong from moderate on a threshold the rule set carries. That is arithmetic
  once two numbers are known - where the codon falls and how long the protein is - and both
  are published: VEP reports the protein position, Ensembl reports the translation length.

[What it deliberately does not supply]
  The same `protein_region` record is read by two other gates, and neither is arithmetic.
  NF04 asks whether a critical functional region is disrupted; NF06 asks whether the region
  is biologically relevant. Those are judgments about what the lost residues do, which no
  coordinate answers, so this leaves them unset and PVS1 stops at NF06 asking for them - with
  the measurement already in hand, so the curator answers one question rather than three.

[The measurement]
  lost_residues counts the codon itself and everything downstream: total - position + 1. For
  a frameshift that is an approximation - translation continues into a new frame until it
  meets a stop - but it is the same approximation the fraction rule is stated in terms of,
  and the alternative is to claim a precision the input does not have.
"""

from __future__ import annotations

import hashlib

from acmg_pipeline.providers.nmd import TRUNCATING, NmdPredictionProvider

METHOD = "ensembl_translation_length_protein_loss"


class ProteinRegionProvider:
    """Builds `protein_region` evidence carrying the protein-loss fraction's two inputs."""

    name = "Ensembl translation length"

    def __init__(self, client, release):
        if not release:
            raise ValueError("Ensembl release must be recorded")
        self.client = client
        self.release = release
        # The VEP request is the one the NMD rule already makes, so this is a cache hit.
        self._consequences = NmdPredictionProvider(client, release)

    def _translation_length(self, ensembl_transcript):
        response = self.client.fetch(
            f"https://rest.ensembl.org/lookup/id/{ensembl_transcript}?expand=1",
            response_format="json", dataset_version=self.release)
        translation = (response["body"] or {}).get("Translation") or {}
        length = translation.get("length")
        if not isinstance(length, int) or length <= 0:
            return None, None, response
        return length, translation.get("id"), response

    def get_protein_region(self, variant, gene, transcript, hgvsc, protein_start):
        """One automated record when both numbers are known, else none."""
        if not (gene and transcript and hgvsc) or not isinstance(protein_start, int):
            return []
        if protein_start < 1:
            return []
        _response, candidates = self._consequences._consequences(hgvsc, gene)
        chosen = self._consequences._for_transcript(candidates, transcript)
        if not chosen:
            return []
        ensembl_transcripts = {item.get("transcript_id") for item in chosen}
        if len(ensembl_transcripts) != 1:
            return []
        ensembl_transcript = ensembl_transcripts.pop()
        if not ensembl_transcript:
            return []
        total, protein_id, response = self._translation_length(ensembl_transcript)
        if total is None or protein_start > total:
            return []
        lost = total - protein_start + 1
        digest = hashlib.sha256(
            f"{response['body_sha256']}:{protein_start}".encode("utf-8")).hexdigest()
        return [{
            "category": "protein_region",
            "variant_key": variant.key,
            "evidence_id": f"urn:sha256:{digest}:protein-region:{transcript}:{variant.key}",
            "source": self.name,
            "source_version": self.release,
            "retrieved_at": response["retrieved_at"],
            "quality_status": "PASS",
            "transcript": transcript,
            "gene": gene,
            "protein_id": protein_id,
            "lost_residues": lost,
            "total_protein_length": total,
            # Left unset on purpose - see the module docstring. NF04 and NF06 still ask.
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": self.release,
            "ensembl_transcript": ensembl_transcript,
            "protein_start": protein_start,
            "policy_note": (
                "Measurement only. lost_residues is total - protein_start + 1, counting the "
                "termination codon and everything downstream. Whether the lost region is "
                "critical (NF04) or biologically relevant (NF06) is not derived here."
            ),
            "response_sha256": response["body_sha256"],
        }]


__all__ = ["METHOD", "TRUNCATING", "ProteinRegionProvider"]
