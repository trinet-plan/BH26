"""Nonsense-mediated decay prediction for PVS1's NF02 gate, from exon numbering.

[The rule]
  ClinGen's PVS1 specification (Abou Tayoun et al. 2018) treats a premature termination
  codon as escaping NMD when it lies in the last exon, or within the last 50 nucleotides of
  the penultimate exon. Everything upstream of that is predicted to undergo NMD.

[What this can and cannot decide]
  Exon numbering alone settles the upstream case: a PTC in exon 2 of 34 is nowhere near the
  final junction, whichever way the last exons are sized. It does not settle the boundary
  case - "within the last 50 nucleotides of the penultimate exon" needs that exon's length,
  which the numbering does not carry. So a variant in either of the last two exons produces
  no record, PVS1 reports the NMD prediction as unresolved, and a curator decides. Guessing
  there would be guessing at exactly the distinction the rule exists to draw.

[Why it asks Ensembl again instead of reusing the annotation]
  The annotation this pipeline already has does not carry exon numbers: the VEP request that
  produces it does not ask for them, and widening that request would change the URL the
  committed offline caches are keyed on - tests/fixtures/ensembl-cache replays those exact
  requests, and an offline demo run asserts no network was used. This provider makes its own
  request with `numbers=1`, so the existing cache stays valid and this one is additive.

[Transcript matching]
  VEP answers about Ensembl transcripts, and the variant is evaluated against a RefSeq
  accession, so the consequence is matched on MANE Select where VEP reports it and otherwise
  on the gene symbol plus a truncating consequence. When more than one candidate disagrees
  about the exon, no record is emitted rather than one of them being picked.
"""

from __future__ import annotations

import hashlib
from urllib.parse import quote

RULE_SOURCE = "ClinGen PVS1 2018"
METHOD = "clingen_pvs1_2018_nmd_exon_rule"
TRUNCATING = {"stop_gained", "frameshift_variant"}
# The rule's own boundary: the last exon, and the last 50 nt of the one before it.
UNDECIDABLE_FROM_NUMBERING = 2


def parse_exon(value):
    """"8/34" -> (8, 34). None when VEP did not report a position."""
    if not isinstance(value, str) or "/" not in value:
        return None
    index, _, total = value.partition("/")
    # An exon span ("8-9/34") means the variant crosses a junction; the rule is about where
    # the PTC falls, so the first exon is the one that matters.
    index = index.split("-")[0]
    try:
        return int(index), int(total)
    except ValueError:
        return None


class NmdPredictionProvider:
    """Builds `nmd_prediction` evidence from VEP exon numbering, for the decidable cases."""

    name = "Ensembl VEP exon numbering"

    def __init__(self, client, release):
        if not release:
            raise ValueError("Ensembl release must be recorded")
        self.client = client
        self.release = release

    def exon_on_transcript(self, gene, transcript, hgvsc):
        """The exon this variant falls in, as VEP numbers it on `transcript`.

        Shared with the MANE provider, which needs it for PVS1's NF03 gate: an exon that VEP
        can number on a transcript is by construction present in that transcript. Returns
        None when no single consequence settles it - the same silence the NMD prediction
        keeps in that case.
        """
        if not (gene and transcript and hgvsc):
            return None
        _response, candidates = self._consequences(hgvsc, gene)
        chosen = self._for_transcript(candidates, transcript)
        if not chosen:
            return None
        exons = {item.get("exon") for item in chosen}
        return exons.pop() if len(exons) == 1 else None

    @staticmethod
    def _for_transcript(candidates, transcript):
        """The consequences VEP ties to `transcript` via MANE Select, else all of them."""
        preferred = [item for item in candidates
                     if (item.get("mane_select") or "").split(".")[0] == transcript.split(".")[0]]
        return preferred or candidates

    def _consequences(self, hgvsc, gene):
        url = ("https://rest.ensembl.org/vep/human/hgvs/"
               f"{quote(hgvsc, safe='')}?hgvs=1&numbers=1&mane=1")
        response = self.client.fetch(url, response_format="json",
                                     dataset_version=self.release)
        body = response["body"]
        if not isinstance(body, list) or not body:
            raise ValueError("Unexpected VEP response")
        candidates = [item for item in body[0].get("transcript_consequences", [])
                      if item.get("gene_symbol") == gene
                      and set(item.get("consequence_terms") or []) & TRUNCATING
                      and parse_exon(item.get("exon"))]
        return response, candidates

    def get_nmd_prediction(self, variant, gene, transcript, hgvsc):
        """One automated record when the exon numbering settles the question, else none."""
        if not (gene and transcript and hgvsc):
            return []
        response, candidates = self._consequences(hgvsc, gene)
        if not candidates:
            return []
        # Prefer the consequence VEP marks as MANE Select for this accession; the variant is
        # evaluated against a RefSeq transcript, and that is the only field tying the two.
        chosen = self._for_transcript(candidates, transcript)
        positions = {parse_exon(item.get("exon")) for item in chosen}
        if len(positions) != 1:
            return []
        index, total = positions.pop()
        if total < 1 or not 1 <= index <= total:
            return []
        if index > total - UNDECIDABLE_FROM_NUMBERING:
            # Last exon, or the penultimate one where the 50-nt boundary decides it.
            return []
        digest = hashlib.sha256(response["body_sha256"].encode("utf-8")).hexdigest()
        item = chosen[0]
        return [{
            "category": "nmd_prediction",
            "variant_key": variant.key,
            "evidence_id": f"urn:sha256:{digest}:nmd:{transcript}:{variant.key}",
            "source": self.name,
            "source_version": self.release,
            "retrieved_at": response["retrieved_at"],
            "quality_status": "PASS",
            "transcript": transcript,
            "gene": gene,
            "predicted": True,
            "rule_source": RULE_SOURCE,
            "exon": item.get("exon"),
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": RULE_SOURCE,
            "ensembl_transcript": item.get("transcript_id"),
            "policy_note": (
                "NMD predicted because the premature termination codon lies upstream of the "
                "final two exons. A variant in the last exon, or in the penultimate one "
                "where the rule's 50-nucleotide boundary applies, yields no record: exon "
                "numbering does not carry the distance that decides it."
            ),
            "response_sha256": response["body_sha256"],
        }]
