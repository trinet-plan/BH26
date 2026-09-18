"""Nonsense-mediated decay prediction for PVS1's NF02 gate, from exon numbering.

[The rule]
  ClinGen's PVS1 specification (Abou Tayoun et al. 2018) treats a premature termination
  codon as escaping NMD when it lies in the last exon, or within the last 50 nucleotides of
  the penultimate exon. Everything upstream of that is predicted to undergo NMD.

[Canonical splice donor/acceptor variants, since 2026-09-17]
  A splice donor/acceptor loss has no exon number of its own from VEP - it reports an
  INTRON number instead (e.g. "27/34"), since the variant itself sits at the exon-intron
  boundary. This module now accepts that case too, applying the identical last-two-
  positions rule to the intron count as a stand-in for the exon count (a gene's intron
  count and exon count differ by exactly one, so this is an approximation, not exact - the
  real position depends on which exon is actually skipped or which intron retained, which
  needs RNA-level data this project does not have). Every such record is flagged
  (position_type: "intron") and acmg_pipeline.criteria.pvs1 adds a review_points caveat
  whenever an intron-derived answer is what actually decided NF02/NF03 - checked against
  MYBPC3 c.2905+1G>A (VEP: intron 27/34 on the MANE transcript; ClinGen: PVS1_very_strong):
  27/34 is nowhere near either end, so this reaches "predicted" the same as the real answer.

[What this can and cannot decide]
  Exon numbering settles both ends of the rule and neither is a guess. A PTC upstream of the
  final two exons is nowhere near the junction whichever way those exons are sized, so NMD is
  predicted. A PTC in the last exon is past the junction entirely, so NMD is escaped - the
  50-nucleotide distance does not enter into it, because that distance is measured within the
  *penultimate* exon.

  The penultimate exon is the one case numbering cannot settle: "within the last 50
  nucleotides" needs that exon's length and the variant's position in it. A PTC there
  produces no record, PVS1 reports the prediction as unresolved, and a curator decides.
  Guessing there would guess at exactly the distinction the rule exists to draw.

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
CANONICAL_SPLICE = {"splice_donor_variant", "splice_acceptor_variant"}


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


def _is_decidable(item):
    """A truncating consequence with a numbered exon, or a canonical splice donor/acceptor
    consequence with a numbered intron - see the module docstring for the splice case."""
    terms = set(item.get("consequence_terms") or [])
    if terms & TRUNCATING:
        return bool(parse_exon(item.get("exon")))
    if terms & CANONICAL_SPLICE:
        return bool(parse_exon(item.get("intron")))
    return False


def _position(item):
    """(raw "N/M" string, "exon"|"intron") for whichever field _is_decidable() accepted."""
    exon = item.get("exon")
    if parse_exon(exon):
        return exon, "exon"
    return item.get("intron"), "intron"


class NmdPredictionProvider:
    """Builds `nmd_prediction` evidence from VEP exon numbering, for the decidable cases."""

    name = "Ensembl VEP exon numbering"

    def __init__(self, client, release):
        if not release:
            raise ValueError("Ensembl release must be recorded")
        self.client = client
        self.release = release

    def exon_on_transcript(self, gene, transcript, hgvsc):
        """The exon (or, for a canonical splice variant, intron) this variant falls in,
        as VEP numbers it on `transcript`.

        Shared with the MANE provider, which needs it for PVS1's NF03 gate: a position VEP
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
        positions = {_position(item)[0] for item in chosen}
        return positions.pop() if len(positions) == 1 else None

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
                      and _is_decidable(item)]
        return response, candidates

    def get_nmd_prediction(self, variant, gene, transcript, hgvsc):
        """One automated record when the exon/intron numbering settles the question, else none."""
        if not (gene and transcript and hgvsc):
            return []
        response, candidates = self._consequences(hgvsc, gene)
        if not candidates:
            return []
        # Prefer the consequence VEP marks as MANE Select for this accession; the variant is
        # evaluated against a RefSeq transcript, and that is the only field tying the two.
        chosen = self._for_transcript(candidates, transcript)
        positions = {_position(item) for item in chosen}
        if len(positions) != 1:
            return []
        (raw_position, position_type), = positions
        index, total = parse_exon(raw_position)
        if total < 1 or not 1 <= index <= total:
            return []
        if index == total - 1:
            # The penultimate exon (or, for a splice variant, intron): the rule's
            # 50-nucleotide boundary decides it, and the numbering does not carry that
            # distance.
            return []
        predicted = index < total
        digest = hashlib.sha256(response["body_sha256"].encode("utf-8")).hexdigest()
        item = chosen[0]
        note = (
            "NMD predicted when the premature termination codon lies upstream of the "
            "final two exons, and escaped when it lies in the last exon - the rule's "
            "50-nucleotide distance is measured inside the penultimate exon, so it does "
            "not bear on either. A PTC in the penultimate exon yields no record: exon "
            "numbering does not carry the distance that decides it."
        ) if position_type == "exon" else (
            "This variant is a canonical splice donor/acceptor loss, not a direct "
            "truncating consequence: VEP reports an intron number, not an exon number, so "
            "the same last-two-positions rule is applied to the intron count as a "
            "flagged approximation (a gene's intron and exon counts differ by exactly "
            "one) - not the exact exon a real skipped-exon or retained-intron outcome "
            "would fall in, which needs RNA-level data this does not have."
        )
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
            "predicted": predicted,
            "rule_source": RULE_SOURCE,
            "exon": raw_position,
            "position_type": position_type,
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": RULE_SOURCE,
            "ensembl_transcript": item.get("transcript_id"),
            "policy_note": note,
            "response_sha256": response["body_sha256"],
        }]
