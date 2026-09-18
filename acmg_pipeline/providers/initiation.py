"""The downstream start codon PVS1's IC02 gate asks about, from the CDS sequence -
plus, since 2026-09-17, a flagged first-pass answer for IC01.

[What it supplies]
  When the initiation codon is lost, PVS1 asks whether translation could restart at an
  in-frame ATG further along. That is a property of the coding sequence and nothing else:
  read the CDS, step three bases at a time, look for ATG. Ensembl publishes the sequence, so
  the answer is derivable and exact.

[IC01 - now a flagged Ensembl-derived prediction, not left unset]
  IC01 asks whether an intact, biologically relevant alternative transcript rescues the
  loss. A naive version of this check - "does the gene have more than one Ensembl
  transcript?" - was tried and rejected: GJB2 alone has 11 transcript models, and every one
  of them shares the exact same start-codon genomic position (confirmed against real
  Ensembl data), so "multiple transcripts exist" would have wrongly read as "an alternative
  start exists" for a gene that, in reality, has none. What this checks instead is whether
  any OTHER protein_coding transcript's own start codon sits at a DIFFERENT genomic
  position (correctly accounting for strand - the start codon is at Translation.end on the
  minus strand, Translation.start on the plus strand) than the transcript actually being
  evaluated. Finding none (as for GJB2's real 11 transcripts) sets
  intact_alternative_transcript=False, matching what the real biology is - GJB2 has no
  rescuing alternative isoform. Finding a genuinely different start position sets it True.
  There is no known real case in this project where the correct answer is True to check
  this against (the two real cases below were decided by a VCEP override of the general
  framework instead), so - more than NF04/NF06/SP01 - this is UNVALIDATED: always disclosed
  in acmg_pipeline.criteria.pvs1's review_points as an Ensembl-derived prediction, never a
  curator's own review.

  IC03 asks whether a pathogenic variant has been reported upstream of that downstream start,
  which decides moderate against supporting. It needs a ClinVar region search bounded by the
  codon found here - see providers/upstream_pathogenic.py.

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
IC01_METHOD = "ensembl_alternative_start_position"
START_CODON = "ATG"
START_LOST = "start_lost"


def _start_codon_genomic_position(transcript_record):
    translation = transcript_record.get("Translation")
    if not isinstance(translation, dict):
        return None
    return translation.get("end") if transcript_record.get("strand") == -1 else translation.get("start")


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

    def _alternative_start(self, gene, ensembl_transcript):
        """{} if unresolvable (honest gap), else IC01 fields to merge - see module docstring."""
        try:
            response = self.client.fetch(
                f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{quote(gene, safe='')}"
                "?expand=1", response_format="json", dataset_version=self.release)
        except ValueError:
            return {}
        body = response["body"]
        transcripts = body.get("Transcript") if isinstance(body, dict) else None
        if not isinstance(transcripts, list):
            return {}
        by_id = {item.get("id"): item for item in transcripts if isinstance(item, dict)}
        evaluated = by_id.get(ensembl_transcript)
        if evaluated is None:
            return {}
        own_position = _start_codon_genomic_position(evaluated)
        if own_position is None:
            return {}
        alternative = next(
            (other_id for other_id, other in by_id.items()
             if other_id != ensembl_transcript and other.get("biotype") == "protein_coding"
             and _start_codon_genomic_position(other) not in (None, own_position)),
            None,
        )
        return {
            "intact_alternative_transcript": alternative is not None,
            "alternative_transcript_method": IC01_METHOD,
            "alternative_transcript_source_version": self.release,
            "alternative_transcript_candidate": alternative,
        }

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
        alternative_start = self._alternative_start(gene, ensembl_transcript)
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
            # IC03 (upstream_pathogenic_evidence) still asks - see providers/upstream_pathogenic.py.
            "downstream_start_codon": codon,
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": self.release,
            "ensembl_transcript": ensembl_transcript,
            "policy_note": (
                "downstream_in_frame_start/downstream_start_codon are read directly from the "
                "coding sequence. intact_alternative_transcript (IC01), when present, is a "
                "flagged Ensembl-derived prediction (alternative_transcript_method), not a "
                "curator's own review and not validated against a known real case - see this "
                "module's own docstring. upstream_pathogenic_evidence (IC03) is not derived here."
            ),
            "response_sha256": response["body_sha256"],
            **alternative_start,
        }]
