"""MANE Select as an automated transcript-relevance signal for PVS1's NF01 gate.

[Why this exists]
  PVS1 asks, before anything else about the variant, whether it sits on a biologically
  relevant transcript. That is `transcript_assessment.relevance`, a curated judgment, and
  without one the tree stops at NF01 - which is where the truncating variants in the
  ground-truth set stopped once ClinGen dosage let them past G01. MANE Select is the
  RefSeq/Ensembl joint pick of the representative transcript for a gene, published and
  versioned, so it answers the same question for the common case.

[Why a non-MANE transcript is unresolved, not irrelevant]
  This repository already recorded the trap: test_data/resolve_erepo_transcripts.py found
  that MANE Select is not always the transcript a variant is actually curated against - for
  TNNT2 the MANE pick is NM_001276345.2 while demo-data and ERepo use a different accession.
  Treating "not MANE" as NOT_RELEVANT would declare PVS1 inapplicable for every such
  variant, and be wrong every time the curation used a legitimate alternative transcript.
  So only a match produces a record; anything else produces none, and PVS1 keeps reporting
  the transcript as unresolved for a curator.

[What it deliberately does not supply]
  `exon_relevance`, which NF03 needs after an NMD prediction, asks whether the affected exon
  is present in the relevant transcripts. That is a comparison of exon structures across
  transcripts, not something the MANE pick answers, so no value is invented for it and NF03
  still stops.

[Version matching]
  Accessions are compared without their version suffix. A transcript's version changes when
  its sequence record is updated; the MANE summary pins one version, and refusing
  NM_000256.3 because the summary now says NM_000256.4 would withhold the signal over a
  difference that is not about which transcript was chosen. The version actually seen and
  the version MANE names are both recorded.
"""

from __future__ import annotations

import csv
import hashlib
import io
import re

SUMMARY_DIRECTORY = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
METHOD = "mane_select_transcript"
_SUMMARY_NAME = re.compile(r"MANE\.GRCh38\.v([\d.]+)\.summary\.txt\.gz")


def accession(transcript):
    """The accession without its version suffix, or None."""
    return transcript.split(".")[0] if transcript else None


def parse(text):
    """Return {gene symbol: RefSeq nucleotide accession} for the MANE Select set."""
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not reader.fieldnames or "symbol" not in reader.fieldnames:
        raise ValueError("MANE summary is missing its header")
    selected = {}
    for row in reader:
        if (row.get("MANE_status") or "").strip() == "MANE Select":
            symbol = (row.get("symbol") or "").strip()
            nucleotide = (row.get("RefSeq_nuc") or "").strip()
            if symbol and nucleotide:
                selected[symbol] = nucleotide
    return selected


class ManeTranscriptProvider:
    """Builds `transcript_assessment` evidence when the evaluated transcript is MANE Select."""

    name = "NCBI MANE Select"

    def __init__(self, client, *, summary_url, release):
        if not summary_url or not release:
            raise ValueError("The MANE summary URL and its release must both be recorded")
        self.client = client
        self.summary_url = summary_url
        self.release = release

    @classmethod
    def from_directory(cls, client, *, directory=SUMMARY_DIRECTORY):
        """Discover the current summary file, so the release is read rather than assumed."""
        listing = client.fetch(directory, response_format="text")["body"]
        names = sorted(set(_SUMMARY_NAME.findall(listing)))
        if not names:
            raise ValueError("No MANE summary file found in the release directory")
        release = names[-1]
        return cls(client, summary_url=f"{directory}MANE.GRCh38.v{release}.summary.txt.gz",
                   release=release)

    def get_transcript_assessment(self, variant, gene, transcript):
        """One automated record when `transcript` is the gene's MANE Select, else none."""
        if not gene or not transcript:
            return []
        response = self.client.fetch(self.summary_url, response_format="text-gz")
        text = response["body"]
        selected = parse(text)
        mane = selected.get(gene)
        if not mane or accession(mane) != accession(transcript):
            return []
        digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        return [{
            "category": "transcript_assessment",
            "variant_key": variant.key,
            "evidence_id": f"urn:sha256:{digest}:mane:{transcript}:{variant.key}",
            "source": self.name,
            "source_version": f"MANE v{self.release}",
            "retrieved_at": response["retrieved_at"],
            "quality_status": "PASS",
            "transcript": transcript,
            "gene": gene,
            "relevance": "RELEVANT",
            # Not supplied on purpose - see the module docstring. NF03 stops here.
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": f"MANE v{self.release}",
            "mane_select_accession": mane,
            "evaluated_accession": transcript,
            "policy_note": (
                "Relevance asserted because the evaluated transcript is this gene's MANE "
                "Select. A transcript that is not MANE Select yields no record rather than "
                "NOT_RELEVANT: curation legitimately uses other transcripts. exon_relevance "
                "is not derived here."
            ),
            "response_sha256": digest,
        }]
