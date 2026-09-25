"""The protein-loss measurement PVS1's NF07 gate weighs, from Ensembl - plus,
since 2026-09-17, a UniProt-derived first-pass answer for NF04/NF06.

[What it supplies]
  When a premature termination codon escapes NMD, PVS1 asks how much of the protein is lost
  and splits strong from moderate on a threshold the rule set carries. That is arithmetic
  once two numbers are known - where the codon falls and how long the protein is - and both
  are published: VEP reports the protein position, Ensembl reports the translation length.

[NF04/NF06 - now a flagged UniProt-derived prediction, not left unset]
  The same `protein_region` record is read by two other gates. NF04 asks whether a critical
  functional region is disrupted; NF06 asks whether the region is biologically relevant.
  Both used to be left for a curator to answer from scratch. They now get a first-pass
  answer from UniProt's own feature annotations for the lost interval
  [protein_start, total_protein_length] - always flagged as a prediction (`review` on the
  criterion result), never presented as a curator's own reviewed judgment:

  - critical_region_disrupted = True only when the lost interval overlaps a UniProt Domain/
    Active site/Coiled coil feature, or a Region feature NOT described as "Disordered" - see
    uniprot_features.CRITICAL_FEATURE_TYPES's own comment for why "Binding site" and "Motif"
    are excluded (single-residue biochemical contacts and short linear motifs are not what
    ACMG PVS1's "critical functional domain" means, and both were confirmed real false-
    positive sources: MYOC/PTEN, 2026-09-24). "Chain" (UniProt's whole-protein-length
    feature) is deliberately excluded too - it always overlaps any interval and would make
    this check meaningless. Left unset (not False) otherwise: no named domain there is not
    proof nothing important is there.
  - region_biologically_relevant defaults to True whenever UniProt features were fetched
    successfully, and only False when the ENTIRE lost interval is covered by "Disordered"-
    described Region feature(s). This default direction (relevant unless proven otherwise)
    was checked against the one real case with a known answer this project has: VHL
    c.610G>T (ClinGen expert panel: PVS1_Moderate, i.e. the panel treated the lost 10
    C-terminal residues as biologically relevant) has no UniProt Domain over that interval
    (its domains stop at residue 192) and no "Disordered" annotation there either - a
    stricter "must overlap a named Domain to count as relevant" reading would have given
    the wrong answer on this exact variant, so it is not the rule used. See
    tests/test_pvs1_automated_gates.py's own VHL fixture for the ground truth this was
    checked against.
  Left unset (None, not False) when UniProt itself could not be resolved (no accession, no
  features, a fetch error) - an honest gap, same as everywhere else in this project.

[The measurement]
  lost_residues counts the codon itself and everything downstream: total - position + 1. For
  a frameshift that is an approximation - translation continues into a new frame until it
  meets a stop - but it is the same approximation the fraction rule is stated in terms of,
  and the alternative is to claim a precision the input does not have.
"""

from __future__ import annotations

import hashlib

from acmg_pipeline.criteria.reference_links import gene_to_uniprot_accession
from acmg_pipeline.providers.nmd import TRUNCATING, NmdPredictionProvider
from acmg_pipeline.providers.uniprot_features import (
    CRITICAL_FEATURE_TYPES, fetch_features, feature_interval, is_disordered, overlaps,
)

METHOD = "ensembl_translation_length_protein_loss"
UNIPROT_METHOD = "uniprot_feature_overlap"


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

    def _uniprot_relevance(self, gene, protein_start, total):
        """A flagged NF04/NF06 first-pass answer from UniProt features, or {} if
        UniProt itself could not be resolved (honest gap - see module docstring).
        """
        accession = gene_to_uniprot_accession(gene)
        if not accession:
            return {}
        features, response = fetch_features(self.client, accession)
        if features is None:
            return {}

        critical = None
        disordered_intervals = []
        for feature in features:
            interval = feature_interval(feature)
            if interval is None or not overlaps(interval[0], interval[1], protein_start, total):
                continue
            feature_type = feature.get("type")
            if feature_type not in CRITICAL_FEATURE_TYPES:
                continue
            if feature_type == "Region" and is_disordered(feature):
                disordered_intervals.append(interval)
                continue
            critical = True

        # region_biologically_relevant defaults True once UniProt answered at all,
        # unless the WHOLE lost interval is covered by disordered-region annotation -
        # see the module docstring for why this default direction, not the reverse.
        covered = set()
        for start, end in disordered_intervals:
            covered.update(range(max(start, protein_start), min(end, total) + 1))
        relevant = not (disordered_intervals and covered >= set(range(protein_start, total + 1)))

        return {
            "critical_region_disrupted": True if critical else None,
            "region_biologically_relevant": relevant,
            "region_relevance_method": UNIPROT_METHOD,
            "region_relevance_source": "UniProt",
            "region_relevance_source_version": accession,
            "region_relevance_evidence_id": (
                f"https://rest.uniprot.org/uniprotkb/{accession}.json#{protein_start}-{total}"
            ),
        }

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
        relevance = self._uniprot_relevance(gene, protein_start, total)
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
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": self.release,
            "ensembl_transcript": ensembl_transcript,
            "protein_start": protein_start,
            "policy_note": (
                "lost_residues is total - protein_start + 1, counting the termination codon "
                "and everything downstream. critical_region_disrupted/"
                "region_biologically_relevant, when present, are a flagged UniProt-derived "
                "first-pass answer (region_relevance_method), not a curator's own review - "
                "see this module's docstring."
            ),
            "response_sha256": response["body_sha256"],
            **relevance,
        }]


__all__ = ["METHOD", "TRUNCATING", "ProteinRegionProvider"]
