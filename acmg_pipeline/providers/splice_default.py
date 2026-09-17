"""SP01/SP02's default first-pass answer for canonical splice donor/acceptor
variants - a literature base-rate default, not a per-variant sequence
computation.

[Why not compute splice_outcome/reading_frame_disrupted from exon length]
  acmg_pipeline.criteria.pvs1's own SP01 investigation (doc/
  external_data_coverage_ja.md, 2026-09-17) tried exactly that - the same
  naive rule AutoPVS1 uses - against the one real splice case this project
  has a known answer for: MYBPC3 c.2905+1G>A (ClinGen expert panel: PVS1
  very_strong). The naive rule computes the skipped exon's length modulo 3,
  finds it in-frame, and would leave PVS1 unmet - the OPPOSITE of the real
  answer. A splice donor/acceptor loss is not reliably a clean exon skip:
  intron retention (usually frameshift) and cryptic splice site use are
  common alternative outcomes sequence length alone cannot distinguish
  (this is exactly why ClinGen's 2023 splicing recommendations lean on RNA
  evidence rather than a sequence-only rule).

[What this provides instead]
  A literature base-rate default for canonical splice donor/acceptor
  variants specifically (this provider is never invoked for anything less
  certain, like a splice region variant further from the junction):
  splice_outcome="OUT_OF_FRAME_PTC", reading_frame_disrupted=True,
  alternative_rescue=False - i.e. assume the empirically more common
  loss-of-function outcome (frameshift/intron retention) rather than guess
  which specific outcome this variant produces. This routes PVS1 through
  _truncating_path(), which reaches MET/very_strong on real truncating
  variants (see tests/test_pvs1_automated_gates.py) - matching MYBPC3
  c.2905+1G>A's real classification for the one case this can be checked
  against.

[Always flagged, never presented as a curator's own review]
  Every field here is a default assumption, not a variant-specific
  measurement or a curator's judgment - acmg_pipeline.criteria.pvs1 adds an
  explicit review_points caveat whenever this provider (not a curator)
  answered SP01/SP02, naming exactly that: alternative_rescue in
  particular has no derivation at all (it would need curator knowledge of
  tissue-specific isoform usage) and is defaulted to the value that lets
  PVS1 continue, precisely the direction this project's establishing
  policy (2026-09-17) requires disclosing rather than silently assuming.
"""

from __future__ import annotations

METHOD = "clingen_splicing_2023_canonical_default"


class SpliceDefaultProvider:
    """One default splice_assessment record per canonical splice variant."""

    name = "BH26 canonical splice default policy"

    def __init__(self, policy_version):
        if not policy_version:
            raise ValueError("policy_version must be recorded")
        self.policy_version = policy_version

    def get_splice_assessment(self, variant, transcript, generated_at):
        if not transcript or not generated_at:
            return []
        return [{
            "category": "splice_assessment",
            "variant_key": variant.key,
            "evidence_id": f"urn:bh26:splice-default:{variant.key}:{transcript}",
            "source": self.name,
            "source_version": self.policy_version,
            "retrieved_at": generated_at,
            "quality_status": "PASS",
            "transcript": transcript,
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": self.policy_version,
            "alternative_rescue": False,
            "splice_outcome": "OUT_OF_FRAME_PTC",
            "reading_frame_disrupted": True,
            "policy_note": (
                "Default assumption for a canonical splice donor/acceptor variant, not a "
                "variant-specific measurement: assumes the empirically more common "
                "frameshift/intron-retention outcome rather than computing an exon-skip "
                "in-frame/out-of-frame prediction from sequence, which was checked against "
                "MYBPC3 c.2905+1G>A (real classification: PVS1_very_strong) and found to "
                "give the wrong answer there. alternative_rescue=False (no tissue-specific "
                "alternative-splicing rescue) is likewise a default, not a finding - there is "
                "no automated way to check it."
            ),
        }]


__all__ = ["METHOD", "SpliceDefaultProvider"]
