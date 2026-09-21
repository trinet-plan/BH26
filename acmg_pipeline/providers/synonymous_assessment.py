"""BP7's `synonymous_assessment` evidence, derived from data already fetched for
every annotated variant - no new HTTP request of its own.

[Why no new fetch is needed]
  VEP's own consequence classification already flags a synonymous variant within 1-3
  bases of an exon/intron boundary as splice_region_variant (Sequence Ontology's own
  definition, which VEP implements exactly), and the VEP REST SpliceAI plugin
  (ensembl.py's VEP_OPTIONS, "SpliceAI=2") is already requested for every variant this
  project annotates, regardless of consequence type - see
  EnsemblIdentityProvider._prediction_evidence(). BP7 only needs to read what map_
  record_with_evidence() already returned alongside the annotation.

[The two conditions this answers - ClinGen SVI 2023, PMID 37352859]
  "A synonymous (silent) variant ... for which SpliceAI predicts no impact on splicing
  (score <0.1)", outside the splice-critical window. No conservation requirement - see
  acmg_pipeline/criteria/bp7.py's own comment for why a prior version of this project's
  BP7 rule required one and no longer does. The position half of that sentence is read
  off VEP's own SO-term classification: splice_region_variant co-occurring with
  synonymous_variant IS the "within the splice-critical window" case, since VEP assigns
  it to exactly the 1-3 exonic / 3-8 intronic bases nearest the boundary.

[contradictory_rna_evidence]
  Always False here: this project has no automated RNA/splicing-assay evidence source.
  Same honest-default convention as providers/splice_default.py's alternative_rescue
  and providers/region_repeat.py's repetitive=False - a curator's own reviewed record
  still wins, since curated_context() prefers any record naming a curator and review
  date over an assessment_method="automated" one only when both exist; a real
  contradiction found later replaces this record rather than fighting it.
"""

from __future__ import annotations

import hashlib

SPLICEAI_NO_IMPACT_THRESHOLD = 0.1
METHOD = "clingen_svi_2023_bp7_spliceai_position"
CALIBRATION_SOURCE = "ClinGen SVI Splicing Subgroup 2023 (PMID:37352859): SpliceAI score <0.1"
POSITION_RULE = "ensembl_vep_splice_region_consequence"


def get_synonymous_assessment(variant, annotation, predictions, policy_version):
    """One `synonymous_assessment` record when a SpliceAI score was resolved for this
    transcript, else none - an honest gap, not a guessed score."""
    if not policy_version:
        raise ValueError("policy_version must be recorded")
    if annotation is None:
        return []
    transcript = annotation.get("transcript")
    consequences = set(annotation.get("consequences") or [])
    if not transcript or "synonymous_variant" not in consequences:
        return []
    splice = next(
        (item for item in predictions
         if item.get("predictor") == "SpliceAI" and item.get("transcript") == transcript),
        None,
    )
    if splice is None or splice.get("score") is None:
        return []
    digest = hashlib.sha256(
        f"{splice.get('evidence_id')}:{variant.key}".encode("utf-8")).hexdigest()
    return [{
        "category": "synonymous_assessment",
        "variant_key": variant.key,
        "evidence_id": f"urn:sha256:{digest}:bp7-synonymous-assessment:{transcript}",
        "source": "BH26 BP7 synonymous assessment (derived from VEP consequence + SpliceAI)",
        "source_version": splice.get("source_version") or splice.get("predictor_version"),
        "retrieved_at": splice["retrieved_at"],
        "quality_status": "PASS",
        "transcript": transcript,
        "assessment_method": "automated",
        "method": METHOD,
        "policy_version": policy_version,
        "outside_splice_critical_region": "splice_region_variant" not in consequences,
        "no_predicted_splice_impact": float(splice["score"]) < SPLICEAI_NO_IMPACT_THRESHOLD,
        "contradictory_rna_evidence": False,
        "splice_prediction_evidence": splice.get("evidence_id") or splice.get("source"),
        "calibration_source": CALIBRATION_SOURCE,
        "position_rule_version": f"{POSITION_RULE}:{policy_version}",
        "spliceai_score": splice["score"],
        "policy_note": (
            "outside_splice_critical_region is read off VEP's own splice_region_variant "
            "co-annotation, no_predicted_splice_impact off the SpliceAI max delta score "
            "already fetched for every variant - see this module's docstring for the "
            "ClinGen SVI 2023 recommendation (PMID:37352859) this implements. "
            "contradictory_rna_evidence is always False here: no automated RNA evidence "
            "source exists; a curator's own reviewed synonymous_assessment record "
            "overrides this one whenever it exists."
        ),
    }]


__all__ = ["METHOD", "SPLICEAI_NO_IMPACT_THRESHOLD", "get_synonymous_assessment"]
