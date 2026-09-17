"""Build gene-disease draft and source-anchored assessment documents."""

import hashlib
import json

from acmg_pipeline.providers.gene_disease_draft import GeneDiseaseDraftProvider


def _canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def build_draft_document(document):
    if not isinstance(document, dict):
        raise ValueError("Gene-disease source document must be an object")
    policy = document.get("policy")
    groups = document.get("groups")
    generated_at = document.get("retrieved_at")
    if not isinstance(groups, list) or not groups:
        raise ValueError("Gene-disease source document requires nonempty groups")
    if not isinstance(generated_at, str) or not generated_at:
        raise ValueError("Gene-disease source document requires retrieved_at")
    provider = GeneDiseaseDraftProvider(policy)
    records = []
    for group in groups:
        if not isinstance(group, dict):
            raise ValueError("Every gene-disease source group must be an object")
        records.extend(provider.build(generated_at=generated_at, **group))
    ids = [item["evidence_id"] for item in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Generated gene-disease draft evidence IDs are not unique")
    return {
        "schema_version": "1.0",
        "source_version": document.get("source_version"),
        "generated_at": generated_at,
        "evidence": records,
    }


def build_assessment_document(document, base_evidence=None):
    """Expand reviewed group decisions to exact-variant evidence records.

    Automated assessments remain explicitly machine-authored. They can exercise the demo
    criteria because their method and policy are versioned, but ``human_signoff`` stays false.
    """
    if not isinstance(document, dict):
        raise ValueError("Gene-disease assessment document must be an object")
    required = ("source_version", "assessed_at", "method", "policy_version")
    if not all(isinstance(document.get(key), str) and document[key] for key in required):
        raise ValueError("Gene-disease assessment provenance is incomplete")
    groups = document.get("groups")
    if not isinstance(groups, list) or not groups:
        raise ValueError("Gene-disease assessment requires nonempty groups")
    records = []
    for group in groups:
        if not isinstance(group, dict):
            raise ValueError("Every gene-disease assessment group must be an object")
        variant_keys = group.get("variant_keys")
        if not isinstance(variant_keys, list) or not variant_keys:
            raise ValueError("Assessment group requires variant_keys")
        required_group = ("gene", "condition", "transcript", "lof_mechanism_established",
                          "pp2_applicable", "bp1_applicable")
        if any(group.get(key) is None for key in required_group):
            raise ValueError("Gene-disease assessment group is incomplete")
        if any(type(group[key]) is not bool for key in
               ("lof_mechanism_established", "pp2_applicable", "bp1_applicable")):
            raise ValueError("Gene-disease assessment applicability and LoF fields must be boolean")
        sources = group.get("sources")
        if not isinstance(sources, list) or not sources:
            raise ValueError("Gene-disease assessment group requires sources")
        decisions = {key: value for key, value in group.items() if key != "variant_keys"}
        fingerprint = hashlib.sha256(_canonical({
            "source_version": document["source_version"],
            "policy_version": document["policy_version"],
            "decisions": decisions,
        }).encode()).hexdigest()
        for variant_key in variant_keys:
            if not isinstance(variant_key, str) or len(variant_key.split(":")) != 5:
                raise ValueError("Assessment variant key must be assembly:chrom:pos:ref:alt")
            records.append({
                "category": "gene_disease",
                "variant_key": variant_key,
                "evidence_id": (
                    f"urn:bh26:gene-disease-assessment:{fingerprint}:"
                    f"{hashlib.sha256(variant_key.encode()).hexdigest()[:16]}"
                ),
                "source": "BH26 source-anchored automated gene-disease review",
                "source_version": document["source_version"],
                "retrieved_at": document["assessed_at"],
                "quality_status": "PASS",
                "assessment_method": "automated",
                "method": document["method"],
                "policy_version": document["policy_version"],
                "human_signoff": False,
                **decisions,
            })
    existing = [] if base_evidence is None else base_evidence.get("evidence")
    if not isinstance(existing, list):
        raise ValueError("Base evidence document must contain an evidence list")
    combined = [*existing, *records]
    ids = [item.get("evidence_id") for item in combined]
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("Combined evidence IDs must be present and unique")
    return {
        "schema_version": "1.0",
        "source_version": document["source_version"],
        "generated_at": document["assessed_at"],
        "evidence": combined,
    }
