"""Reconcile independent identity candidates before accepting genomic coordinates."""

from acmg.core.models import Variant
from acmg.core.reference import normalize


def reconcile(record, candidates, reference):
    """Candidates must come from an identity adapter, not from source classification labels.

    candidate: variant, source, source_version, matched_identifiers, retrieved_at.
    All claimed identifiers are matched exactly against this record, including transcript version.
    An inconsistent source coordinate may be corrected only by corroborated identity evidence.
    """
    result = {**record, "issues": list(record["issues"])}
    result["identity_status"] = "PENDING"
    result["resolution"] = None
    if any(issue in {"ASSEMBLY_CONFLICT", "MISSING_VCF_COLUMN_HEADER", "INVALID_COLUMN_COUNT"}
           or issue.startswith("Duplicate INFO key") for issue in record["issues"]):
        result["issues"].append("INVALID_SOURCE_STRUCTURE")
        return result
    if record["raw_variant"].get("assembly") != "GRCh38":
        result["issues"].append("UNSUPPORTED_ASSEMBLY")
        return result
    identities = {k: v for k, v in record["identity"].items()
                  if k in {"TRANSCRIPT", "HGVSC", "CLNVARIATIONID"}}
    normalized_candidates = {}
    candidate_errors = []
    for candidate in candidates:
        matched = candidate.get("matched_identifiers", {})
        hgvs_matched = all(key in matched for key in ("TRANSCRIPT", "HGVSC"))
        accession_matched = bool(matched.get("CLNVARIATIONID"))
        if not (hgvs_matched or accession_matched):
            candidate_errors.append("CANDIDATE_MISSING_IDENTITY_MATCH")
            continue
        if any(identities.get(key) != value for key, value in matched.items()):
            candidate_errors.append("CANDIDATE_IDENTIFIER_MISMATCH")
            continue
        if not all(candidate.get(key) for key in ("source", "source_version", "retrieved_at")):
            candidate_errors.append("CANDIDATE_MISSING_PROVENANCE")
            continue
        try:
            variant = normalize(Variant(**candidate["variant"]), reference)
        except (ValueError, KeyError, TypeError) as exc:
            candidate_errors.append(f"INVALID_IDENTITY_CANDIDATE: {exc}")
            continue
        normalized_candidates.setdefault(variant.key, {"variant": variant, "evidence": []})[
            "evidence"].append(candidate)
    result["issues"].extend(candidate_errors)
    if candidate_errors:
        return result
    if len(normalized_candidates) > 1:
        result["issues"].append("CONFLICTING_IDENTITY_CANDIDATES")
        return result
    if not normalized_candidates:
        result["issues"].append("IDENTITY_EVIDENCE_UNAVAILABLE")
        return result
    resolved = next(iter(normalized_candidates.values()))
    corroborated = set()
    for candidate in resolved["evidence"]:
        corroborated.update(candidate["matched_identifiers"])
    # Where coding HGVS is supplied, an accession alone cannot settle a HGVS disagreement.
    if identities.get("HGVSC") and not {"TRANSCRIPT", "HGVSC"}.issubset(corroborated):
        result["issues"].append("HGVS_NOT_VERIFIED")
        return result
    original = record.get("parsed_variant")
    variant = resolved["variant"]
    same = original == variant.to_dict()
    result["identity_status"] = "VERIFIED" if same else "CORRECTED"
    result["resolution"] = {
        "variant": variant.to_dict(), "evidence": resolved["evidence"],
        "reason": "Confirmed source coordinates" if same else "Identity mapping and reference validation",
    }
    return result


def evaluation_inputs(records):
    """Allowlist crossing the audit/evaluation boundary: never pass source_record or labels."""
    return [
        {"record_id": record["record_id"], "source": record["source"],
         "variant": record["resolution"]["variant"],
         "identity_provenance": record["resolution"]["evidence"],
         "transcript": record["identity"].get("TRANSCRIPT")}
        for record in records if record["identity_status"] in {"VERIFIED", "CORRECTED"}
    ]
