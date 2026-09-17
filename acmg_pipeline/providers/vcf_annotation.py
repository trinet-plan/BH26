"""Evidence read from an annotated VCF, for the fields that arrive with provenance.

[Why this exists, against a decision that said it could not]
  services/resolve.py records that evidence cannot be read off the INFO column, and for the
  demo VCFs that is exactly right: they carry identity and context, already-decided answers,
  and scores with no version. The reasoning was never about VCF as a format, though - it was
  about what those particular files contain. A VCF written by VEP or bcftools annotate
  carries the same population frequencies and consequences the providers fetch, and its
  headers name the tool and release that produced them. Where they do, the evidence contract
  is satisfiable from the file.

  So this reads a VCF under the same rules a provider is held to, and drops on the floor
  anything that cannot meet them. It does not relax the contract; it declines the fields
  that do not fit it.

[Three refusals, each answering one line of that decision]
  No version, no evidence. Every INFO field must be declared with a Source and Version in
  its ##INFO header, or name a header the file states elsewhere. A field whose release
  cannot be pinned is skipped, and the providers fetch it instead.

  Already-decided answers are never read. CLNSIG, ACMG_CODES and their kin are a curator's
  conclusion, and feeding a conclusion back as its own support is the circularity
  test_clingen_positive.py exists to catch. They are refused by name, and the mapping is a
  whitelist besides, so a field nobody mapped is ignored rather than guessed at.

  Uncalibrated scores stay uncalibrated. A predictor score is emitted with the version its
  header names; computational.py then accepts it only if a calibration matches that version.
  Nothing here makes an unusable score usable.

[The one thing a fetched record has that this cannot]
  A provider queried a source; a VCF was handed to us by whoever submitted the variant. That
  is a real difference in what the evidence is worth, and it does not show in the values. So
  every record carries origin="vcf_input", and the run manifest records the file's checksum.
  A curator reading the result can tell which half of the evidence was independently
  retrieved.
"""

from __future__ import annotations

import hashlib
import re

CATEGORY_POPULATION = "population"
CATEGORY_ANNOTATION = "annotation"
CATEGORY_COMPUTATIONAL = "computational"

# A curator's conclusion about this variant. Never evidence for it, whatever a mapping says.
REFUSED = frozenset({
    "CLNSIG", "CLNSIGCONF", "CLNREVSTAT", "CLNVARIATIONID", "CLNDN",
    "ACMG_CODES", "ACMG_CLASS", "DISEASE_ASSOCIATION", "NOTE",
})

_SOURCE = re.compile(r'Source="([^"]*)"')
_VERSION = re.compile(r'Version="([^"]*)"')


def field_provenance(definition, meta):
    """(source, version) for one INFO field, or None when it cannot be pinned.

    The VCF spec lets an ##INFO line carry Source and Version, and an annotating tool that
    means its output to be usable writes them. Where the line is silent, the file's own
    ##source and a dated header are accepted as a fallback - a file that says what produced
    it and when is pinned, even if each field does not repeat it.
    """
    description = (definition.description or "") if definition else ""
    source = _SOURCE.search(description)
    version = _VERSION.search(description)
    if source and version:
        return source.group(1).strip(), version.group(1).strip()
    file_source = (meta.get("source") or "").strip()
    file_version = (meta.get("fileDate") or meta.get("reference") or "").strip()
    if file_source and file_version:
        return file_source, file_version
    return None


class VcfEvidenceProvider:
    """Builds normalized evidence from INFO fields a mapping names and the file can pin."""

    name = "Annotated VCF INFO"

    def __init__(self, mapping):
        """`mapping` is {INFO id: {"category": ..., "field": ...}} from configuration.

        Policy lives in the configuration, not here: which INFO field means which evidence
        field is a property of the annotation pipeline that wrote the file.
        """
        if not isinstance(mapping, dict):
            raise ValueError("VCF INFO mapping must be an object")
        refused = REFUSED & set(mapping)
        if refused:
            raise ValueError(
                f"These INFO fields carry a curator's conclusion and cannot be mapped to "
                f"evidence: {sorted(refused)}")
        self.mapping = mapping

    def get_evidence(self, variant, parsed):
        """Every record the file supports, grouped and provenanced per category."""
        info = parsed.record.info
        digest = None
        grouped: dict[str, dict] = {}
        skipped = []
        for info_id, rule in self.mapping.items():
            if info_id not in info:
                continue
            category = rule.get("category")
            field = rule.get("field")
            if category not in {CATEGORY_POPULATION, CATEGORY_ANNOTATION,
                                CATEGORY_COMPUTATIONAL} or not field:
                skipped.append({"info": info_id, "reason": "UNSUPPORTED_MAPPING"})
                continue
            provenance = field_provenance(parsed.info_defs.get(info_id), parsed.meta)
            if provenance is None:
                # No release to pin it to, so a provider fetches it instead.
                skipped.append({"info": info_id, "reason": "NO_SOURCE_VERSION"})
                continue
            source, version = provenance
            key = (category, rule.get("group") or source)
            record = grouped.get(key)
            if record is None:
                if digest is None:
                    digest = hashlib.sha256(
                        "\n".join(f"{k}={v}" for k, v in sorted(info.items()))
                        .encode("utf-8")).hexdigest()
                record = grouped[key] = {
                    "category": category,
                    "variant_key": variant.key,
                    "evidence_id": (f"urn:sha256:{digest}:vcf:{category}:"
                                    f"{_slug(key[1])}:{variant.key}"),
                    "source": source,
                    "source_version": version,
                    "retrieved_at": (parsed.meta.get("fileDate") or "").strip() or None,
                    "quality_status": "PASS",
                    # A provider queried a source; this was handed to us with the variant.
                    "origin": "vcf_input",
                    "assessment_method": "automated",
                    "method": "annotated_vcf_info",
                    "policy_version": version,
                    "info_fields": [],
                    "response_sha256": digest,
                }
            record[field] = info[info_id]
            record["info_fields"].append(info_id)
        return list(grouped.values()), skipped


def _slug(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "vcf"
