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

  (VEP's own VCF output packs everything into one CSQ field rather than writing flat ones;
  both shapes are handled - see below.)

[Three refusals, each answering one line of that decision]
  No version, no evidence - unless a deployment declares one. Every INFO field must be
  pinned to a release: by Source and Version on its own ##INFO line, or by the file's
  ##source together with a date. A field that neither pins is skipped and the providers
  fetch it instead.

  A real file may still carry none of that - an annotator that writes "annotated with VEP"
  and no release. Rather than either refusing such a file outright or quietly inventing a
  version for it, a deployment can declare one: which source it ran, at which version, who
  says so and why. The declaration names the ##source it applies to, so it cannot attach
  itself to a different file, and every record made under it carries version_status
  "DECLARED" with the declaration beside it. The same shape the computational calibrations
  use for a predictor whose release Ensembl does not report.

  Already-decided answers are never read. CLNSIG, ACMG_CODES and their kin are a curator's
  conclusion, and feeding a conclusion back as its own support is the circularity
  test_clingen_positive.py exists to catch. They are refused by name, and the mapping is a
  whitelist besides, so a field nobody mapped is ignored rather than guessed at.

  Uncalibrated scores stay uncalibrated. A predictor score is emitted with the version its
  header names; computational.py then accepts it only if a calibration matches that version.
  Nothing here makes an unusable score usable.

[Two shapes of INFO, and why the compound one is safe to read]
  bcftools annotate writes one field per value - gnomAD_AF, gnomAD_AC - and a mapping names
  each. VEP and SnpEff instead pack every subfield into one value, CSQ or ANN, pipe
  separated, and declare the order in the ##INFO header's own "Format:" clause. That header
  is what makes the compound form readable without guessing: the file states which position
  holds the consequence and which holds the transcript, so nothing is inferred from the data.

  A compound field also carries one entry per transcript, and picking among them is the
  question MANE Select answers elsewhere. It is not answered here. Only the entry whose
  feature matches the transcript being evaluated is read; no match, more than one match, or
  no transcript to match against all yield nothing.

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
from urllib.parse import unquote

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
# VEP writes "... Format: Allele|Consequence|...", SnpEff "Functional annotations: 'Allele |
# Annotation | ...'". Both state the order; this finds whichever clause the file used.
_FORMAT = re.compile(r"(?:Format|Functional annotations)\s*:\s*(.+)$", re.IGNORECASE)


def subfield_order(definition):
    """The subfield names a compound INFO field declares, in order, or None."""
    match = _FORMAT.search((definition.description or "") if definition else "")
    if not match:
        return None
    names = [name.strip().strip("'\"").strip() for name in match.group(1).split("|")]
    return [name for name in names if name] or None


def compound_entries(value, order):
    """One dict per transcript entry, keyed by the declared subfield names."""
    entries = []
    for chunk in str(value).split(","):
        parts = chunk.split("|")
        if len(parts) < 2:
            continue
        entries.append({name: parts[index].strip() if index < len(parts) else ""
                        for index, name in enumerate(order)})
    return entries


def _accession(value):
    return str(value).split(".")[0] if value else None


def declared_version(declaration, meta, group):
    """(version, declaration) a deployment states for `group`, or (None, None).

    A declaration is only usable when it names the ##source it applies to, says who declared
    it and why, and states a version for this particular group. Anything less would be a
    version appearing from nowhere.
    """
    if not isinstance(declaration, dict):
        return None, None
    required = ("matches_source", "declared_by", "justification")
    if any(not declaration.get(field) for field in required):
        return None, None
    file_source = (meta.get("source") or "").strip()
    if declaration["matches_source"] not in file_source:
        return None, None
    version = (declaration.get("declared_versions") or {}).get(group)
    if not version:
        return None, None
    return version, {
        "matches_source": declaration["matches_source"],
        "declared_by": declaration["declared_by"],
        "justification": declaration["justification"],
        "declared_at": declaration.get("declared_at"),
        "file_source": file_source,
    }


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

    def __init__(self, mapping, version_declaration=None):
        """`mapping` is {INFO id: {"category": ..., "field": ...}} from configuration.

        Policy lives in the configuration, not here: which INFO field means which evidence
        field is a property of the annotation pipeline that wrote the file.
        `version_declaration` is the deployment's statement of the releases behind a file
        whose headers do not name them - see declared_version().
        """
        if not isinstance(mapping, dict):
            raise ValueError("VCF INFO mapping must be an object")
        refused = REFUSED & set(mapping)
        if refused:
            raise ValueError(
                f"These INFO fields carry a curator's conclusion and cannot be mapped to "
                f"evidence: {sorted(refused)}")
        self.mapping = mapping
        self.version_declaration = version_declaration

    def get_evidence(self, variant, parsed, transcript=None):
        """Every record the file supports, grouped and provenanced per category.

        `transcript` is the accession being evaluated. A compound field carries one entry per
        transcript, and choosing among them is not this provider's question, so without it
        such a field is skipped.
        """
        info = parsed.record.info
        digest = None
        grouped: dict[str, dict] = {}
        skipped = []
        for info_id, rule in self.mapping.items():
            if info_id not in info:
                continue
            category = rule.get("category")
            compound = rule.get("compound")
            field = rule.get("field")
            if category not in {CATEGORY_POPULATION, CATEGORY_ANNOTATION,
                                CATEGORY_COMPUTATIONAL} or not (field or compound):
                skipped.append({"info": info_id, "reason": "UNSUPPORTED_MAPPING"})
                continue
            group = rule.get("group")
            declaration = None
            provenance = field_provenance(parsed.info_defs.get(info_id), parsed.meta)
            if provenance is None:
                version, declaration = declared_version(
                    self.version_declaration, parsed.meta, group)
                if not version:
                    # No release to pin it to, so a provider fetches it instead.
                    skipped.append({"info": info_id, "reason": "NO_SOURCE_VERSION"})
                    continue
                source = group
            else:
                source, version = provenance
            key = (category, group or source)
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
                if declaration:
                    # The version did not come from the file; the result says so.
                    record["version_status"] = "DECLARED"
                    record["version_declaration"] = declaration
            if compound:
                values, reason = self._compound_values(
                    info[info_id], parsed.info_defs.get(info_id), compound, transcript)
                if reason:
                    skipped.append({"info": info_id, "reason": reason})
                    if not record["info_fields"]:
                        del grouped[key]
                    continue
                record.update(values)
                record["compound_format"] = "declared_in_header"
            else:
                record[field] = _decode(info[info_id])
            record["info_fields"].append(info_id)
        return list(grouped.values()), skipped

    @staticmethod
    def _compound_values(value, definition, compound, transcript):
        """(evidence fields, skip reason) for the entry matching `transcript`."""
        order = subfield_order(definition)
        if not order:
            # Without the header's Format clause the positions would be guesswork.
            return None, "NO_DECLARED_FORMAT"
        feature = compound.get("transcript_subfield")
        fields = compound.get("fields") or {}
        if not feature or not fields:
            return None, "UNSUPPORTED_MAPPING"
        if not transcript:
            return None, "NO_TRANSCRIPT_TO_MATCH"
        entries = [entry for entry in compound_entries(value, order)
                   if _accession(entry.get(feature)) == _accession(transcript)]
        if len(entries) != 1:
            # No entry for this transcript, or several - choosing is not this question.
            return None, "NO_SINGLE_TRANSCRIPT_ENTRY"
        entry = entries[0]
        separator = compound.get("list_separator") or "&"
        lists = set(compound.get("list_fields") or [])
        values = {}
        for subfield, target in fields.items():
            raw = _decode(entry.get(subfield, ""))
            if raw == "":
                continue
            values[target] = ([part for part in raw.split(separator) if part]
                              if target in lists else raw)
        return values, None


def _decode(value):
    """VCF percent-encodes the characters its own delimiters use, so "%2C" is a comma."""
    return unquote(value) if isinstance(value, str) else value


def _slug(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "vcf"
