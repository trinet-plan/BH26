"""
vcf_record.py
Classes for the embedded VCF text found in a
democase/case*_api_input_*.json file's "vcf" field (e.g.
democase/case1_api_input_case1-noise2.json): that field is a *string*
containing header meta-lines + one column-header line + exactly 1 data
line (per this project's own "1 variant per request" contract - doc/
docker_api_deployment_plan_v1_ja.md section 0-1), not structured JSON on
its own, so parse_vcf() turns it into:
  - `meta`: every ##key=value header line (fileformat, source, reference,
    the free-text ##triage_expansion/##caveat provenance notes, etc.)
  - `info_defs`: every ##INFO=<ID=...,Number=...,Type=...,Description=...>
    declaration, used to type-cast the record's INFO values (Float for
    AM_PATHOGENICITY, etc.) instead of leaving everything as raw strings
  - `record`: the single VariantRecord, with core VCF columns
    (chrom/pos/id/ref/alt/qual/filter) plus a typed `info` dict

parse_vcf() raises ValueError if the text doesn't contain exactly 1 data
line - there is no other caller of this module that reads multi-row VCF
text, so it isn't built to handle that case.

This only has to read this already-annotated demo text back into Python
objects - producing/annotating the VCF itself (real VEP output, plugin
config, etc.) is a separate, external concern this module doesn't touch.

[Undeclared INFO keys - a real gap in the demo data itself]
  Cross-checking the "declared via ##INFO" list against what actually shows
  up in data rows found GNOMAD_AF and TAVTIGIAN_POINTS used in
  case1_api_input_case1-noise2.json / case1-var1.json / case1-var2.json
  without a matching ##INFO=<ID=GNOMAD_AF,...> / ##INFO=<ID=TAVTIGIAN_
  POINTS,...> declaration anywhere in that same file's header. Rather than
  raising on this, an undeclared key falls back to best-effort numeric
  sniffing (see _coerce_unknown) so GNOMAD_AF=0.0130 still comes back as a
  float instead of the raw string "0.0130".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

_CORE_COLUMNS = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")


@dataclass
class InfoFieldDef:
    """One ##INFO=<ID=...,Number=...,Type=...,Description="..."> declaration."""
    id: str
    number: str  # "1", "0" (Flag), "." (variable/list), or a small integer as a string
    type: str  # String / Integer / Float / Character / Flag
    description: str = ""


@dataclass
class VariantRecord:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: Optional[float]
    filter: str
    info: dict[str, Any] = field(default_factory=dict)


@dataclass
class ParsedVcf:
    meta: dict[str, str]
    info_defs: dict[str, InfoFieldDef]
    record: VariantRecord


def _parse_info_def_line(line: str) -> InfoFieldDef:
    """Parses a ##INFO=<ID=...,Number=...,Type=...,Description="..."> line.

    The angle-bracket body is a comma-separated key=value list, except
    Description's value is a quoted string that may itself contain commas -
    so Description (always declared last per the VCF spec, and last in
    every line this project's demo data actually uses) is split off first
    by locating `Description="` and taking everything up to the final `"`,
    and only the remainder is split on plain commas.
    """
    body = line.split("=<", 1)[1].rstrip(">\n")
    desc_marker = 'Description="'
    desc_idx = body.find(desc_marker)
    description = ""
    if desc_idx != -1:
        description = body[desc_idx + len(desc_marker):].rstrip('"')
        body = body[:desc_idx].rstrip(",")

    fields: dict[str, str] = {}
    for part in body.split(","):
        if "=" not in part:
            continue
        key, _, value = part.partition("=")
        fields[key.strip()] = value.strip()

    return InfoFieldDef(
        id=fields["ID"],
        number=fields.get("Number", "1"),
        type=fields.get("Type", "String"),
        description=description,
    )


def _coerce_unknown(raw: str) -> Any:
    """Best-effort typing for an INFO key with no ##INFO declaration (see
    module docstring: GNOMAD_AF/TAVTIGIAN_POINTS in the current demo data).
    Tries int, then float, then falls back to the raw string."""
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        return raw


def _cast_scalar(raw: str, type_: str) -> Any:
    if type_ == "Integer":
        return int(raw)
    if type_ == "Float":
        return float(raw)
    return raw  # String / Character, and anything unrecognized


def _cast_info_value(raw: str, info_def: Optional[InfoFieldDef]) -> Any:
    if info_def is None:
        return _coerce_unknown(raw)
    if info_def.type == "Flag":
        return True
    if info_def.number == "1":
        return _cast_scalar(raw, info_def.type)
    # Number="." or a fixed count >1: a comma-separated list of values.
    # NOTE/ACMG_CODES both declare Number="." - see the module docstring on
    # why splitting on comma is safe for this project's demo data (free-text
    # fields were written without embedded commas specifically to keep this
    # simple split-based parsing correct).
    return [_cast_scalar(v, info_def.type) for v in raw.split(",")]


def _parse_info_string(raw_info: str, info_defs: dict[str, InfoFieldDef]) -> dict[str, Any]:
    info: dict[str, Any] = {}
    for entry in raw_info.split(";"):
        if not entry:
            continue
        if "=" not in entry:
            # Bare flag (e.g. "SOMEFLAG" with no "="): present -> True.
            info[entry] = True
            continue
        key, _, value = entry.partition("=")
        info[key] = _cast_info_value(value, info_defs.get(key))
    return info


def parse_vcf(text: str) -> ParsedVcf:
    """Parses one embedded mini-VCF string (the "vcf" field of a
    democase/case*_api_input_*.json file) into meta/info_defs/record.
    Raises ValueError if no #CHROM column-header line is found (a truncated
    or non-VCF input), or if the data line count isn't exactly 1 (this
    project's api_input.json contract - see the module docstring)."""
    meta: dict[str, str] = {}
    info_defs: dict[str, InfoFieldDef] = {}
    records: list[VariantRecord] = []
    saw_column_header = False

    for line in text.splitlines():
        if not line:
            continue
        if line.startswith("##INFO=<"):
            info_def = _parse_info_def_line(line)
            info_defs[info_def.id] = info_def
        elif line.startswith("##"):
            key, _, value = line[2:].partition("=")
            meta[key] = value
        elif line.startswith("#CHROM"):
            saw_column_header = True
        else:
            fields = line.split("\t")
            if len(fields) < len(_CORE_COLUMNS) + 1:
                raise ValueError(f"Expected {len(_CORE_COLUMNS) + 1} tab-separated columns, got {len(fields)}: {line!r}")
            chrom, pos, vid, ref, alt, qual, filter_, raw_info = fields[:8]
            records.append(VariantRecord(
                chrom=chrom,
                pos=int(pos),
                id=vid,
                ref=ref,
                alt=alt,
                qual=None if qual == "." else float(qual),
                filter=filter_,
                info=_parse_info_string(raw_info, info_defs),
            ))

    if not saw_column_header:
        raise ValueError("No #CHROM column-header line found; this does not look like a VCF")
    if len(records) != 1:
        raise ValueError(f"Expected exactly 1 VCF data line, got {len(records)}")

    return ParsedVcf(meta=meta, info_defs=info_defs, record=records[0])
