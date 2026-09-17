"""Convert VCF files into the shared criterion input classes.

This adapter belongs to the integration boundary.  Criterion modules receive
``VariantRecord`` and ``ClinicalNoteExtraction`` only; they do not parse VCF.
"""

from __future__ import annotations

import re
from pathlib import Path

from acmg_pipeline.clinical_note import (
    ClinicalNoteExtraction,
    Proband,
    ProbandGenotype,
)
from acmg_pipeline.vcf_record import InfoFieldDef, ParsedVcf, VariantRecord


_INFO_DEF = re.compile(r"^##INFO=<(?P<body>.*)>$")


def _definition_fields(body: str) -> dict[str, str]:
    """Parse comma-separated VCF structured metadata, respecting quoted commas."""
    fields = []
    start = 0
    quoted = False
    escaped = False
    for index, character in enumerate(body):
        if escaped:
            escaped = False
        elif character == "\\" and quoted:
            escaped = True
        elif character == '"':
            quoted = not quoted
        elif character == "," and not quoted:
            fields.append(body[start:index])
            start = index + 1
    if quoted:
        raise ValueError("Unterminated quote in INFO definition")
    fields.append(body[start:])
    result: dict[str, str] = {}
    for field in fields:
        key, separator, value = field.partition("=")
        if not separator or not key or key in result:
            raise ValueError(f"Invalid INFO definition field: {field}")
        if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
            value = value[1:-1].replace(r'\"', '"').replace(r"\\", "\\")
        result[key] = value
    return result


def _convert_scalar(value: str, value_type: str):
    if value == ".":
        return None
    if value_type == "Integer":
        return int(value)
    if value_type == "Float":
        return float(value)
    if value_type == "Flag":
        raise ValueError("Flag INFO values must not use '='")
    if value_type == "Character" and len(value) != 1:
        raise ValueError(f"Character INFO value must have length 1: {value}")
    return value


def _convert_info_value(value: str | bool, definition: InfoFieldDef | None):
    if definition is None:
        return value
    if definition.type == "Flag":
        if value is not True:
            raise ValueError(f"Flag INFO field {definition.id} must not have a value")
        return True
    if value is True:
        raise ValueError(f"INFO field {definition.id} requires a value")
    values = value.split(",") if definition.number != "1" else [value]
    converted = [_convert_scalar(item, definition.type) for item in values]
    return converted[0] if definition.number == "1" else converted


def _parse_info(text: str, definitions: dict[str, InfoFieldDef]) -> dict:
    result = {}
    if text == ".":
        return result
    for item in text.split(";"):
        key, separator, value = item.partition("=")
        if not key or key in result:
            raise ValueError(f"Invalid or duplicate INFO field: {key}")
        raw_value: str | bool = value if separator else True
        result[key] = _convert_info_value(raw_value, definitions.get(key))
    return result


def _select_alt_info(info: dict, definitions: dict[str, InfoFieldDef], alt_index: int, alt_count: int):
    """Select allele-specific INFO values when a multi-ALT row is decomposed."""
    selected = dict(info)
    for key, value in info.items():
        definition = definitions.get(key)
        if definition is None or not isinstance(value, list):
            continue
        if definition.number == "A":
            if len(value) != alt_count:
                raise ValueError(
                    f"INFO/{key} Number=A has {len(value)} values for {alt_count} ALT alleles"
                )
            selected[key] = value[alt_index]
        elif definition.number == "R":
            if len(value) != alt_count + 1:
                raise ValueError(
                    f"INFO/{key} Number=R has {len(value)} values for {alt_count} ALT alleles"
                )
            selected[key] = [value[0], value[alt_index + 1]]
    return selected


def parse_vcf(path, *, assembly: str | None = None) -> list[ParsedVcf]:
    """Parse every VCF row and decompose every ALT into a shared ``VariantRecord``.

    The adapter handles the eight fixed VCF columns. Sample/FORMAT columns are
    intentionally outside ``VariantRecord`` and must be represented in the
    separately supplied ``ClinicalNoteExtraction`` when clinically relevant.
    """
    path = Path(path)
    meta: dict[str, str] = {}
    info_defs: dict[str, InfoFieldDef] = {}
    data_lines: list[tuple[int, str]] = []
    header_seen = False

    for line_number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line:
            continue
        match = _INFO_DEF.match(line)
        if match:
            fields = _definition_fields(match.group("body"))
            required = {"ID", "Number", "Type"}
            if not required.issubset(fields):
                raise ValueError(f"Incomplete INFO definition at {path}:{line_number}")
            definition = InfoFieldDef(
                id=fields["ID"],
                number=fields["Number"],
                type=fields["Type"],
                description=fields.get("Description", ""),
            )
            if definition.id in info_defs:
                raise ValueError(f"Duplicate INFO definition {definition.id} at {path}:{line_number}")
            info_defs[definition.id] = definition
            continue
        if line.startswith("##"):
            key, separator, value = line[2:].partition("=")
            if separator:
                meta[key] = value
            continue
        if line.startswith("#CHROM\t"):
            header_seen = True
            continue
        if line.startswith("#"):
            continue
        data_lines.append((line_number, line))

    if not header_seen:
        raise ValueError(f"VCF column header is missing: {path}")

    declared_assembly = meta.get("reference")
    if assembly and declared_assembly and assembly != declared_assembly:
        raise ValueError(
            f"Requested assembly {assembly} conflicts with VCF reference {declared_assembly}"
        )
    resolved_assembly = assembly or declared_assembly
    parsed: list[ParsedVcf] = []
    for line_number, line in data_lines:
        fields = line.split("\t")
        if len(fields) < 8:
            raise ValueError(f"VCF record has fewer than 8 columns at {path}:{line_number}")
        try:
            position = int(fields[1])
        except ValueError as exc:
            raise ValueError(f"Invalid POS at {path}:{line_number}: {fields[1]}") from exc
        alts = fields[4].split(",")
        info = _parse_info(fields[7], info_defs)
        for alt_index, alt in enumerate(alts):
            record_info = _select_alt_info(info, info_defs, alt_index, len(alts))
            if resolved_assembly:
                record_info["ASSEMBLY"] = resolved_assembly
            parsed.append(ParsedVcf(
                meta=dict(meta),
                info_defs=dict(info_defs),
                record=VariantRecord(
                    chrom=fields[0],
                    pos=position,
                    id=fields[2],
                    ref=fields[3],
                    alt=alt,
                    qual=fields[5],
                    filter=fields[6],
                    info=record_info,
                ),
            ))
    return parsed


def _clinical_note_from_vcf(record: VariantRecord) -> ClinicalNoteExtraction:
    """Create only the clinical fields explicitly represented by the demo VCF."""
    zygosity = record.info.get("ZYGOSITY")
    if not isinstance(zygosity, str):
        return ClinicalNoteExtraction()
    return ClinicalNoteExtraction(
        proband=Proband(genotype=ProbandGenotype(variant_status=True, zygosity=zygosity))
    )


def inputs_from_vcf(
    path,
    *,
    clinical_note: ClinicalNoteExtraction | None = None,
    assembly: str | None = None,
) -> list[tuple[VariantRecord, ClinicalNoteExtraction]]:
    """Convert a VCF into the two-object input accepted by every criterion.

    ``clinical_note`` applies to all records in this VCF (the current demo uses
    one case per file).  When omitted, only INFO/ZYGOSITY is mapped; no
    phenotype, family history, or de-novo fact is inferred.
    """
    if clinical_note is not None and not isinstance(clinical_note, ClinicalNoteExtraction):
        raise TypeError("clinical_note must be a ClinicalNoteExtraction")
    return [
        (item.record, clinical_note or _clinical_note_from_vcf(item.record))
        for item in parse_vcf(path, assembly=assembly)
    ]
