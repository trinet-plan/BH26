"""Helpers for reading the shared criterion input objects."""

from __future__ import annotations

from typing import Any

from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def info_value(variant: VariantRecord, key: str, default: Any = "") -> Any:
    """Read an INFO field case-insensitively."""
    wanted = key.upper()
    for name, value in variant.info.items():
        if name.upper() == wanted:
            return value
    return default


def variant_identity(variant: VariantRecord) -> tuple[str, str, str, list[str]]:
    """Return ``gene, hgvsc, hgvsp, equivalents`` from a VariantRecord."""
    gene = str(info_value(variant, "GENE"))
    hgvsc = str(info_value(variant, "HGVSC"))
    hgvsp = str(info_value(variant, "HGVSP", "N/A") or "N/A")
    raw_equivalents = info_value(variant, "EQUIVALENTS", [])
    if isinstance(raw_equivalents, str):
        equivalents = [raw_equivalents] if raw_equivalents else []
    else:
        equivalents = [str(value) for value in raw_equivalents]
    return gene, hgvsc, hgvsp, equivalents


def empty_clinical_note() -> ClinicalNoteExtraction:
    """Explicit empty note for literature-only callers."""
    return ClinicalNoteExtraction()
