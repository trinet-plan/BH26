"""Shared VCF record types used at criterion boundaries."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class InfoFieldDef:
    """One ``##INFO=<...>`` declaration from a VCF header."""

    id: str
    number: str
    type: str
    description: str = ""


@dataclass
class VariantRecord:
    """One normalized VCF data row."""

    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: str
    filter: str
    info: dict[str, Any] = field(default_factory=dict)


@dataclass
class ParsedVcf:
    meta: dict[str, str]
    info_defs: dict[str, InfoFieldDef]
    record: VariantRecord
