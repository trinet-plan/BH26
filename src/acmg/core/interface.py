"""Adapters between shared team inputs and the internal evidence model."""

from dataclasses import asdict

from acmg.core.models import Variant
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


CONTEXT_FIELDS = (
    "record_id",
    "source",
    "transcript",
    "identity_provenance",
    "condition",
    "condition_label",
    "inheritance",
    "disease_frequency_threshold",
    "ba1_exception_assessment",
)


def _info_value(variant: VariantRecord, key: str, default=None):
    wanted = key.casefold()
    return next(
        (value for name, value in variant.info.items() if name.casefold() == wanted),
        default,
    )


def inputs_from_prepared_record(record: dict) -> tuple[VariantRecord, ClinicalNoteExtraction]:
    """Turn one validated prepared-JSON record into the shared input classes."""
    raw_variant = record["variant"]
    info = {
        "ASSEMBLY": raw_variant["assembly"],
        **{key: record[key] for key in CONTEXT_FIELDS if key in record},
    }
    clinical_data = record.get("clinical_note") or {}
    if not isinstance(clinical_data, dict):
        raise ValueError("clinical_note must be an object")
    variant = VariantRecord(
        chrom=str(raw_variant["chrom"]),
        pos=raw_variant["pos"],
        id=str(record.get("record_id", "")),
        ref=raw_variant["ref"],
        alt=raw_variant["alt"],
        qual=str(record.get("qual", ".")),
        filter=str(record.get("filter", "PASS")),
        info=info,
    )
    return variant, ClinicalNoteExtraction.from_json(clinical_data)


def criterion_input(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
) -> dict:
    """Build the allowlisted internal view consumed by existing decision logic."""
    if not isinstance(variant, VariantRecord):
        raise TypeError("variant must be a VariantRecord")
    if not isinstance(clinical_note, ClinicalNoteExtraction):
        raise TypeError("clinical_note must be a ClinicalNoteExtraction")
    normalized = Variant(
        # The shared VariantRecord contract has no assembly member. This
        # evaluator is GRCh38-only, so an omitted INFO/ASSEMBLY means GRCh38;
        # any explicitly supplied different assembly is still rejected.
        assembly=str(_info_value(variant, "ASSEMBLY", "GRCh38")),
        chrom=variant.chrom,
        pos=variant.pos,
        ref=variant.ref,
        alt=variant.alt,
    )
    data = {"variant": normalized.to_dict(), "clinical_note": asdict(clinical_note)}
    for key in CONTEXT_FIELDS:
        value = _info_value(variant, key)
        if value is not None:
            data[key] = value
    if "inheritance" not in data and clinical_note.family.inheritance_pattern:
        data["inheritance"] = clinical_note.family.inheritance_pattern
    return data
