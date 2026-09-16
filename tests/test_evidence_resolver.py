"""The evidence path that replaced the `normalized_evidence` argument.

Covers the two rules that made provider resolution the only acceptable
source (see acmg/services/resolve.py): every record must arrive with real
provenance, and a variant absent from a provider must not become a
zero-frequency observation.
"""

import os

import pytest

os.environ.setdefault("VLLM_BASE_URL", "http://127.0.0.1:8000/v1")
os.environ.setdefault("VLLM_API_KEY", "test-only")

from acmg.core.models import Variant
from acmg.services.resolve import ProviderEvidenceResolver, StaticEvidenceResolver
from acmg_pipeline.pipeline import _automated_variant, _identity_from_info
from acmg_pipeline.vcf_record import VariantRecord

VARIANT = Variant(assembly="GRCh38", chrom="14", pos=23883114, ref="G", alt="A")


def _vcf_record(**extra_info) -> VariantRecord:
    return VariantRecord(
        chrom="14", pos=23883114, id="case-1", ref="G", alt="A",
        qual=".", filter="PASS",
        info={
            "ASSEMBLY": "GRCh38",
            "GENE": "MYH7",
            "TRANSCRIPT": "NM_000257.4",
            "HGVSC": "c.1594T>C",
            "HGVSP": "p.Ser532Pro",
            **extra_info,
        },
    )


def test_identity_is_an_allowlist_not_the_whole_info_column():
    """CLNSIG/ACMG_CODES are conclusions; nothing may reach them via identity."""
    identity = _identity_from_info(_vcf_record(
        CLNSIG="Likely_Pathogenic",
        ACMG_CODES="PVS1_VeryStrong,PM2_Moderate",
        AM_PATHOGENICITY=0.98,
        DISEASE_ASSOCIATION="hypertrophic_cardiomyopathy",
    ))
    assert identity == {
        "GENE": "MYH7",
        "TRANSCRIPT": "NM_000257.4",
        "HGVSC": "c.1594T>C",
        "HGVSP": "p.Ser532Pro",
    }
    assert "CLNSIG" not in identity
    assert "ACMG_CODES" not in identity
    assert "AM_PATHOGENICITY" not in identity


def test_identity_keeps_the_clinvar_accession_when_present():
    identity = _identity_from_info(_vcf_record(CLNVARIATIONID="VCV000042013"))
    assert identity["CLNVARIATIONID"] == "VCV000042013"


def test_identity_drops_empty_and_is_case_insensitive():
    identity = _identity_from_info(_vcf_record(**{"gene": "", "clnvariationid": "VCV000000001"}))
    assert identity["GENE"] == "MYH7"  # the populated spelling wins, "" is dropped
    assert identity["CLNVARIATIONID"] == "VCV000000001"


def test_automated_variant_defaults_to_grch38_and_rejects_others():
    assert _automated_variant(_vcf_record()).key == "GRCh38:14:23883114:G:A"
    record = _vcf_record()
    record.info.pop("ASSEMBLY")
    assert _automated_variant(record).key == "GRCh38:14:23883114:G:A"
    record.info["ASSEMBLY"] = "GRCh37"
    with pytest.raises(ValueError, match="GRCh38"):
        _automated_variant(record)


def _resolver(**providers) -> ProviderEvidenceResolver:
    """A resolver whose HTTP clients are never used because every provider is stubbed."""
    resolver = ProviderEvidenceResolver.__new__(ProviderEvidenceResolver)
    resolver._client = object()
    resolver._external = object()
    resolver._ensembl_release = "116"
    resolver._gnomad_release = "4.1.1"
    resolver._clinvar_release = "2026-09-15"
    resolver._ensembl = object()
    for name, value in providers.items():
        setattr(resolver, name, value)
    return resolver


def test_absent_from_gnomad_contributes_no_population_record(monkeypatch):
    """None means 'not observed', which must never become AF=0."""
    import acmg.services.resolve as resolve

    class _Absent:
        name = "gnomAD"

        def __init__(self, *args, **kwargs):
            pass

        def get_frequency(self, variant):
            return None

    monkeypatch.setattr(resolve, "GnomadProvider", _Absent)
    resolver = _resolver()
    resolved = resolve.ResolvedEvidence()
    resolver._add_population(VARIANT, resolved)

    assert resolved.records == []
    assert resolved.failures == []  # absence is not a failure


def test_provider_error_is_recorded_not_silently_empty(monkeypatch):
    """An unreachable provider must be distinguishable from 'no such evidence'."""
    import acmg.services.resolve as resolve

    class _Broken:
        name = "gnomAD"

        def __init__(self, *args, **kwargs):
            pass

        def get_frequency(self, variant):
            raise ValueError("gnomAD query failed")

    monkeypatch.setattr(resolve, "GnomadProvider", _Broken)
    resolver = _resolver()
    resolved = resolve.ResolvedEvidence()
    resolver._add_population(VARIANT, resolved)

    assert resolved.records == []
    assert resolved.failures == [{"provider": "gnomAD", "error": "gnomAD query failed"}]


def test_annotation_for_a_different_variant_is_rejected():
    """A transcript HGVS resolving elsewhere must not annotate this variant."""
    import acmg.services.resolve as resolve

    class _Elsewhere:
        def map_record_with_evidence(self, record):
            return None, {"variant_key": "GRCh38:14:99999999:C:T",
                          "transcript": "NM_000257.4",
                          "consequences": ["missense_variant"]}, []

    resolver = _resolver()
    resolver._ensembl = _Elsewhere()
    resolved = resolve.ResolvedEvidence()
    annotation, predictions = resolver._annotate(
        {"TRANSCRIPT": "NM_000257.4", "HGVSC": "c.1594T>C"}, VARIANT, resolved
    )

    assert annotation is None and predictions == []
    assert resolved.failures[0]["provider"] == "Ensembl"
    assert "GRCh38:14:23883114:G:A" in resolved.failures[0]["error"]


def test_static_resolver_round_trips_records_and_rejects_bad_input():
    records = [{"category": "population", "evidence_id": "https://example.test/x"}]
    assert StaticEvidenceResolver(records).resolve({}, VARIANT).records == records
    with pytest.raises(TypeError, match="list of dictionaries"):
        StaticEvidenceResolver({"not": "a list"})
