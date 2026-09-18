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

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.services.resolve import (
    ProviderEvidenceResolver,
    ResolvedEvidence,
    StaticEvidenceResolver,
)
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
    resolver._population_sources = None
    resolver._population = None
    resolver._clinvar_release = "2026-09-15"
    resolver._ensembl = object()
    for name, value in providers.items():
        setattr(resolver, name, value)
    return resolver


def test_absent_from_togovar_contributes_no_population_record():
    """None means 'not observed', which must never become AF=0."""
    import acmg_pipeline.services.resolve as resolve

    class _Absent:
        name = "TogoVar"

        def __init__(self, *args, **kwargs):
            pass

        def get_frequency(self, variant):
            return None

    resolver = _resolver(_population=[_Absent()])
    resolved = resolve.ResolvedEvidence()
    resolver._add_population(VARIANT, resolved)

    assert resolved.records == []
    assert resolved.failures == []  # absence is not a failure


def test_provider_error_is_recorded_not_silently_empty():
    """An unreachable provider must be distinguishable from 'no such evidence'."""
    import acmg_pipeline.services.resolve as resolve

    class _Broken:
        name = "TogoVar"

        def __init__(self, *args, **kwargs):
            pass

        def get_frequency(self, variant):
            raise ValueError("TogoVar query failed")

    resolver = _resolver(_population=[_Broken()])
    resolved = resolve.ResolvedEvidence()
    resolver._add_population(VARIANT, resolved)

    assert resolved.records == []
    assert resolved.failures == [{"provider": "TogoVar", "error": "TogoVar query failed"}]


def test_independent_population_databases_are_aggregated_with_failure_isolation():
    observation = {
        "category": "population", "variant_key": VARIANT.key,
        "evidence_id": "https://other.example/evidence/1",
        "source": "OtherFrequencyDB", "source_version": "2026.1",
        "retrieved_at": "2026-09-17T00:00:00Z", "population": "global",
        "AC": 1, "AN": 10000, "AF": "0.0001", "callable": True,
        "quality_status": "PASS", "filters": [], "variant_flags": [],
    }

    class _Good:
        name = "OtherFrequencyDB"

        def get_frequency(self, variant):
            return [observation]

    class _Broken:
        name = "UnavailableFrequencyDB"

        def get_frequency(self, variant):
            raise ValueError("service unavailable")

    resolver = _resolver(_population=[_Good(), _Broken()])
    resolved = ResolvedEvidence()
    resolver._add_population(VARIANT, resolved)

    assert resolved.records == [observation]
    assert resolved.failures == [{
        "provider": "UnavailableFrequencyDB", "error": "service unavailable",
    }]


def test_annotation_for_a_different_variant_is_rejected():
    """A transcript HGVS resolving elsewhere must not annotate this variant."""
    import acmg_pipeline.services.resolve as resolve

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


# ============================================================================
# _apply_curated_context() - config/curated-context.json's BA1 exception list,
# wired into the live pipeline for the first time on 2026-09-17 (previously
# only automated_cli.py's separate batch path ever called load_context()/
# apply_context(), so BA1's exception check always came back "unavailable"
# here regardless of what curated-context.json already had resolved).
# ============================================================================

import json

from acmg_pipeline.pipeline import _apply_curated_context


def _write_curated_context(tmp_path, *, exceptions=None, complete=True):
    document = {
        "context_version": "test-1",
        "source": "test",
        "ba1_exceptions": {
            "source": "test list", "source_version": "1", "reviewed_at": "2026-09-17",
            "entry_method": "manual_transcription", "complete": complete,
            "variants": exceptions or [],
        },
        "records": {}, "record_contexts": {},
    }
    path = tmp_path / "curated-context.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def test_a_known_exception_resolves_to_is_exception_true(tmp_path):
    """HFE c.845G>A (p.Cys282Tyr) is a real entry in the project's committed
    config/curated-context.json - this uses a synthetic copy of the same shape."""
    path = _write_curated_context(tmp_path, exceptions=[
        {"variant_key": "GRCh38:6:26092913:G:A", "gene": "HFE", "hgvs_c": "c.845G>A",
         "caid": "CA113795", "resolved_by": "test"},
    ])
    variant = VariantRecord(chrom="6", pos=26092913, id="", ref="G", alt="A",
                            qual="", filter="", info={"GENE": "HFE"})
    _apply_curated_context(variant, {"curated_context_path": path})
    assert variant.info["ba1_exception_assessment"]["is_exception"] is True


def test_a_variant_absent_from_a_complete_list_resolves_to_false_not_unavailable(tmp_path):
    path = _write_curated_context(tmp_path, exceptions=[])
    variant = VariantRecord(chrom="1", pos=100, id="", ref="A", alt="G",
                            qual="", filter="", info={"GENE": "TEST"})
    _apply_curated_context(variant, {"curated_context_path": path})
    assert variant.info["ba1_exception_assessment"]["is_exception"] is False


def test_an_incomplete_list_resolves_nothing(tmp_path):
    """An incomplete list cannot say a variant is absent from it - see
    acmg_pipeline.automated_core.context.apply_context()'s own docstring."""
    path = _write_curated_context(tmp_path, exceptions=[], complete=False)
    variant = VariantRecord(chrom="1", pos=100, id="", ref="A", alt="G",
                            qual="", filter="", info={"GENE": "TEST"})
    _apply_curated_context(variant, {"curated_context_path": path})
    assert "ba1_exception_assessment" not in variant.info


def test_no_configured_path_is_a_no_op():
    variant = VariantRecord(chrom="1", pos=100, id="", ref="A", alt="G",
                            qual="", filter="", info={"GENE": "TEST"})
    _apply_curated_context(variant, {})
    assert variant.info == {"GENE": "TEST"}
