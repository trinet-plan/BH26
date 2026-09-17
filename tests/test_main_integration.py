import json
import os
import asyncio
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest
from ga4gh.va_spec.base.core import EvidenceLine
from jsonschema import Draft202012Validator

os.environ.setdefault("VLLM_BASE_URL", "http://127.0.0.1:8000/v1")
os.environ.setdefault("VLLM_API_KEY", "test-only")

from acmg_pipeline.services.resolve import StaticEvidenceResolver
from acmg_pipeline.classification import ALL_ACMG_CODES
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.common import MatchStatus, PaperContribution, VariantMatchingResult
from acmg_pipeline.criteria import segregation
from acmg_pipeline.export import build_evidence_line
from acmg_pipeline.pipeline import evaluate_variant_evidence_lines
from acmg_pipeline.vcf_record import VariantRecord


def test_integrated_entrypoint_returns_one_ordered_line_per_acmg_code():
    variant = VariantRecord(
        chrom="1",
        pos=100,
        id="case-1",
        ref="A",
        alt="G",
        qual=".",
        filter="PASS",
        info={"ASSEMBLY": "GRCh38", "GENE": "GENE1", "HGVSC": "c.1A>G"},
    )
    clinical_note = ClinicalNoteExtraction()
    config = json.loads((Path("config") / "demo-rules.json").read_text(encoding="utf-8"))

    with (
        patch(
            "acmg_pipeline.pipeline.judge_variant_from_shared_input",
            new=AsyncMock(return_value={}),
        ),
        patch(
            "acmg_pipeline.export.reference_links.reference_url_for_criterion",
            side_effect=lambda code, _variant: f"https://example.test/{code}",
        ),
        patch.dict("acmg_pipeline.export._CODE_CURATOR_INFO_FETCHERS", {}, clear=True),
    ):
        lines = asyncio.run(
            evaluate_variant_evidence_lines(
                variant,
                clinical_note,
                automated_config=config,
                normalized_evidence=[],
                mcp=object(),
                erepo_client=object(),
            )
        )

    assert len(lines) == 28
    assert tuple(line["specifiedBy"]["methodType"] for line in lines) == ALL_ACMG_CODES
    assert len({line["id"] for line in lines}) == 28

    schema_validator = Draft202012Validator(EvidenceLine.model_json_schema())
    for line in lines:
        EvidenceLine.model_validate(line)
        schema_validator.validate(line)
        extensions = {item["name"]: item["value"] for item in line["extensions"]}
        code = line["specifiedBy"]["methodType"]
        assert extensions["bh26AssessmentDetails"]["criterion"] == code
        assert extensions["referenceLink"] == f"https://example.test/{code}"
        assert "reportedIn" not in line


def test_integrated_entrypoint_requires_shared_input_classes():
    with pytest.raises(TypeError, match="VariantRecord"):
        asyncio.run(
            evaluate_variant_evidence_lines(
                {},
                ClinicalNoteExtraction(),
                automated_config={},
                normalized_evidence=[],
                mcp=object(),
                erepo_client=object(),
            )
        )


def test_pp1_reference_is_attached_to_real_literature_line():
    variant = VariantRecord(
        chrom="14", pos=23883114, id="case-pp1", ref="G", alt="A",
        qual=".", filter="PASS", info={"GENE": "MYH7", "HGVSC": "c.2155C>T"},
    )
    judgment = segregation.SegregationJudgment(
        variant_matching=VariantMatchingResult(match_status=MatchStatus.MATCHED),
        families=[segregation.FamilySegregationData(
            family_id="F1", affected_with_variant=3, affected_without_variant=0,
            unaffected_with_variant=0, unaffected_without_variant=2,
        )],
        overall_evidence=segregation.SegregationEvidence(
            direction=segregation.SegregationDirection.PP1,
            rationale="Variant segregates with disease.",
        ),
    )
    result = segregation.finalize(judgment, pmid="12345678")
    aggregated = segregation.aggregate_multi_paper_results([
        PaperContribution("12345678", result)
    ])

    with patch(
        "acmg_pipeline.export.reference_links.reference_url_for_criterion",
        return_value="https://example.test/PP1",
    ):
        line = build_evidence_line(
            aggregated, "MYH7", "c.2155C>T", "PP1", variant=variant,
        )

    extensions = {item["name"]: item["value"] for item in line["extensions"]}
    assert extensions["referenceLink"] == "https://example.test/PP1"
    assert line["reportedIn"][0]["pmid"] == "12345678"


def test_reference_lookup_failure_does_not_remove_criterion_line():
    variant = VariantRecord(
        chrom="14", pos=23883114, id="case-pm1", ref="G", alt="A",
        qual=".", filter="PASS", info={"GENE": "MYH7", "HGVSC": "c.2155C>T"},
    )
    with (
        patch(
            "acmg_pipeline.export.reference_links.reference_url_for_criterion",
            side_effect=RuntimeError("TogoID unavailable"),
        ),
        patch.dict("acmg_pipeline.export._CODE_CURATOR_INFO_FETCHERS", {}, clear=True),
    ):
        from acmg_pipeline.export import build_workflow_evidence_line

        line = build_workflow_evidence_line(
            "PM1", variant, status="unknown", description="Review required."
        )

    assert line["specifiedBy"]["methodType"] == "PM1"
    assert not any(item["name"] == "referenceLink" for item in line["extensions"])
