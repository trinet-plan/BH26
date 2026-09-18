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
from acmg_pipeline.automated_core.models import CriterionResult
from acmg_pipeline.classification import ALL_ACMG_CODES
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.common import MatchStatus, PaperContribution, VariantMatchingResult
from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.criteria import segregation
from acmg_pipeline.export import build_automated_evidence_line, build_evidence_line
from acmg_pipeline.pipeline import evaluate_variant_evidence_lines, judge_variant_from_shared_input
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
                evidence_resolver=StaticEvidenceResolver([]),
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
        assert extensions["referenceLink"] == f"https://example.test/{code}"
        assert "reportedIn" not in line


def test_integrated_entrypoint_requires_shared_input_classes():
    with pytest.raises(TypeError, match="VariantRecord"):
        asyncio.run(
            evaluate_variant_evidence_lines(
                {},
                ClinicalNoteExtraction(),
                automated_config={},
                evidence_resolver=StaticEvidenceResolver([]),
                mcp=object(),
                erepo_client=object(),
            )
        )


def test_literature_search_uses_clinical_note_diagnosis_not_vcf_disease():
    variant = VariantRecord(
        chrom="1", pos=100, id="case-1", ref="A", alt="G", qual=".", filter="PASS",
        info={
            "GENE": "GENE1", "HGVSC": "c.1A>G", "HGVSP": "p.Lys1Arg",
            "DISEASE_ASSOCIATION": "wrong_vcf_disease",
        },
    )
    clinical_note = ClinicalNoteExtraction(
        diagnosis="right clinical diagnosis", condition_id="MONDO:0005045")
    mcp = object()
    erepo = object()

    with patch(
        "acmg_pipeline.pipeline.resolve_pmids_for_variant",
        new=AsyncMock(return_value=([], "pubmed_search")),
    ) as resolve_pmids:
        result = asyncio.run(judge_variant_from_shared_input(
            variant, clinical_note, mcp, erepo, criteria=("PS3",),
        ))

    assert result == {}
    resolve_pmids.assert_awaited_once_with(
        mcp, erepo, "GENE1", "c.1A>G", "p.Lys1Arg", "right clinical diagnosis")


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


def test_unknown_automated_line_keeps_review_points_as_curator_hints():
    variant = VariantRecord(
        chrom="1", pos=100, id="case-pm1", ref="A", alt="G",
        qual=".", filter="PASS", info={"GENE": "GENE1", "HGVSC": "c.1A>G"},
    )
    result = CriterionResult(
        "PM1", CriterionStatus.UNKNOWN,
        {"assembly": "GRCh38", "chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
        "Reviewed region unavailable",
        missing_inputs=["region"],
        review_points=["Curate a disease-relevant functional region"],
    )
    with (
        patch(
            "acmg_pipeline.export.reference_links.reference_url_for_criterion",
            return_value=None,
        ),
        patch.dict("acmg_pipeline.export._CODE_CURATOR_INFO_FETCHERS", {}, clear=True),
    ):
        line = build_automated_evidence_line(result, variant)

    extensions = {item["name"]: item["value"] for item in line["extensions"]}
    assert "reviewPoints" not in extensions["bh26AssessmentDetails"]
    # An implemented-but-inconclusive line reports not_met, not unknown (see
    # export.build_automated_evidence_line()'s note, 2026-09-18) - the real
    # review_points-derived hint is joined by a second one disclosing that.
    assert extensions["bh26AssessmentDetails"]["status"] == "not_met"
    assert {
        "severity": "caution",
        "category": "review",
        "message": "Curate a disease-relevant functional region",
    } in extensions["curatorHints"]
    assert any(h["category"] == "unevaluated" for h in extensions["curatorHints"])
