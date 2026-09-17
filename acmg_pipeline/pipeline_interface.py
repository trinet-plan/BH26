"""
pipeline_interface.py
API層とパイプライン担当者(別の人)との唯一の境界(doc/docker_api_deployment_plan_v1_ja.md
4章)。run_pipeline() のシグネチャと PipelineOutput の形だけがAPI層との契約であり、
中身の実装(LLM/PubMed MCP呼び出し、ERepoゲート、基準ごとの判定)はここでは行わない。

[現状の中身について]
実際の判定ロジック(acmg_pipeline/pipeline.py の judge_variant() 等)がまだこの関数に
接続されていないため、暫定的に全28コードをNOT_EVALUATEDとして返す
(acmg_pipeline.criteria.stubs.all_stub_evidence() がすでに表現している
「未評価」をそのまま反映するだけで、判定を捏造しない)。接続する担当者は、
実装済み5コード(classification.IMPLEMENTED_CODES)について実際のCriterionEvidenceを
作り、all_stub_evidence() の代わりに classify() へ渡すよう置き換えればよい。
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from acmg_pipeline.classification import CriterionEvidence, Strength, classify
from acmg_pipeline.clinical_note import extract_clinical_note
from acmg_pipeline.constants import ALL_ACMG_CODES, CriterionStatus
from acmg_pipeline.gate import ERepoClient
from acmg_pipeline.pipeline import evaluate_variant_evidence_lines
from acmg_pipeline.vcf_record import VariantRecord, parse_vcf


@dataclass
class PipelineOutput:
    classification: Any
    evidence_lines: dict[str, dict]

def load_automated_config() -> dict:
    """Load server-owned rule settings; callers cannot override them."""
    path = Path(__file__).resolve().parents[1] / "config" / "demo-rules.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _assessment(line: dict) -> dict:
    for extension in line.get("extensions", []):
        if extension.get("name") == "bh26AssessmentDetails":
            return extension.get("value") or {}
    return {}


def _evidence_from_line(code: str, line: dict) -> CriterionEvidence:
    details = _assessment(line)
    status = CriterionStatus(details.get("status", CriterionStatus.UNKNOWN.value))
    if status != CriterionStatus.MET:
        return CriterionEvidence(code=code, status=status, source="integrated_pipeline")
    strength = line.get("strengthOfEvidenceProvided", {}).get("primaryCoding", {}).get("code")
    return CriterionEvidence(code=code, status=status, strength=Strength(strength), source="integrated_pipeline")


async def run_pipeline(
    variant: VariantRecord,
    clinical_note: str,
    *,
    mcp,
    erepo_client: ERepoClient,
    vcep_name: str | None = None,
    full_text_cache=None,
) -> PipelineOutput:
    """Evaluate 19 implemented and 9 UNKNOWN stub criteria in ACMG order."""
    extraction = extract_clinical_note(clinical_note)
    lines = await evaluate_variant_evidence_lines(
        variant, extraction,
        automated_config=load_automated_config(),
        mcp=mcp, erepo_client=erepo_client, vcep_name=vcep_name,
        full_text_cache=full_text_cache,
    )
    evidence_lines = {code: line for code, line in zip(ALL_ACMG_CODES, lines)}
    return PipelineOutput(
        classification=classify([_evidence_from_line(code, evidence_lines[code]) for code in ALL_ACMG_CODES]),
        evidence_lines=evidence_lines,
    )


def parse_request_vcf(vcf: str) -> VariantRecord:
    return parse_vcf(vcf).record
