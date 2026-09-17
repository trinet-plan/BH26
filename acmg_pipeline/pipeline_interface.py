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
from typing import Any, Optional

from acmg_pipeline.classification import CriterionEvidence, Strength, classify
from acmg_pipeline.clinical_note import extract_clinical_note
from acmg_pipeline.constants import ALL_ACMG_CODES, LITERATURE_CODES, CriterionStatus
from acmg_pipeline.gate import ERepoClient
from acmg_pipeline.pipeline import evaluate_selected_criteria, evaluate_variant_evidence_lines
from acmg_pipeline.vcf_record import VariantRecord, parse_vcf


@dataclass
class PipelineOutput:
    classification: Any
    evidence_lines: dict[str, dict]
    # clinical_noteから抽出された自由文の診断名。VA-Spec Statementの
    # objectConditionに使う(acmg_pipeline/va_spec_statement.py参照)。
    # extract_clinical_note()が未実装のため現状は常にNone。
    diagnosis: Optional[str] = None

def load_automated_config() -> dict:
    """Load server-owned rule settings; callers cannot override them."""
    path = Path(__file__).resolve().parents[1] / "config" / "demo-rules.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _assessment(line: dict) -> dict:
    for extension in line.get("extensions", []):
        if extension.get("name") == "bh26AssessmentDetails":
            return extension.get("value") or {}
    return {}


# The inverse of acmg_pipeline.automated_va_spec.STRENGTHS (this project's internal
# Strength enum -> GA4GH's own ACMG coding string, e.g. "very_strong" -> "very strong",
# "stand_alone" -> "standalone"). Written out explicitly rather than inverting that dict:
# STRENGTHS has two internal keys ("stand_alone" and "standalone") mapping to the same
# GA4GH string, so a naive {v: k for k, v in STRENGTHS.items()} silently picks whichever
# entry iterates last - found 2026-09-17 when a real PVS1 "very strong" line (the GA4GH
# space-separated form) reached Strength("very strong") directly and raised ValueError,
# since acmg_pipeline.classification.Strength.VERY_STRONG's own value is "very_strong".
_STRENGTH_FROM_GA4GH_CODE = {
    "standalone": Strength.STAND_ALONE,
    "very strong": Strength.VERY_STRONG,
    "strong": Strength.STRONG,
    "moderate": Strength.MODERATE,
    "supporting": Strength.SUPPORTING,
}


def _evidence_from_line(code: str, line: dict) -> CriterionEvidence:
    details = _assessment(line)
    status = CriterionStatus(details.get("status", CriterionStatus.UNKNOWN.value))
    if status != CriterionStatus.MET:
        return CriterionEvidence(code=code, status=status, source="integrated_pipeline")
    ga4gh_code = line.get("strengthOfEvidenceProvided", {}).get("primaryCoding", {}).get("code")
    strength = _STRENGTH_FROM_GA4GH_CODE.get(ga4gh_code)
    if strength is None:
        raise ValueError(f"Unrecognized GA4GH ACMG strength code for {code}: {ga4gh_code!r}")
    return CriterionEvidence(code=code, status=status, strength=strength, source="integrated_pipeline")


async def run_pipeline(
    variant: VariantRecord,
    clinical_note: str,
    *,
    mcp,
    erepo_client: ERepoClient,
    vcep_name: str | None = None,
    full_text_cache=None,
    llm_cache=None,
) -> PipelineOutput:
    """Evaluate 19 implemented and 9 UNKNOWN stub criteria in ACMG order."""
    extraction = extract_clinical_note(clinical_note)
    lines = await evaluate_variant_evidence_lines(
        variant, extraction,
        automated_config=load_automated_config(),
        mcp=mcp, erepo_client=erepo_client, vcep_name=vcep_name,
        full_text_cache=full_text_cache, llm_cache=llm_cache,
    )
    evidence_lines = {code: line for code, line in zip(ALL_ACMG_CODES, lines)}
    return PipelineOutput(
        classification=classify([_evidence_from_line(code, evidence_lines[code]) for code in ALL_ACMG_CODES]),
        evidence_lines=evidence_lines,
        diagnosis=extraction.diagnosis,
    )


def needs_literature_workflow(criteria: tuple[str, ...]) -> bool:
    """True if any of `criteria` requires the LLM/PubMed literature workflow.

    The API layer uses this to route a /v1/get_evidence_line_by_target_criteria
    request: False means it
    can answer synchronously (run_selected_criteria() with mcp/erepo_client
    left None); True means it must go through the async job pattern.
    """
    return bool(set(criteria) & LITERATURE_CODES)


async def run_selected_criteria(
    variant: VariantRecord,
    clinical_note: str,
    criteria: tuple[str, ...],
    *,
    mcp=None,
    erepo_client: ERepoClient | None = None,
    vcep_name: str | None = None,
    full_text_cache=None,
    llm_cache=None,
) -> dict[str, dict]:
    """Evaluate only `criteria`, returning {code: VA-Spec EvidenceLine}.

    The API-facing counterpart of run_pipeline(): where run_pipeline()
    always evaluates and classifies all 28 codes, this answers "does this
    variant meet criterion X" for an arbitrary caller-chosen subset, without
    running (or classifying against) the other 27.
    """
    extraction = extract_clinical_note(clinical_note)
    return await evaluate_selected_criteria(
        variant, extraction, criteria,
        automated_config=load_automated_config(),
        mcp=mcp, erepo_client=erepo_client,
        vcep_name=vcep_name, full_text_cache=full_text_cache, llm_cache=llm_cache,
    )


def parse_request_vcf(vcf: str) -> VariantRecord:
    return parse_vcf(vcf).record
