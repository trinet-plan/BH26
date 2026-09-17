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

from dataclasses import dataclass

from acmg_pipeline.api_input import ApiCaseInput
from acmg_pipeline.classification import ClassificationResult, classify
from acmg_pipeline.criteria.stubs import all_stub_evidence


@dataclass
class PipelineOutput:
    classification: ClassificationResult
    evidence_lines: dict[str, dict]

# TODO: mock関数
def run_pipeline(case: ApiCaseInput) -> PipelineOutput:
    case.parse_vcf()  # 1バリアント契約の検証(不正な場合の ValueError はそのまま呼び出し元に伝播する)
    return PipelineOutput(classification=classify(all_stub_evidence()), evidence_lines={})
