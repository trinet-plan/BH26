"""
va_spec_statement.py
doc/docker_api_deployment_plan_v1_ja.md 5章の変換ロジック。
PipelineOutput(pipeline_interface.py)から、GET /v1/variant/{job_id} が返す
va_specフィールドの中身(1つのJSON文書)を組み立てる。

acmg_pipeline/export.py が最初hand-builtなdict -> 後に実ga4gh.va_specモデルへ
移行した前例(export.pyのモジュールdocstring参照)に倣い、まずは手組みのdictとして
実装する。EvidenceLine単位の詳細な変換は既にexport.build_evidence_line()が担って
いるので、ここでは PipelineOutput.evidence_lines(その出力を基準コード別に集めたdict)
と ClassificationResult をひとつのJSON文書にまとめるだけ。
"""

from __future__ import annotations

from datetime import datetime, timezone

from acmg_pipeline.pipeline_interface import PipelineOutput


def _safe_id_fragment(hgvsc: str) -> str:
    # acmg_pipeline/export.py の build_evidence_line() と同じ変換規則(idに使える文字だけにする)
    return hgvsc.replace(">", "_").replace(".", "_").replace("+", "p").replace("-", "m")


def build_variant_statement(output: PipelineOutput, gene: str, hgvsc: str, hgvsp: str = "") -> dict:
    result = output.classification
    evidence_codes = set(output.evidence_lines)

    return {
        "id": f"stmt:{gene}_{_safe_id_fragment(hgvsc)}",
        "subjectVariant": {"gene": gene, "hgvsc": hgvsc, "hgvsp": hgvsp},
        "classification": {
            "category": result.category.value,
            "score": result.score,
            "ba1Override": result.ba1_override,
        },
        "hasEvidenceLines": list(output.evidence_lines.values()),
        "criteriaNotEvaluated": result.not_evaluated_codes,
        # classification.met/not_met のうち、詳細なEvidenceLineが無いコード(=stub、
        # またはLLM以外の手段で判定されたコード)だけの簡易表現(5章参照)。
        "criteriaMet": [
            {"code": e.code, "strength": e.strength.value}
            for e in result.met if e.code not in evidence_codes
        ],
        "criteriaNotMet": [
            {"code": e.code} for e in result.not_met if e.code not in evidence_codes
        ],
        "generatedAt": datetime.now(timezone.utc).isoformat(),
    }
