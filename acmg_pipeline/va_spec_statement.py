"""
va_spec_statement.py
doc/docker_api_deployment_plan_v1_ja.md 5章の変換ロジック。
GET /v1/classify_criteria/{job_id} が返す va_spec フィールドの中身(1つのJSON文書)を
組み立てる。PipelineOutput(pipeline_interface.py)から、GA4GH VA-Spec
`acmg_2015` プロファイルの `VariantPathogenicityStatement` を構築する。

[v2, 2026-09-17: 手組みdictから実 ga4gh.va_spec モデルへ移行]
  acmg_pipeline/export.py がEvidenceLine側で先に実施した移行(そのモジュールdocstring
  参照)に倣い、ここでも手組みdictをやめ、`ga4gh.va_spec.acmg_2015.models.
  VariantPathogenicityStatement` を実際に構築して `.model_dump()` する。
  設計判断の経緯は z_tmp_va_spec_statement_decisions.md 参照:
    - `proposition.subjectVariant` は VRS Allele化せず、`iriReference`として
      hgvsc文字列をそのまま使う(subjectVariantの型は
      `MolecularVariation | CategoricalVariant | iriReference` のUnionで、
      iriReferenceの実体はただの文字列のため)。
    - `proposition.objectCondition` は clinical_note から抽出された診断名
      (PipelineOutput.diagnosis)を `MappableConcept.name` にそのまま載せる
      (MedGen/OMIM等のコード化は行わない)。抽出未実装の間はClinVar自身の
      "not_provided" 慣習にならい "not provided" をプレースホルダーとする。
    - `strength` は `definitive` / `likely` の2値のみ(EvidenceLine側の
      supporting/moderate/strong/very_strongとは別語彙)。
  VA-Spec本体が持たないフィールド(Tavtigianスコア以外の内訳: 実装外基準の
  met/not_met/not_evaluatedの一覧)は `extensions` の `bh26ClassificationSummary`
  に退避する(実モデルは `BaseModelForbidExtra` で未知フィールドを許さないため)。
"""

from __future__ import annotations

from datetime import datetime, timezone

from ga4gh.core.models import Coding, Extension, MappableConcept, iriReference
from ga4gh.va_spec.acmg_2015.models import METHOD, AcmgClassification, VariantPathogenicityStatement
from ga4gh.va_spec.base.core import Contribution, Direction, VariantPathogenicityProposition
from ga4gh.va_spec.base.enums import System

from acmg_pipeline.classification import ClassificationCategory
from acmg_pipeline.export import PIPELINE_AGENT
from acmg_pipeline.pipeline_interface import PipelineOutput

# ClassificationResult.category(Title Case) -> (direction, strengthコード, AcmgClassification)。
# strengthコードがNone のときは strength フィールド自体を省略する
# (Uncertain Significance は neutral 方向で確信度も主張しない)。
_CLASSIFICATION_TO_STATEMENT = {
    ClassificationCategory.PATHOGENIC: (Direction.SUPPORTS, "definitive", AcmgClassification.PATHOGENIC),
    ClassificationCategory.LIKELY_PATHOGENIC: (Direction.SUPPORTS, "likely", AcmgClassification.LIKELY_PATHOGENIC),
    ClassificationCategory.UNCERTAIN_SIGNIFICANCE: (Direction.NEUTRAL, None, AcmgClassification.UNCERTAIN_SIGNIFICANCE),
    ClassificationCategory.LIKELY_BENIGN: (Direction.DISPUTES, "likely", AcmgClassification.LIKELY_BENIGN),
    ClassificationCategory.BENIGN: (Direction.DISPUTES, "definitive", AcmgClassification.BENIGN),
}


def _safe_id_fragment(text: str) -> str:
    # acmg_pipeline/export.py の evidence_line_id() と同じ変換規則(idに使える文字だけにする)。
    # subjectVariant(HGVS文字列そのもの)には適用しない - HGVS構文を壊してしまうため。
    return text.replace(">", "_").replace(".", "_").replace("+", "p").replace("-", "m").replace(":", "_")


def build_variant_statement(output: PipelineOutput, gene: str, hgvsc: str, hgvsp: str = "") -> dict:
    result = output.classification
    evidence_codes = set(output.evidence_lines)
    direction, strength_code, classification_value = _CLASSIFICATION_TO_STATEMENT[result.category]

    variant_key = hgvsc or f"gene:{gene or 'unknown-gene'}"

    proposition = VariantPathogenicityProposition(
        subjectVariant=iriReference(root=variant_key),
        objectCondition=MappableConcept(
            name=output.diagnosis or "not provided", conceptType="Disease"
        ),
        geneContextQualifier=MappableConcept(name=gene) if gene else None,
    )

    strength = (
        MappableConcept(primaryCoding=Coding(code=strength_code, system=System.ACMG.value))
        if strength_code
        else None
    )
    classification = MappableConcept(
        primaryCoding=Coding(code=classification_value.value, system=System.ACMG.value)
    )

    summary_extension = Extension(
        name="bh26ClassificationSummary",
        value={
            "ba1Override": result.ba1_override,
            "criteriaNotEvaluated": result.not_evaluated_codes,
            # classification.met/not_met のうち、詳細なEvidenceLineが無いコード
            # (=stub、またはLLM以外の手段で判定されたコード)だけの簡易表現。
            "criteriaMet": [
                {"code": e.code, "strength": e.strength.value}
                for e in result.met
                if e.code not in evidence_codes
            ],
            "criteriaNotMet": [
                {"code": e.code} for e in result.not_met if e.code not in evidence_codes
            ],
        },
    )

    statement = VariantPathogenicityStatement(
        id=f"stmt:{gene}_{_safe_id_fragment(variant_key)}",
        proposition=proposition,
        direction=direction,
        strength=strength,
        classification=classification,
        score=result.score,
        specifiedBy=METHOD,
        contributions=[
            Contribution(
                contributor=PIPELINE_AGENT,
                activityType="evaluated",
                date=datetime.now(timezone.utc),
            )
        ],
        hasEvidenceLines=list(output.evidence_lines.values()) or None,
        extensions=[summary_extension],
    )
    return statement.model_dump(mode="json", exclude_none=True)
