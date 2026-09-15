# PM1設計 — hotspotとcritical functional domainの2ルート

`docs/PLAN.md` 工程4の詳細化。Phase 1〜3は実装済み、Phase 4は次工程。

## 定義への対応

ACMG PM1は "mutational hot spot **and/or** critical and well-established functional
domain, without benign variation" であり、成立経路は2つある。共通条件はbenign variationの
欠如だけで、pathogenic variantの集積はhotspotルートの条件である。

```
                         PM1
           ┌──────────────┴──────────────┐
     Hotspot route                Critical-domain route
  mutational_hotspot           critical_functional_region
           AND                           AND
    benign_depletion               benign_depletion
           └──────────────┬──────────────┘
                         MET
```

初期実装は3項目すべてを要求していたため、機能的重要性が独立に確立されたactive siteでも
報告症例が少なければNOT_METとなり、定義より厳しかった。またhotspotは
`pathogenic_enrichment` + `benign_depletion` から導かれるため、`critical_region` に
hotspotを含めたうえで両方を要求すると同じEvidenceを二重に要求する構造になっていた。

## Evidence契約

`category="region"`、variantごとに1件。共通のtransport要件（`variant_key`、`evidence_id`、
`source`、`source_version`、`retrieved_at`、`quality_status="PASS"`、transcript一致）に加え、
`region_type` でルートを宣言する。protein座標は `protein_id` 一致と
`start <= protein_start <= protein_end <= end` の完全包含を要求する（部分重複はMANUAL_REVIEW）。

由来は2種類あり、`assessment_method` で区別する。

| | 人手 | 自動 |
| --- | --- | --- |
| 受理条件 | `curator` + `reviewed_at` | `assessment_method="automated"` + `method` + `policy_version` |
| 使えるルート | 両方 | **hotspotのみ** |

自動Evidenceでcritical functional domainを主張した場合はNOT_EVALUATED
（`missing: curated_criticality`）。報告密度から機能的重要性を導く循環を構造的に禁止する。

### hotspotルート

`pathogenic_count` / `benign_count`（非負int）と、configと一致する `method` /
`policy_version` を要求する。判定はEvidence側のbooleanではなくcountsとpolicyから導出し、
Evidenceが矛盾するbooleanを宣言していればMANUAL_REVIEW。

- `benign_count > max_benign` → **NOT_MET**（窓内のbenign variationは実証的な反証）
- `pathogenic_count < min_pathogenic` → **NOT_EVALUATED**（報告不足であり、hotspotでないことの証明ではない）
- 両方満たす → **MET**（moderate）

よく検査される遺伝子ほど閾値に到達しやすいため、この2つを同じ結論にまとめない。

### critical-domainルート

`critical_functional_region` と `benign_depletion` のbooleanを要求する。
`pathogenic_enrichment` は要求しない。

## 疾患コンテキスト

PM1は領域そのものに関するprotein-levelの主張なので、`condition` なしでも評価する。
PS1と同じく、判定と疾患関連性を分離して記録する。

- `condition` あり かつ region評価のcondition一致 → `condition_assessment: MATCHED`
- それ以外 → `condition_assessment: NOT_EVALUATED`。MET時は確認事項を併記する

`condition` を指定した場合、別疾患のregion評価は `EvidenceService` の段階で除外され、
流用されない。

## 閾値policy

閾値はコードに持たず、rules JSONの `PM1.hotspot` にversion付きで置く。
`window_aa` / `min_pathogenic` / `max_benign` / `method` / `policy_source` /
`policy_version` が揃わなければNOT_EVALUATED。

`config/demo-rules.json` の `PM1-hotspot-v1`（±5 aa、P/LP>=3、B/LB=0）はfastVEP型の
hotspot proxyであり、ACMG/ClinGenの普遍基準ではない。ClinGen VCEPはPM1をgene/disease固有に
定義し、明確なhotspotを定義できない遺伝子ではPM1を使用しない仕様もある。

## Provider

`ClinVarHotspotProvider`（`src/acmg/providers/clinvar.py`）。遺伝子単位のesearch
（missense限定）とesummaryのバッチ取得で、残基±N aaに入る変異をP/LPとB/LBに数える。

- 数えたVCVのaccession・protein_change・分類・conditionを保存し、後のrelease で再検証できる
- conflicting分類はどちらにも数えず、監査用に残す
- 単一残基置換（例 R248W）以外は位置づけ不能として除外する
- 検索が打ち切られた場合（`count > retmax`）はdensity Evidenceを出さず、
  `quality_status="INCOMPLETE_SEARCH"` の `region_search` だけを残す
- 送信するのは遺伝子記号のみ。応答は内容ハッシュ付きキャッシュに固定し `--offline` で再生する

既知の制約として、ClinVar esummaryの `protein_change` は遺伝子レベルの表記であり、
こちらの `protein_id`/transcriptに紐づかない。遺伝子一致は確認しているが、別isoform由来の
座標が混じる可能性は残るため、`protein_change_source` に出典を明示している。

## 実行

```powershell
$env:PYTHONPATH = 'src'
.venv/Scripts/python.exe -m acmg prepare-demo-online --input-dir demo-data --cache-dir tests/fixtures/ensembl-cache --evidence-cache-dir tests/fixtures/external-cache --output-dir work/pm1 --ensembl-release 116 --with-gnomad --with-clinvar --clinvar-release 2026-09-15 --with-pm1-hotspot --rules config/demo-rules.json --offline
```

固定キャッシュでは10件のhotspot Evidenceが生成され、PM1はMET 2件、NOT_MET 6件、
NOT_EVALUATED 20件（報告密度不足10件、非missenseでregion Evidenceなし10件）。
MET 2件はconditionが未指定のため確認事項付き。

## 検証

- 2ルートの成立・不成立、critical domainがpathogenic enrichmentを要求しないこと
- benign報告ありのNOT_METと、報告不足のNOT_EVALUATEDを取り違えないこと
- policy未設定・policy不一致・counts欠落・counts矛盾の分岐
- 自動Evidenceがcritical domainを主張できないこと
- 別疾患のregion評価を流用しないこと
- 固定キャッシュによる全28件のオフライン再現と、VA-Spec Evidence Lineの検証

## 未完了（Phase 4）

1. gene/disease-specific対応。強度可変（特定residueはPM1、周辺domainはPM1_supporting）は
   region Evidenceに `strength` を持たせればVA-Spec出力まで通る。
2. PM1を使用しない遺伝子の表明。policyに除外遺伝子を置き、NOT_EVALUATEDではなく
   NOT_APPLICABLEを返して「Evidenceがない」と区別する。
3. critical functional domainのreviewed Evidence入力経路（InterPro/Pfam境界の取り込みを含む）。
   domain境界だけではPM1は成立しないため、criticalityの判断と出典が別途必要。
4. `condition` の入力経路。PM1のcondition_assessmentをMATCHEDにできるようにする。
