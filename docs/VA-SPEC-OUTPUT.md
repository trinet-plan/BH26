# VA-Spec 1.0.1 JSON出力仕様

## 適用範囲

本プロジェクトは、GA4GH VA-Spec stable 1.0.1の
`VariantPathogenicityEvidenceLine`（ACMG 2015 profile）を出力対象とする。
VA-Spec 1.0.1が参照するGKS-Coreは1.0.0である。

公式例の全体構造は、最終的な `Statement` が `proposition` と
`hasEvidenceLines` を持ち、各 `EvidenceLine` が根拠となる `StudyResult` 等を
`hasEvidenceItems` に持つ形である。

```text
Statement
├─ proposition: VariantPathogenicityProposition
└─ hasEvidenceLines[]: EvidenceLine
   └─ hasEvidenceItems[]: StudyResult または IRI参照
      ├─ sourceDataSet: DataSet
      ├─ cohort: StudyGroup
      └─ specifiedBy: Method
         └─ reportedIn: Document
```

現行入力には疾患・表現型を表すconditionと最終病原性分類がない。そのため、値を推測して
`Statement` や `VariantPathogenicityProposition` を作らず、独立して妥当な
`EvidenceLine` を成果物とする。gnomAD集団頻度Evidenceは、公式例にならって
`CohortAlleleFrequencyStudyResult` としてEvidenceLine内に埋め込む。

## 出力ファイル

| パス | 内容 | VA-Spec成果物か |
| --- | --- | --- |
| `results.json` | 全基準の内部評価結果。未評価、対象外、要確認、非推奨も保持 | いいえ |
| `summary.tsv` | 内部評価結果の表形式要約 | いいえ |
| `evidence-lines.json` | record、variant、Evidence参照をまとめたBH26監査用envelope | いいえ |
| `va-spec-1.0.1/*.json` | 1ファイルにつき1個のACMG EvidenceLine | はい |
| `run-manifest.json` | schema ID、ローカルprofileのSHA-256、生成ファイル一覧 | いいえ |

`evidence-lines.json` は転送・監査を容易にする独自コンテナであり、それ自体を
VA-Spec objectとは称さない。相互運用時は `va-spec-1.0.1` 配下のJSONを使用する。

## EvidenceLineのフィールド対応

| VA-Specフィールド | 出力規則 |
| --- | --- |
| `id` | variant、criterion、outcomeから決定的に生成した `urn:bh26:evidence-line:...` |
| `type` | 常に `EvidenceLine` |
| `name` | criterionとGRCh38 variantを含む表示名 |
| `description` | 内部判定の要約。判定根拠の監査用説明 |
| `specifiedBy.type` | `Method` |
| `specifiedBy.methodType` | `PM2`、`BA1`等のACMG criterion code |
| `specifiedBy.reportedIn` | ACMG/AMP 2015文献（DOI、PMID、URL） |
| `hasEvidenceItems` | 構造化できるEvidenceは埋込みobject、それ以外はIRI参照 |
| `directionOfEvidenceProvided` | 病原性側MET=`supports`、良性側MET=`disputes`、NOT_MET=`neutral` |
| `strengthOfEvidenceProvided` | METの場合のみ。ACMG強度のMappableConcept |
| `evidenceOutcome` | criterionの成立・不成立を表すMappableConcept |
| `targetProposition` | 現行は省略。疾患conditionを取得後に追加する |

`MappableConcept` には `primaryCoding.system`、`primaryCoding.code` と、必要に応じて
`name` を出力する。GKS-Core 1.0.0の定義に存在しない `type: MappableConcept` は付けない。

NOT_METの内部directionは `none` だが、VA-Spec 1.0.1のmachine-readable JSON Schemaで
許可される列挙値は `supports`、`neutral`、`disputes` であるため、JSON出力時だけ
`neutral` に正規化する。NOT_METには強度を付けない。

## gnomAD集団頻度StudyResult

人口頻度Evidenceは次のように対応づける。

| StudyResultフィールド | gnomAD Evidence |
| --- | --- |
| `id` | キャッシュ内の一意な `evidence_id` |
| `type` | `CohortAlleleFrequencyStudyResult` |
| `focusAllele` | 正規化済みGRCh38 variantから生成した安定URN |
| `focusAlleleFrequency` | `AF` |
| `focusAlleleCount` | `AC` |
| `locusAlleleCount` | `AN` |
| `sourceDataSet` | gnomAD releaseを持つ `DataSet` |
| `cohort` | population名を持つ `StudyGroup` |
| `specifiedBy` | 頻度計算方法とgnomADヘルプを表す `Method` / `Document` |
| `qualityMeasures` | PASS状態、callable、filter、variant flag、取得日時 |

AFは0以上1以下、ACは0以上、ANは1以上でなければ出力を拒否する。source、release、populationが
欠ける場合も、不完全なStudyResultを作らずエラーにする。gnomADに未登録という応答はAF=0へ
変換しないため、このStudyResultにもならない。

## 出力例

以下は形を示す短縮例である。URNと値は実行対象から決定される。

```json
{
  "id": "urn:bh26:evidence-line:<sha256>",
  "type": "EvidenceLine",
  "name": "PM2 assessment for GRCh38:1:2:C:T",
  "description": "Population frequency does not satisfy PM2",
  "specifiedBy": {
    "type": "Method",
    "name": "ACMG/AMP 2015 with ClinGen General Guidance",
    "methodType": "PM2",
    "reportedIn": {
      "type": "Document",
      "name": "Richards et al., 2015, Genet Med.",
      "doi": "10.1038/gim.2015.30",
      "pmid": "25741868",
      "urls": ["https://pubmed.ncbi.nlm.nih.gov/25741868/"]
    }
  },
  "hasEvidenceItems": [
    {
      "id": "https://example.org/evidence/frequency",
      "type": "CohortAlleleFrequencyStudyResult",
      "focusAllele": "urn:bh26:variant:<sha256>",
      "focusAlleleFrequency": 0.00001,
      "focusAlleleCount": 1,
      "locusAlleleCount": 100000,
      "sourceDataSet": {
        "type": "DataSet",
        "name": "gnomAD v4.1.1",
        "version": "4.1.1"
      },
      "cohort": {
        "type": "StudyGroup",
        "name": "exome:global"
      }
    }
  ],
  "directionOfEvidenceProvided": "neutral",
  "evidenceOutcome": {
    "name": "ACMG 2015 PM2 criterion not met",
    "primaryCoding": {
      "system": "ACMG Guidelines, 2015",
      "code": "PM2_not_met"
    }
  }
}
```

## 検証

各EvidenceLineはファイル書込み前に次の順で検証する。

1. `ga4gh.va-spec` の `VariantPathogenicityEvidenceLine` modelで基本構造を正規化する。
2. VA-Spec 1.0.1/GKS-Core 1.0.0に固定した自己完結型の制限schemaで検証する。
3. criterion、outcome、direction、strengthのACMG cross-field整合性を検証する。
4. run manifestに公式canonical schema IDと制限schemaのSHA-256を記録する。

成立時の出力は、次の専用合成fixtureでE2E検証する。

| record | 入力 | 成立結果 | VA-Spec方向・強度・outcome |
| --- | --- | --- | --- |
| `fixture:pm2-met` | AC=0、AN=10000、AF=0 | PM2 MET | `supports` / `supporting` / `PM2_supporting` |
| `fixture:ba1-met` | AC=501、AN=10000、AF=0.0501、BA1例外でない | BA1 MET | `disputes` / `standalone` / `BA1` |

fixtureは `tests/fixtures/population-met-{prepared,evidence,rules}.json` に分離し、評価入力、
人口頻度Evidence、ルール設定の境界も実運用と同じにする。テストでは内部statusだけでなく、
生成された個別VA-Spec JSONを再読込し、direction、strength、outcome、AF、AC、AN、
BA1例外確認EvidenceのIRIまで確認する。

加えて、demo-dataのcase2-var1と固定済みgnomAD 4.1.1応答を使う回帰テストを行う。
観測された最大集団AFは約0.000013602であり、`case2-pm2-rules.json` のテスト限定閾値
0.000014ではPM2 METになる。これは実データ経路の動作確認用であり、疾患別に承認された
臨床閾値ではない。case2-var2のgnomAD未登録応答にはANとcallabilityがないため、AF=0の
StudyResultを作らずNOT_EVALUATEDとする。

公式schema IDは次である。

```text
https://w3id.org/ga4gh/schema/va-spec/1.0.1/acmg-2015/json/VariantPathogenicityEvidenceLine
```

## 現時点で意図的に省略する情報

- `targetProposition` と最上位 `Statement`: 疾患conditionが未入力。
- 最終classification: 本工程は基準別Evidence評価であり、最終5段階分類は対象外。
- 完全なVRS/CatVRS Allele object: 現在は正規化variantの決定的IRI参照を使う。
- curator、contribution、最終判定日時: キュレーションworkflowが未実装。

完全な公式Statement例へ拡張するには、少なくとも疾患ID・疾患名、最終classification、
分類方法と版、VRS/CatVRS variant ID、curator/contribution provenanceが必要になる。

## 参照

- [VA-Spec 1.0.1 Examples](https://va-spec.ga4gh.org/en/stable/examples/)
- [ACMG Variant Pathogenicity Statement with Evidence](https://va-spec.ga4gh.org/en/stable/examples/acmg-variant-pathogenicity-statement-with-evidence.html)
- [VariantPathogenicityEvidenceLine schema](https://w3id.org/ga4gh/schema/va-spec/1.0.1/acmg-2015/json/VariantPathogenicityEvidenceLine)
- [CohortAlleleFrequencyStudyResult schema](https://w3id.org/ga4gh/schema/va-spec/1.0.1/base/json/CohortAlleleFrequencyStudyResult)
