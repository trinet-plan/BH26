# 実装の俯瞰

更新: 2026-09-16。数値は固定キャッシュによるオフライン実行の結果で、下記コマンドで再現できる。

| 文書 | 役割 |
| --- | --- |
| この文書 | 全体像・検証範囲・16基準の状況・判定方針・残件 |
| [計画](PLAN.md) | 当初計画と不変条件、対象範囲 |
| [進捗](PROGRESS.md) | 工程順の実装記録 |
| [PM1設計](design/PM1-PLAN.md) | PM1の2ルート設計とhotspot policy |
| [PVS1設計](design/PVS1-PLAN.md) | PVS1 General decision treeとEvidence契約 |
| [ClinGen正例リファレンス](CLINGEN-POSITIVE-REFERENCES.md) | 14 criterionのClinGen正例と再現方法 |
| [MET一覧](MET-RESULTS.md) | MET 52件の変異・強度・根拠数値 |
| [予測Evidence仕様](PREDICTION-EVIDENCE.md) | 蛋白注釈・予測値の取得と校正policy |
| [VA-Spec出力](VA-SPEC-OUTPUT.md) | JSONのフィールド対応と制約 |

```powershell
$env:PYTHONPATH = 'src'
.venv/Scripts/python.exe -m acmg prepare-demo-online --input-dir demo-data --cache-dir tests/fixtures/ensembl-cache --evidence-cache-dir tests/fixtures/external-cache --output-dir work/full/prepared --ensembl-release 116 --with-gnomad --with-clinvar --clinvar-release 2026-09-15 --with-dbnsfp --with-pm1-hotspot --rules config/demo-rules.json --offline
.venv/Scripts/python.exe -m acmg evaluate --input work/full/prepared/variants.json --evidence work/full/prepared/evidence.json --config config/demo-rules.json --context config/curated-context.json --criteria all --offline --output-dir work/full/evaluated
```

## パイプライン

```
audit-demo  →  prepare-demo-online  →  evaluate
（原本監査）    （同定・補正・Evidence取得）   （16基準評価・VA-Spec出力）
                                              ← --config  閾値policy（版付き）
                                              ← --context 人のキュレーション
```

demo-dataの4症例28 ALTレコードを原本非改変で監査し、VERIFIED 24・CORRECTED 4・PENDING 0で
同定する。閾値と人の判断はコードに持たず、それぞれ`--config`と`--context`から入る。

## データ源

| Provider | 取得内容 | Evidence件数 |
| --- | --- | ---: |
| Ensembl VEP（release 116） | version付きRefSeqの蛋白注釈 | annotation 17 |
| MyVariant.info dbNSFP 4.8a | REVEL・AlphaMissense（版を固定） | computational 22 |
| Ensembl VEP | SpliceAI・保存性（版は宣言依存） | computational 42 |
| gnomAD 4.1.1 | 集団別頻度 | population 293 |
| ClinVar E-utilities | VCV同定・PS1 exact・PM5 residue・PM1密度 | 61 |

全応答は内容ハッシュ付きキャッシュに固定し、`--offline`でSHA-256一致再現する。
人のキュレーションは`config/curated-context.json`から入る（現在はClinGen SVI BA1例外リスト
9変異、手入力）。

## ClinGen正例スイート

PP5/BP6を除く実装対象14 criterionは、ClinGen Evidence Repositoryの公開VCEP解釈から選んだ
13変異で最低1件ずつMETを再現する。ERepoの`metCodes`は期待値manifestにだけ置き、evaluatorには
prepared variantと独立に正規化したEvidenceだけを渡すため、専門家パネルの結論をそのまま根拠に
する自己参照はない。VA-Spec 1.0.1 Evidence Lineの生成・schema検証までE2Eで確認する。

```powershell
$env:PYTHONPATH = 'src'
.venv/Scripts/python.exe -m unittest tests.test_clingen_positive -v
```

2026-09-16取得のERepo 13,265件を全走査した結果、14 criterionすべてに現行MET例が存在した。
BP3も23件あり、過去の限定556件走査で0件だった記録を更新した。個々のUUID、座標、適用コード、
General Guidance evaluatorとの強度差は[正例一覧](CLINGEN-POSITIVE-REFERENCES.md)に固定している。

## 検証範囲（demo-data + fixture）

demo-dataと`tests/fixtures/`の4スイートを合わせた49レコード・758判定が現在の検証範囲。
fixtureのEvidenceは各スイートの対象criterionをMETにする分だけを積むため、
NOT_EVALUATEDが多いのは設計どおりで、網羅性の不足を意味しない。

| スイート | レコード | 判定数 | 出力検証 |
| --- | ---: | ---: | --- |
| demo-data（4症例） | 28 | 448 | VA-Spec 1.0.1 schema |
| clingen-positive | 13 | 273 | VA-Spec 1.0.1 schema |
| clingen-gate | 5 | 80 | 内部評価 |
| synthetic | 1 | 16 | 内部評価 |
| population-met | 2 | 32 | 内部評価 |
| 合計 | 49 | 758 | |

### criterionごとのMET到達

| code | demo MET | fixture MET | 到達 |
| --- | ---: | ---: | --- |
| PVS1 | 0 | 1 | fixtureのみ |
| PS1 | 0 | 1 | fixtureのみ |
| PM1 | 2 | 1 | 両方 |
| PM2 | 0 | 4 | fixtureのみ |
| PM4 | 0 | 1 | fixtureのみ |
| PM5 | 1 | 1 | 両方 |
| PP2 | 0 | 1 | fixtureのみ |
| PP3 | 6 | 2 | 両方 |
| PP5 | - | - | DEPRECATED |
| BA1 | 11 | 2 | 両方 |
| BS1 | 0 | 1 | fixtureのみ |
| BP1 | 0 | 1 | fixtureのみ |
| BP3 | 0 | 1 | fixtureのみ |
| BP4 | 12 | 2 | 両方 |
| BP6 | - | - | DEPRECATED |
| BP7 | 0 | 1 | fixtureのみ |

実装対象14 criterionすべてがMET終端に到達する。demo-dataだけで到達するのは
PM1・PM5・PP3・BA1・BP4の5件で、残り9件はfixture側でのみ到達する。
clingen-gateのMET 2件（PP3・BP4）はゲート確認の副産物として付随的に成立したもの。

### status合算

| status | demo | fixture | 合計 |
| --- | ---: | ---: | ---: |
| MET | 32 | 20 | 52 |
| NOT_MET | 108 | 14 | 122 |
| NOT_APPLICABLE | 147 | 113 | 260 |
| NOT_EVALUATED | 101 | 145 | 246 |
| MANUAL_REVIEW | 4 | 2 | 6 |
| DEPRECATED | 56 | 16 | 72 |
| 計 | 448 | 310 | 758 |

## MET一覧

MET 52件（demo-data 32・fixture 20）の変異・強度・根拠数値は
[MET一覧と成立理由](MET-RESULTS.md)に分離した。demo-dataのMETはPM1・PP3・BP4・BA1・PM5の
5 criterionに集中し、fixture 20件が残る9 criterionを覆う。

## 16基準の状況

### 判定に到達している（7基準）

| 基準 | MET | NOT_MET | 根拠 |
| --- | ---: | ---: | --- |
| BP4 | 12 | 13 | REVEL + SpliceAI |
| BA1 | 11 | 15 | gnomAD + BA1例外リスト |
| PP3 | 6 | 19 | REVEL + SpliceAI |
| PM1 | 2 | 0 | ClinVar missense density |
| PM5 | 1 | 17 | ClinVar residue検索 |
| PS1 | 0 | 18 | ClinVar exact検索 |
| PM2 | 0 | 26 | gnomAD |

### 構造的に対象外を確定できる（4基準、147件）

PM4（28）・BP3（28）・PVS1（25）・BP7（21）。転写産物注釈から確定し、再評価を要しない。
demo-dataにin-frame indelとstop lossが無いためPM4とBP3のゲートは実データで開かない。
ただしこれはdemo-dataの構成によるもので、判定ロジックが未検証という意味ではない。
ClinGen Evidence Repository由来の5変異（`tests/fixtures/clingen-gate-*.json`）で、
PVS1・PM4・BP7が対象consequenceで開き、それ以外では閉じることを確認している。
さらに`clingen-positive-*` fixtureでは、PM4（`clingen-positive:pm4-otc`、moderate）と
BP3（`clingen-positive:bp3-foxg1`、supporting）がゲート通過後の
`protein_length_change`・`nonfunctional_repeat`・`repetitive`・`functional_importance`
まで走ってMET終端に到達する。BP3はgate fixtureを持たず正例fixtureのみで確認している。
PVS1 evaluator自体はnonsense/frameshift、canonical/実証splice LoF、start-lossのGeneral
decision treeと4段階強度に対応済み。デモのPVS1候補3件はgene_disease等の独立Evidenceが
未投入のためNOT_EVALUATEDであり、判定ロジック未実装を意味しない。

### Evidence待ち（4基準、101件）

| 待っているもの | 基準 | 件数 |
| --- | --- | ---: |
| gene_disease（人手キュレーション） | PP2・BP1・PVS1 | 39 |
| 疾患と疾患別頻度閾値 | BS1 | 28 |
| ClinVar報告密度・region評価 | PM1 | 26 |
| synonymous_assessment・gnomAD未登録 | BP7・PM2・BA1 | 8 |

### 非推奨（2基準）

PP5・BP6を各28件DEPRECATED。外部assertionを点数化せず一次Evidenceを取りに行く方針による。

### 集計（demo-dataのみ。fixtureを含む合算は[検証範囲](#検証範囲demo-data--fixture)）

```
MET 32 / NOT_MET 108 / NOT_APPLICABLE 147 / NOT_EVALUATED 101 / MANUAL_REVIEW 4 / DEPRECATED 56
VA-Spec Evidence Line 140件（MET + NOT_MET）
```

## 判定方針

| 原則 | 具体 |
| --- | --- |
| ラベルを評価入力に渡さない | CLNSIG・ACMG_CODES・NOTEは監査専用 |
| 「報告がない」と「成立しない」を分ける | gnomAD未登録・密度不足・検索未完了はNOT_EVALUATED |
| 自己参照の禁止 | PM1密度とcomparatorから評価対象自身を除外 |
| 循環の禁止 | 自動Evidenceでcritical functional domainを主張できない |
| 版の確定を要求 | 校正できない予測値は不使用。例外は出典付き宣言のみ |
| 機序を加算しない | PP3はOR（最強区間）、BP4はAND（全機序が良性側） |
| 疾患文脈の分離 | conditionは任意（BS1のみ必須）、condition_assessmentに記録 |
| 由来の申告 | curatedリストはentry_method必須。手入力はそう記録する |
| 座標の厳密性 | isoform表記の取り違えをefetch照合で排除 |
| 検索scopeの区別 | exact検索をPM5の不成立根拠に流用しない |

## 残件（優先順）

1. gene_disease（39件）。ClinGen CSpec RegistryにMYH7・MYBPC3・TNNI3の仕様がある。
2. BS1の疾患と閾値（28件）。入力経路は完成済みで、閾値のキュレーションのみ。
3. BA1例外リストの第三者照合。手入力のため原典Table 1との再確認が未実施。
4. SpliceAIの版が宣言依存。版を確定できる配布元からの取得が望ましい。
5. PM1のgene/disease対応、強度可変、critical domainのreviewed Evidence入力経路。
6. BP7の位置policy。PVS1はNMD・MANE/transcript・RNA・protein region Evidenceの自動取得拡張。
7. CI未設定。

## 対象外

最終5段階分類、スコア加点、ランキング、Web UI、CSpecの自動適用、残り12基準、GRCh37、CNV/SV。
`run-manifest.json`は`final_classification: OUT_OF_SCOPE`を明示する。
