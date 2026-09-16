# 実装の俯瞰

更新: 2026-09-16。数値は固定キャッシュによるオフライン実行の結果で、下記コマンドで再現できる。

| 文書 | 役割 |
| --- | --- |
| この文書 | 全体像・検証範囲・MET一覧・16基準の状況・判定方針・残件 |
| [計画](PLAN.md) | 当初計画と不変条件、対象範囲 |
| [進捗](PROGRESS.md) | 工程順の実装記録 |
| [PM1設計](PM1-PLAN.md) | PM1の2ルート設計とhotspot policy |
| [PVS1設計](PVS1-PLAN.md) | PVS1 General decision treeとEvidence契約 |
| [ClinGen正例リファレンス](CLINGEN-POSITIVE-REFERENCES.md) | 14 criterionのClinGen正例と再現方法 |
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

## MET一覧と成立理由

判定に使った数値・校正区間・除外理由は`results.json`の`provenance`と`evidence`に全件残る。
以下は上記コマンドで再現できるMET 52件の内訳。

### demo-data（32件）

PM1（2件）: ClinVar missense density、route=mutational_hotspot、窓±5aa。閾値は
`config/demo-rules.json`の`PM1.hotspot`から入り、評価対象の変異自身は密度に数えない。

| record | gene | 変異 | 強度 | P/LP | benign |
| --- | --- | --- | --- | ---: | ---: |
| case1:26:1 | KCNJ5 | p.Arg155Gln | moderate | 6 | 0 |
| case3:24:1 | MYH7 | p.Arg719Trp | moderate | 9 | 0 |

PP3（6件）: REVEL（dbNSFP 4.8a、Pejaver et al. 2022校正）。6件とも蛋白機序で成立し、
splicing機序は区間に届かない。PP3はORなので最強の区間を採る。

| record | gene | 変異 | 強度 | REVEL | SpliceAI |
| --- | --- | --- | --- | ---: | ---: |
| case1:26:1 | KCNJ5 | p.Arg155Gln | strong | 0.934 | 0 |
| case3:24:1 | MYH7 | p.Arg719Trp | moderate | 0.814 | 0.02 |
| case1:30:1 / case4:30:1 | TNNI3 | p.Pro82Ser | supporting | 0.76 | 0 |
| case3:25:1 | MYBPC3 | p.Glu334Lys | supporting | 0.757 | 0.02 |
| case4:26:1 | DSG2 | p.Phe531Cys | supporting | 0.744 | 0 |

BP4（12件）: REVEL + SpliceAIのAND。synonymous 4件はdbNSFPがnonsynonymous SNV用で
REVELを持たないため、splicing機序のみで成立する。

| gene / 変異 | 件数 | 強度 | REVEL | SpliceAI |
| --- | ---: | --- | ---: | ---: |
| MYBPC3 p.Ser236Gly | 2 | moderate | 0.173 | 0.01 |
| MYBPC3 p.Arg382Trp | 4 | supporting | 0.265 | 0.07 |
| MYBPC3 p.Val158Met | 2 | supporting | 0.228 | 0.01 |
| LMNA p.Ala287= | 2 | supporting | - | 0 |
| TNNI3 p.Arg68= | 1 | supporting | - | 0.01 |
| TNNI3 p.Glu66= | 1 | supporting | - | 0.02 |

BA1（11件）: gnomAD 4.1.1。全件がClinGen SVI例外リスト9変異に非該当。AN不足の集団は
`rejected_observations`に`INSUFFICIENT_AN`として理由が残る。下2行はglobalが5%未満でも
特定集団で超える例。

| record | 変異 | exome:global AF | 最大集団 |
| --- | --- | ---: | --- |
| case1:29:1 / case3:27:1 | 11:47348490T>C | 0.1202 | nfe 0.1331 |
| case2:28:1 / case3:28:1 | 11:47350047C>T | 0.0895 | fin 0.1047 |
| case1:31:1 / case4:27:1 | 1:156135237T>C | 0.0814 | afr 0.4586 |
| case2:30:1 | 19:55156279C>A | 0.0501 | fin 0.0809 |
| case2:26:1 / case3:30:1 / case4:28:1 | 1:201361238G>C | 0.0088 | eas 0.1884 |
| case4:31:1 | 19:55156285C>T | 0.0031 | afr 0.1096 |

PM5（1件）: case3:24:1 MYH7 p.Arg719Trp、moderate。ClinVar residue検索で同一残基の
別アミノ酸変化にP/LP報告。`condition_assessment: NOT_EVALUATED`のため確認事項を併記する。

### fixture（20件）

clingen-positive 15件が14 criterionすべてを覆う。demo-dataで未評価のPVS1・PP2・BP1・BS1が
ここで成立するのは、fixtureがgene_disease・region・disease_frequency_thresholdを
最初から積んでいるため。demo側に足りないのはロジックではなくキュレーションであることを示す。

| criterion | 強度 | gene | 変異 | ERepo metCode | 成立理由 |
| --- | --- | --- | --- | --- | --- |
| PVS1 | very_strong | PAH | 12:102852850GA>G | PVS1 | relevant transcript/exonでNMD予測 |
| PS1 | strong | GCK | 7:44149809C>A | PS1 | exact検索で同一アミノ酸変化のP報告 |
| PM1 | moderate | PAX6 | 11:31802793C>T | PM1 | critical_functional_domain（curated） |
| PM2 | supporting | GUCY2D | 17:8009531T>C | PM2_Supporting | 全観測がrarity policyを満たす |
| PM4 | moderate | OTC | X:38369845GAAG>G | PM4 | in-frame deletion + reviewed region |
| PM5 | moderate | OTC | X:38367331C>T | PM5 | residue検索で同一残基の別変化にP報告 |
| PP2 | supporting | GCK | 7:44150049A>G | PP2 | gene_level機序 + variant spectrum |
| PP3 | strong | PAX6 | 11:31802793C>T | PP3_Moderate | REVEL 0.967 |
| BA1 | stand_alone | ITGB3 | 17:47283530T>C | BA1 | AF>5%、例外リスト非該当 |
| BS1 | strong | MYH7 | 14:23420189C>T | BS1 | max_credible_AF 0.0002（MONDO:0005045、AD） |
| BP1 | supporting | BRCA2 | 13:32333210G>C | BP1_Strong | gene_level機序 + variant spectrum |
| BP3 | supporting | FOXG1 | 14:28767509C>CGCCGCC | BP3 | 反復領域のin-frame挿入、機能的重要性なし |
| BP4 | moderate | SLC6A8 | X:153688650G>A | BP4 | REVEL 0.079 |
| BP7 | supporting | PAX6 | 11:31790720G>A | BP7 | 位置・splice予測・保存性が全て良性側 |
| PM2 | supporting | MYH7 | 14:23420189C>T | - | bs1-myh7で付随的に成立 |

PP3（strong / ERepo `PP3_Moderate`）とBP1（supporting / ERepo `BP1_Strong`）はVCEPが
criterion-specificに調整した強度との差で、CSpec自動適用を対象外としているため意図的に残す。
個々のUUIDと差の根拠は[正例一覧](CLINGEN-POSITIVE-REFERENCES.md)に固定している。

clingen-gate 2件（PP3 supporting: MYBPC3 c.405A>G、SpliceAI 0.27 / BP4 supporting:
MYH7 c.3036C>T、SpliceAI 0）はBP7ゲート確認用synonymousの副産物。REVELを持たないため
splicing機序のみで成立する。synthetic・population-met 3件（PM2 2・BA1 1）は
`policy_source: synthetic-test-only`の配線確認用で、臨床閾値ではない。

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
