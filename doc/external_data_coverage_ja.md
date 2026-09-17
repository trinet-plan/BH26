# 外部データ取得の対応状況と未対応事項

この文書は、criterion の**判定実装の有無**と、判定に必要な**外部データを取得できるか**を分けて記録する。`UNKNOWN` は「陰性」ではなく、必要な根拠を安全に取得・確認できなかった状態である。

## 現在の外部データ経路

`acmg_pipeline.services.resolve.ProviderEvidenceResolver` が、VCFの座標を起点に次のプロバイダから normalised evidence を組み立てる。VCF INFO は監査用入力であり、版・assembly・transcript・取得元が明示された根拠へ自動変換しない。

| データ | 現在の経路 | 利用するcriterionの例 |
|---|---|---|
| consequence、HGVS、transcript、予測値 | Ensembl VEP | PVS1、PM4、BP3、BP7、PP3、BP4 |
| REVEL、SpliceAI等の追加予測値 | Ensembl VEP / dbNSFP由来の注釈 | PP3、BP4、PS1、PM5 |
| AF、AC、AN、FILTER、集団別頻度 | gnomAD GraphQL | PM2、BA1、BS1 |
| individual ClinVar record、同一残基比較、hotspot検索 | NCBI ClinVar | PS1、PM1、PM5、PP5、BP6 |
| gene--disease mechanism、頻度閾値、領域ポリシー | サーバー設定のreview済みcontext | PVS1、PM1、PP2、BP1、BS1 |

## 現在も未対応、またはreview済みデータ投入が必要な項目

### 1. 実装済みcriterionでも根拠データが不足しうる部分

判定コードがあっても、以下は変異単体の外部検索だけでは確定できない。review済み・版付きデータがない場合は `UNKNOWN` とする。

| 項目 | 影響するcriterion | 不足時の扱い |
|---|---|---|
| 遺伝子--疾患のLoF / missense機序 | PVS1、PP2、BP1 | `UNKNOWN`。遺伝子だけで一般化しない |
| 疾患別頻度閾値 | BS1（および頻度評価） | `UNKNOWN`。単一の汎用閾値で置換しない |
| hotspot・critical domain・良性変異スペクトラム | PM1 | `UNKNOWN`。近傍の病的登録数だけでは決めない |
| 対象transcript、NMD、開始コドン・exon評価 | PVS1 | `UNKNOWN` または限定的な強度。VEP consequenceだけでPVS1を確定しない |
| predictorの固定版 | PP3、BP4、PS1、PM5 | predictor versionは `UNKNOWN` と明記したまま、取得日時・Ensembl release・レスポンスhash・許可ポリシーを残して使用できる。固定版が必須の運用では `UNKNOWN` 判定にする |
| ClinVar assertionの詳細、review status、独立性 | PS1、PM5、PP5、BP6 | `UNKNOWN`。表示ラベルだけで支持根拠にしない |

### 2. 現在はstubとして残すcriterion

次の9基準は実判定を行わず、全28本のVA-Spec EvidenceLineには `UNKNOWN` として出力する。

| criterion | 実装に必要な主な追加情報 |
|---|---|
| PS2 / PM6 | 親子検体、血縁確認、de novo確認度 |
| PM3 | phase、相手allele、疾患・遺伝形式、症例数 |
| PP1 / BS4 | 家系内segregation、表現型、浸透率 |
| PP4 | 疾患特異的表現型・検査所見 |
| BS2 | 健常成人の観察、疾患の年齢依存性・浸透率 |
| BP2 | cis/trans情報、もう一方の変異の病原性 |
| BP5 | 代替となる分子診断と表現型の整合 |

## TogoVarへ置換する場合の境界

TogoVarはGRCh38の座標・REF・ALTから、VEP transcript/HGVS/consequence、gnomAD等の頻度、ClinVar/MGeNDラベル、SIFT/PolyPhen/AlphaMissenseを取得できる。従って次は置換候補である。

| 判定群 | TogoVarのみでの見込み | 残る不足情報 |
|---|---|---|
| PM2 / BA1 / BS1 | 基本頻度は取得可能 | 集団粒度、callability、release/version、疾患別閾値 |
| PM4 / BP3 / BP7 | consequenceとHGVSで概ね可能 | 採用transcriptの固定 |
| PP3 / BP4 | 部分的 | REVEL、SpliceAI、各predictorの固定版 |
| PVS1 / PS1 / PM1 / PM5 / PP2 / BP1 / PP5 / BP6 | 完全置換不可 | 上表の機序、比較、領域、詳細ClinVarレビュー情報 |

## 安全な取得優先順位

1. VCFに assembly・transcript・データ版・取得元を伴う値があれば、それをnormalised evidenceとして保存する。
2. 欠ける変異単体データはTogoVarを `GRCh38:chrom-pos-ref-alt` で照会する。
3. TogoVarにない、またはprovenanceが判定要件を満たさないデータは補完・推測しない。
4. そのcriterionを `UNKNOWN` とし、missing inputs と理由をVA-Specの `bh26AssessmentDetails` に残す。

TogoVarだけを唯一の外部経路とする場合、既存のEnsembl・gnomAD・ClinVar直接呼出しへ暗黙にfallbackしない。必要なら、どのcriterionで別のreview済みデータを許可するかを設定で明示する。

## キュレーターが関わる箇所

自動処理は候補データの取得と、承認済みポリシーへの機械的な照合までを担う。疾患特異的な知識の採否、根拠の独立性、患者・家系情報の真偽はキュレーターが確認する。

### 事前にreview・版管理するデータ

| キュレーターの作業 | 影響するcriterion | 成果物 |
|---|---|---|
| gene--disease mechanismを確認する（LoF、missense、truncatingの疾患機序） | PVS1、PP2、BP1 | 疾患・遺伝形式・根拠文献・review日・版を伴うcontext record |
| 疾患別の最大許容頻度を確認する | BS1、BA1、PM2 | 疾患別・遺伝形式別の閾値、対象集団、計算法、review日 |
| hotspot / critical domainの採用範囲を確認する | PM1 | transcript/protein範囲、良性変異評価、適用ポリシー、根拠文献 |
| predictorの利用可否を承認する | PP3、BP4、PS1、PM5 | predictor、版または`UNKNOWN`許可、適用範囲、閾値、理由 |
| ClinVar assertionの採用基準を確認する | PS1、PM5、PP5、BP6 | accepted review status、対象疾患一致条件、競合時の扱い |

### 変異ごとに確認する事項

| 確認対象 | 主なcriterion | キュレーターが判断する内容 |
|---|---|---|
| transcriptと変異の同一性 | 全criterion、特にPVS1 | VCF座標、REF/ALT、assembly、HGVS、採用transcriptが同じ変異を指すこと |
| PVS1のdecision tree | PVS1 | NMD、重要exon、開始コドン、spliceの影響、LoF機序への適合 |
| ClinVar比較候補 | PS1、PM5、PP5、BP6 | 同一残基/別変異の妥当性、対象疾患、review status、独立性、競合分類 |
| 領域・集積の妥当性 | PM1 | hotspot/domainが対象疾患に適用可能か、良性変異の反証がないか |
| 頻度の解釈 | PM2、BA1、BS1 | ancestry、品質FILTER、callability、疾患別閾値との適合、集団間の矛盾 |
| version不明の予測値 | PP3、BP4、PS1、PM5 | `predictor_version=UNKNOWN`を補助根拠として採用してよいか。採用しても固定版と誤認しないこと |
| 文献判定 | PS3、BS3、PS4 | LLMが抽出した実験・症例対照の要約と原著を照合し、ACMG strengthを承認すること |
| 患者・家系固有の根拠 | PS2、PM3、PM6、PP1、PP4、BS2、BS4、BP2、BP5 | phase、de novo、segregation、表現型、健常観察、代替診断を原資料から確認すること |

### 出力上の扱い

キュレーター確認が未了、または必要な根拠が不足する場合、自動処理は病原性・良性の方向へ推測しない。criterionは `UNKNOWN` とし、理由、missing inputs、review points を各EvidenceLineの `bh26AssessmentDetails` に残す。キュレーターが確認済みのcontext/policyを更新した後に、同じ入力を再評価する。
