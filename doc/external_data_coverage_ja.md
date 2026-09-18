# 外部データ取得の対応状況と未対応事項

この文書は、criterion の**判定実装の有無**と、判定に必要な**外部データを取得できるか**を分けて記録する。`UNKNOWN` は「陰性」ではなく、必要な根拠を安全に取得・確認できなかった状態である。

> **この文書はコードより遅れることがある。** 2026-09-18 の照合で、本文書の記述と実装が
> 4箇所食い違っていた（PM2の総不在の扱い、BS1のDRAFT閾値でのMET可否、PM2の不一致件数、
> PVS1のdemo-data到達状況）。**課題や精度を判断する前に、必ず再実行して裏を取ること。**
>
> ```bash
> ./.venv/Scripts/python.exe test_automated_criteria_ground_truth.py   # キュレーター一致率
> ./.venv/Scripts/python.exe -m pytest tests/ -q                       # 回帰
> ```
>
> システムのPythonには依存が入っていない。`.venv/Scripts/python.exe` を使うこと。
> `work/` 配下の古いラン成果物も同様で、`NOT_EVALUATED` / `NOT_APPLICABLE` /
> `MANUAL_REVIEW` を含むsummary.tsvは旧語彙（現行は `met` / `not_met` / `unknown` の3値）。

## 設計思想

**機械的に判断できないところは、その理由と判断材料を書いて人に渡す。**

推測で埋めない。既定値で黙って通さない。「データが無い」を「陰性」に読み替えない。
出力には、どこで止まったか・何が足りないか・何を仮定したかが残り、キュレーターがそれを
見て判断できる状態にする。

この原則が具体的にどう現れているか:

- `UNKNOWN` は陰性ではなく「決められなかった」。`missing_inputs` に何が足りないかを、
  `review_points` に何を確認すべきかを列挙する
- 自動導出した値は `assessment_method: "automated"` と `method` / `policy_version` を
  名乗り、キュレート済み判断と取り違えられないようにする（[1-4](#1-4-pvs1決定木を公開データでどこまで通したか)のPVS1ゲート、PP3/BP4の較正、PM1のhotspot）
- 仮定を置いた場合は仮定したことを出力に残す（BS1の統計量・比較演算子、
  `unknown_version_policy` による予測器バージョン不明の受け入れ）
- 判断が臨床的なものは、暫定値を入れる場合も `PLACEHOLDER` と明記する
  （`BS1.default_max_credible_af`、`PM2.max_af`）

以下の「未対応事項」は、この原則に照らして**まだ人に渡せていない**か、**人に渡す形は
できているが判断そのものが未了**かを区別して記録する。

## 決定記録: 未承認の疾患別頻度閾値（PM2 / BS1）

**決定日: 2026-09-17。** 疾患別の閾値が curator / VCEP により承認されていない場合、
公開情報から導出した値を最終的な ACMG/AMP 判定の根拠にしない。必要ならキュレーター向けの
`DRAFT` 候補として出力する。

**実装との差異（2026-09-18 に照合）。** この決定記録は「`MET` にせず `UNKNOWN` とする」と
書いているが、`acmg_pipeline/criteria/bs1.py` の実装はそうなっていない。未承認（`DRAFT`）の
閾値を超えた観測に対して `MET / strong` を返し、`draft_threshold_candidate` と review_points
で「承認待ちの予測であって最終判断ではない」と明示する形になっている。デモ28変異では
**23件が `threshold_scope: "default"` / `policy_status: "DRAFT"` のまま BS1 = MET** になる。
どちらが正しい方針かは未決。実装を決定記録に合わせるか、決定記録を実装に合わせるかを
決める必要がある（[1-2](#1-2-bs1の既定閾値が継承様式を見ていない未解決)も参照）。

- **PM2** は ClinGen SVI の勧告に従い、成立し得ても強度は `PM2_Supporting` とする。
  公開DBの未登録は、座位coverage・callabilityが確認できる場合だけ absence の根拠になり得る。

  **2026-09-18 更新。** 「TogoVarからレコードが返らないことは `UNKNOWN` とする」と書いて
  あったが、現在の実装は違う。`criteria/common.py` の `population_context()` は
  `treat_total_absence_as_evidence` を持ち、**PM2だけがこれを渡す**。照会した全ソースが
  正常に解決したうえで何も返さなかった場合（`NO_OBSERVATION`、取得エラーではない）に
  限り、その不在自体を弱い希少性の根拠として `MET / supporting` とし、「座位coverageは
  未確認」という review point を付ける。ソースに到達できなかった場合は従来どおり
  `UNKNOWN`。BA1/BS1（「頻度が高すぎる」を見る）は同じ沈黙から何も言えないため、この
  扱いはPM2に限定されている。`tests/test_pm2_absence_policy.py` が6本で固定している。
- **BS1** は、疾患有病率、遺伝形式、浸透率、遺伝的／アレル不均一性から導く
  maximum credible AF を用いる。公開された計算法やVCEPの閾値は、候補作成・レビューの
  出典には使えるが、他疾患へそのまま移植しない。
- `DRAFT` 候補には、計算法、入力値、有病率、遺伝形式、浸透率、最大遺伝子寄与・
  最大アレル寄与、対象集団、DB/release、AC/AN/AF（BS1では使用統計量と比較演算子）を残す。
  `policyStatus: "DRAFT"`、レビュー理由、出典URL・版を必須とする。
- 承認済みの疾患別仕様は、対象疾患・遺伝子・遺伝形式・DB版・例外変異を照合でき、
  `reviewed_at` を記録したものだけを runtime の `disease_frequency_threshold` として用いる。

公開根拠:

- ClinGen SVI, *Recommendation for Absence/Rarity (PM2)* v1.0 (2020-09-04):
  <https://clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf>
- Whiffin et al., *Using high-resolution variant frequencies to empower clinical genome
  interpretation* (2017): <https://pubmed.ncbi.nlm.nih.gov/28518168/>
- ClinGen Criteria Specification Registry: <https://cspec.clinicalgenome.org/>

## 現在の外部データ経路

`acmg_pipeline.services.resolve.ProviderEvidenceResolver` が、VCFの座標を起点に次のプロバイダから normalised evidence を組み立てる。VCF INFO は監査用入力であり、版・assembly・transcript・取得元が明示された根拠へ自動変換しない。

| データ | 現在の経路 | 利用するcriterionの例 |
|---|---|---|
| consequence、HGVS、transcript、予測値 | Ensembl VEP | PVS1、PM4、BP3、BP7、PP3、BP4 |
| REVEL、SpliceAI等の追加予測値 | Ensembl VEP / dbNSFP由来の注釈 | PP3、BP4、PS1、PM5 |
| AF、AC、AN、FILTER、データセット別頻度 | TogoVar API 0.9.1（GRCh38） | PM2、BA1、BS1 |
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

### 1-2. BS1の既定閾値が継承様式を見ていない（未解決）

疾患別閾値がない場合、BS1は `config` の既定閾値へフォールバックする。このとき
`inheritance` は参照されない。疾患別閾値の照合には使うが、既定閾値経路では無視される。

`config/bs1_thresholds_draft.json` の転記が示すとおり、継承様式は本来この判定に効く軸で
ある。ClinGen Hearing Loss VCEPはAR `0.003` / AD `0.0002` と**15倍**の差を付けている。
単一の既定値をどちらにも当てると、一方は緩すぎ、他方は厳しすぎる。

実例: demo-dataの `case4-var1`（DSG2 c.1592T>G、ARVC、**ホモ接合**症例）は、東アジア集団
でのみ AC=21/AN=39,696（FAF 0.00037）と観測され、既定閾値 0.0001 を超えて BS1 = MET に
なる。他の17集団はすべて AC=0 で、ToMMo・HGVDでも同様に日本人集団で 1e-4 台。頻度の
観測自体は堅い。一方でキュレーターはBS1を挙げておらず（`PS3_Strong` / `PS4_Strong` のみ）、
劣性/複合的な機序を想定していれば、ヘテロでの集団頻度はBS1の根拠として弱くなる。

現状の出力は `threshold_scope: "default"` と review_points で「疾患別頻度で確認せよ」と
明示するため、キュレーターは弾ける。ただし継承様式を見ない既定判定である事実は
summaryに出ない。

**2026-09-18 追記: 影響範囲は1件ではなく23件。** デモ28変異を再評価すると、BS1は
**23件 MET / 3件 not_met / 2件 unknown**。MET 23件はすべて `threshold_scope: "default"`、
`policy_status: "DRAFT"`、`draft_threshold_candidate` 付きである。既定閾値
`1e-4` ＋ `faf95` ＋ 継承様式非参照の組み合わせで、デモデータのほぼ全件が閾値を超える。
`case4-var1` は「キュレーターが挙げていないのにengineがMETにする」例として下に挙げていたが、
実際には例外ではなく多数派の1つだった。

想定される対応（いずれも臨床判断が要る、未決）:

- 既定閾値経路でも `inheritance` を必須にし、不明なら `UNKNOWN`
- AR / AD 別の既定値を持つ
- 既定閾値でのBS1をMETにせず、review必須の別状態として出す

### 1-3. PM2がキュレーター判断と一致しない（2026-09-18 解決済み）

> **この節は解決済みの記録である。** 現在のPM2はデモ4症例のキュレーター判断 **5/5 で一致**
> する。以下は経緯として残す。全体の一致率は 8/15 → 11/15 → 13/15 → **14/15** と推移した。

`test_automated_criteria_ground_truth.py` の照合で、キュレーターがMETとした15件のうち
7件をengineが確認できなかった。**うち5件がPM2**で、原因は別々の2つだった。

#### (a) 閾値 `max_af = 0` が厳しすぎた（3件）→ 暫定値を設定して解消

```
case1-var2  curator: MET -> engine: not_met  [highest AF 0.0000232 exceeds cutoff (0)]
case2-var1  curator: MET -> engine: not_met  [highest AF 0.0000136 exceeds cutoff (0)]
case3-var1  curator: MET -> engine: not_met  [highest AF 0.0000299 exceeds cutoff (0)]
```

`config/demo-rules.json` の `PM2.max_af = 0` は、ACMG/AMP 2015の文言
"Absent from controls" をそのまま実装したもの。しかしgnomAD規模のデータでは、実在する
希少変異はほぼ必ず数アレル観測される。キュレーターは **AF 1e-5台の変異にPM2を付与**して
おり、文字どおりのゼロは現行データに対して機能しない。

なお `max_af` に既定値はない（未設定は `UNKNOWN`）。`0` は明示的な設定値であって、
設定漏れではなかった。

**暫定対応（2026-09-17）**: `PM2.max_af = 0.00005` を設定。デモ28変異の最大観測AFは
キュレーターがPM2を付けた3件が `1.36e-5`〜`2.99e-5`、次に高い変異が `7.73e-4` で、
その間に26倍のギャップがある。この区間ならどこに置いても同じ結果になる。
`5e-5` を選んだのは `BS1.default_max_credible_af`（`1e-4`）より低く保つため。
同値にすると「最大AF ≤ 閾値ならPM2」「いずれかのAF > 閾値ならBS1」が完全な補集合になり、
すべての変異がPM2かBS1のどちらかに必ず該当してしまう。ACMGが想定するはずの中間帯
（どちらも成立しない領域）を残している。ただしデモデータでは `5e-5`〜`1e-4` に該当する
変異がなく、中間帯は 0/28 で、この設定の効果は現データでは観測されない。

値は臨床レビュー未実施。`threshold_source` にPLACEHOLDERと明記してある。

#### (b) 頻度データセットに登録がない変異が `UNKNOWN` になっていた（2件）→ 実装して解消

```
case1-var1  curator: MET -> engine: unknown  [No reliable population observation]
case2-var2  curator: MET -> engine: unknown  [No reliable population observation]
            provider_failures: [{'provider': 'TogoVar', 'reason': 'NO_OBSERVATION'}]
```

`acmg_pipeline/services/resolve.py` の原則「未登録をAF=0として扱わない」による。
BA1・BS1（「頻度が高すぎる」を見る）では、未登録から何も言えないのは正しい。

**しかしPM2では、未登録そのものが求めている根拠**である。同じ原則を3基準に一律適用した
結果、PM2だけ意味が反転している。

区別に必要なのは **その座位のcoverage**。「未登録かつ十分にcallable」なら absent from
controls だが、「未登録かつcoverage不明」は判断不能。現在のTogoVar APIレスポンスは座位
coverageを返さない。`usable_observations()` は `AC == 0` の観測に `callable: true` を要求して
いるが、レコード自体が返らない場合はcallability情報が存在しない。

**対応（実装済み）**: `population_context()` に `treat_total_absence_as_evidence` を足し、
**PM2だけがこれを渡す**。照会した全ソースが正常に解決し、いずれも `NO_OBSERVATION` を
返した場合に限り、不在そのものを弱い希少性の根拠として `MET / supporting` とする。
座位coverageが取れない事実は消さず、「総不在をAF観測ではなく弱い希少性の根拠として
扱った」ことを review point に残す。取得エラーがあれば従来どおり `UNKNOWN`（探索が
不完全なだけで、不在が確認できたわけではない）。BA1/BS1には渡さない。

上の「同じ原則を3基準に一律適用した結果、PM2だけ意味が反転している」が問題の記述で、
**基準ごとに沈黙の意味が違う**というのが解になった。座位coverageを取得する経路
（下表の3）は依然として未実装で、これがあればreview pointを外して確定的な absence に
できる。

#### 対応の方向

| # | 内容 | 状態 |
|---|---|---|
| 1 | ~~`PM2.max_af` に非ゼロの閾値を設定する~~ | 暫定 `5e-5` を設定済み。正式値は臨床レビュー待ち |
| 2 | BS1と同じく `frequency_statistic` / `comparison` を対で持たせる | 未実装。PM2で保守的なのは信頼区間の**上限**（「稀少と言い切れるか」は最悪ケースで判断する）で、BS1のFAF（下限）とは逆側 |
| 3 | coverageを取得できるreview済み経路を追加し、「未登録 + 十分なcoverage」を absent as evidence として扱う | 未実装。現在は総不在を弱い根拠として扱い、coverage未確認をreview pointで渡している。TogoVar APIだけでは不足 |

1が3件、`treat_total_absence_as_evidence` が2件の不一致を解消し、PM2は **5/5 一致**に
なった。全体の一致率は 8/15 → 11/15 → 13/15 と推移した。

### 1-4. PVS1決定木を公開データでどこまで通したか

PVS1の決定木（545行、4経路17ノード）は、実データに対して長く動いていなかった。ERepoで
`PVS1=met` とされた変異でも全件が `G01`（LoF機序）で停止し、変異型を読む `V01` にすら
到達していなかった。必要な値がすべてキュレーター判断で、まだ存在しなかったため。

#### 通した経路と、埋めたゲート

| ゲート | 埋めた方法 | 種別 |
|---|---|---|
| `G01`/`G02` LoF機序 | ClinGen Dosage の haploinsufficiency スコア（3 → 確立） | 公開データの参照 |
| `NF01`/`IC00` transcript関連性 | MANE Select 一致 | 公開データの参照 |
| `NF02` NMD予測 | ClinGen PVS1 2018 の規則 + VEP エクソン番号 | 規則の適用 |
| `NF03` exon関連性 | `NF01` の帰結 | 論理的帰結 |
| `NF07` タンパク損失率 | VEP のタンパク位置 + Ensembl の翻訳長 | 測定 |
| `IC02` 下流in-frame開始コドン | Ensembl CDS配列を3塩基ずつ読む | 測定 |
| `IC03` 上流の病的変異 | ClinVar 遺伝子検索（ミスセンス限定を外す） | 公開データの参照 |

いずれも `assessment_method: "automated"` を名乗り、`method` と `policy_version` を持つ。
キュレート済みレコードがあればそちらが優先される。

#### 到達状況（2026-09-17時点）

| 経路 | 到達点 | 残る判断 | 実変異 |
|---|---|---|---|
| truncating | **`MET / very_strong`** | なし | `RUNX1 c.601C>T`、`MYBPC3 c.278delA`、`c.836del`（いずれもパネル判断と強度まで一致） |
| 領域（NMD逃れ） | `NF06` で停止 | 領域の生物学的関連性 | `VHL c.610G>T`、`c.634G>T`（損失率4.7%は算出済み） |
| 開始コドン喪失 | `IC01` で停止 | 無傷の代替transcriptの有無 | `GJB2 c.2T>C`（codon 34）、`HNF4A c.3G>A`（codon 71） |
| スプライス | `SP01` で停止 | 自動化しない判断（後述） | `MYBPC3 c.2905+1G>A` 他25件 |

**残る3経路の停止点は、いずれも「生物学的に関連か」型の判断1つだけ**になっている。
キュレーターは材料が揃った状態で1問に答えればよい。

**この表はERepo由来の変異についてのもので、demo-dataの28変異ではない。** demo-dataを
評価するとPVS1は28件すべて `unknown` になるが、内訳は **25件が「truncating /
canonical-splice / start-loss / RNA確認のいずれでもない」= そもそもPVS1の対象外**、
残り3件が `transcript_assessment` 2件と `splice_assessment` 1件での停止である。機序
ゲート（G01）で止まっているわけではない。demo-dataのPVS1 unknown件数を「機序データ
不足」と読むのは誤り。

#### 機械的に決めない箇所（意図的に未解決のまま返す）

| 状況 | 返り値 | 理由 |
|---|---|---|
| haploinsufficiency = 30（常染色体劣性） | レコードを出さない | 「ハプロ不全でない」であって「LoF機序でない」ではない。PAHは30だがPKU VCEPは `PAH c.806delT` にPVS1を適用している |
| 評価対象transcriptがMANE Selectでない | レコードを出さない | キュレーションは正当に別transcriptを使う（TNNT2の実例、`test_data/resolve_erepo_transcripts.py`） |
| PTCが**最終前**エクソン内 | レコードを出さない | 「末端50nt以内」の距離はエクソン番号に含まれない。**最終**エクソンは最終ジャンクションを越えているので距離は無関係で、`predicted=False` と言い切れる |
| エクソン番号が取れない | `exon_relevance` を付けない | 番号が無いことは「エクソンが存在しない」を意味しない |
| CDSが3の倍数でない/短すぎる | レコードを出さない | 問い自体が成立しない。下流ATGが**無い**ことは答えとして記録する |
| `IC03` で未配置の病的変異が残る | レコードを出さない | ClinVarは単一残基置換にしか位置を持たない。未配置のナンセンス変異が上流の証拠かもしれないので「無い」とは言えない |
| `IC03` の検索が打ち切られた | レコードを出さない | どちらの答えも支持できない |

`IC03` は **Yes と No で基準が違う**。配置確定した病的置換が1件でも上流にあれば Yes、
No には加えて「検索が返した全件が配置できた」ことを要求する。

#### 設計上の制約を2回踏んだ — 1カテゴリ1レコード

`_select_context_record()` は transcript ごとに1レコードしか受け付けず、2件あると
conflict を返す。これを2度踏んだ。

- `NF03`（`transcript_assessment`）— エクソン取得を共有ヘルパーに切り出し、MANE providerの
  レコードに渡す形にした
- `IC03`（`initiation_assessment`）— 別レコードで出したところ実際に
  `MANUAL_REVIEW / Resolve initiation_assessment assessments` になった。フィールドを
  `upstream_` 接頭辞付きでIC02のレコードにマージする形に変更

**1つのゲート群に対して provider は1つ**。複数の情報源が同じカテゴリを埋める場合は、
レコードを増やさずフィールドをマージする。

#### 引いてくるか、予測するか — 今回引いた線

自動化したゲートと、しなかったゲートの違いはここにある。

| | 内容 | 例 |
|---|---|---|
| **引いてくる** | 既に決着している事実を公開データから参照する | ClinGen Dosage のスコア、MANE Select、ClinVar の報告 |
| **規則を適用する** | 公開された規則を、入力が揃っている範囲で適用する | ClinGen PVS1 2018 のNMD規則（境界は適用しない） |
| **測定する** | 座標や配列から一意に決まる値を計算する | タンパク損失率、下流ATGの位置 |
| **予測する** | 何が起きるかを推定する | スプライス供与部位喪失の帰結 → **自動化しない** |
| **判断する** | 生物学的・臨床的な意味づけ | 領域の関連性、代替transcriptの有無 → **人に渡す** |

上3つは自動化し、下2つは人に渡した。`G01`〜`NF03`・`NF07`・`IC02`・`IC03` は上3つに
収まるが、`SP01` は「予測する」に入るため線の向こう側にある。

`IC01`・`NF04`・`NF06` はいずれも「生物学的に関連か」を問う判断で、材料を揃えたうえで
未解決のまま返している。

#### `SP01`（スプライス経路）は自動化しない — 実測による判断

`SP01` は `alternative_rescue` / `splice_outcome` / `reading_frame_disrupted` を要求する。
このうち `splice_outcome` は「エクソンスキップが読み枠を壊すか」なので、エクソン長から
計算できるように見える。AutoPVS1もこの計算を行う。

**実測すると、この素朴な規則は手元の唯一の実スプライス症例で専門家判断と食い違う。**

```
MYBPC3 c.2905+1G>A   ClinGen判定: PVS1 very_strong
  VEP:      intron 27/34  (ENST00000545968 = MANE NM_000256.3)
  Ensembl:  exon 27 の長さ = 168 bp  ->  168 % 3 == 0  ->  IN_FRAME
```

素朴な規則なら `IN_FRAME` → `reading_frame_disrupted = False` → `_region_path` へ分岐し、
`very_strong` にならない。スプライス供与部位の喪失は単純なエクソンスキップとは限らず、
イントロン保持（多くはframeshift）や隠れスプライス部位の使用など複数の帰結があり、配列
だけでは決まらない。ClinGenのスプライス仕様(2023)がRNA解析を重視するのはこのため。

`alternative_rescue`（代替スプライシングがLoFを回避するか）はさらに導出元がない。組織
特異的アイソフォームの知識を要するキュレーター判断で、`False` を既定にすればPVS1は進むが、
それは過剰判定の方向であり、上表の保守的な選択と逆になる。

**`G01`〜`NF03` は「既に決まっていることを引いてくる」だったが、`SP01` は「予測する」。**
自動化するなら次のいずれかで、いずれも方針決定が要る（未決）:

- SpliceAI の delta score を使う（取得済みだが `calibration_eligible: False`、モデル版不明）
- RNA解析データを入力として受け取る（PVS1は `rna` 経路を実装済み: `lof_effect_confirmed`）
- `alternative_rescue` の扱いを設定でポリシーとして明示する

現状 `MYBPC3 c.2905+1G>A` は `SP01` で停止し、`missing: ['splice_assessment']` を返す。
キュレーターには「何が足りないか」が出るので、判断の材料は渡せている。

### 1-5. BP4が「棄権した予測器」を反証として数えていた（2026-09-18 解決済み）

較正済み予測器のうち1つでも良性を支持しない場合、BP4はNOT_METになっていた。根拠は
ClinGen SVI / Walker et al. 2023 の「低いタンパクスコアは壊れたスプライス部位について
何も言っていない」で、これ自体は正しい。問題は**「影響なしを示した」と「何も示さなかった」
を同一視していた**こと。

使用中の2較正はどちらも2つのバンドの間に明示的な空白帯を持つ。

| 較正 | BP4バンド | 空白帯 | PP3バンド |
|---|---|---|---|
| REVEL（Pejaver et al. 2022） | ≤ 0.290 | **0.290〜0.644** | ≥ 0.644 |
| SpliceAI（Walker et al. 2023） | ≤ 0.1 | **0.1〜0.2** | ≥ 0.2 |

空白帯は較正が「どちらとも言えない」と明示している区間で、反証ではない。これを反証として
数えた結果、`MYH7 c.4472C>G`（REVEL 0.398、SpliceAI 0）はキュレーターのBP4に対して
NOT_METを返していた。タンパク予測器が棄権しただけなのに、反対したものとして数えられていた。

**対応（実装済み）**: 各スコアを自分側のバンドと**反対側の criterion のバンドの両方**に
当て、反対側のバンドに入るものだけがブロックする。結果には「何を支持しているか」を
名指しする（`Another predicted mechanism supports PP3 rather than a benign effect
(SpliceAI 0.9)`）。棄権した mechanism は良性を支持もしないので、握りつぶさず
`provenance.uninformative_mechanisms` と review point に残す。良性判定が実際より広く
読まれないようにするため。PP3は変更なし（1つのmechanismで成立する仕様）。

デモ28変異でBP4のMETは12→16、キュレーター一致率は13/15→**14/15**。

残る不一致 `TNNT2 c.801+41C>G` は適用できる予測器がSpliceAI 0.13の1つだけで、空白帯の
ど真ん中にある。較正が判定を拒否している値からBP4を主張はできないため、これはNOT_METの
まま**正しい**。較正を緩めない限り一致しないし、緩めるべきではない。

**未決**: この `not_met` は `unknown` の方が正確かもしれない。「唯一の予測器が棄権した」は
陰性ではなく gap であり、本文書の「データが無い ≠ 陰性」という原則に当たる。

### 1-6. `gene_disease` は1カテゴリに2つの生産者が同居している

同じ `gene_disease` カテゴリに、粒度も想定消費者も違う2種類のレコードが入る。

| | 生産者A | 生産者B |
|---|---|---|
| 出所 | ClinGen Dosage / Gene2Phenotype（`--with-clingen-dosage` / `--with-gene2phenotype`） | `build-gene-disease-evidence`（`config/gene-disease-source-data.json`） |
| `transcript` | 無し（遺伝子レベルの主張に transcript は無い） | 有り |
| 機序フィールド | `lof_mechanism_established` のみ | ＋ `missense_mechanism_established`、`spectrum_review_complete`、`low_benign_missense_variation`、`predominantly_truncating`、`pp2_applicable`、`bp1_applicable` |
| 答えられる問い | PVS1のG01ゲート | PVS1 ＋ PP2/BP1 |

生産者Aは設計上PP2/BP1に答えられない（ミスセンスが機序かを記録していない）。ところが
選別は `curated_context()` の transcript 厳密一致が担っており、`transcript` を持たない
生産者Aを「**別の** transcript でreviewされた」として弾いていた。存在しないtranscript
不整合の解決をキュレーターに指示することになる。PVS1側の `_candidates()` は最初から
transcript無しを許容していたので、2つのマッチャが偶然だけで食い違っていた。

**対応（実装済み）**: `curated_context()` は `required_fields`（その criterion の判定が
実際に読むフィールド）で選別し、transcript は PVS1 と同じ規則にした。拒否理由は
「必要なフィールドを持たない別種の assessment」と、不足フィールド名を名指しする。
`answered_by` は、単独で criterion に答えるフィールドを別扱いにする（VCEPが
`pp2_applicable: false` と記録したレコードは機序所見を意図的に持たないため、機序
フィールドを要求すると専門パネルの判断を捨ててしまう）。

判定結果は変わらない。生産者Bを投入していれば以前から PP2/BP1 は判定に到達していた
（`work/gene-disease-20260916` が実例）。変わったのは、投入していないときにキュレーターへ
何を伝えるか。

**未解決**: 生産者Bの投入が標準の run レシピに載っていない。`config/gene-disease-source-data.json`
単体では `build_assessment_document()` が provenance 不足で失敗し、
`config/gene-disease-review-decisions.json` との組み合わせ手順が残っていない。これが
現状 PP2/BP1 が判定に到達するかどうかを決めている唯一の要因である。

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

## TogoVar経由で取得する場合の境界

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

頻度プロバイダーは設定の `population_sources` で選択する。現在登録されているオンラインプロバイダーはTogoVarで、その内部ソースとして `gnomad` と `tommo` のみを有効にしている。前者はTogoVarレスポンスの `gnomad_exomes` / `gnomad_genomes`、後者は `tommo` に展開される。選択した各データセットは、重複する可能性があるため統合せず別々の観測として保持する。TogoVar内で指定可能なグループは `gnomad`、`tommo`、`jga`、`ncbn`、`gem_j`、`hgvd`。APIが各上流データセットの版を返さないため、`upstream_dataset_version` は `NOT_PROVIDED_BY_TOGOVAR_API` と明示し、TogoVar API版、取得日時、キャッシュ本文hashを監査情報として保持する。

TogoVarとは別の頻度DBを追加する場合は、共通のnormalised evidence形式を返すアダプターを実装し、`population_registry.py` にfactoryを明示登録する。設定ファイルから任意のPythonコードをimportする方式は採用しない。`gnomad.py` は既存キャッシュの参考・互換資産として残すが、このregistryおよび通常の実行経路からは利用しない。

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

公開済みのClinGen/VCEP BS1仕様は [`config/bs1_thresholds_draft.json`](../config/bs1_thresholds_draft.json) に
review用DRAFTとして転記する。DRAFTはruntimeへ自動投入しない。各行について対象疾患・遺伝形式・対象遺伝子・gnomAD release・
AF/FAF・境界比較・最低AN/AC・例外変異をキュレーターが確認し、`APPROVED`、`reviewed_at`、版付きの根拠を記録してから
`disease_frequency_threshold`として利用する。

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
