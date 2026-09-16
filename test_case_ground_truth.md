# PP4・Layer1(16基準)・デモケース全基準ground truthテストケース収集

**version: v1**

作成:2026-09-15、2026-09-16更新
位置づけ:`ps3_bs3_ps4_implementation_v10.md`(PS3/BS3/PS4/PP1/BS4の実装設計そのもの)から切り出した別文書。こちらは**判定ロジックの実装ではなく、28基準全部を横断するground truthテストケースデータの収集経緯**を記録する。PP4とLayer1の16基準(`acmg_pipeline/criteria/stubs.py`の`AUTOMATED_RULE_BASED_OTHER_TEAM`からPP4を除いたもの)の判定ロジック自体は別メンバーの担当であり、本プロジェクト側の役割は(a)そのメンバーが実装を検証するためのフィクスチャ、(b)`classification.classify()`を将来フル28基準で試す際のテストデータ、としてground truthを整備することに限定している(実装設計との切り分けの経緯は`ps3_bs3_ps4_implementation_v10.md`セクション15-13を参照)。

成果物:`test_data/full_criteria_ground_truth.py`(データ本体)、`test_data/erepo_full_requery_2026-09-15.json`(ERepo生スナップショット)、`test_full_criteria_ground_truth.py`(健全性チェック+classify()比較レポート)。

---

## 1. PP4・Layer1(16基準)・democase変異体の全基準ground truthデータセット整備(2026-09-15)

**背景**:PP4とLayer1の16基準について、当初は判定ロジックの実装まで検討したが、**判定ロジックの実装自体は別メンバーの担当**であることが確認できたため、本プロジェクト側のスコープは**28基準全部を横断するground truthテストケースの整備のみ**に絞った。あわせて、democase配下の4変異体(MYH7 c.2155C>T, MYBPC3 c.2905+1G>A, MYBPC3 c.1000G>A, DSG2 c.1592T>G)もテストケースに追加した。

### 1-1. データソースA:既存22変異体のERepoフル再照会

`pipeline.py`の`test_cases`に既にある22種類のユニーク変異体について、`ERepoClient.lookup()`をそのまま使い、**PS3/BS3/PS4/PP1/BS4だけでなく評価済みの全コード**(`evidenceCodeStatus`の中身をフィルタせず全部)を再取得した。追加のクライアント実装は不要(既存の`ERepoClient`をそのまま再利用)。

### 1-2. データソースB:「ClinGenから逆引き」によるPP4実例の追加収集(ユーザー提案)

democaseのKCNJ5の一例(ERepoに載らない人間キュレーションのみの例)だけに頼るのではなく、**ERepo APIに`hgvs`を指定せず`gene`のみで問い合わせると、その遺伝子の全curated変異体が返る**ことを確認し、これを使ってPP4を実際に使っている変異体を広く探索した。候補26遺伝子(既存22変異体の遺伝子 + ATM, VHL, ERCC6, ERCC8, HNF1A)をスキャンした結果:

- **PP4は想定より広く使われていた**:HNF1A(24件)・RUNX1(21件)・TP53(24件)・CDH1(25件)・LDLR(25件)・DYSF(23件)・USH2A(21件)・GJB2(11件)等、多くのVCEPで採用例あり(gene単位のヒット数、strength内訳はPP4/PP4_Moderate/PP4_Strongが混在)
- 特に**HNF1A**はPP4使用例が最多——ClinGen SVI 2023年ガイダンス論文(`ps3_bs3_ps4_implementation_v10.md`セクション15-13-2で言及)がPP4の具体例として挙げているMODY確率計算機ベースの仕様を持つVCEPで、想定通りの結果
- ここから遺伝子多様性のため追加で**HNF1A×3・ATM×1**を新規変異体としてフル再照会し、データセットに加えた(既存22変異体の中にもPTEN・TP53・MTOR・CDH1・DYSF・GUCY2D・LDLR・MSH2・RUNX1・SCN2A・USH2Aが元々PP4を含んでおり、既存変異体の再照会だけでも自然にPP4データの大半が集まった)

### 1-3. データソースC:democase 4変異体+関連変異(手動転記)

`democase/real_cases_groundtruth_integrated_v6_ja.md`から、4つの標的変異体と、教材的価値の高い関連変異(**KCNJ5 c.464G>A のPP4=Reject——表現型ミスマッチの実例**、MYBPC3 c.278delA、MYBPC3 c.836del、MYH7 c.3382G>A(ノイズ、BS1))を手動で転記した。出典の確からしさ(Tier A/B/C、democase文書独自の階層)をそのままフィールドに残している。DSG2 c.1592T>Gは「単一の正解なし」という文書自体の結論を尊重し、`variant_outcome=None`のまま記録した(無理に単一の正解をでっち上げない)。

### 1-4. 成果物

- `test_data/full_criteria_ground_truth.py`:`GroundTruthEntry`データクラス(判定ロジックは一切含まない、純粋なデータ)。ERepo由来357件 + democase由来23件 = **計380件、34変異体、28基準中26基準をカバー**(PP5・BP6のみ0件——VCEPが自前の一次エビデンスで判定するため他ソースへの参照系基準はそもそも引用されにくいという傾向で、収集漏れではないと判断)
- `test_data/erepo_full_requery_2026-09-15.json`:ERepo由来357件の生スナップショット(再現手順はdocstringに記載)
- `test_full_criteria_ground_truth.py`:(1)データセット自体の健全性チェック(コードの妥当性・strength整合性・重複無し等、10件)、(2)`classification.classify()`に本プロジェクト実装済みの5基準分だけを与えて、各変異体の実際の分類(ERepoの`outcome`またはdemocase文書の最終分類)とどれだけ一致するかを**レポートする**(pass/fail判定はしない)ステップ、の2部構成
- **classify()比較の結果**:実際の分類が既知で、かつ5基準のうち最低1つでも評価済みの29変異体のうち、**一致したのは7件のみ**。残り22件は「5/28基準しか見ていないため過小評価/過大評価になる」という、事前に予想されていた通りの結果(ユーザー自身「あんまり合わないだろうけど」とコメント済み)。この結果は、`classification.py`のロジックが誤っているのではなく、**他23基準(Layer1の16基準・患者記録系6基準・PP4)が未実装であることの定量的な裏付け**として位置づける

既存227件の回帰テストに影響なし(全件パス確認済み、新規10件と合わせて計237件)。

## 2. PP5・BP6の実例収集を断念(2026-09-16)

1節時点でPP5・BP6のみground truthが0件だった件について、対象遺伝子を26→**115のClinGen VCEP関連遺伝子**(心筋症・難聴・RASopathy・遺伝性腫瘍・LGMD・てんかん・不整脈・網膜疾患・血液凝固異常等、主要VCEPを広くカバー)まで拡大して再スキャンしたが、**PP5・BP6のヒットは0件のまま**だった。

これにより、1-4節で立てた仮説(VCEPは自前の一次エビデンスで判定するため、他情報源への参照系基準であるPP5/BP6はそもそも`evidenceCodeStatus`にほぼ現れない)が、より広いサンプルでも裏付けられたと判断する。

**方針決定**:これ以上ERepoを探索してもPP5・BP6の実例は見つかりにくいと判断し、収集を断念した。代替手段(ClinVar提出者側のACMG内訳記載、CGBenchデータセットの再確認)も検討したが、**「実データで裏付けられたground truthのみ使う」という本プロジェクトの一貫した方針**(捏造しない)に照らし、無理に別ソースを探すより「28基準中26基準はカバーできたが、PP5・BP6の2基準は実データでの収集ができなかった」という結果を正直に記録することを選んだ。`test_data/full_criteria_ground_truth.py`・`test_full_criteria_ground_truth.py`はコード変更なし(元々この2基準を0件として正しく扱っていたため)。

## 3. データセットの倍増(2026-09-16、batch 2)

現状把握のために基準ごとの件数を整理したところ、ユーザーから「今の倍くらいに増やせるか」という依頼があった。

**方法**:2節のPP5/BP6探索で使った115遺伝子の候補リストを再利用し、今度は各遺伝子クエリが返す**変異体を全部**(PP5/BP6の有無だけでなく)取得した。ERepoの`gene`のみクエリは1回あたり最大25変異体を返す仕様(offset/start/skip/page/fromいずれのパラメータを試してもページングは効かないことを確認済み)。既存34変異体と未重複の候補が約1,435変異体見つかり、そのうち**まだ登場していない30遺伝子から、evidenceCodes件数が最も多い変異体を1件ずつ選出**(遺伝子多様性を優先し、既存遺伝子から追加で拾うより新規遺伝子を優先)。

**結果**:30変異体・648件を追加し、**計1028件・64変異体**に拡大。`test_data/erepo_full_requery_2026-09-16_batch2.json`として保存し、`full_criteria_ground_truth.py`の`_load_erepo_entries()`が両スナップショット(9/15分・9/16分)を結合して読み込むよう変更した。

新規追加した30遺伝子:MYO7A, SLC26A4, CDH23, MLH1, APC, CAPN3, SGCA, SGCB, SGCG, SGCD, PIK3CA, SCN1A, FBN1, KCNQ1, TECTA, OTOF, HRAS, TPM1, SCN8A, ACTC1, FOXG1, RPE65, MYL3, BMPR2, RPGR, TNNI3, MYL2, NRAS, TNNT2, KRAS(難聴・LGMD・チャネロパチー・遺伝性腫瘍・心筋症・RASopathy等、幅広いVCEPカテゴリをカバー)。

既存227件の回帰テスト・`test_full_criteria_ground_truth.py`の健全性チェック10件、いずれも全件パス確認済み(計237件)。`classify()`比較は59変異体中16件一致(母数・一致数ともに増えたが、比率は元の7/29と概ね同水準)。

## 4. 現状サマリー(2026-09-16、拡張後)

### 4-1. 回帰テスト(判定ロジック単体)

| ファイル | 件数 |
|---|---|
| `test_ps3_bs3_judgment.py` | 34 |
| `test_ps3_bs3_ps4_gate.py` | 18 |
| `test_ps3_bs3_ps4_gate_full.py` | 151 |
| `test_classification.py` | 24 |
| `test_full_criteria_ground_truth.py` | 10 |
| **合計** | **237**(全件パス) |

### 4-2. LLMパイプラインのground truth(`pipeline.py`、実装済み5基準)

22変異体・25件(PS3=15, BS3=4, PS4=2, PP1=2, BS4=2)。ClinGen ERepoでクロス検証済み。

### 4-3. 全基準ground truthデータセット(`test_data/full_criteria_ground_truth.py`)

| | |
|---|---|
| 総件数 | 1028件 |
| ユニーク変異体 | 64種 |
| データソース | ERepo 1005件 / democase手動転記 23件 |
| Tier | A 1014件 / B 10件 / C 4件 |

**基準グループ別カバレッジ**:

| グループ | カバレッジ |
|---|---|
| 実装済み5基準(PS3/BS3/PS4/PP1/BS4) | 5/5 |
| PP4 | 1/1 |
| Layer1(他チーム担当16基準) | 14/16(**PP5・BP6のみ0件、収集断念済み。2節参照**) |
| 患者臨床記録系(対象外6基準:PS2/PM3/PM6/BS2/BP2/BP5) | 6/6(※ERepoに元々含まれていたため副産物的にカバー) |

### 4-4. 基準ごとの件数(全28基準)

**病原性側**

| 基準 | 件数 | 基準 | 件数 |
|---|---|---|---|
| PVS1 | 31 | PM4 | 33 |
| PS1 | 40 | PM5 | 43 |
| PS2 | 41 | PM6 | 36 |
| PS3 | 54 | PP1 | 49 |
| PS4 | 52 | PP2 | 34 |
| PM1 | 46 | PP3 | 57 |
| PM2 | 59 | PP4 | 39 |
| PM3 | 31 | PP5 | 0 |

**良性側**

| 基準 | 件数 | 基準 | 件数 |
|---|---|---|---|
| BA1 | 42 | BP1 | 28 |
| BS1 | 44 | BP2 | 31 |
| BS2 | 31 | BP3 | 31 |
| BS3 | 38 | BP4 | 40 |
| BS4 | 34 | BP5 | 32 |
| | | BP6 | 0 |
| | | BP7 | 32 |

PM2(59件)・PP3(57件)・PS3(54件)が最多。全基準ほぼ均等に30〜60件のレンジに収まっており、PP5・BP6のみ0件という構図は変わらない。

### 4-5. classify()の実力チェック

59変異体で実際の分類(ERepo/democase記載)と比較し、本プロジェクト実装済み5基準だけでは**16/59件のみ一致**(比率は拡張前の7/29とほぼ同水準——残り23基準が未評価であることによる過小/過大評価が、データ量が増えても一貫して観測される)。
