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
