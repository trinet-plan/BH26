# PP1 / BS4 / PP4 の定義とシステム処理

作成: 2026-09-18
出典: Richards et al. 2015, Genet Med. (DOI:10.1038/gim.2015.30, PMID:25741868) / Biesecker et al. 2024, Am. J. Hum. Genet. 111:24-38 (DOI:10.1016/j.ajhg.2023.11.009)
位置づけ: `acmg_pipeline/criteria/pp1_bs4_pp4_engine.py` / `acmg_pipeline/criteria/pp4_pp1_bs4.py` / `acmg_pipeline/pp4_literature_search.py` の実装が何をしているかの参考資料。

---

## 1. ACMG/AMP 2015 原典での定義

| 基準 | 方向 | 定義 |
|---|---|---|
| **PP1** | 病的支持（Supporting） | 疾患の原因であることが確立している遺伝子において、罹患している複数の家族員で変異が共分離する |
| **BS4** | 良性支持（Strong） | 家系内で変異と疾患の共分離が見られない |
| **PP4** | 病的支持（Supporting） | 患者の表現型または家族歴が、単一の遺伝的病因を持つ疾患に非常に特異的である |

2015年時点では、これらをどの程度の強度で使うかの詳細な運用基準は示されていなかった。

## 2. ClinGen 2024ガイダンス（Biesecker et al.）による定量化

3基準を「同一の遺伝子座に対する一連の証拠」として統合的に扱う、ベイズ点数方式を提示した論文。このプロジェクトはこれを「PP1/BS4/PP4に関する現行のACMG/ClinGen実務ガイドライン」として採用している。

| 表 | 内容 |
|---|---|
| **Table 2** | 検出率（diagnostic yield）→ PP4のベイズ点数への変換表 |
| **Table 3** | 家系内の共分離者数・遺伝形式ごとのPP1/BS4点数 |
| **Table 4** | PP4+PP1合算点数 → 最終的な証拠強度（Supporting〜Strong）への変換表 |

いずれの点数も、1アレルあたり上限+5.0でキャップされる。

用語について: 「診断的yield」は一般的な日本語表現ではないため、正式な統計量の名称としては英語表記 **diagnostic yield** をそのまま用いる。日本語で言い換える場合は「検出率」「陽性率」を使う。

---

## 3. システム上の処理（`pp1_bs4_pp4_engine.evaluate()`）

### 入力
- `variant`：VCFの`GENE`アノテーションから遺伝子名を取得
- `clinical_note`：LLMが自由文から抽出した`diagnosis`（診断名）と`family.relatives`（家族の罹患・変異保有状況）

### PP4の処理

```
① clinical_note.diagnosis を検索クエリにして、その場でPubMed文献検索を実行
② 該当論文からLLMが検出率（diagnostic yield）を抽出
③ Table2で点数化
④ 見つからなければ正直にUNKNOWN（推測しない）
```

事前キュレーションDBは持たず、毎回ライブ検索する（PS3/BS3/PS4の文献ワークフローと同じ方式）。

#### 文献からの検出率抽出の詳細（`pp4_literature_search.search_diagnostic_yield()`）

1. 遺伝子名＋診断名でPubMed検索し、候補論文のPMIDを取得
2. 候補論文の本文（PMC収載分）またはアブストラクトを取得
3. 論文本文をLLMに読ませ、以下をJSON形式で抽出させる

   | フィールド | 内容 |
   |---|---|
   | `yield_percent_stated` | 論文に明記された数値（推測・計算は禁止） |
   | `is_overall_yield` | **検査対象者全体**に対する割合か、**既に何らかの変異が確定した患者群内**での内訳かの判定 |
   | `denominator_description` | 分母の説明 |
   | `sample_size` | 母数（分母） |
   | `quote` | 根拠となった原文引用 |
   | `no_yield_statistic_found` | 統計量が論文に存在しない場合はtrue |

4. `is_overall_yield: true` かつ数値が妥当な最初の候補のみ採用。それ以外は次候補へ、全て不採用ならUNKNOWN

**分母の罠**: 検出率は「疾患Xで検査した患者全員」を分母にすべきだが、論文によっては「既に何らかの変異が確定した患者群」を分母にした内訳（例：「MYH7は遺伝学的に解決済み症例の33%を占める」）を報告している場合があり、これを誤って検出率として使うと実際の値より大きく過大評価する（実験で3倍以上の差を確認）。`is_overall_yield`の確認はこれを防ぐための必須チェック。

**既知の弱点**: 分子（変異保有者の実数）は独立したフィールドとして抽出しておらず、LLMが計算済みの`yield_percent_stated`をそのまま採用している。分子も抽出してこちらで計算し直すクロスチェックは未実装。

### PP1/BS4の処理

```
family.relativesの各家族員について：
  ・罹患状況＋変異保有状況の組み合わせ
  ・遺伝形式（AD/AR/X連鎖）
  を見て、Table3に基づき加点/BS4判定
```

**PP4の文献検索とは完全に独立**して計算される。PP4の検索が失敗（診断名なし、または検出率が見つからない）してもPP1/BS4は家族歴データがあれば評価される。

### 合算・出力

```
PP4点数 ＋ PP1点数使用分 → Table4で強度に変換 → CriterionEvidence → VA-Spec EvidenceLine
```

出力に必ず付く注記：
- PP4がヒットした場合：「文献検索による自動抽出、未確認、要人間確認」という開示を`source`に必ず含める（PS3/BS3の`strengthEstimationMethod`と同じ考え方）
- 検出率の母数が20人未満の場合：判定はブロックせず、注意書きを追加

### 設計方針（事前キュレーションDBを持たない理由）

当初は`config/pp4_reference_records.json`という手動キュレーション登録簿（DRAFT/APPROVED方式）で設計したが、以下の理由から廃止し、常時ライブ検索方式に統一した。

- このツールに「AIは下書きを出すだけ、人間が承認するまで使わない」という原則は無く、PS3/BS3/PS4も含め、LLMの判定はそのままCriterionEvidenceとして提示し、開示情報で「未確認である」ことを示す方式で統一されている
- 人間が事前に入力しないと動かない登録簿は、「前例の無いケースの判定を助ける」というツールの目的そのものと矛盾する
- 検索結果を登録簿へ自動キャッシュする中間案も検討したが、PS3/BS3/PS4の非永続キャッシュ方針と一貫性が無く、PP4の発生頻度自体が低いため都度検索のコストは許容範囲と判断し不採用とした

---

## 4. 精度検証（2026-09-18、参照元：PubMed）

`pp4_literature_search`のLLM抽出・検索ロジックについて、実際のPubMed MCP＋ローカルLLMを使い、正解が既知の事例で検証を行った。**本節（4-1・4-2）の検証は、本番システムが実際に使っているPubMed（PubMed MCP経由）を参照元とした結果**。GeneReviewsとの比較は5章を参照。

### 4-1. 抽出ロジックの精度検証（分母判定は正しく機能するか）

正解を事前に独立して確認した4件のPMIDで、`_judge_paper_for_yield()`の判定を照合した。

| PMID | 遺伝子/表現型 | 期待される判定 | 実際の判定 | 結果 |
|---|---|---|---|---|
| 12442267（Kalatzis et al. 2002） | CTNS / cystinosis | 採用（95.8%, n=108、全体に対する割合） | 採用（95.8%, n=108、完全一致） | ✅一致 |
| 41998504 | MYH7 / HCM | 却下（解決済み症例内の内訳） | 却下（`is_overall_yield: false`） | ✅一致 |
| 28000701 | GJB2 / 難聴 | 却下（陽性67例内の13.4%） | 「統計なし」判定 | ⚠️不一致（原因判明、下記） |
| 41488457（Wang et al. 2025） | MYBPC3 / 新生児期致死性HCM | 「統計なし」（単一症例報告） | 「統計なし」 | ✅一致 |

**4件中3件が期待通り。** 唯一の不一致（GJB2）は、LLMの誤判定ではなく、**論文の全文取得が失敗し（既知の別不具合）、アブストラクトにフォールバックした結果、そもそも13.4%という数字自体がテキストに含まれていなかった**ことが原因。LLMは見ていない数字を捏造せず、正しく「統計なし」と回答していた。→ **分母判定ロジック自体の誤りは、今回の検証では0件。**

### 4-2. 検索カバレッジ調査（実際にどれだけ見つかるか）

実際の`search_diagnostic_yield()`（検索・本文取得・LLM抽出をフェイク無しで通し実行）を、多様な8組の遺伝子×表現型で実行。

| 遺伝子 | 表現型 | 結果 | 出典 |
|---|---|---|---|
| MYH7 | 肥大型心筋症 | 見つからず | - |
| FBN1 | マルファン症候群 | **見つかった** | 全文（93%） |
| DMD | デュシェンヌ型筋ジストロフィー | **見つかった** | アブストラクトのみ（73%, n=249） |
| SMN1 | 脊髄性筋萎縮症 | 見つからず | - |
| BRCA1 | 遺伝性乳がん卵巣がん | **見つかった** | アブストラクトのみ（12%, n=101） |
| PKD1 | 常染色体優性多発性嚉腎 | 見つからず | - |
| GJB2 | 難聴 | 見つからず | - |
| CTNS | シスチン症 | 見つからず（下記参照） | - |

**8件中3件（37.5%）で確定的な検出率を取得。** うち2/3件は本文が読めず**アブストラクトのみで取得**できており、アブストラクトだけでも十分な情報が載っている論文が一定数存在することを確認。

**注目すべき発見**：CTNS×シスチン症は4-1で「95.8%の正解論文（PMID:12442267）が存在する」ことを確認済みだが、本調査の検索では**その論文自体が候補に出てこなかった**（別の新しい4論文がヒットし、いずれも検出率の記載なし）。これはLLM抽出の誤りではなく、**PubMed検索のクエリ設計・候補選定（リコール）の問題**であることが判明。

---

## 5. GeneReviewsとの比較検証（2026-09-18、参照元：GeneReviews）

4-2と同じ8組の遺伝子×表現型について、**PubMedではなくGeneReviews（NCBI Bookshelf）を参照元**として、同じ分母チェック基準（検査対象者全体に対する割合か／既に変異が確定した患者群内での内訳か）で確認した。

注記：GeneReviewsは本番の`pp4_literature_search.py`からは参照できない（PubMed MCP経由の`get_full_text_article`はPMC収載論文のみが対象で、NCBI Bookshelfは対象外）。本検証はWebブラウジングツールで個別に確認したものであり、本番システムの検索結果ではない。

| 遺伝子 | 表現型 | GeneReviewsでの結果 | PubMed（4-2の結果） |
|---|---|---|---|
| MYH7 | 肥大型心筋症 | ❌却下（33%は解決済み症例内訳、30%は遺伝子非特定） | 見つからず |
| **FBN1** | **マルファン症候群** | ✅**採用**（93%、n=93） | ✅採用（93%、全文） |
| DMD | デュシェンヌ型筋ジストロフィー | ❌却下（テスト手法内の内訳のみ） | 見つかった（73%, n=249） |
| SMN1 | 脊髄性筋萎縮症 | ❌却下（検査感度・解決済み症例内訳のみ） | 見つからず |
| BRCA1 | 遺伝性乳がん卵巣がん | ❌却下（解決済み症例内でのBRCA1/BRCA2内訳のみ） | 見つかった（12%, n=101） |
| PKD1 | 常染色体優性多発性嚉腎 | ❌却下（解決済み症例内での内訳のみ） | 見つからず |
| GJB2 | 難聴 | ❌却下（検査感度の記載のみ） | 見つからず |
| CTNS | シスチン症 | ❌却下（手法別の内訳のみ、95.8%の記載は無い） | 見つからず |

**GeneReviews: 8件中1件（12.5%）が採用可能。PubMed: 8件中3件（37.5%）。**

### 分かったこと

1. **GeneReviewsの方がPubMedより成績が悪かった**（今回の8件では）。ClinGen論文が「GeneReviewsの変異yield表が良い出典」と推奨しているのとは逆の結果。
2. 原因：**GeneReviewsの標準的な"Table 1"形式自体が、ほぼ常に「解決済み症例内での遺伝子別内訳」を報告する形式**になっており、そのままでは検査対象者全体に対する検出率として使えない（MYH7・PKD1・CTNSで確認）。
3. 唯一採用できたFBN1も、**表形式の数字ではなく、本文の文章中に埋め込まれた記述**（"86 (93%) of 93"）だった。しかもこの数字はPubMedで見つけたのと**完全に同じ93%**——GeneReviewsも結局は元の研究論文（Loeys et al. 2004）を引用しているだけで、一次情報源としては同じところに行き着く。
4. 検証の過程で、質問文（プロンプト）の言い回しが曖昧だと検証者（LLM）自身も分母を誤判定しうることを実際に確認した（FBN1の初回検証で誤って却下→言い回しを精緻化して再検証し訂正）。分母の定義（「臨床的にその表現型と診断されテストを受けた人全員」か「既に何らかの遺伝子で変異確定済みの人」か）を明確に区別する指示の重要性を裏付ける結果となった。

### 結論

GeneReviewsを本番システムに追加接続する優先度は、当初想定より低い。それよりも、4-2で見つかった**PubMed検索クエリ自体の改善（正解論文を候補に挙げられていない問題）**の方が、投資対効果が高いと考えられる。

### 4-3. まとめ

| 観点 | 状況 |
|---|---|
| 精度（見つけた時に正しいか） | 良好。分母判定ロジックの誤りは今回の検証では検出されず |
| 再現率（そもそも見つかるか） | 8件中3件（37.5%）。原因は主に①本文取得できない論文が多い、②検索で正解論文自体を拾えていない、の2種類 |

次の改善余地：検索クエリの工夫（②への対応）、複数候補を試す際のクエリバリエーション追加。

---

## 6. ガイドラインとの乖離点の洗い出しと対応（2026-09-18）

`pp4_pp1_bs4.py`/`pp1_bs4_pp4_engine.py`をBiesecker et al. 2024の本文・Table 3・Table 4と改めて突き合わせ、3件の乖離を発見した。

### 6-1. X連鎖劣性（XLR）で保因者女性が一切カウントされていなかった（修正済み）

`_score_family_segregation()`はXLRのPP1/BS4を`father`/`son`/`brother`等の「既知の男性関係」を持つrelativeにしか与えておらず、Table 3脚注e「Additional segregations can be counted for obligate heterozygous females」および論文自身のFigure 6のワークド例（非罹患の保因者女性に+1.0点）と矛盾していた。加えて、AD/AR-unaffected行だけに付く脚注aの「親を除外する」ルールが、脚注aの付いていないXLR行にも誤って適用されていた。

**対応**: 非罹患・variant保有・既知の女性関係を持つrelativeをTable 3脚注e通りPP1(+1.0)としてカウントするよう修正し、親除外ルールをXLR以外に限定した（[pp4_pp1_bs4.py](../acmg_pipeline/criteria/pp4_pp1_bs4.py)）。

### 6-2. `locus_model`が常に"heterogeneous"固定だった（yield閾値による近似で対応済み）

PP4の文献検索由来のreferenceは常に`locus_model="heterogeneous"`に固定されており、論文の中心的主張である「locus homogeneity + 診断的yield>90%の場合はPP1をPP4に追加してはいけない（二重計上になる）」というルール（`evaluator.py`の`high_yield_homogeneous`判定）が本番では絶対に発火しない状態だった。

**論文の厳密な定義との差**: 論文の"Summary of this heuristic approach"のstep 2は「表現型がsubstantially specificで、かつ診断的yieldが90%を超えている」ことを求めており、さらにstep 2.bは「locus homogeneityだが診断的yieldが低い（＝他に原因遺伝子があるかもしれない）」というケースを明示的に想定している。つまりlocus homogeneityは診断的yieldの数字だけから機械的に導出できる値ではなく、本来は「他にこの表現型の原因になる遺伝子が知られていない」という、yieldとは別の文献的事実（CTNS・FBN1の例で論文が別途明記している）が必要。

**採用した対応（ユーザー承認、2026-09-18）**: 厳密な「本当にこの遺伝子だけが原因か」の確認には新しいLLM文献判定（新たなハルシネーション面）か登録簿（このプロジェクトが明示的に廃止した設計）が必要になるため、簡易的な近似として、**取得したyieldが90%を超えていたら`locus_model="homogeneous"`とみなす**方式を採用した（[pp1_bs4_pp4_engine.py](../acmg_pipeline/criteria/pp1_bs4_pp4_engine.py)の`_build_reference_from_literature()`）。

```python
locus_model = "homogeneous" if result.yield_fraction > 0.90 else "heterogeneous"
```

**この近似の実害範囲は限定的**: `DIAGNOSTIC_YIELD_POINT_TABLE`はyield 81.6%の時点で既にPP4単独で+5.0点（1アレル上限）に達するため、この近似が発火する90%超のケースでは、PP4だけで既にcapに達している。したがって近似が誤っていても、PP1+PP4の合計点・最終分類は変わらず、変わるのは「PP1が個別にMET表示されるか、NOT_MET表示されるか」という表示上の正確さのみ。

### 6-3. PP1/BS4評価に必要なパラメータが本番呼び出しで常にNone（未対応、今後の検討課題）

`pp1_bs4_pp4_engine.evaluate()`は`evaluate_locus_evidence()`を呼ぶ際、`inheritance_mode`/`ar_case_mode`/`fully_penetrant`/`low_phenocopy`を全て`None`固定で渡している。

- `inheritance_mode`：`clinical_note.family.inheritance_pattern`にフォールバックするが、`clinical_extraction.py`のプロンプトは「ノートに明記されている場合のみ埋める」方針のため、実際の症例ノートではほぼ常にNoneのまま（デモ4症例全てで確認）。
- `ar_case_mode`：常にNoneのため、AR（常染色体劣性）症例のPP1は本番で一切スコアされない。
- `fully_penetrant`/`low_phenocopy`：常にNoneのため、**BS4（非分離）は本番では実質発火しない**。

**検討した対応案（未実装）**:
- `inheritance_mode`：今回のマージで入った`gene2phenotype.py`（G2P、PVS1の遺伝子-疾患知識ベース）が遺伝子×MONDO疾患ペアごとに継承様式を持っているため、新たに実装した`diagnosis`→MONDO変換と組み合わせて機械的に取得できる可能性がある。ただしG2Pの`x_linked`はdominant/recessiveを区別しないため、`X_LINKED_RECESSIVE`への単純マッピングは近似になる。
- `ar_case_mode`：`clinical_extraction.py`が既に抽出している`proband.genotype.zygosity`（"homozygous"/"compound_heterozygous"等）から直接導出可能。新しいデータソース不要の低リスクな修正。
- `fully_penetrant`/`low_phenocopy`：既存のデータソースが無く、LLM文献判定（新規ハルシネーション面）か登録簿（廃止済み設計）のどちらかが必要になる。優先度は低いと判断し保留。

---

## 7. 診断名のMONDO変換（2026-09-18、新規追加）

`diagnosis`（自由文の診断名）をHPO変換と同じ仕組みでMONDO疾患IDに変換する`resolve_diagnosis_mondo()`を`acmg_pipeline/hpo_mondo_extraction.py`（旧`hpo_extraction.py`をリネーム）に追加した。TogoMCPの`mondo`データベース（HPOとは別データベース）に対し、SPARQL検索→LLM最終判定という同じパターンで実装。結果は`ClinicalNoteExtraction.mondo_id`に格納される。

### 検証結果

デモ4症例の診断名で、実際のTogoMCP+LLMを3回ずつ（計12回）実行し、EBI OLS4（`https://www.ebi.ac.uk/ols4`、TogoMCPとは独立した公式API）で確認した正解データと比較した。

| 診断名 | 正解（EBI OLS4） | 結果（3回とも） |
|---|---|---|
| hypertrophic cardiomyopathy | MONDO:0005045 | MONDO:0005045 ✅ |
| hypertrophic cardiomyopathy (HCM) complicated by a left ventricular apical aneurysm | MONDO:0005045 | MONDO:0005045 ✅ |
| arrhythmogenic right ventricular cardiomyopathy | MONDO:0016587 | MONDO:0016587 ✅ |

12回全てが正解データと一致し、実行ごとのブレも無かった。

### 開発中に見つかった2件のバグ（判定精度ではなく、検索クエリ自体の構文エラー）

検証の過程で、LLMが生成するSPARQLクエリ自体が失敗するケースが2件見つかった。いずれもLLMの意味判定（ハルシネーション）ではなく、SPARQL構文レベルの問題だった。

1. **`PREFIX bif:`宣言がVirtuosoの予約語と衝突**：プロンプトで`bif:contains`に言及したところ、LLMが（未使用でも）`PREFIX bif: <...>`を宣言し、Virtuosoエンドポイントがクエリ全体をHTTP 400で拒否した。→ `bif:contains`への言及を削除し、標準の`CONTAINS(LCASE(?label), ...)`のみを指示するよう修正。
2. **`ORDER BY STRLEN(?label) ASC`はSQL構文であり有効なSPARQLではない**：正しくは`ORDER BY ASC(STRLEN(?label))`。プロンプトで「意図」を英語で説明するだけでは誤った構文を生成したため、正確な構文を直接明記するよう修正。

「ラベルが短い順に並べる」という並び替え自体は、"hypertrophic cardiomyopathy"のような一般名が"hypertrophic cardiomyopathy 4"のような遺伝子座別サブタイプ名（20件以上存在）に候補リストの`LIMIT`枠から押し出されてしまう問題（4-2で見つかったHPO検索の`myocardial necrosis`と同型の問題）への対策として追加した。
