# PS3 / BS3 / PS4 / PP1 / BS4 実装設計書

**version: v10**(v9から、以下を追加実装・確定。詳細はセクション15を参照:
1. `ERepoClient`が実際にHTTP接続できることを確認し、本番実装を完了(9-4節の「ネットワーク制約」は本セッションで解消)
2. 出力形式をGA4GH VA-Spec `EvidenceLine`(JSON)に確定(打ち合わせ決定、5節を置き換え)、`va_spec_export.py`として実装
3. curator向けhintの重複バグ・曖昧な文言を修正
4. 共通インフラを`evidence_common.py`に汎用化し、`JudgmentEngine`パターンで新基準を追加しやすい構造に整理
5. 文献読解で判定可能な残り2基準、**PS4(症例対照)・PP1/BS4(家系内分離)**を新規実装
6. ground truthの再監査で、8-9節の36件監査対象外にも同種のCGBenchラベル付けバグが実在することを確認・修正
実装ファイルは`ps3_bs3_ps4_gate.py`・`ps3_bs3_judgment.py`・`ps4_judgment.py`・`segregation_judgment.py`・`evidence_common.py`・`va_spec_export.py`・`ps3_bs3_llm_pipeline.py`および付随テスト一式として別途提供)

作成:2026-09-14(v9)、2026-09-15更新(v10)
位置づけ:Layer3のうち「文献読解が必要」な基準について、AcmGENTIC(Saadat & Fellay, arXiv:2604.00075)・CGBench(Queen, Zhang, Zou, NeurIPS 2025, arXiv:2510.11985)の実証結果と、Anthropic公式PubMedコネクタでの実地検証を踏まえた実装設計として出発し(セクション0〜5)、その後CGBench実データセットを用いた大規模検証(セクション6〜11)、実装(セクション13〜14)、そしてVA-Spec出力確定・基準拡張(セクション15、v10)へと発展した、一連の技術調査の全記録。当初はPS3/BS3/PS4の3基準を対象としていたが、v10でPP1/BS4を追加し、doc/BH26_participant_briefing_v3_en.mdのLayer3分類のうち「実際に論文を読んで判定する」5基準(PS3, BS3, PS4, PP1, BS4)全てを対象とするに至った(15-5節参照)。

---

## 0. 現実的な期待値の調整(CGBenchの知見、重要)

実装設計に入る前に、既存研究が示す「どこまで難しいか」を正直に共有する。

- **CGBenchのVCI E-Score(tertiary code、PS3/PS4含む全基準横断)では、8モデル全てが不正解だった症例が74.6%に達した**。最良モデル(o4-mini)でもPrecision@5は0.420
- **PMC全文が公開されている引用論文は、VCIで約20%、GCIで約23%のみ**——今回私たちがPubMedコネクタで確認した実感(全PubMedの1/6程度がPMC収録)と一致する、構造的な制約
- **「推論モデルは細かいタスクに強いが、高レベルの解釈は非推論モデルの方が良い」**——単純に「推論モデルが常に優れている」とは言えない、タスクの粒度依存性がある
- **CGBench自身の限界として、補足資料(supplementary)・図表を読めないことが明記**されている。AcmGENTICのdirect modeも同様の限界を持つ

**この設計の位置づけ**:上記を踏まえ、**「完全自動化」ではなく「人間キュレーターの一次スクリーニングを高速化する下書き生成」**として設計する。`not_clear`を恥じずに多用する設計にする。

## 1. 設計方針(AcmGENTIC・CGBenchの知見、および今回の全検証を踏まえた4つの制約)

1. **方向性(damaging/normal)のみを狙う。強度(Very Strong〜Supporting)は狙わない**——AcmGENTICで方向性は96%(推論モデル)だが、強度は33〜36%しか出ない
2. **変異の同定を明示的なゲート(必須の中間ステップ)にする**——AcmGENTICでも26〜30%が変異同定に失敗。ここを飛ばして直接方向性判定に進むと誤帰属のリスクがある
3. **判定不能な場合は`not_clear`(判定保留)を正式な出力値として認める**——無理に断定せず、人間キュレーターに戻す設計。CGBenchの「74.6%不正解」を踏まえると、この設計は必須というより前提
4. **【v8で基本方針に格上げ】文献読解パイプラインを起動する前に、対象変異が既にClinGenでキュレーション済みかどうかを必ず確認する**——詳細は次節。今回発見した問題(MYH7 R719WのPS4騒動、met_statusの系統的バイアスに見えた誤差、RASopathyアッセイ表とRIT1の食い違い)は、いずれも「ClinGenに既に答えがあるのに、それを見ずに独自に再導出しようとした」ことに起因していた

## 2. 変異照合ゲートと検証/本番の二重モード設計(v8、基本方針)

### 2-1. 背景

4-1節(旧)の「手順0」(PS4向けに「ClinGen VCEPレビュー済みの変異なら、まずERepoの完全なエントリに直接当たる」)は、当初PS4専用の応急処置だった。しかし7〜11節の全検証を通じて、**この発想はPS4に限らずPS3/BS3を含む全基準に適用すべき、パイプライン全体の基本設計**であることが明らかになった。既にキュレーション済みの変異について、AIが論文を読んで独自に再導出する必要は本来ない。

### 2-2. 変異照合ゲート(全基準共通、最優先)

```
0. 変異照合ゲート(PS3/BS3/PS4のいずれの判定にも先立って実行)
   対象変異をClinGen Allele Registry(CAID)またはERepo API
   (erepo.clinicalgenome.org/evrepo/api/classifications?gene=<遺伝子>&hgvs=<変異>)
   で照合する
   → 既にVCEPキュレーション済みなら、そのMet/Not Metコードをそのまま採用し、
     3〜4節の「文献から独自導出」フローには進まない
   → 未収載なら、初めて3〜4節のフローに進む
```

**照合には正規化されたIDを使うこと**。7-2節で、生のテキスト一致による照合が座標1塩基ずれ・括弧表記ゆれ等で失敗する例を複数確認した。ClinGen Allele Registry(CAID)はこの正規化を既に解決済みのため、自前でHGVS表記のゆれを吸収するロジックを組むより、CAID経由の照合を優先すべき。

### 2-3. バリデーション・パラドックスとその解決(二重モード)

このゲートには自明な問題がある——**ゲートが有効な状態では、キュレーション済みの変異でロジックが一度も実行されないため、ロジック自体の精度検証ができない**。かといって、精度を検証したいなら「正解が分かっている変異」、すなわちキュレーション済みの変異でロジックを走らせるしかない。「本当に価値があるのは未収載の新規変異への適用」だが、そこには比較対象の正解が存在せず検証もできない、というパラドックスになる。

これは機械学習における「学習・評価用データ」と「本番の未知データ」の区別と同じ構造であり、**モードを分離することで解決する**:

| モード | 対象 | ゲートの扱い | 目的 |
|---|---|---|---|
| **検証モード**(ブラインド評価) | 既にClinGen収載済みの変異を意図的に選ぶ | あえて無効化し、正解を見せずに論文本文だけをロジックに渡す。判定が終わってから初めてClinGenの公式コードと突き合わせる | ロジックの精度を定量的に測る |
| **本番モード** | 実際のキュレーション対象(未知・既知混在) | 常に有効。既知の変異はClinGenの答えをそのまま使い、未収載の変異だけロジックを起動する | 実務での効率化 |

**7〜10節で行った検証は、まさにこの「検証モード」そのものだった。** CGBenchの289行から抽出した既知変異について、met_statusを見ずに(あるいは見た後で答え合わせのためだけに)本文から方向性を判定し、事後にClinGenの公式コードと突き合わせる、という手順を踏んでいた。ゲートを本番パイプラインに組み込むことと、この検証手法とは矛盾しない——**ゲートは本番モードの設計であり、検証はゲートをあえて無効にして行う別の作業**、という整理になる。

### 2-4. BH26のデモへの含意

この二重モード設計は、デモの見せ方としても有効である。「既知の変異についても、あえてAIに(答えを見せずに)独自に判定させ、公式記録と何%一致するか」を、7〜10節で実際に行った形式でライブで見せれば、ロジックの信頼性を定量的に示しつつ、本番では効率的にゲートを使う、という両立を視覚的に伝えられる。CGBenchの「74.6%不正解」という数字の中身を初めて具体化した今回の調査結果そのものが、このデモの土台になる。同時に、「未収載の新規変異にだけAIが挑戦し、既知の変異は正直にClinGenに譲る」という設計は、藤原さんの論文が主張する「過信を避け、限界を明示する」という姿勢とも一貫する。

## 3. PS3 / BS3 共通パイプライン

### 2-1. 候補論文の特定(優先順位付き)

```
1. ClinVar/ClinGen Evidence Repositoryが既にPMIDを引用している場合
   → そのPMIDを最優先で使用(高信頼度・低コスト。Case2-var1のPS4、Case3-var1のPS3で実証済み)
2. 引用が無い場合、PubMed:search_articles で「遺伝子名+変異名(HGVSp)+動物モデル名」を検索
   → 「変異の考察内容を言い換えたキーワード」では0件になることを実証済み。
     タイトル・変異表記に近い語で検索すること
3. 見つからない場合 → insufficient_data(文献が存在しない/見つけられない)
```

### 2-2. 変異同定ゲート(必須、スキップ禁止)

論文がターゲット変異を実際にテストしているかを、以下の当量表記(equivalents set)と照合する:

| 照合方法 | 例 |
|---|---|
| タンパク3文字表記 | `Arg719Trp` |
| タンパク1文字表記 | `R719W` |
| 遺伝子内の略称 | `RW`(論文によっては変異名を略記号で通す場合がある。今回の実例で実際に必要だった) |
| rsID | あれば最優先 |
| HGVSc | `c.2155C>T` |

**`match_status`は`matched` / `heuristic` / `single_variant_study` / `unsuccessful`の4値**。`unsuccessful`なら以降の処理を行わず`not_clear`で確定する。

### 2-3. 全文取得

```python
# 1. PMIDからPMC IDへの変換(必要な場合)
PubMed:convert_article_ids(ids=[pmid], id_type="pmid")
# 2. 全文取得
PubMed:get_full_text_article(pmc_ids=[pmc_id])
```

PMC収録が無い場合(CGBenchの実測でも約20%のみPMC収録)、Abstractのみで判定を試み、`confidence`を`low`に落とす。補足資料(supplementary)・図表中のデータはこの設計では読めない(AcmGENTIC・CGBench共通の限界)。

### 2-4. 実験抽出プロンプト(構造化出力)

```
あなたはACMG/AMP PS3/BS3の機能的証拠評価を支援する臨床遺伝キュレーターです。

対象変異: {gene} {hgvsc} ({hgvsp}, 別名: {equivalents})
論文全文: {full_text}

以下の手順で判定してください:
1. 変異同定:この論文は対象変異を直接テストしているか?
   match_status: matched / heuristic / single_variant_study / unsuccessful
   ※unsuccessfulの場合は以降を全てnullとし、ここで終了すること

2. 実験抽出(matched/heuristic/single_variant_studyの場合のみ):
   各実験について: assay_type, experimental_system(細胞株/マウス/患者由来等),
   readout, comparator(対照), result_direction
   (functionally_abnormal/functionally_normal/intermediate/mixed/unclear)

3. モデル系統の注記:
   in vivo動物モデル(特にマウスノックイン)の場合、model_system_caveatとして
   「ヒト組織由来の直接的証拠ではない」旨を明記すること。人間キュレーターは
   マウスモデルの評価に確信度が低いと繰り返し表明している(Case3・Case4の実例)

4. 総合判定:
   direction: PS3(damaging) / BS3(normal) / not_clear
   strength_hint: not_clear固定(強度判定はこのシステムでは狙わない)
   rationale: 2-3文の根拠

JSON形式で出力してください。
```

### 2-5. 実例による検証(Case3-var1、実際に実行した結果)

対象:MYH7 c.2155C>T (p.Arg719Trp)、PMID:24829265(PubMedコネクタで実際に全文取得済み)

```json
{
  "variant_matching": {
    "match_status": "matched",
    "match_type": "protein_notation",
    "confidence": "high",
    "notes": "論文の主対象は別の変異(V606M)だが、対象変異(R719W)は同一論文内で
              新たに作成・特性解析された比較用マウス系統として明示的に報告されている"
  },
  "experiments": [{
    "assay_type": "knock-in mouse model, echocardiography + histology",
    "experimental_system": "Mouse (Arg719Trp knocked into endogenous alpha-MHC gene)",
    "result_direction": "functionally_abnormal",
    "key_findings": [
      "RW/+マウス:26週齢で進行性の求心性肥大(壁厚がWTより20%以上増加、p<0.05)、心筋線維乱れ、間質線維化",
      "RW/RW(ホモ接合)マウス:生後9日以内に100%死亡、ANP/BNP著明上昇"
    ],
    "model_system_caveat": "マウスノックインモデル。血行動態低下→交感神経/RAAS活性化→心筋リモデリングという機序はマウス表現型からの推論であり、ヒト組織からの直接的証拠ではない"
  }],
  "overall_evidence": {
    "direction": "PS3",
    "strength_hint": "not_clear",
    "rationale": "ホモ接合マウスは新生児期に重度の心肥大マーカーを伴い死亡、ヘテロ接合マウスはヒトのRW/+表現型(肥大・線維化・線維乱れ)を再現。PS3(有害な影響)を支持する"
  }
}
```

**この結果は、既存のGround Truth(Tier A、PS3、確信度に留保あり)と完全に一致**した。パイプラインが実際に機能することを確認済み。**ただし、CGBenchの「74.6%不正解」という基準率を踏まえると、この1件の成功だけで一般化はできない**点は留意する。

---

## 4. PS4 専用パイプライン

PS3/BS3と異なり、PS4は「症例対照/コホート研究」という特定の研究デザインを扱う。

### 3-0. 既存の実装・検証状況(訂正)

- **AcmGENTIC**:PS3/BS3専用とタイトル・スコープで明記。PS4は対象外
- **CGBench**:VCIタスク(E-Score/E-Ver)は引用PMIDがあれば全基準を評価対象にしており、**PS4が含まれている可能性は高いが、個別の精度は論文中で報告されていない**(全tertiary codeの集計値のみ)。GCI(実験的証拠抽出)の13カテゴリは全て機能実験系で、症例対照(PS4)に相当するカテゴリは存在しない
- **結論:PS4単独の検証実績を持つ既存研究は見つかっていない。** 以下の設計は、PS3/BS3の手法をPS4に類推適用した、未検証のアイデアという位置づけを維持する

### 3-1. 候補論文の特定

```
0. 【2026-09-14追加、最優先】ClinGen VCEPレビュー済みの変異の場合、
   まずClinGen Evidence Repository(erepo.clinicalgenome.org)の
   完全なエントリページに直接当たる。基準ごとの引用論文リストが
   全て確認できるため、二次資料(人間キュレーション等)が抜き出した
   単一のPMIDだけを信用しない
1. ClinVarの引用PMIDを最優先(Case2-var1で実証済み:VCV000042666がPMID:27532257を引用)
2. 見つからない場合、PubMed:search_articlesで「遺伝子名+変異名+"case-control" OR "cohort" OR "families"」を検索
3. 見つからない場合 → insufficient_data
```

**この手順0が必要になった経緯**:当初、Case3-var1のPS4を「PMID:19645038単独では8人という主張を支持しない」としてRejectに訂正したが、これは二次資料(人間キュレーションが抜き出した1個のPMID)だけを検証した結果だった。ClinGen ERepoの一次情報源を直接確認すると、**PS4は実際には9本以上の論文を束ねた「20人以上で報告」という根拠で正式に成立しており、「8人」はPP1(別基準)の、別の論文群からの数字**だったことが判明した(詳細は`real_cases_groundtruth_integrated_v8`参照)。**二次資料の1点だけを検証して「不成立」と判定するのは、一次資料に当たっていない分、過剰な訂正になるリスクがある**。

### 3-2. 抽出プロンプト(PS4専用)

```
あなたはACMG/AMP PS4の症例対照研究評価を支援する臨床遺伝キュレーターです。

対象変異: {gene} {hgvsc} ({hgvsp})
論文全文: {full_text}

1. 変異同定:(PS3/BS3と同じゲート)

2. 研究デザインの確認:
   study_design: case_control / family_cohort / case_series / not_case_control
   ※not_case_controlの場合はPS4に使えないためnot_clearで終了

3. 症例対照データの抽出:
   affected_carriers: 罹患者における当該変異の保有者数
   unaffected_carriers: 非罹患者(対照群)における保有者数
   total_affected: 罹患者の総数
   total_unaffected: 対照群の総数
   (数値が明記されていない場合はnullとし、記述のみをrationaleに残す)

4. 総合判定:
   direction: PS4(有意に増加) / not_clear
   strength_hint: not_clear固定
```

### 3-3. 実例(Case3-var1、2026-09-14に再検証した正しい内容)

**当初の記録(誤り)**:「ClinVar引用なし(ClinGen Evidence Repositoryが直接評価)、PS4はPMID:19645038を独自に引用。抽出結果:常染色体優性HCMの中国人非血縁6家系を調査し、罹患者8人中8人で変異検出→PS4」

**実際に確認した内容**:ClinGen ERepoの完全なエントリ(https://erepo.clinicalgenome.org/evrepo/ui/interpretation/7b17a8bd-b169-46a1-8efa-b7388c4ecdef)を直接取得すると、PS4は**9本以上の論文(PMID:9829907, 8282798, 9822100, 12974739, 22429680, 23816408, 12707239, 19645038, 27532257, SHaRe consortium/30297972)を束ねた「20人以上の罹患者で報告」**という根拠だった。PMID:19645038はその中の1本に過ぎず、単独で読むと「6家系中1家系のみ連鎖、Arg719Trp単独ではなくAla26Val+Arg719Trpの二重変異として発見」という、より限定的な内容だった。**「8人」という数字は、実際にはPP1(家族内共分離)の基準で、別の論文群(PMID:9829907等)から来ていた**。

**教訓**:PS4のように**複数論文の集約で成立する基準**の場合、パイプラインが「1本の論文を読んで判定する」設計だけでは、ERepoが実際に束ねている根拠の全体像を再現できない。**複数PMIDが列挙されている場合は、その全リストを候補論文として扱い、個々の論文が集合としてどう積み上がっているかを見る設計に拡張する必要がある**(現状の3-1・4-1節の設計は、この「複数論文の集約」パターンを想定していなかった)。

---

## 5. VA-Specへの統合

上記パイプラインの出力を、そのままEvidence Lineへマッピングする:

| パイプライン出力 | VA-Spec field |
|---|---|
| `overall_evidence.direction`("PS3"/"BS3"/"PS4") | `specifiedBy` |
| `overall_evidence.strength_hint` | `strengthOfEvidenceProvided`(`not_clear`の場合は人間キュレーターに委ね、Evidence Lineを`draft`状態のまま保留する運用を推奨) |
| `overall_evidence.rationale` + `model_system_caveat` | `description` |
| `paper.pmid`/`paper.doi` | `reportedIn`(Documentとして) |

**`not_clear`は、Evidence Lineを作らずに「Evidence Gap」として3-3-2節のEvidence Gap Analysisに戻すのが一貫した設計**だと考える(このパイプライン自体が判定不能を返した場合、それは「証拠が無い」ではなく「証拠はあるが確信を持って判定できない」という異なる状態のため、区別して扱う)。

---

## 6. 未実装・今後の課題

- 上記は設計とCase3-var1での1件の実地検証のみ。Case2-var1(PS4)・Case3-var2/Case4-var1(PS3、確信度留保あり)は2026-09-14に追加検証済み(Case2-var1は「遺伝子レベルの統計研究であり変異specificではない」という別の注意点が判明。詳細は`real_cases_groundtruth_integrated_v8`参照)
- AcmGENTIC・CGBench共通の限界として、補足資料(supplementary)・図表は今回のスコープでも扱わない。本文PDFのみを対象とする
- `insufficient_data`と`not_clear`の区別、およびその後の人間への戻し方(UIが無い今回の設計でどう表現するか)は未確定
- **CGBenchの「74.6%不正解」を踏まえ、Day1での実演は「うまくいく例」だけでなく「うまくいかない例」も正直に見せる方が、藤原さんの論文の主張(過信を避け、限界を明示する)と一貫する**
- **【2026-09-14追加】PS4のように複数論文の集約で成立する基準への対応が未設計**。現在の3-1・4-1節の「候補論文の特定」は単一のPMIDを起点にする設計だが、ClinGen ERepoのように複数PMIDが列挙される場合、それら全てを個別に処理した上で「集合としてどう積み上がるか」を判断する仕組みが必要。今回はこれが無かったために、単一PMIDの検証結果だけで「基準不成立」と誤って結論づけてしまった
## 7. CGBenchデータを用いたパイロット検証(新規)

### 6-1. データソース

`VCI/clingen_vci_pubmed_fulltext_dedup_pmid.csv`(289行、ClinGen VCIキュレーション済み変異単位のACMG基準・met/not_met正解ラベル・論文全文が埋め込み済み)。うちPS3系70件・BS3系32件、計102件を対象とした。

**注意**:同データセットの`GCI/evidence_tables/experimental_evidence/evidence_cleaned_fulltext.csv`は、変異単位のPS3/BS3データではなく、**遺伝子-疾患関連性(Gene Validity)キュレーション用のデータ**であり、目的が異なる。混同しないこと。

### 6-2. 変異同定ゲートの改良

v3設計の変異同定ゲート(3-2節)を以下の点で改良し、実データに適用した:

- タンパク表記の空白入りバリエーション対応(`c.1012G > A`)
- 括弧付き表記対応(`p.(Gly338Arg)`)
- indel系変異の座標±2塩基の許容照合(ヒューリスティック)

結果、マッチ率は素朴な完全一致(初期パイロット10件で50%)から、**58/102件(56.9%)**に改善した(PS3系41/70、BS3系17/32)。

**発見**:MT-ND6の欠失変異で、論文本文の表記(`m.14512_14513del`)とGround Truthの表記(`m.14513_14514del`)が1塩基ずれているケースを、座標許容照合で発見した。境界表記の慣習差と見られ、厳密一致では原理的に検出できない。

残り44件は`unsuccessful`(変異同定ゲート失敗)。主因は、飽和変異アッセイ論文(例:MSH2/PTEN/CDH1の大規模機能アッセイ)で結果が本文でなく補足表・図にのみ存在するケースと見られる。v3の「補足資料・図表は扱わない」という既知の限界が、実データでも高頻度で顕在化した。

### 6-3. 方向性判定の系統的バイアス(最重要発見)

ゲートを通過した58件(重複を除き57件)について、v3の抽出プロンプト設計に沿って本文の記述から方向性判定を行い、Ground Truthのmet_statusと照合した。

| 結果 | 件数 | 割合 |
|---|---|---|
| Ground Truthと一致 | 17 | 30% |
| **不一致** | 20 | 35% |
| 判定保留(not_clear相当、安全側) | 20 | 35% |

「本文に明確な効果が書かれている」37件に限ると、**一致率は45.9%**——コイントス以下だった。

**さらに重要な点**:不一致20件のうち**19件が同じ方向のエラー**だった。論文の生データは明確に「有害」を示しているのに、公式にはnot_metだったケースを`met`と誤判定する、**過剰陽性の系統的バイアス**である(逆方向の誤り、すなわちGTがmetなのに誤ってnot_metと判定したケースは1件のみ)。

アッセイ種別で見ると、このバイアスには明確な勾配があった:

| カテゴリ | 不一致率 |
|---|---|
| ミトコンドリアsingle-fiber分離解析/cybrid研究 | **6/6 = 100%** |
| RASopathy系遺伝子の過剰発現アッセイ(PIK3CA/PIK3R2/RIT1等) | 4/7 = 57% |
| その他 | 10/24 = 42% |

**解釈**:論文自体は科学的に妥当なことを主張している。しかしClinGen VCEPが基準を正式に適用するには、論文の結論とは別に、そのVCEPが事前に承認したアッセイ手法か、十分な再現数があるか、といった**手続き的・制度的な条件**を満たす必要がある。文献のみを読むパイプラインには、この情報が原理的に見えない。

### 6-4. 改善案:VCEP仕様書(SOP)の活用(検証済み・部分的に有効)

CGBenchデータセットには、ClinGen各VCEPが公開している**遺伝子/疾患別のACMG基準仕様書(SOP)が実データとして含まれている**(`VCI/parsing_csr_criteria/version_csv_individual/`、約150ファイル)。これを実際に参照し、改善の余地を検証した。

**RIT1(RASopathy VCEP)の例**:
```
PS3, Moderate: Two or more different approved assays.
PS3, Supporting: One approved assay.
```
→ 「承認済みアッセイが何本必要か」という再現数の要件が明文化されている。パイプラインにこの情報を持たせることは技術的に可能。

**ミトコンドリアVCEPの例**:
```
PS3, Supporting: Functional validation is present in cybrid studies or single fiber analysis
```
→ **single fiber解析は明示的に承認されたアッセイ手法だった**。7-3節で立てた「アッセイ手法自体が非承認だから」という仮説は、少なくともミトコンドリアに関しては誤りだったことが判明した。

**結論(重要)**:
1. SOPをパイプラインに読み込ませ、「このコードには何本のアッセイが必要か」という構造的ルールを反映することは、技術的に可能かつ有効な改善である
2. **しかしそれだけでは今回の不一致は説明しきれない**。ミトコンドリアの仕様書は該当アッセイ手法を明示的に許可しているにもかかわらず、実データでは6件全滅だった。仕様書に明文化されていない、より細かい審査基準(統計的有意性の厳密な閾値、他基準との重複排除、パネル内部の追加審査等)が存在すると考えられる
3. RIT1についても、「承認済みアッセイのリスト」自体はこのCSVには含まれておらず、VCEPの原論文を別途参照する必要がある

### 6-5. その他の発見

- **GUCY2Dの2件(R588W・P575L)で矛盾を発見**:コメント欄が「閾値を満たした」と明記しているにもかかわらず`met_status`が`not_met`になっている。CGBench側の抽出ミスの可能性があり、要個別検証。

---

## 8. ミトコンドリア6件全滅の深掘り結果(重要な訂正)

7-3節で「6件全滅」としたミトコンドリアのPS3-Supporting不一致について、`summary`列(ClinGenキュレーターの評価要約そのもの)を精査した。

**発見**:6件全てで、`summary`列が明示的に基準名(例:「PS3_supporting」)を記載し、変異全体の`assertion`もPathogenic/Likely Pathogenicであるにもかかわらず、`met_status`列だけが`not_met`になっていた。特にLDLR c.1444G>C(idx 8146)は`summary`列に**英語でそのまま「So, PS3_Moderate was met.」と明記**されており、`met_status`列との矛盾は疑いようがない。

この矛盾パターンを、PS3/BS3系102件全体に対して機械的に検索したところ、**summary列を持つnot_met行62件中33件(53%)で同じ矛盾を確認した**。逆方向(met行なのにsummaryがnot metを示唆する)はsummaryを持つmet行27件中1件のみで、完全に一方向のバグだった。

### 8-1. 結論の訂正

7-3節で「文献を信じた結果、公式判定と系統的に食い違う」と結論づけたが、**この解釈は部分的にしか正しくなかった**。不一致20件のうち少なくとも13件(65%)は、パイプラインの誤判定ではなく、**CGBenchデータセットの`met_status`列自体に含まれるラベル付けバグ**で説明できる。真に「論文は明快だが公式判定は違う制度的理由による」と言えるのは、summary列が空欄で矛盾を確認できない少数のケース(例:PIK3R2/PIK3CAの3件、いずれもsummary欄が空)に絞られる。

**実務上の教訓**:CGBenchのようなキュレーション済みベンチマークを「正解ラベル」として無条件に信頼してはならない。少なくとも今回のCSVでは、met_status列を使う前に、同じ行のsummary/comments列と矛盾しないかを機械的にクロスチェックする前処理が必須。BH26のデモで「ベンチマークとの比較」を出す場合、この前処理をしないと、パイプラインの側ではなくベンチマーク側の欠陥を「AIの誤り」として誤って報告するリスクがある。

## 9. ClinGen ERepo公式APIによる全数検証(v5、確定結果)

7-3節で機械的に検出した「summary列が基準名を肯定的に記載しているのに`met_status`がnot_met」の疑わしい行、計36件(本体33件+参考3件)全てを、ClinGen ERepo公式API(`erepo.clinicalgenome.org/evrepo/api/classifications?gene=<遺伝子>&hgvs=<変異>`)で手動照合した。この環境からはネットワーク許可リストの制約でAPIに直接アクセスできなかったため、ユーザーが別環境から取得したJSONレスポンスを都度貼り付ける形で、36件全件の検証を完了した。

### 9-1. 最終結果

| 結果 | 件数 | 割合 |
|---|---|---|
| ✅ CGBenchのmet_status(not_met)は正しかった | 10件 | 28% |
| ⚠️ **本物のmet_statusバグ(公式にはMet)** | **26件** | **72%** |

**この矛盾候補カテゴリでは、`met_status`列は単独では信頼性がほぼ半分以下だった。**

### 9-2. 内訳詳細

**✅ 正しかった10件**(summary列に除外理由の説明があったもの、またはPS4/PS2/PM1等の強い基準だけで分類に到達していたもの):
USH2A c.12295-3T>A(idx 985、PVS1優先)/ MYOC c.731G>T・c.1412A>G・c.1037G>C(idx 3162, 3184, 4479、閾値未達)/ RUNX1 c.496C>T(idx 3859、PVS1優先)/ HNF1A c.347C>T(idx 7053、MDEP基準未達)/ GUCY2D c.3271C>T(idx 9547、PVS1経路でBS3不要)/ PIK3R2 c.1117G>A・PIK3CA c.2740G>A・c.2176G>A(idx 2747, 2759, 2761、PS4/PS2/PM1で分類済み)

**⚠️ 本物の矛盾26件**(遺伝子・VCEP問わず広範囲):
MT-ND1 m.3890G>A(3786)/ DICER1 c.2642T>C(3964)/ LDLR c.-136C>G・c.1003G>T・c.1444G>C・c.1739C>T(4046, 5072, 8146, 8149)/ MT-TL1 m.3291T>C(4263)/ MT-TN m.5690A>G(4870)/ MT-TS1 m.7497G>A(5109)/ HNF1A c.66C>G(5545)/ MT-ND6 m.14513_14514del・m.14487T>C(6600, 4852)/ GCK c.447C>A(6612)/ RAG1 c.2690G>A(7035)/ VHL c.273C>A(7223)/ MT-TK m.8340G>A(7456)/ MT-ND6 m.14597A>G(7655)/ RIT1 c.268A>G(8336)/ DYSF c.953T>A・c.1906G>A・c.1717C>T(9142, 9144, 9155)/ SGCB c.452C>G(9200)/ GUCY2D c.1762C>T・c.1724C>T(9523, 9551)/ TP53 c.329G>A(1587)/ MYOC c.898G>A(3208)

矛盾は最終的に**13遺伝子・13以上のVCEP**にまたがっており、特定の遺伝子・VCEPに限定されない、**CGBenchデータセットの`met_status`列全体に及ぶ系統的なラベル付けバグ**であることが確定的に裏付けられた。

### 9-3. 興味深い個別知見

- **同一論文・同一遺伝子内でも結果が割れるケースがあった**:MYOC変異3件(3162, 3184は同一論文でPS3/BS3ともに正しくnot_met、898G>Aのみ本物の矛盾)。met_status列の破損はランダムに近く、遺伝子や論文単位で予測できるものではない
- **HNF1Aも同様に割れた**:c.66C>G(5545)は矛盾、c.347C>T(7053)は正しい、という対比により、「summary欄の記述内容」こそが有効な予測因子であり、遺伝子は予測因子にならないことが再確認された
- **PI3K経路3変異(PIK3R2/PIK3CA×2)は3/3正しかった**:summary欄が空欄のケースは、既に強い基準(PS4/PS2/PM1)で分類が完結しており、PS3/BS3が公式に評価されていない、という制度的に妥当なケースだった
- **GUCY2Dのペア(9523/9551)は当初の最有力候補通り矛盾、もう1件(9547、ナンセンス変異のレスキュー実験)は正しかった**——ナンセンス変異はPVS1経路で分類されるため、BS3は別途要求されないという制度的理由

### 9-4. この環境のネットワーク制約について

`erepo.clinicalgenome.org`・`erepo.genome.network`はいずれもこのセッションのbash_toolネットワーク許可リストに含まれておらず、`curl`は一貫して403(`host_not_allowed`)を返した。web_fetchツールも「検索結果に出てきたURLしか開けない」制約があり、任意のAPIクエリを直接構築して叩くことはできなかった。今回はユーザーが別環境(制限のないマシン)から`gene=<遺伝子>&hgvs=<変異>`形式のクエリを36回実行し、その都度JSONレスポンスを本セッションに貼り付けてもらうことで全数検証を完了した。将来的にこのドメインがネットワーク許可リストに追加されれば、同様の検証をAIエージェント単体で自動化できる。

## 10. ラベル修正後の再検証(v6、最終決着)

9節で確定した26件の修正を`met_status`列に実際に適用し、7-3節の「テキストの方向性判定 vs Ground Truth」分析をやり直した。

### 11-1. 劇的な改善

| | 一致 | 不一致 | 判定保留 | 明確な信号があった件数中の一致率 |
|---|---|---|---|---|
| **ラベル修正前**(v4時点) | 17件 | 20件 | 20件 | **45.9%**(コイントス以下) |
| **ラベル修正後**(v6) | 29件 | 8件 | 20件 | **78.4%** |

**7-3節で「AIパイプラインは有害方向に系統的に偏る」と結論づけた現象の大部分(20件中12件)は、パイプラインの欠陥ではなく、CGBenchデータセット自体のラベル付けバグが原因だった。** ラベルを修正しただけで、一致率はコイントス以下から実用レベル(約8割)まで跳ね上がった。

### 11-2. 修正後もなお残る不一致(8件、こちらが「本物の限界」)

| idx | 遺伝子 | 内容 |
|---|---|---|
| 2747, 2759, 2761 | PIK3R2/PIK3CA | ✅ClinGen ERepoで直接確認済み・真に正しいnot_met。PS4/PS2/PM1で既に分類完結しており、機能データが強くても制度上PS3は評価されない(9-3節参照) |
| 1023 | KCNQ4 | 「null control」として使われた既知病的変異。それ自体の再評価対象ではない |
| 1843 | ITGB3 | 未検証(APIで裏取りしていない)。VCEP承認済みアッセイでない可能性 |
| 4950 | PALB2 | 方向は合っているが、強度コード不一致(BS3 vs BS3-Supporting相当) |
| 8403 | MSH2 | 本文中の別表(スプライシング解析)を誤って参照。実際の機能データは図表のみに存在 |
| **1587** | **TP53** | ⚠️**修正後も残った本物の逆転ケース**。本文は「Functionally abnormal」と明記(有害を示唆)だが、公式にはBS3_Supporting(良性)がMet。同一論文内に複数アッセイがあり、ClinGenが採用したのは別のアッセイ結果だったと推測される |

**つまり、ラベル修正後に残った不一致は、ほぼ全て「論文を読むだけでは分からない制度的・技術的理由」による、正真正銘のパイプラインの限界だった。** これは7-3節で当初想定していた仮説(「文献のみを読むAIは、制度的な事情を知らないので判定を誤る」)そのものであり、データバグというノイズを取り除いたことで、初めて本来検証したかった仮説を正しく評価できたことになる。

### 11-3. 結論

- **修正済みデータ(`clingen_vci_pubmed_fulltext_dedup_pmid_CORRECTED.csv`)を今後の評価・デモの基準として使うべき**
- パイプラインの真の弱点は、当初考えていたより小さい:方向性判定ロジック自体は、明確な信号がある場合78%の精度がある。残る限界は、(1)複数アッセイが混在する論文での参照先の取り違え、(2)制度上の適格性(承認済みアッセイか、既に他基準で分類完結しているか)の2種類に整理できる
- BH26のデモでは、「AIは有害方向に偏る」という当初の(誤った)結論ではなく、「**ベンチマークデータの品質検証が、AIパイプライン自体の検証と同じくらい重要である**」というメタな教訓を主軸に据える方が、より正確かつ価値のあるメッセージになる

## 11. RASopathy VCEP承認済みアッセイリスト(v7、新規入手)

11節で「未取得」としていたRASopathy VCEPの承認済み機能アッセイリスト(Wilcox et al. 2025, Genet Med Open, Supplementary Table 2)を、ユーザーが論文の補足PDF(オープンアクセス、CC BY 4.0)をアップロードしたことで入手できた。この環境のネットワーク制約(PMC/ScienceDirect/GIM Openいずれもボット検知でブロック)のため、直接のダウンロードは失敗していたが、ユーザー自身がブラウザから取得したPDFを渡すことで解決した。

### 13-1. 承認済みアッセイ一覧

| アッセイ | 測定内容 | PS3判定基準 | 適用可能な遺伝子 |
|---|---|---|---|
| RAS Activation Assay | RAF1/RBDと免疫沈降したRAS結合量 | RAS/RAF1またはRAS/RBD複合体の増加 | MRAS, HRAS, KRAS, SOS1, SOS2, NRAS, LZTR1 |
| MEK Activation Assay | リン酸化MEK/非リン酸化MEK比(基礎値・RTK刺激後) | リン酸化の増加 | MRAS, BRAF, HRAS, KRAS, MAP2K1, MAP2K2, PTPN11, RAF1, SOS1, SOS2, NRAS, LZTR1 |
| ERK Activation Assay | リン酸化ERK/非リン酸化ERK比(基礎値・刺激後) | リン酸化の増加 | MRAS, BRAF, HRAS, KRAS, MAP2K1, MAP2K2, PTPN11, RAF1, SOS1, SOS2, NRAS, LZTR1 |
| SHP-2 Phosphatase Activity | リン酸化/脱リン酸化SHP2比 | 脱リン酸化の増加 | PTPN11 |
| BRAF Kinase Activity | MEK/ERKをリン酸化するキナーゼ活性 | リン酸化の増加 | BRAF |
| RAF1 Kinase Activity | 同上 | リン酸化の増加 | RAF1 |
| LZTR1 Stability/Localization | LZTR1タンパク量・局在 | タンパク量減少または局在異常 | LZTR1 |

本文には「BS3はRASopathyでは適用不可」と明記されている(生理条件・刺激・細胞株依存で野生型的な結果を示すP変異が既知のため)。これは7-3節でRASopathy系にBS3が一件も観測されなかったことと整合する。

### 13-2. 新たな食い違い:RIT1がリストに含まれていない

この表の「適用可能な遺伝子」列にRIT1が一つも含まれていない。しかし9節でClinGen ERepo公式APIにより確定させたidx 8336(RIT1 c.268A>G、ERK1/2リン酸化アッセイ、HEK293T細胞)は、**PS3_Supporting: Met**だった。使用された手法(ERKリン酸化測定)はこの表の「ERK Activation Assay」の説明と完全に一致するにもかかわらず、遺伝子リストにRIT1が明記されていない。

考えられる説明(未検証、いずれも推測の域を出ない):
- この論文の88件のpilot遺伝子(RIT1含む15件)の評価が、本表の版に反映される前に固定された可能性
- 単純な表の更新漏れ
- RIT1固有の追加要件が本文中の別箇所(未確認)に記載されている可能性

**教訓**:VCEPの公式仕様書(SOP)であっても、実際の個別変異のキュレーション結果(ERepo)と完全に同期しているとは限らない。仕様書とERepoの間で食い違いがあった場合、最終的な正解は常にERepo(実際にキュレーターが適用した結果)であり、仕様書は「参考」に留める必要がある——これはCGBenchのmet_status列で学んだ教訓(9〜10節)と同じ構造の問題である。

## 13. 実装:変異照合ゲート・判定ロジック・LLM統合パイプライン(v9、新規)

セクション2〜3の設計を、実際に動くPythonコードに落とし込んだ。実装ファイルは以下の3本(+テスト一式):

| ファイル | 役割 |
|---|---|
| `ps3_bs3_ps4_gate.py` | 変異照合ゲート本体(PS3/BS3/PS4共通)。ERepo参照、production/validation二重モード、承認済みアッセイ照合(PS3/BS3専用) |
| `ps3_bs3_judgment.py` | 判定ロジック本体。プロンプトテンプレート、構造化出力スキーマ、人間キュレーター向け補助情報の生成、および後述の自動セーフティネット4種 |
| `ps3_bs3_llm_pipeline.py` | 上記2つを実際のLLM(vLLM上のgemma-4)・PubMed MCPに接続する統合スクリプト |

### 13-1. 変異照合ゲートの実装(ps3_bs3_ps4_gate.py)

- 表記ゆれ吸収(タンパク1/3文字表記、括弧有無、c./m.表記の空白ゆれ)+座標許容照合(±2塩基)
- `ERepoClient`(本番用スタブ)/ `OfflineERepoClient`(JSON形式でのテスト用)/ `CsvBackedERepoClient`(CGBench形式CSVを直接ソースにする、バルク運用にも使える実装)の3系統
- `run_gate()`でproduction/validationの二重モードを実装。production は既知変異でパイプラインを止め(`should_run_pipeline=False`)、validation は正解を`_hidden_ground_truth`に隠したまま強制的に後段を走らせる
- 承認済みアッセイ照合(`check_approved_assay`)にRASopathy VCEP表(11節)を搭載。**実装中に見つけた実バグ**:アッセイ名との一致を最優先せず、説明文のキーワード部分一致だけに頼ると、"ERK Activation Assay"が"MEK Activation Assay"と(共通する一般語"stimulation"等を介して)誤って一致することが判明。アッセイ名一致を1st passにし、キーワード一致は2nd passのフォールバックにする2段階方式に修正
- CGBench修正済みCSV(289行、PS3/BS3系102行)全件に対する回帰テストで、6〜9節の結果(matched 57/heuristic 1/unsuccessful 44、および26件の修正・10件の正解)を完全に再現することを確認

### 13-2. 判定ロジックの実装(ps3_bs3_judgment.py)

設計書2-4節のプロンプトを構造化出力スキーマ(`PS3BS3Judgment`)とともに実装。当初は日本語プロンプトだったが、**実運用のLLM対話は英語で行うという方針に合わせ、プロンプト本体・コード内コメントを全て英語に切り替えた**(このドキュメントおよびユーザーとの会話は引き続き日本語)。

後処理として、以下を実装:
- 承認済みアッセイ照合結果を、判定結果とは独立の「キュレーター向け補助情報」(`CuratorHint`)として出力
- 「単一研究のみに基づく判定か」「論文内に複数アッセイが混在していないか」の自己申告をLLMに求め、該当する場合は注意喚起

### 13-3. gemma-4での実地テスト:4つの自動セーフティネットが生まれた経緯

このパイプラインを、実際にvLLM上のgemma-4(ユーザー環境)に接続し、MYH7・PTEN・TP53の3変異(いずれもClaudeが事前に手動でブラインド判定していたものと同一)で繰り返しテストした。その過程で、**gemma-4が「有害方向へ過剰に断定する」という同一のバイアスを、毎回異なる形で表出させる**ことが分かり、都度、実際の失敗データに基づいて自動検出ロジックを追加していった。

| # | 発見された失敗パターン | 実例 | 追加した検出 |
|---|---|---|---|
| 1 | enum値に説明を括弧で付け足す(`"PS3(damaging)"`) | MYH7・PTEN | プロンプトの曖昧さ修正+パース側の防御的正規化(`_clean_enum_token`) |
| 2 | 実験間で結果が矛盾しているのに、無視して断定 | PTEN(VAMP-seq=正常 vs Akt活性化=異常、両方あるのにPS3と断定) | `detect_experiment_conflict`:カテゴリフィールド同士の矛盾を検出し`not_clear`に強制上書き |
| 3 | 根拠文は正常寄りと書いているのに、カテゴリラベルは結論に合わせて改ざん | PTEN(rationaleは"WT-like"と明言、しかし全実験が`functionally_abnormal`) | `detect_narrative_categorical_mismatch`:自由記述とカテゴリの食い違いを検出 |
| 4 | 具体的なデータが本文に無いのに、断定する(実験0件) | TP53(1回目) | `detect_empty_experiments_with_definitive_direction`:実験0件での断定を矛盾として検出 |
| 5 | 構造化フィールドが文字列に潰れる(`"variant_matching": "matched"`) | MYH7・PTEN(3回目実行) | `PS3BS3Judgment.from_json`に防御的フォールバック追加+プロンプトにリテラルJSON例を追加 |
| 6 | 具体的な測定値が無いまま、論文全体の傾向から一般化して断定 | TP53(2〜4回目、実験はあるが数値が無い) | `detect_no_quantitative_evidence`:測定値らしい数値パターン(小数・%・p値・単位付き数値)の不在を検出 |

**パターン6を追う過程で、検出ロジック自体に2つのバグを作り込み、都度実データで発見・修正した**:
- 単純な正規表現一致では、"lack WT-like activity"(否定文、実際は有害を示す)を正常寄りのシグナルと誤読する → 否定語チェック(`_has_uncontradicted_normal_finding`)を追加
- 「数字が1つでもあるか」という基準では、"p53"・"Nutlin-3"・"A549"のような遺伝子/薬剤/細胞株名の数字を拾ってしまい、ほぼ全てのテキストが「数値あり」になる → 小数・%・p値・単位付き数値に限定した正規表現に絞り込み

### 13-4. 最終的な実行結果(3変異、複数回)

| 変異 | LLMの生判定 | 発動したセーフティネット | 採用判定 | Ground Truth | 結果 |
|---|---|---|---|---|---|
| MYH7 c.1594T>C | PS3 | (発動なし、素直に正しい) | PS3 | PS3=met | ✅的中 |
| PTEN c.112C>T | PS3 | 実験間矛盾検出 | not_clear | PS3-Moderate=not_met | ✅誤りを回避 |
| TP53 c.875A>G | PS3 | 数値データ不在検出 | not_clear | BS3=met | ✅誤りを回避 |

3変異とも、**最終的に「誤った断定」はゼロ**になった。ただし三度とも安全側(not_clear)に倒れており、「LLM単体で正しく断定できる」ケースはMYH7の1件のみ。セーフティネットに頼らずに正しく判定できる比率(=プロンプト自体の一次的な精度)は、この3件のサンプルでは高いとは言えない。

### 13-5. テスト規模

`ps3_bs3_ps4_gate.py`・`ps3_bs3_judgment.py`合わせて、実データに基づく回帰テストが200件超(ゲート単体16〜17件、ゲート全件回帰151件、判定ロジック26件)。うち大半が、今回の一連の実行で実際に発生した失敗を再現する形の実例ベーステストになっている。

## 14. 現時点での正直な評価と今後の論点(v9時点)

- **後処理によるセーフティネットは「誤った断定」の防止には有効だが、「正しく断定できる」ケースを増やす方向には寄与しない**。3件中2件がnot_clearに倒れた現状は、安全だが実用上の歩留まりとしては物足りない。プロンプト自体の精度向上(特にTP53のような大規模飽和変異スクリーニング論文への対処)が次の課題として残る
- **モグラ叩き的な対策の追加は、13-3節で見た通り、対策自体が新しいバグを生むリスクと表裏一体**。パターン6の検出ロジックだけで2回作り直しが必要だった。今後、新しいセーフティネットを追加する際は、必ず実データでの回帰テストとセットにすることを徹底する
- **サンプル数が3件(MYH7・PTEN・TP53)のみ**。統計的な結論を出すには全く足りない。次に取り組むべきは、CGBenchの未使用プール(3757行、8-9節参照)から新規に10〜20件程度を選び、同じ手順でgemma-4に投げて、セーフティネット発動率・的中率を定量化すること
- **TP53のような「本文に個別データが無い大規模スクリーニング論文」への対処は、後処理の検出だけでなく、変異同定ゲート(match_status)自体の判定基準をプロンプト側でより厳格にする方向の改善も検討の余地がある**(例:「個別の変異名・数値が本文中で直接言及されていない場合はunsuccessfulとする」という指示を、現在の「含まれていればmatched可」という緩い基準に対して明確に追加する)

- **met_status列のバグの発生源は未特定**(CGBenchのデータ収集・整形スクリプト側の問題か、ClinGen ERepo APIの応答形式変化への追従漏れか等)。ClinGen側から見て、なぜsummary列(記述テキスト)とmet_status列(構造化フラグ)の間でこれほど広範な不整合が生じているのか、CGBench開発チームへの問い合わせが必要
- **RASopathy承認済みアッセイ表とERepo実測結果の食い違い(11-2節、RIT1)は未解明**。表の版差か更新漏れかは未検証
- 変異同定ゲートの残り44件(補足表・図のみに存在するケース)への対応は引き続き未設計
- **今回確定した「met_status列は単独では信頼性が低い」という知見を踏まえ、CGBenchデータを今後BH26のデモや評価に使う場合は、summary/comments列との整合性チェックを前処理として必須にすべき**(7-3節で提案した通り)
- idx 1843(ITGB3)は今回API未検証のまま残っている

---

## 15. v10更新まとめ(2026-09-15、別クライテリア実装者向け)

v9以降、ユーザーとの対話の中で実装した内容をまとめる。**他の基準(criterion)を実装する同僚向けに、特に15-4節(共通インフラの使い方)と15-5節(スコープの線引き)を重点的に読んでほしい。**

### 15-1. ERepoClient:実HTTP実装が完了(9-4節の制約は解消)

9-4節では「`erepo.clinicalgenome.org`はこのセッションのネットワーク許可リストに含まれておらず403だった」と記録していたが、**本セッションの環境では普通に到達可能だった**(`requests.get`で200 OK)。`ps3_bs3_ps4_gate.py`の`ERepoClient.lookup()`を実装済み:

```python
def lookup(self, gene: str, hgvs: str) -> ERepoLookupResult:
    resp = requests.get(
        f"{self.BASE_URL}/evrepo/api/classifications",
        params={"gene": gene, "hgvs": hgvs},
        timeout=self.TIMEOUT_SEC,
    )
    resp.raise_for_status()
    data = resp.json()
    return OfflineERepoClient([data]).lookup(gene, hgvs)  # パースはOfflineERepoClientに委譲、重複実装しない
```

実データ(`RUNX1 c.601C>T`)で`evidenceLinks`4件・`PS4:Met`を確認済み(design doc本文が既に触れていた実例と一致)。`ERepoLookupResult`に`evidence_pmids`フィールド(`evidenceLinks`から抽出したPMIDリスト)を追加した。

**環境依存の注意**:ネットワーク到達性はセッション・実行環境によって変わりうる。新しい環境で動かない場合は、まずこのURLへの`curl`到達性を再確認すること。

### 15-2. VA-Spec出力の確定(5節を置き換え、`va_spec_export.py`)

**打ち合わせにより、出力形式はGA4GH VA-Spec `EvidenceLine`(JSON)に確定した。** 5節の素朴なマッピング表は、複数論文集約(13節で実装した`AggregatedJudgment`)より前の設計であり、実装と対応していなかったため置き換える。

**v2追記(2026-09-15、同日中に再修正)**:当初`va_spec_export.py`は公式ドキュメントとサンプルyamlを読んで**dictを手組み**していた(公式スキーマライブラリを未導入・未使用)。ユーザーからの確認をきっかけに`pip install ga4gh.va_spec`(実際のPydanticスキーマ実装、`ga4gh.va_spec.acmg_2015`/`ga4gh.va_spec.base.core`)を導入し、**実際のPydanticモデル(`EvidenceLine`/`Document`/`Method`/`Contribution`/`Agent`/`MappableConcept`/`Coding`/`Extension`)を構築して`.model_dump(mode="json", exclude_none=True)`で出力する実装に書き換えた**。この切り替えで、手組み実装の**2箇所の実際の誤り**が見つかった:
- `Contribution.activityType`は**素のstr**であり、`MappableConcept`(`{"name": ...}`)ではなかった
- `Contribution.contributor`は任意のdictではなく**`Agent`モデル**だった

それ以外(Direction enumの値が`supports`/`neutral`/`disputes`であること、`System.ACMG`の値が`"ACMG Guidelines, 2015"`と完全一致すること、Document/Method/MappableConcept/Coding/Extensionの各フィールド名)は手組み実装と一致しており、公式ドキュメントの読み込みが概ね正確だったことも確認できた。なお`StrengthOfEvidenceProvided`の正式な語彙(`standalone`/`very strong`/`strong`/`moderate`/`supporting`)に合わせ、strength tierの表記を`Supporting`/`Moderate`/`Strong`(独自の大文字表記)から`supporting`/`moderate`/`strong`(公式語彙、小文字)に変更した。

**実行環境への影響**:`ga4gh.va_spec`(および依存の`ga4gh.vrs`・`ga4gh.cat_vrs`・`bioutils`等)のpipインストールが必要(`requirements.txt`等の依存管理ファイルはこのプロジェクトに存在しないため、`requests`同様その都度手動インストールが必要)。

**粒度**:1つの(変異, criterion)ペアにつき1つの`EvidenceLine`(集約後)を作り、実際に読んだ論文それぞれを`hasEvidenceItems`内の入れ子`EvidenceLine`として持たせる(サンプルyaml`va_spec_example.yaml`のEvidenceLine002→003の入れ子パターンと同型)。

| 実装側のデータ | VA-Spec側のフィールド |
|---|---|
| `AggregatedJudgment.aggregated_direction` | `directionOfEvidenceProvided`(`supports`/`disputes`/`neutral`) |
| `PaperContribution`(論文ごと)×N | `hasEvidenceItems`内の入れ子`EvidenceLine` |
| `FinalResult.judgment.overall_evidence.rationale` | `description`(自由記述、`InformationEntity`継承フィールド) |
| `FinalResult.curator_hints`(severity付き複数件) | `extensions`(構造化データの置き場、severityを保持できる) |
| PMID | `reportedIn`(`Document`型、`pmid`/`urls`/`doi`/`title`フィールドを持つ) |
| 独立して同方向を示した論文数 | `strengthOfEvidenceProvided`/`evidenceOutcome`(下記参照、ヒューリスティック) |

**`not_clear`の扱い(重要な設計変更)**:5節は「`not_clear`はEvidence Lineを作らずEvidence Gapに戻す」としていたが、**VA-Spec仕様を実際に確認したところ`directionOfEvidenceProvided`には`neutral`という値が公式に用意されており**(`supports, disputes, neutral`の3値、[公式ドキュメント](https://va-spec.ga4gh.org/en/latest/core-information-model/entities/information-entities/evidence-line.html)で確認)、これが`not_clear`の受け皿として適切だと判断した。よって**`not_clear`でもEvidence Line自体は作り**、`directionOfEvidenceProvided: neutral`、`reportedIn`(調べた論文は残す)、`description`(なぜ確信が持てなかったか)を設定し、strength系フィールドだけ省略する。

**強度(strength)の扱い(ヒューリスティック、要注意)**:このパイプラインは設計上「方向のみを判定し、強度は対象外」(`overall_evidence.strength_hint`は常に固定文字列で、プロンプト側でLLMに実質判断させていない)。しかしユーザーから「精度が低くてもいいから入れたい」との要望があり、**唯一実際に計算しているシグナル(独立して同方向を示した論文数)から代理指標を作った**:

| 独立して合意した論文数 | strength |
|---|---|
| 1本 | Supporting |
| 2本 | Moderate |
| 3本以上 | Strong |

ACMGコードは**「テスト対象のcriterion」ではなく「実際に確立された方向」でラベル付けする**(`direction.value`を使う)。例:PS3を評価していてBS3方向の証拠が出た場合は`PS3_Moderate`ではなく`BS3_Moderate`と正しく表記する(実装時に見つけて直したバグ)。全てのstrength付きEvidence Lineには、これがVCEP非校正のヒューリスティックである旨を明記した`extensions`(`strengthEstimationMethod`)を必ず添付する。

`va_spec_export.py`はcriterion非依存の作りで、`OverallDirection`のような特定のEnum型を一切importしない。`direction.value == "not_clear"`という文字列比較で判定するため、PS4/PP1/BS4を含むどの基準の判定結果でもそのまま使える(15-4節参照)。

### 15-3. curator hintの品質バグ(2種類、実データ実行で発覚)

ユーザーから「同じメッセージが複数出て何を指すか分からない」との指摘を受け、`ps3_bs3_judgment.py`の`generate_curator_hints()`を精査した結果、**同じ構造のバグを2箇所で発見**:

- `NO_SPEC_AVAILABLE`(VCEPに承認済みアッセイリストが無い)・`BS3_NOT_APPLICABLE`(そのVCEPではBS3自体が使えない)の判定は、**VCEP+criterionだけで決まり実験内容に依存しない**にもかかわらず、実験ごとのループの中でチェックしていたため、実験数だけ同一メッセージが複製されていた(実例:PTEN c.278A>Gの1論文で同一文言が×10)。→ ループの外で一度だけ判定するよう修正。
- `model_system_caveat`(実験ごとのモデル系限界の注記)も、どの実験のcaveatか無記名で、かつ複数実験が同一文言を持つ場合に重複していた。→ アッセイ名を明記し、同一文言は重複排除。

**教訓(新基準実装者へ)**:「実験ごとのループの中でVCEP/criterion単位の判定を行っていないか」は、新しいcuratorヒントを追加するたびに確認すべき典型的なバグパターン。

また、「This paper」「This judgment」のような曖昧な指示語は、hintがVA-Specの`extensions`に埋め込まれて元の文脈から切り離されて読まれる可能性があるため、**全てのhint生成関数にオプションの`pmid`引数を追加し、本文中にPMIDを埋め込むよう修正**した(`generate_curator_hints(..., pmid=None)`、`finalize(..., pmid=None)`)。新基準のhintを書くときも、この慣習(pmidを引数で受け取り、"This paper (PMID:xxx)"のように埋め込む)を踏襲すること。

### 15-4. 共通インフラの汎用化:`evidence_common.py` / `JudgmentEngine`パターン

**別クライテリアを実装する同僚が最も参照すべき節。** PS4・PP1/BS4を追加するにあたり、PS3/BS3専用だった以下のクラス・関数を`evidence_common.py`に切り出し、criterion非依存にした:

- `MatchStatus` / `VariantMatchingResult`(論文が対象変異を実際に扱っているかの判定。全criterionで共通)
- `CuratorHint`(severity + message)
- `FinalResult`(`judgment: Any`, `curator_hints`, `effective_direction: Any`)
- `PaperContribution`(pmid + FinalResult)
- `AggregatedJudgment`(複数論文の集約結果)
- `aggregate_multi_paper_results(contributions, not_clear)`:複数論文の結果を1つの変異レベル判定に集約する汎用ロジック。**`not_clear`センチネルを呼び出し側が明示的に渡す**設計にすることで、`OverallDirection.NOT_CLEAR`・`CaseControlDirection.NOT_CLEAR`・`SegregationDirection.NOT_CLEAR`など、どのEnumでも同じ実装を使い回せる(全てのdirection Enumが`NOT_CLEAR = "not_clear"`という文字列を持つ、という共通の作法に依存している)

**新しいcriterionモジュールを追加する手順**(`ps4_judgment.py`・`segregation_judgment.py`が実例):

1. 自分のcriterion専用の構造化出力スキーマを定義する(`@dataclass`)。`variant_matching: VariantMatchingResult`は必ず持つ。`overall_evidence`は`direction`(自分のEnum)・`strength_hint`(固定"not_clear")・`rationale`を持つ`XxxEvidence`にする(この命名を揃えておくと、パイプライン側の`judgment.overall_evidence.direction.value`/`.rationale`アクセスがそのまま使い回せる)
2. 自分のdirection Enumを定義する:必ず`NOT_CLEAR = "not_clear"`を含める。ACMGの実際のコード名(例:`"PS4"`, `"PP1"`, `"BS4"`)をそのまま値にする
3. `build_prompt(gene, hgvsc, hgvsp, equivalents, full_text) -> str`を書く
4. `from_json(data: dict)`(自分のJudgmentクラスの`@staticmethod`)を書く
5. `finalize(judgment, pmid=None) -> FinalResult`を書き、自分のcriterionに特有の矛盾検出(セーフティネット)をここに入れる
6. `aggregate_multi_paper_results(contributions)`は、`evidence_common`の汎用実装への薄いラッパーにする:
   ```python
   def aggregate_multi_paper_results(contributions):
       return _generic_aggregate_multi_paper_results(contributions, not_clear=MyDirection.NOT_CLEAR)
   ```
7. `ps3_bs3_llm_pipeline.py`に`JudgmentEngine`インスタンスを1つ追加し、`ENGINE_BY_CRITERION`に自分のcriterion文字列を登録する:
   ```python
   MY_ENGINE = JudgmentEngine(
       name="MY_CRITERION",
       build_prompt=my_mod.build_prompt,
       from_json=my_mod.MyJudgment.from_json,
       finalize=lambda judgment, gene, vcep_name, criterion, pmid: my_mod.finalize(judgment, pmid=pmid),
       aggregate=my_mod.aggregate_multi_paper_results,
   )
   ENGINE_BY_CRITERION["MY_CRITERION"] = MY_ENGINE
   ```
   これで`judge_variant`/`judge_single_paper`/`va_spec_export.build_evidence_line`/`score_direction`は一切変更せずに新criterionが動く。
8. (任意だが推奨)自分のJudgmentクラスに`structured_evidence_items(self) -> list[dict]`を実装する。詳細は15-4-1節参照

`va_spec_export.py`側も同じ理由で、特定のEnumをimportせず`direction.value == "not_clear"`という文字列比較のみで動くようにしてある(15-2節末尾参照)。新criterionを追加してもVA-Spec出力側は無改修で対応できることを、PS4/PP1/BS4追加時に実際に確認済み。

#### 15-4-1. curator UI向け:`structured_evidence_items()`規約(v10追加)

実際のcurator向けUI要件資料(`doc/recs for expert board.docx`)を確認したところ、evidence cardのモックアップ(PM1の例)は「criterionバッジ + strengthバッジ + 定義文」だけでなく、**個別の判定項目ごとにチェックマークと具体的な数値を見せるチェックリスト**という構造を持っていた(例:「✓ Exonic hotspot: 15 pathogenic...が21bp領域内に」)。`description`の自由記述文だけではこの種のUIを描画できないため、各criterionモジュールのJudgmentクラスに、以下の規約に従う**任意のメソッド**を実装できるようにした:

```python
def structured_evidence_items(self) -> list[dict]:
    """各要素は {"label": str, "checked": bool, "detail": str} の形。"""
    ...
```

- `label`:チェックリストの見出し(例:`"VAMP-seq: functionally_normal"`、`"Family F1"`)
- `checked`:その項目単体で見て「陽性/条件を満たす」かどうかの真偽値。**論文全体の採用された方向とは独立**(例:ある家系の分離データ単体が矛盾を示しているかどうかは、そのproduct全体のPP1/BS4判定結果と別に判定できる)
- `detail`:具体的な数値・所見を含む説明文

`va_spec_export.py`はこのメソッドを`getattr(judgment, "structured_evidence_items", None)`で**ダックタイピングにより存在すれば呼び出す**だけで、特定のJudgmentクラスをimportしない(`_is_not_clear()`と同じ設計方針)。実装しなければ単に`structuredEvidenceItems`extensionが省略されるだけなので、後方互換も保たれる。出力は各論文の入れ子`EvidenceLine`の`extensions`に`{"name": "structuredEvidenceItems", "value": [...]}`として追加される。

現在の実装例:
- PS3/BS3(`ps3_bs3_judgment.py`):実験(`ExperimentExtraction`)ごとに1項目、`checked`は「functionally_abnormal/functionally_normalという確定的な結果か(intermediate/mixed/unclearでないか)」
- PS4(`ps4_judgment.py`):罹患者側カウント・対照側カウント・OR/p値・study_designの4種類の項目
- PP1/BS4(`segregation_judgment.py`):家系(`FamilySegregationData`)ごとに1項目、`checked`は「その家系単体で非分離シグナルが無いか」

新しいcriterionを実装する場合も、このメソッドを実装しておくとcuratorのUI側でチェックリスト表示にそのまま使える。

### 15-5. スコープの確定:「文献で判定できる基準」だけを対象にする

ユーザーの指摘で、`doc/BH26_participant_briefing_v3_en.md`(66行目)にLayer3(手動レビュー要)の分類が既にあることが判明した:

> PS2, PS3, PS4, PM3, PM6, PP1, BS2, BS3, BS4, BP2, BP5

このうち、**実際にPubMed論文を読んで判定できるのはPS3, BS3, PS4, PP1, BS4の5つ**で、残り(PS2, PM3, PM6, BS2, BP2, BP5)は基本的に**患者本人の臨床・遺伝子検査記録**(親子鑑定、phasing検査、本人の家系図、健常キャリア記録等)に依存するもので、論文を読んでも出てこない性質のものと判断した。よって**このLLMパイプラインのスコープはPS3/BS3/PS4/PP1/BS4の5基準に確定**し、他の6基準は対象外とする。

PP1とBS4は同一の家系分離データに基づく方向ペア(PS3/BS3と同じ構造)なので、`segregation_judgment.py`1モジュールで両方をカバーしている。PS4には対になる「否定側」のACMGコードが存在しない(BS2は関連はするが別の判定データソース)ため、`ps4_judgment.py`の`CaseControlDirection`は`PS4`/`not_clear`のみ。

### 15-6. PS4専用モジュール(`ps4_judgment.py`)

4節の草案プロンプト(症例対照データの抽出)を実装に落とし込んだ。構造化スキーマ:

- `study_design`: `case_control` / `family_cohort` / `case_series` / `not_case_control`
- `case_control_data`: `affected_carriers` / `unaffected_carriers` / `total_affected` / `total_unaffected` / `odds_ratio` / `p_value`(全てnull許容)
- `overall_evidence.direction`: `PS4` / `not_clear`

**実データで見つかった設計バグ(重要)**:当初`case_control_data`は生の人数(affected_carriers/total_affected)しか持たせていなかったが、実際にgemma-4を走らせたところ(GJB2 c.109G>A、PMID:31160754)、論文は**オッズ比とp値(OR=20, 95%CI 17-24, p<0.0001)のみで有意性を示しており、生の人数は本文中に無かった**。当初のセーフティネット(`detect_definitive_without_numbers`)は「生の人数が無ければnot_clearに強制」という設計だったため、この**科学的に妥当で強い証拠を誤ってnot_clearに落としていた**。症例対照研究がOR/p値のみで有意性を報告するのは通常の実務であり、`odds_ratio`/`p_value`フィールドを追加し、セーフティネットを「生の人数もOR/p値も両方無い場合のみnot_clearに強制」に修正した。

**教訓**:PS3/BS3のセーフティネット(13節)は「安易にPS3/BS3と断定しすぎる」方向の誤りを防ぐために作られたが、**セーフティネット自体が「安全側に倒しすぎて正しい証拠を捨てる」新しい誤りを生みうる**ことを、PS4で実際に確認した。新しいセーフティネットを追加するときは、必ず実データで両方向(過剰検出・過小検出)を確認すること。

**実行結果**(GJB2 c.109G>A、1論文、修正後):LLM判定=PS4、ground truth(ERepo実測)=PS4 met → 一致。

### 15-7. PP1/BS4専用モジュール(`segregation_judgment.py`)

新規設計(design docに事前の草案なし)。構造化スキーマ:

- `families`: 家系ごとのリスト。各要素は`family_id` / `affected_with_variant` / `affected_without_variant`(**BS4の鍵となるシグナル**) / `unaffected_with_variant`(**同じくBS4の鍵**) / `unaffected_without_variant` / `informative_meioses` / `notes`
- `overall_evidence.direction`: `PP1`(完全共分離)/ `BS4`(明確な非分離)/ `not_clear`

**セーフティネット(`detect_segregation_direction_conflict`)**:LLMが「PP1(完全共分離)」と結論したのに、抽出した家系データ自体が非分離シグナル(`affected_without_variant>0`または`unaffected_with_variant>0`)を含む場合、矛盾として`not_clear`に強制する(PS3/BS3の`detect_experiment_conflict`と同じ発想)。フィクスチャテストで動作確認済み。

**実行結果**(GJB2 c.109G>A、同一の1論文をPP1・BS4それぞれの視点で評価):
- PP1として評価 → LLM判定=PP1、ground truth(ERepo実測)=PP1_Strong met → 一致
- BS4として評価 → LLM判定=PP1(=BS4は不成立を意味)、ground truth=BS4 not_met → 一致

同じ論文・同じ家系データから、PP1とBS4という別々のcriterionを評価しても矛盾のない結果が得られることを確認した。

### 15-8. Ground truthの再監査(重要な教訓の再発)

9-10節で「CGBenchのmet_status列には系統的なラベル付けバグがある」ことを確認し36件を修正済みだったが、**その監査対象に入っていなかった変異でも同じ種類のバグが残っていた**ことを、今回`main()`の全テストケース(25件)をライブERepo APIと機械的に突き合わせる形で発見した:

| 変異 | 修正前のground_truth | ライブERepoの実際 |
|---|---|---|
| PTEN c.112C>T | PS3-Moderate = not_met | PS3_Moderate = **met** |
| PTEN c.278A>G | PS3-Supporting = not_met | PS3_Supporting = **met** |
| MTOR c.4447T>C | PS3-Supporting = not_met | PS3_Supporting = **met** |
| MYOC c.1187_1188insCCCAGA | PS3-Moderate = not_met | PS3_Moderate = **met** |

いずれも、ERepo連携より前の初期(design doc記載の「genuinely blind」検証)に手打ちされた変異で、後から追加した変異群(12件+PS4/PP1/BS4の6件)のように個別のライブAPI裏取りをしていなかったもの。**この4件のうち`PTEN c.278A>G`は、後述15-9節の19件実行で唯一の「不一致」として報告していたケースであり、実際にはground truthの誤りで、LLMの判定(PS3)は正しかった**ことが判明した。

**恒久的なルール(以降の全作業に適用)**:CGBenchのCSV・過去の記録にあるground truthを**そのまま信用せず、テストケースとして使う前に必ずライブERepo APIの`evidence_code_status`と突き合わせる**。今回は以下のワンライナーで全件を機械的に検証できた:

```python
r = ERepoClient().lookup(gene, hgvsc)
live_status = r.evidence_code_status.get(code)  # code: "PS3_Moderate"のようにアンダースコア表記
assert live_status.value == expected_status
```
(CSVは`PS3-Moderate`のようにハイフン表記、ERepoのAPIレスポンスは`PS3_Moderate`のようにアンダースコア表記なので変換が必要な点に注意)

### 15-9. 実行結果まとめ(v10時点)

- **PS3/BS3**:19件のground truth全件をライブERepoで確認済み(15-8節)。実際にgemma-4で実行した結果(修正前のground truthベース):match=3, mismatch=1(実際はground truth誤りでmatch), reserved(not_clear)=15。判定確定分の一致率は実質4/4=100%(母数が小さく統計的意味は限定的)
- **PS4/PP1/BS4**:GJB2 c.109G>Aの1論文で3criterion(PS4/PP1/BS4)とも一致を確認。ERepoで実データ検証済みの候補は他に5件用意済み(RUNX1 c.601C>T、USH2A c.5581G>A/c.8559-2A>G、RUNX1 c.442_449del)だが、フル実行はまだ行っていない
- 既存の回帰テスト(`test_ps3_bs3_judgment.py`・`test_ps3_bs3_ps4_gate.py`・`test_ps3_bs3_ps4_gate_full.py`、計203件)は本セッションの全変更を通して壊れていない

### 15-10. 関連研究メモ

**"Benchmarking generative AI tools for literature retrieval and summarization in genomic variant interpretation"**([bioRxiv 2025.09.29.679212](https://www.biorxiv.org/content/10.1101/2025.09.29.679212v1.full))

5つの生成AIツール(ChatGPT, MistralAI, VarChat, Perplexity, ScholarAI)を、40件の変異(germline 20・somatic 20)に対する文献ベース要約タスクで比較したベンチマーク論文。6名の専門家によるブラインド評価(要約精度・ハルシネーション耐性・可読性・文献網羅性・所要時間の5指標)。VarChatが総合最高評価、GPT-4oが僅差2位。文献引用特化のはずのScholarAI/Perplexityは相対的に低評価。参照文献が10件未満だとどのツールも性能が明確に低下する点、専門家によるバリデーション抜きでの研究利用は推奨できないとする結論が特徴。

**このプロジェクトとの関係**:ACMG/AMP基準(PS3/BS3/PS4/PP1/BS4)自体の判定精度(direction/strength)は評価対象外で、あくまで文献要約タスク一般のベンチマーク。ただし「参照文献が少ないと性能が落ちる」「ハルシネーションリスクがあり専門家の検証が必須」という結論は、本設計で既に採用している設計方針(証拠不十分な場合は`not_clear`/`neutral`に倒す、curator_hintsで人間の確認を促す、13節のセーフティネット群)と方向性が一致する外部エビデンスとして参照できる。

### 15-11. パッケージ化、および分類(classification)モジュールの新設(v10追記、2026-09-15)

**背景**:別チームとのインターフェース(API仕様)がまだ固まらない見込みのため、優先順位を「classificationモジュール(複数criteriaの証拠を統合して最終分類を出す)を先に作る」に変更した。理由は、これが自チーム5criteriaだけでも今すぐ作成・検証でき、他チームのcriteria追加時にも無改修で拡張できる基盤だから(15-4節のJudgmentEngineパターンと同じ思想)。あわせて、ファイル数が増えてきたためパッケージ構成に整理した。

#### 15-11-1. パッケージ構成

フラットに置かれていた7つのモジュールを`acmg_pipeline/`パッケージ配下に再配置した:

```
acmg_pipeline/
  __init__.py
  common.py          (旧 evidence_common.py)
  gate.py            (旧 ps3_bs3_ps4_gate.py)
  export.py          (旧 va_spec_export.py)
  pipeline.py        (旧 ps3_bs3_llm_pipeline.py)
  classification.py  (新規)
  criteria/
    __init__.py
    ps3_bs3.py       (旧 ps3_bs3_judgment.py)
    ps4.py           (旧 ps4_judgment.py)
    segregation.py   (旧 segregation_judgment.py、PP1/BS4)
    stubs.py         (新規)
    registry.py      (新規)
```

内部import文を全てパッケージパス(`from acmg_pipeline.gate import ...`等)に更新し、既存テスト3ファイル(計203件)も新パスに追随させた上で全件パス確認済み。`common.py`には`export.py`が個別に持っていた`_is_not_clear`/`_strength_tier`ヘルパーを`is_not_clear`/`strength_tier_from_paper_count`として公開関数化し、後述の`classification.py`と共有するようにした(同じヒューリスティックを2箇所に別々に持たせない)。

#### 15-11-2. `classification.py`:Tavtigian点数法による最終分類算出

複数criteriaの証拠(met/not_met/not_evaluated + strength)を1つの最終分類(Pathogenic〜Benignの5段階)に統合する。手法は**Tavtigian et al. 2018の点数法**(Richards et al. 2015のACMG/AMPカテゴリカル併合規則を数値化した、公表済み・標準的な拡張手法。本設計書や`doc/recs for expert board.docx`のTable 2/3画像で既に言及済み)。各metの基準がstrengthに応じて点数を持つ(Supporting=1・Moderate=2・Strong=4・Very Strong=8、benign側は同じ絶対値の負数)。合計点数を閾値(Pathogenic≥10・Likely Pathogenic 6-9・VUS 0-5・Likely Benign -1〜-6・Benign≤-7)に当てはめて分類する。BA1がmetの場合は、点数計算を経由せず直接Benignとする特例(Richards et al. 2015のカテゴリカル規則そのもの、Tavtigianの点数拡張とは別)も実装した。

主要な型:
- `CriterionEvidence(code, status, strength, source)`:1基準1変異分の証拠。`status`は`gate.py`の既存`CriterionStatus`(MET/NOT_MET/UNKNOWN)をそのまま再利用し、新しい列挙型を増やさなかった
- `ClassificationResult(category, score, met, not_met, not_evaluated_codes, ba1_override)`
- `classify(evidence: list[CriterionEvidence]) -> ClassificationResult`:本体ロジック
- `from_aggregated_judgment(aggregated, code) -> CriterionEvidence`:このプロジェクトのLLMパイプライン出力(`AggregatedJudgment`)を証拠に変換するヘルパー。**評価対象のcode(例:"BS4")と、実際に確立された方向(`aggregated_direction.value`、例:"PP1")が食い違う場合は正しくNOT_METとする**(export.pyのstrength表記バグ修正時と同じ教訓を再適用)

実データでの確認(GJB2 c.109G>A、以前の実行結果を再現したfixtureで検証):PS4=met(Supporting)・PP1=met(Supporting)・BS4=not_met(方向がPP1だったため)→合計スコア2 →Uncertain Significance。現状25/28コードが未評価であることも正しく可視化される。

#### 15-11-3. `criteria/stubs.py`・`criteria/registry.py`:未実装23基準のスタブ化

`classification.classify()`が「一度も評価していない基準」を「評価した上でnot_met」と取り違えないよう、**未実装の23基準(ALL_ACMG_CODES − IMPLEMENTED_CODES)全てに対してNOT_EVALUATED(`CriterionStatus.UNKNOWN`)を返す共通スタブ**を用意した。23個の個別ファイルではなく1ファイルに集約し、以下の2種類に分類して理由をコメントで明記している(いずれも`doc/BH26_participant_briefing_v3_en.md`・15-5節のスコープ整理と対応):

- **文献では判定できず、患者本人の臨床・遺伝子検査記録が必要**(PS2, PM3, PM6, BS2, BP2, BP5、計6件)
- **ルール・lookupベースで、他チームの自動判定モジュールが担当する想定**(PVS1, PS1, PM1, PM2, PM4, PM5, PP2〜PP5, BA1, BS1, BP1, BP3, BP4, BP6, BP7、計17件)

`registry.py`の`get_criterion_evidence(code, real_evidence)`が、実装済みcodeなら渡された実証拠を、未実装codeならスタブを返す統一インターフェースを提供する。これにより、呼び出し側は`ALL_ACMG_CODES`を単純にループするだけで28基準全部のCriterionEvidenceを組み立てられる。

新しくcriterionを実装する同僚は、該当コードを`IMPLEMENTED_CODES`に追加し、`stubs.py`の対象から自然に外れる形で移行できる(スタブ→実装への切り替えに`classification.py`側の変更は不要)。

回帰テストとして`test_classification.py`(24件、閾値境界・BA1特例・重複拒否・`from_aggregated_judgment`の方向不一致処理・28基準ループを検証)を追加し、既存203件と合わせて計227件全てパス。

### 15-12. `pipeline.py`への`classify()`組み込み(v10追記、2026-09-15)

`main()`のテストケース実行ループに、変異体単位の最終分類算出を追加した。狙いは、15-11節で作った`classification.classify()`が実際のパイプライン出力(LLM判定→`AggregatedJudgment`)から本当に最終分類まで一気通貫で出せることを確認すること。

**実装**:各テストケースの`judge_variant()`実行後、`from_aggregated_judgment(aggregated, case["criterion"])`で`CriterionEvidence`を作り、`(gene, hgvsc)`をキーにした辞書に蓄積する。全テストケース終了後、変異体ごとにその辞書の中身(実際に評価したcriterionの証拠) + `stubs.all_stub_evidence()`(未実装23基準のNOT_EVALUATEDプレースホルダ)を`classify()`に渡し、分類結果(category・score・met/not_met一覧)を表示する。

**設計判断(`registry.get_criterion_evidence()`を経由しない理由)**:`registry.py`の統一インターフェースは「実装済み5基準は必ず実証拠が要る」という契約(15-11-3節)だが、`main()`の`test_cases`は変異体ごとに1〜3criterionしか評価しない設計になっている(例:MYH7 c.1594T>CはPS3のみ、RUNX1 c.601C>TはPS4/PP1/BS4の3つ)。そのため`registry`をそのまま使うと「このデモでは評価しなかった実装済みcriterion」に対してValueErrorが飛んでしまう。今回は`registry`を経由せず、実際に評価した証拠+スタブ証拠だけを直接`classify()`に渡す方式にした。`classify()`は元々「渡されなかったcodeは`not_evaluated_codes`に落とす」設計になっているため、渡し忘れた実装済みcriterionと真に未実装(スタブ)のcriterionを`classify()`の出力上でも正しく区別できる(前者は`not_evaluated_codes`に残り、後者は既にUNKNOWN状態の証拠として渡されているので`not_evaluated_codes`には現れない)。出力ではこの「デモ内で未評価だった実装済みcriterion」を`not run in this demo`として明示し、「このプロジェクトが対応していない基準」と混同しないようにしている。

**動作確認**(MCP/vLLM接続なしで、実際のtest_casesと同じ入力組み合わせをfixtureで再現):
- RUNX1 c.601C>T(PS4=met/Supporting、PP1=met/Supporting、BS4=not_met——test_casesの3エントリと同じ組み合わせ)→ スコア2 → Uncertain Significance、`not run in this demo: BS3, PS3`
- MYH7 c.1594T>C(PS3のみ評価)→ スコア0 → Uncertain Significance、`not run in this demo: BS3, BS4, PP1, PS4`

いずれも期待通り。既存227件の回帰テストにも影響なし(全件パス確認済み)。実際のLLM/PubMed MCPを介した通し実行(`python3 -m acmg_pipeline.pipeline`)は、ネットワーク到達性がサンドボックス外前提のため未実施——このロジック自体の妥当性はfixtureベースの確認で担保している。

### 15-13. PP4の実装可否に関する検討ログ(未実装・議論のみ、2026-09-15)

現時点でPP4は`stubs.py`で`AUTOMATED_RULE_BASED_OTHER_TEAM`(他チームのルール/lookupベース自動判定モジュール担当)に分類しているが、この分類が実態に合っているか疑問が出たため調査した。**本節は設計上の論点整理であり、コード変更は未実施**(`stubs.py`の分類は現状のまま)。

#### 15-13-1. Richards et al. 2015原文でのPP4の定義

`doc/gim201530.pdf`(Richards et al. 2015、ACMG/AMP標準ガイドライン原著論文)Table 3およびPP4詳細解説(本文418ページ)を確認した。PP4は以下の4条件を要求する:

1. その遺伝子の検査の臨床感度が高い(陽性の患者が多い)
2. 患者が他疾患と紛らわしくない明確な症候群を呈している
3. 遺伝子に良性バリアントが多くない(gnomAD等で確認可能)
4. 家族歴が遺伝形式と矛盾しない

一見すると(1)(3)は遺伝子レベルの一般統計(文献・データベースから取得可能)だが、(2)(4)は**患者個人の表現型・家族歴**に強く依存しており、当初は「論文を読んでバリアント自体を判定する」という本パイプラインの型(PS3/BS3/PS4/PP1/BS4)に乗らない、むしろ患者カルテが必要な基準(15-5節のPS2/PM3/PM6/BS2/BP2/BP5グループ)に近いのではないか、という仮説が立った。

#### 15-13-2. Web調査:ClinGen SVIによる公式の定量化ガイダンス(2023年)

Web検索の結果、2015年のRichards論文の定性的な記述は既に更新されており、ClinGen SVI(Sequence Variant Interpretation)ワーキンググループが2023年に**PP4とPP1/BS4(家系内分離)を統合したベイズ点数法の正式ガイダンス**を*American Journal of Human Genetics*に発表していることが判明した([PMC10806742](https://pmc.ncbi.nlm.nih.gov/articles/PMC10806742/)、[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0002929723004044)、[PubMed](https://pubmed.ncbi.nlm.nih.gov/38103548/))。

要点:
- PP4とPP1/BS4は**分離不可能な一体の証拠**として扱う(同じ遺伝子座に対する重複証拠のため)
- **診断陽性率(diagnostic yield)**——「その表現型の患者の何%が当該遺伝子に病的バリアントを持つか」という、GeneReviews等の文献から得られる**遺伝子レベルの統計**——をルックアップ表(Table 2)で点数化する(例:陽性率70%→+4.0点)
- 遺伝子座の均一性(locus heterogeneity)を確認し、単一遺伝子で説明できる疾患(陽性率>90%)ならPP1(分離)は別途加算しない/複数遺伝子が原因になりうる疾患ならPP1点数(Table 3、家系分離データ)を追加する
- PP4+PP1の合計点数は**1アレルあたり+5.0点で上限キャップ**(連鎖不平衡で別の真の病的バリアントを拾っている可能性への配慮)
- 最終的にTable 4で点数をSupporting/Moderate/Strongにマッピングする

この発見により、当初の仮説(PP4=患者カルテが必要な基準)は修正が必要と判明した。実際には**PP4は文献由来の遺伝子レベル統計(診断陽性率)**であり、かつ**本プロジェクトが既に実装済みの`segregation.py`(PP1/BS4、家系分離データ)と表裏一体**という、想定より本パイプラインの型に近い基準だった。ただし、`classification.py`が前提とする「基準ごとに独立したMET/NOT_MET+strength」というカテゴリカルモデルとは相性が悪く、PP4+PP1をペアで点数キャップする専用ロジック(15-11-2節のTavtigian点数法をそのまま流用できない)が必要になる。

#### 15-13-3. 入力設計の検討:HPOは不要という結論

PP4実装の入力として構造化表現型オントロジー(HPO)が必要かを検討した。結論は**基本実装では不要**:ClinGen SVIのTable 2は「Cockayne症候群」のような**命名された疾患単位**でキーされており、HPO粒度の症状の組み合わせまでは要求しない。ERepoから既に取得しているMONDOコード(疾患名)を陽性率ルックアップのキーとしてそのまま使える。HPOが必要になるのは、単一の命名疾患に当てはまらない稀な症状の組み合わせで原因遺伝子候補が絞り込めないケース(ACGS Best Practice Guidelinesが言及する例外的シナリオ)のみで、これは将来の拡張であり初期スコープには含めない。

#### 15-13-4. 入力設計の検討:Clinical noteから拾うべき項目

PP4/PP1統合ロジックの各ステップに対応させて、Clinical noteから抽出すべき項目を5つに整理した:

1. **確定診断名/表現型**(疾患名)——Table 2の診断陽性率ルックアップのキー。ERepoのMONDOコードで代用可能な場合が多い
2. **遺伝形式の手がかりとなる家族歴**(発端者の同胞・親・祖父母等の罹患状況と血縁関係)——既存`segregation.py`の`FamilySegregationData`(affected_with_variant等)にそのままマッピング可能
3. **他遺伝子の検査結果**(パネル/エクソーム検査で他候補遺伝子が陰性だったか)——locus heterogeneity判定と、「他座位が分離しない→実質的に単一遺伝子相当の陽性率に上がる」という特殊ルールに必要
4. **検査手法**(単一遺伝子/パネル/エクソームのいずれか)——文献の陽性率と検査条件を揃えるために必要(揃わないと点数を割り当てられない)
5. **症候群としての特異性を示す具体的な臨床所見**(例:Gorlin症候群の基底細胞癌+掌蹠陥凹+歯原性角化嚢胞)——他疾患と紛らわしくない明確な症候群かどうかの判定材料

逆に、詳細なHPO粒度の全症状リストや、家族歴以外の既往歴・検査値等は、PP4/PP1のロジックが直接使わないため不要と判断した。

#### 15-13-5. 現時点の結論

- **未着手**:`ClinicalNoteData`のような入力スキーマ設計、PP4+PP1統合点数ロジックの実装、`stubs.py`での分類見直し(現状は`AUTOMATED_RULE_BASED_OTHER_TEAM`のまま)、いずれもコード変更なし
- **次にやるとすれば**:(a)`segregation.py`をPP4対応に拡張し、Table 2/3/4相当のルックアップ+キャップロジックを追加、(b)患者表現型・家族歴・他遺伝子検査結果を受け取る新しい入力スキーマを設計、(c)`classification.py`側もPP4+PP1のペア証拠を特別扱いできるよう拡張(現行の「基準ごと独立加算」モデルからの逸脱が必要)、の3点が主な作業になる見込み

### 15-14. PP4・Layer1(16基準)・デモケースのground truthテストケース収集(v10追記、2026-09-16、別文書に分離)

PP4とLayer1の16基準(判定ロジックは別メンバー担当)、およびdemocase 4変異体について、28基準を横断するground truthテストケースデータセットを整備した。この作業は実装設計そのものではなくテストデータ収集の記録であり、性質が異なるため別文書 **`test_case_ground_truth.md`** に切り出した。データセット自体(`test_data/full_criteria_ground_truth.py`、380件・34変異体・28基準中26基準をカバー)とその収集経緯、PP5・BP6が実データで収集できなかった件の記録は同文書を参照。
