# BH26 ExpertBoard(Pan-Asian Variant Review Network)参加者向け説明資料

**version: v3**

---

## 1. このプロジェクトについて

ExpertBoardは、遺伝性疾患のバリアント解釈(ACMG/AMP分類)を、AIが証拠(Evidence)集めを支援し、専門家が最終判断を下す形で半自動化する取り組みです。MedHackathon Asia 2026(シンガポール)での成果を引き継ぎ、今回のBH26で発展させます。

BH26での成果は、Rare Disease and Orphan Drugs Journal誌の特集号 *"Artificial Intelligence Applications in Rare Disease Research and Care"* にPerspectiveとして投稿する計画です。

### 基本原則

1. システムは候補となるEvidence(証拠)を準備する。専門家がその価値を判断する
2. Evidenceは再現可能・追跡可能でなければならない
3. LLMは病原性の最終分類を直接確定しない
4. LLMはEvidenceの発見・構造化・不足分析を積極的に支援してよい
5. オントロジー/DB識別子(HPO ID等)は、LLMの自由生成ではなく、権威あるソースに基づいて確定する
6. Evidenceの正確さは、参照する集団頻度データの多様性に左右される。欧米中心のデータだけでは、アジア系患者に対して誤ったEvidence(頻度に基づく証拠)を生み出しうる

### 一番伝えたいこと

このプロジェクトの主張は「Pan-Asianの頻度データを集めた」ことではありません。**「半自動でEvidenceを集める仕組み自体を、アジア独自で作る」**ことが本題です。Pan-Asianの頻度データによる実演は、その主張を支える説得力ある一例という位置づけです。

---

## 2. 何を作るか

### 2-1. 全体の流れ

```
Case(複数の候補Variant、正解データ付き)
    ↓
Evidence Preparation(基準ごとのモジュール)
    ↓
生成されたEvidence(Criterion + Strength + Source)
    ↓
Classificationモジュール(決定的ルールで算出)
    ↓
比較ツール(正解データとの突合)
    ↓
一致率・不一致箇所のレポート
```

最終classificationよりも先に、**Evidence Record(証拠候補の記録)**を中心データとして扱います。各ツール・LLMは最終判定ではなく「Candidate Evidence」を生成します。今回は、**人がその場でAccept/Reject/Modifyする画面は作りません**。あらかじめ用意した正解データと自動的に突き合わせて検証します。

### 2-2. Evidence Record(最小スキーマ、たたき台)

- **Variant**:どのバリアントに対する証拠か
- **Criterion候補**:例 PM2
- **Strength候補**:例 Supporting / Moderate / Strong / Very Strong
- **Evidence source**:どのツール・データから来たかに加え、実際の値・根拠も記載(例:「gnomAD AF=0.00008」「REVEL=0.97」「ClinVar VCV000012345, reviewed_by_expert_panel」)

これはたたき台です。Day1にチーム全員で見直し、現地で確定します。

### 2-3. Evidence Preparation(5層構造)

以下は固定区分ではなく、BH26時点の暫定実装方式です。Automated/Semi-automated/Manualの境界は「実際どこまで自動化できるか、やってみないと分からない」ため、これを検証すること自体が今回のタスクの一部です。

| 層 | 内容 | 対応するACMGコード | 実装 |
|---|---|---|---|
| 0. トリアージ | 数千バリアント→候補への絞り込み | — | `vep_parser.py`(拡張予定) |
| 1. Automated Evidence | 構造化データから機械的に判定 | PVS1, PS1, PM1, PM2, PM4, PM5, PP2, PP3, PP5, BA1, BS1, BP1, BP3, BP4, BP6, BP7 | InterVar/AutoPVS1等既存の決定的ツール |
| 2. Semi-automated Evidence | 表現型を使う判定 | PP4 | PubCaseFinder(入力は自然文の臨床記述) |
| 3. Evidence Gap / Manual Review Required | 現行データでは判断できない点を明示 | PS2, PS3, PS4, PM3, PM6, PP1, BS2, BS3, BS4, BP2, BP5 | ルール+LLM(Evidence Gap Analysis) |
| 4. Documentation | Evidence・判断結果の自然文整理 | — | ローカルLLM |

**Strengthの根拠**:独自に決めず、既存の較正済み基準に従います(PVS1はAutoPVS1の決定木、PM2はClinGen SVI推奨でSupportingに格下げ、PP3/BP4はPejaver et al. 2022の較正済み閾値、他はACMG/AMP 2015のデフォルト)。

**境界事例**:PM3(トランス配置確認)とBS2(健常者観察)は、trio VCFや集団健康データベースがあれば自動化の余地がありますが、現状はEvidence Gap層に置いています。当日確認したい論点です。

### 2-4. LLMの役割:Evidence Copilot

LLMは最終判定は行いません。以下4つの機能を担う"傘"です。

| 機能 | 内容 |
|---|---|
| 自然文要約 | 既に生成されたEvidenceを踏まえた自然文要約。新しいEvidenceは生成しない |
| Evidence Gap Analysis | 「今あるEvidenceと、判断に必要なEvidenceの差分」を提示。次に何をすべきかまで提示できると良い |
| 判定案生成とClassification評価 | ①Evidenceだけを見た盲検の判定案 → ②Classificationモジュールの算出結果を見せて評価コメント、の2段階 |
| HPO抽出 | 自由記述の臨床記述からHPO用語を抽出(PP4の入力として必要) |

---

## 3. 決定論 vs LLM比較実証(今回のPriority、3本柱)

「AI(LLM)がどこまでEvidenceをサポートできるか」を検証するため、以下3つを優先的に実施します。決定論的な手法・ツールと、LLMに同じ入力を与えて結果を比較する構成です。

1. **PVS1**:AutoPVS1の決定木 vs LLM。AutoPVS1には文献化された明確なロジックがあるため、LLMがそれを暗黙に再現できるかを測れます
2. **Evidence Gap Analysis**:単純な欠損チェック(null check) vs LLMの自由記述。LLMの出力が単なる言い換えでなく、実質的な付加価値を持つかを検証します
3. **HPO Grounding**:LLMによる素朴なHPO抽出(Critical path) vs TogoMCP経由のontology lookupでグラウンディングした抽出(時間が許せば実施)

**おまけ**:生のVCF・臨床記述をそのままLLMに渡し、Evidence整理からClassificationまで一気に出させた結果も、他のモジュールが完成した後の副産物として記録します。ただしDemo Caseの規模では統計的な説得力は薄いため、参考結果として扱います。

---

## 4. Demo Case(実在の症例をベースに作成済み)

4症例を用意しています。いずれも実在の症例報告・実データベースに基づいており、正解データも用意済みです。

| Case | 疾患・出典 | パターン |
|---|---|---|
| **1** | HCM、Han et al. 2026(実症例、中国人患者)DOI: [10.3389/fcvm.2026.1841777](https://doi.org/10.3389/fcvm.2026.1841777) | 表現型ミスマッチ:MYBPC3(真の原因候補)とKCNJ5(偶発的所見、HCMとは無関係な遺伝子)が同時に見つかる |
| **2** | 新生児期致死性HCM、Wang et al. 2025(実症例)DOI: [10.3389/fcvm.2025.1726463](https://doi.org/10.3389/fcvm.2025.1726463) | 複合ヘテロ接合:2つの変異が両方とも原因として必要(1つに絞り込む型ではない) |
| **3** | HCM(変異・頻度は実在、患者の物語は作成)MYBPC3側の出典DOI: [10.1161/CIRCGEN.122.003736](https://doi.org/10.1161/CIRCGEN.122.003736)(MYH7側はClinVar VCV000014104、DOIなし) | Pan-Asian:MYH7(ClinGen公式検証済みPathogenic)とMYBPC3 p.Glu334Lys(シンガポールSG10Kコホートで頻度が判明し判定が変わった実例) |
| **4** | ARVC、Human Genome Variation誌(実症例、日本人患者)DOI: [10.1038/s41439-022-00206-9](https://doi.org/10.1038/s41439-022-00206-9) | Pan-Asian:DSG2 p.Phe531Cys、日本人集団(ToMMo)でgnomAD全体より最大約14倍高頻度 |

各Caseには、AlphaMissense・AlphaGenomeによる独立したアノテーション(ミスセンス病原性予測・スプライシング影響予測)も付与済みです。

**正解データの信頼度には階層があります**:

- **Tier A**:ClinGen Expert Panelが基準ごとに正式検証済み(Case 3のMYH7のみ)
- **Tier B**:論文著者が根拠を明記
- **Tier C**:実在する不一致そのものが正解(判定が割れている、という事実自体が正解)

これは隠すのではなく、比較ツールの出力にも反映します。「証拠の出所・確からしさを保持する」という原則そのものです。

**全体としての合格ラインはCase 1の一致率検証**です。Case 2〜4は余力に応じて取り組みます。

---

## 5. 今回作らないもの

以下は今回のBH26のスコープには含めません。理由は、仕様を決め切るのが難しく仮決めしても現地の状況とずれるリスクがあること、参加者の関心がEvidence Preparation側に強いことです。

- 人がAccept/Reject/Modifyを行うレビュー画面
- ロールベースのBoard承認機構
- Keycloak/RBAC

かわりに、比較ツールによる正解データとの自動突合が、今回の「検証」に相当します。

---

## 6. モジュール分担とインターフェース

Keycloak/RBAC/レビュー画面がスコープ外のため、チームは**Evidence Recordという共通スキーマを介して繋がる、複数の独立したモジュール**として作業を分担します。

| モジュール | 入力 | 出力 |
|---|---|---|
| 0層:トリアージ | VCF+アノテーション | 絞り込まれた候補バリアントリスト |
| 1〜3層:Evidence生成 | 候補バリアント+アノテーション | Evidence Record群 |
| Evidence Gap Analysis | Evidence Record群 | 不足Evidenceのレポート |
| Classificationモジュール | Evidence Record群 | Variantごとのclassification |
| 比較ツール | Evidence Record群・Classification・正解データ | 一致率・不一致レポート |

**運用ルール**:

1. Day1 PM冒頭でEvidence Recordのスキーマを凍結し、全員がそれに沿って実装する
2. 各モジュールは独立して実装・テストできる(ダミーデータで単体テスト可能、他モジュールの完成を待たない)
3. 統合チェックポイントを最終日前に前倒しする(9/17を目安)

---

## 7. チーム編成

| 担当 | 人 | 内容 |
|---|---|---|
| 0層(トリアージ)の修正・拡張 | **Francisさん** | `vep_parser.py`のバグ修正、Pan-Asian対応。Ruthさんも共同でこの層の検証を担当 |
| 1層(Automated Evidence) | 当日決定 | PVS1, PS1-4, PM1-6等 |
| 2層(Semi-automated Evidence) | 当日決定 | PP4、PubCaseFinder連携・HPO抽出 |
| 3層(Evidence Gap / Manual) | 当日決定 | 手動確認が必要な基準のフラグ付け |
| Evidence Gap Analysis | 当日決定 | 複数基準を横断するため1人が通しで担当 |
| Classificationモジュール + 比較ツール | 当日決定 | 分類算出ロジックと正解データ突合ツール |
| 全体統括 | 高橋 | ファシリテーション、仕様調整、統合 |

---

## 8. スケジュール

| 日付 | 内容 |
|---|---|
| 9/13(日) | オープニングワークショップ |
| **9/14(月)** | AM:原則への合意→Evidence Recordスキーマ・Demo Case確認・担当割り振り / PM:各自モジュール実装開始 |
| **9/15(火)** | 実装継続。判断が必要な事項はこの日までに確定 |
| 9/16(水) | 実装継続 |
| 9/17(木) | モジュール統合、Evidence生成→Classification→比較ツールの統合、Demo Case実行 |
| 9/18(金) | AM:最終デバッグ / PM(Writethon):技術デッキ・サマリー作成 |
| 9/19(土) | Wrap-up発表 |

---

## 9. 当日決めること(現地での議論を歓迎します)

- Evidence Recordスキーマの詳細
- 各層・各基準の担当分担
- Pan-Asian実演の具体的な方法(トリアージ層で行うか、Automated Evidence層で行うか)
- 定足数・境界事例(PM3, BS2)の扱い
- `vep_parser.py`のアップデート方針(Francisさんと相談)

---

*詳細な設計根拠・引用文献・Ground Truthの内訳は、別途配布するDemo Case資料一式を参照してください。*
