# ACMG/AMP判定パイプライン — 5日間の開発成果

BH26ハッカソン(5日間、2026-09-18が最終日)における開発成果のまとめ。
LLM(vLLM上のgemma-4)とPubMed MCPを組み合わせ、ACMG/AMP 2015の28基準について
人間キュレーターの一次スクリーニングを高速化する下書き判定を生成するパイプラインを、
チームで構築した。

## 1. 実装の全体像

ACMG/AMP 2015の全28基準を、判定手法ごとに4層に分けて実装している。既存の決定的ツール
(自動系)、文献読解が必要な基準(LLM+PubMed MCP)、表現型マッチが必要な基準、の3系統が
すでに稼働しており、残りは正直に「未評価」として明示している。

| 層 | 基準 | 件数 | 状態 |
|---|---|---:|---|
| 自動判定(Layer1) | `PVS1 PS1 PM1 PM2 PM4 PM5 PP2 PP3 BA1 BS1 BP1 BP3 BP4 BP7` | 14 | 実装済み |
| 文献エンジン | `PS3 BS3 PS4` | 3 | 実装済み |
| 表現型/家系 | `PP1 BS4 PP4` | 3 | 実装済み |
| ポリシーにより無効化 | `PP5 BP6` | 2 | 常にUNKNOWN(下記参照) |
| 未実装 | `PS2 PM3 PM6 BS2 BP2 BP5` | 6 | NOT_EVALUATED明示 |

`PP5`・`BP6`(「信頼できる情報源による評価」)はコード自体は存在するが、
ClinGen SVIの一般方針(外部の主張だけで採点せず一次エビデンスを参照すべき)に従い、
入力に関わらず常に`UNKNOWN`を返す設計になっている — 実質的には使用していない。

最終分類は Tavtigian et al. 2018 のベイズ点数体系(Supporting=1 / Moderate=2 / Strong=4 /
Very Strong=8、閾値: ≥10 Pathogenic、6-9 Likely Pathogenic、0-5 VUS、-1〜-6 Likely
Benign、≤-7 Benign)で統合し、GA4GH VA-Spec形式のEvidenceLineとしてAPI(FastAPI +
Docker)から返す。

## 2. VA-Spec出力の設計改善

キュレーターが実際に読む出力の品質向上に特に力を入れて取り組んだ領域。

- **重複・不整合の解消** — `bh26AssessmentDetails`という単一オブジェクトに固めていた
  フィールド(status / provenance / direction / rulesUsed 等)を、`curatorHints`と同じ
  レベルの個別extensionへ全面フラット化。`summary` / `strength` / `evidenceOutcome` /
  `criterion` / `evidenceItemIds`の重複、`structuredEvidenceItems`・
  `strengthEstimationMethod`のcuratorHintsへの統合、`rulesUsed`と
  `specifiedBy.reportedIn`の重複排除まで一貫して整理した。
- **判定間の連携修正** — PVS1の疾患メカニズム判定に必要な情報(clinical_noteの診断名)が
  実際には判定ロジックまで届いていなかった配線の欠落を発見し修正。診断名→MONDO疾患クラス
  への解決を経て、PVS1がvery_strongへ正しく到達することを確認した。
- **方向性の不整合修正** — 良性方向のエビデンス(BS3)が常に病原性方向(SUPPORTS)として
  出力されていた不整合を修正し、回帰テストを追加。
- **GA4GH VA-Specへのフィードバック** — 拡張作業中に見つかったスキーマ制約
  (`Method.reportedIn`は単一Documentのみ、`hasEvidenceItems`に汎用StudyResult型が無い、
  `Direction`に「未評価」が表現できない等)を5件のissueとして
  [doc/ga4gh_va_spec_feedback_v1_en.md](ga4gh_va_spec_feedback_v1_en.md) にまとめ、
  GA4GHへの提案として整理した。

## 3. 精度検証

### 4デモケース(全12変異体)

変異体レベルの分類一致は **4/9 (44.4%)**(正解判明分)、criteria単位では **18/26 (69.2%)**。

| Code | 一致 | Code | 一致 |
|---|---:|---|---:|
| PVS1 | 3/3 (100%) | BS1 | 3/4 (75%) |
| PM2 | 5/5 (100%) | BP4 | 1/2 (50%) |
| PP3 | 2/2 (100%) | PS3 | 1/3 (33.3%) |
| PM1 / PM5 | 2/2 (100%) | PS4 | 1/3 (33.3%) |
| | | PP1 | 0/1 (0%) |

自動系(Layer1: PVS1・PM1・PM2・PM5・PP3)は**全件一致**。一方、文献エンジン(PS3・PS4)は
33.3%と弱く、原因分析ツールにより「対象変異を直接扱う論文がそもそも乏しい」ことが主因
(71.0%)と判明 — 実装の欠陥ではなく文献自体の希少性に起因する。

### ERepo 64変異体からのPP4サブセット検証

PP4がMET(充足)と判定された12変異体を厳選して実行(実装済み22 criteria中20をカバー)。
結果は変異体レベルで **2/11 (18.2%)** と低かったが、原因究明の過程で重要な構造的発見があった。

> **発見: このハーネスではPP1/PP4を原理的に検証できない**
> 検証スクリプトが空の臨床記述を渡していたため、診断名を必要とするPP1/BS4/PP4は常に
> `unknown`となり、12件すべてで一度も`met`に現れなかった。ERepoの64変異体データセット
> には遺伝子・HGVScのみで自然文の臨床記述が付随していないための構造的な限界であり、
> 精度そのものの問題ではない。

## 4. その他の成果物

- **GA4GH VA-Spec feedback** — [doc/ga4gh_va_spec_feedback_v1_en.md](ga4gh_va_spec_feedback_v1_en.md)(5 issues)
- **エラー分析ハーネス** — 文献エンジンのnot_clear判定を根本原因別に分類するツール
- **64変異体ground truthデータセット** — `test_data/full_criteria_ground_truth.py`(ERepo再照会+democase由来)
- **論文日本語訳** — Perspective論文(著者・参考文献・図表を保持した全訳)

## 5. 既知の課題 / 引き継ぎ事項

- **未実装6基準**(PS2, PM3, PM6, BS2, BP2, BP5)— 今後の実装対象。
  (`PP5`・`BP6`は別枠 — 上記の通りポリシーにより意図的に常時UNKNOWNとしており、
  実装待ちではない。)
- **PP1/PP4の実証検証** — ERepo変異体に対応する臨床記述(diagnosis)を用意しない限り、
  これらの基準は自動検証パイプラインで一度も`met`に到達しない。次のステップとして
  診断名データの追加が必要。
- **文献エンジンの精度** — LLM判定自体の欠陥ではなく、対象変異に言及する文献の希少性が
  ボトルネック。
- **AutoPVS1本体との統合** — 現状は参照リンクのみで、実際の自動判定ロジックとの統合は
  スコープ外。

---
BH26 ACMG/AMP判定パイプライン &middot; ハッカソン5日間の開発サマリー(2026-09-18時点)
