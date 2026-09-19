# ACMG/AMP判定パイプライン — 6日間の開発成果

BH26ハッカソン(6日間、2026-09-19が最終日)における開発成果のまとめ。
[前日版](hackathon_summary_2026-09-18_ja.md)を踏襲し、64変異体フル検証・4デモケース
検証の結果と、それに続くPVS1メカニズムゲートの改善を反映した最新版。

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
| ポリシーにより無効化 | `PP5 BP6` | 2 | 常にUNKNOWN |
| 未実装 | `PS2 PM3 PM6 BS2 BP2 BP5` | 6 | NOT_EVALUATED明示 |

最終分類は Tavtigian et al. 2018 のベイズ点数体系(Supporting=1 / Moderate=2 / Strong=4 /
Very Strong=8、閾値: ≥10 Pathogenic、6-9 Likely Pathogenic、0-5 VUS、-1〜-6 Likely
Benign、≤-7 Benign)で統合し、GA4GH VA-Spec形式のEvidenceLineとしてAPI(FastAPI +
Docker)から返す。

## 2. VA-Spec出力の設計改善

キュレーターが実際に読む出力の品質向上に力を入れた領域。

- **重複・不整合の解消** — `bh26AssessmentDetails`という単一オブジェクトに固めていた
  フィールドを、`curatorHints`と同じレベルの個別extensionへ全面フラット化。
- **方向性の不整合修正** — 良性方向のエビデンス(BS3)が常に病原性方向として出力
  されていた不整合を修正。
- **GA4GH VA-Specへのフィードバック** — 拡張作業中に見つかったスキーマ制約を5件の
  issueとして[doc/ga4gh_va_spec_feedback_v1_en.md](ga4gh_va_spec_feedback_v1_en.md)に
  まとめ、GA4GH本家に[Issue #448](https://github.com/ga4gh/va-spec/issues/448)として
  実際に提出。GA4GH側メンテナ(korikuzma氏)から、次期`1.1.0-ballot.2026-09`で
  issue #2(汎用エビデンス項目型`DataItem`の追加)が既に対応されていることを確認した
  一方、issue #1(`Method.reportedIn`の複数ドキュメント対応)・#3(`Direction`の
  「未評価」表現)は未対応、`EvidenceLine`クラス自体が`Statement`に統合される
  大きな構造変更が進行中であることも判明した。

## 3. 診断名(condition)解決の拡張

PVS1の疾患メカニズムゲートやPP1/BS4/PP4は、患者の診断名がMONDO疾患クラスに解決
されていないと機能しない。これまでの「ERepoのcondition欄を使う」に加え、以下を追加した。

- **MONDO解決のリトライ** — `resolve_diagnosis_mondo()`の最終候補選択はLLM呼び出しで
  あり、ボーダーラインな診断名(例:「hypertrophic cardiomyopathy with apical
  ventricular aneurysm」)で約20%の確率でサンプリングブレによりNOT_FOUNDを返すことを
  実測で確認。3回→5回のリトライに変更し、実測10/10成功まで改善。
- **検証スクリプト自体の配線修正** — `run_integrated_validation_demo.py`が本番API
  パス(`pipeline_interface.py`)と異なり、`resolve_diagnosis_mondo()`を一度も呼んで
  いなかった配線漏れを発見・修正。これによりMYBPC3 c.278delA等でPVS1が正しく
  very_strongに到達するようになった。
- **文献ベースのcondition検索(新規)** — ERepoに存在しない変異体(例: MYBPC3の
  democase由来3変異体)向けに、`acmg_pipeline/condition_from_literature_search.py`を
  新規実装。PubMed検索+LLM判定で「この変異体が引き起こす疾患」を論文から抽出し、
  `resolve_diagnosis_mondo()`に渡してMONDO解決する。略語(例:「HCM」)ではなく
  正式病名で抽出するようプロンプトを調整し、MYBPC3 c.278delAで実証確認。

## 4. PVS1メカニズムゲートの拡張 — 「不明」と「矛盾」を区別する

64変異体フル検証でPVS1のMET率が0%だった原因を調査し、2件の意図的なポリシー変更を行った。

- **関連疾患(parent/child)への適用** — ケースの診断名とキュレーション済み疾患がMONDO
  階層上で親子関係にあるだけ(完全一致ではない)の場合、これまでは常にMANUAL_REVIEWで
  停止していた。全レコードがメカニズム確立に同意している場合はPVS1を適用し、
  「関連疾患経由である」旨をcuratorHintsに開示する方式に変更。RUNX1 c.601C>Tで
  UNKNOWN→MET(very_strong)への改善を実証。
- **遺伝形式が不明な場合への適用** — ケース側の遺伝形式が未記録(「不明」)なだけの
  場合と、明示的に矛盾する場合を区別。前者かつキュレーション済み遺伝形式が一意で
  メカニズム確立に同意している場合はPVS1を適用し、「唯一知られている遺伝形式と仮定した」
  旨を開示。RPE65 c.495+1dupではメカニズムゲート自体は通過するようになったが、
  別の未解決課題(NMD予測データ欠如)によりまだMETには到達していない。

いずれも「情報が無い」ことと「情報が矛盾する」ことを区別し、前者はcuratorHints開示付きで
適用、後者は従来通り保留するという一貫した設計方針。

## 5. 精度検証

### 4デモケース(全4ケース、本番APIパス相当のスクリプトで実行)

正解が判明している3件中 **1/3 (33.3%)** が変異体レベルで一致。criteria単位では
**9/15 (60.0%)**。

| Code | 一致 | Code | 一致 |
|---|---:|---|---:|
| PVS1 | 2/2 (100%) | PS3 | 0/2 (0%) |
| PM1 | 1/1 (100%) | PP1 | 0/1 (0%) |
| PM2 | 3/3 (100%) | PS2 | 0/1 (0%、未実装code) |
| PM5 | 1/1 (100%) | | |
| PP3 | 2/2 (100%) | | |
| PS4 | 1/3 (33.3%) | | |

自動系(PVS1・PM1・PM2・PM5・PP3)は全件一致。文献エンジン(PS3・PS4)とPP1は引き続き弱く、
「対象変異を直接扱う論文が乏しい」という従来からの傾向と一致する。

### ERepo 64変異体フル検証(初の全件実行)

前日までは12件・16件のサブセットに限定していたが、今回**64変異体全件を一晩(約3.5時間、
想定より大幅に速く完了)かけて実行**した。

- **59/64件を評価**(5件スキップ: 座標未解決4件[PTEN, USH2A, TECTA, TPM1] +
  正解データ欠如1件[DSG2])
- **変異体レベル一致: 19/59 (32.2%)**
- **criteria単位(MET→MET再現率): 108/224 (48.2%)**

| Code | 一致 | Code | 一致 |
|---|---:|---|---:|
| BP4 | 8/8 (100%) | PP4 | 2/11 (18.2%) |
| PM2 | 37/38 (97.4%) | PS4 | 4/23 (17.4%) |
| PP3 | 30/31 (96.8%) | PP2 | 1/6 (16.7%) |
| PM1 | 12/15 (80.0%) | PS3 | 0/21 (0%) |
| BS1 | 4/5 (80.0%) | PVS1 | 0/5 (0%)※ |
| PM5 | 4/6 (66.7%) | PS2/PM3/PM6/BP7/BA1/BS2等 | 0%(大半は未実装code) |
| PP1 | 5/15 (33.3%) | | |
| BS3 | 1/5 (20.0%) | | |

※PVS1が0/5だった原因を調査し、上記セクション4の2件のポリシー変更で改善(RUNX1は
実証確認済み、RPE65は部分改善)。フル64件セットへの再反映(再実行)はまだ行っていない。

## 6. その他の成果物

- **GA4GH VA-Spec feedback** — [Issue #448](https://github.com/ga4gh/va-spec/issues/448)として提出済み
- **`condition_from_literature_search.py`(新規)** — ERepoに存在しない変異体向けの文献ベースcondition解決
- **PP1/BS4向けライブ文献検索(`pp1_segregation_search.py`)** — 家系共分離データが無い場合の文献検索。実証済みで一部変異体でPP1が実際にMETに到達
- **64変異体ground truthデータセット** — `test_data/full_criteria_ground_truth.py`
- **論文日本語訳** — Perspective論文の全訳

## 7. 既知の課題 / 引き継ぎ事項

- **未実装6基準**(PS2, PM3, PM6, BS2, BP2, BP5)— 今後の実装対象。
- **RPE65のNMD予測データ欠如** — PVS1のメカニズムゲートは通過するようになったが、
  変異体レベルの決定木でNMD予測情報が取得できずMETに至らない。原因を特定済み:
  Ensembl VEPに`numbers=1`で直接問い合わせて比較したところ、単純な
  splice_donor_variant(例: MYBPC3 c.2905+1G>A)ではintron番号が正しく返るのに対し、
  RPE65 c.495+1dupはエクソン-イントロン境界そのものに重複挿入が生じる
  frameshift_variant + splice_region_variantという複合的な帰結のため、VEPがexon/intron
  番号をどちらも返さない。`acmg_pipeline/providers/nmd.py`は「番号が取れない時は推測せず
  レコード無しとして正直に報告する」設計であり、これは仕様通りの動作。恒久対応(CDS位置
  とトランスクリプト構造から自前でexon番号を計算する等、別のアノテーション経路の追加)は
  相応の実装コストがかかるため、**今回は見送り、課題として引き継ぐ**。現状の
  「キュレーター判断に委ねる(MANUAL_REVIEW/not_met)」動作を維持する。
- **64変異体フル検証への修正反映** — 今回のPVS1関連修正(セクション4)を反映した
  再実行はまだ行っていない。再実行すればMET率のさらなる改善が見込まれる。
- **文献エンジンの精度** — LLM判定自体の欠陥ではなく、対象変異に言及する文献の希少性が
  ボトルネック。
- **AutoPVS1本体との統合** — 現状は参照リンクのみで、実際の自動判定ロジックとの統合は
  スコープ外。

---
BH26 ACMG/AMP判定パイプライン &middot; ハッカソン6日間の開発サマリー(2026-09-19時点)
