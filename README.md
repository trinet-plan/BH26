# BH26 — ACMG/AMP PS3/BS3/PS4/PP1/BS4 判定パイプライン

LLM(vLLM上のgemma-4)とPubMed MCPを組み合わせ、変異のACMG/AMP分類基準のうち
「文献読解が必要な5基準」(PS3, BS3, PS4, PP1, BS4)について、人間キュレーターの
一次スクリーニングを高速化する下書き判定を生成するパイプラインです。

設計方針・検証結果の詳細は [`ps3_bs3_ps4_implementation_v10.md`](ps3_bs3_ps4_implementation_v10.md)
を参照してください。

## ディレクトリ構成

```
acmg_pipeline/            判定パイプライン本体
  classification.py         ACMG/AMPカテゴリの最終分類ロジック
  common.py                 共通データ型(PaperContribution, FinalResult 等)
  export.py                 GA4GH VA-Spec EvidenceLine形式での出力
  gate.py                   ClinGen ERepoの既存キュレーション確認ゲート
  pipeline.py                PubMed MCP + LLMを繋ぐメイン実行スクリプト
  criteria/                  基準ごとの判定ロジック(PS3/BS3, PS4, PP1/BS4)

democase/                 デモ用の臨床ノート・VCF・正解データ
doc/                       設計・参加者向け資料
va_spec_output/           パイプライン実行結果(VA-Spec JSON)
logs/                      実行ログ(git管理対象外)
ref_impl/                  参考実装アーカイブ

test_*.py                  スタンドアロンのテストスクリプト(pytest不要)
demo_ps3_bs3_judgment.py   LLM API接続なしのエンドツーエンドデモ
mcp_sample_multi.py        TogoMCP + PubMed MCP 接続サンプル
```

## セットアップ

### 前提

- Python 3.10+ 推奨(開発環境は Python 3.14.5)
- vLLMサーバー、およびPubMed MCPサーバーへのネットワーク到達性
  (`python -m acmg_pipeline.pipeline` / `mcp_sample_multi.py` の実行時のみ必要)

### 1. リポジトリを取得

```bash
git clone <このリポジトリのURL>
cd <リポジトリ名>
```

### 2. 仮想環境を作成(推奨)

```bash
python3 -m venv .venv
source .venv/bin/activate      # Windowsは .venv\Scripts\activate
```

### 3. 依存パッケージをインストール

`requirements.txt` は未整備のため、コード内のimportに基づき手動でインストールしてください。

```bash
pip install mcp openai requests ga4gh.va_spec ga4gh.core openpyxl
```

- `ga4gh.va_spec` / `ga4gh.core` — `acmg_pipeline/export.py` がVA-Spec形式での出力に使用
- `openpyxl` — `democase/annotation_alphamissense_alphagenome_v1.xlsx` を読む場合に使用

### 4. 環境変数を設定

vLLMサーバーのURL・APIキーは `.env` で管理します(リポジトリには含まれません)。

```bash
cp .env.example .env
```

`.env` を開き、`VLLM_BASE_URL` / `VLLM_API_KEY` を実際の値に書き換えてください
(値の入手方法はプロジェクト管理者に確認してください)。`.env` が無い、または
値が空の場合、`mcp_sample_multi.py` / `acmg_pipeline/pipeline.py` は起動時に
`RuntimeError` を送出します。

## 動作確認

### ネットワーク不要(ロジックのみ)

```bash
python3 test_classification.py
python3 test_ps3_bs3_judgment.py
python3 test_ps3_bs3_ps4_gate.py
python3 test_ps3_bs3_ps4_gate_full.py
python3 demo_ps3_bs3_judgment.py
```

いずれもpytest不要のスタンドアロンスクリプトで、末尾に `N passed, M failed` の
ように結果が表示されます(`test_classification.py` は失敗があると exit code 1)。

### LLM/MCP経由の実行(`.env` 設定 + ネットワーク到達性が必要)

```bash
# TogoMCP + PubMed MCPに接続し、固定の質問に対してLLMがツールを呼び出しながら
# 回答する一連の流れを実行。ログは logs/run_YYYYmmdd_HHMMSS.log に保存。
python3 mcp_sample_multi.py

# PS3/BS3・PS4・PP1/BS4の各判定エンジンで、コード内定義済みのテストケース
# (MYH7・PTEN等)を対象に、PubMed MCPから論文全文取得 → LLM判定 →
# ACMG/AMP分類までを一気通貫で実行。
# 結果は logs/ps3bs3_run_YYYYmmdd_HHMMSS.log と va_spec_output/ 以下のJSONに出力。
python3 -m acmg_pipeline.pipeline
```

CLI引数は用意されていません。対象の遺伝子/変異を変えたい場合は
`acmg_pipeline/pipeline.py` の `main()` 内 `test_cases`(505行目付近)を
直接編集してください。

## ライセンス・注意事項

- `.env` にはAPIキーが含まれるため、絶対にコミットしないでください
  (`.gitignore` で除外済み)。
- `logs/` 以下の実行ログ、`__pycache__/` はgit管理対象外です。
