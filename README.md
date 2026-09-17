# BH26 — ACMG/AMP 統合EvidenceLine生成パイプライン

LLM(vLLM上のgemma-4)とPubMed MCPを組み合わせ、変異のACMG/AMP分類基準のうち
「文献読解が必要な5基準」(PS3, BS3, PS4, PP1, BS4)と、ルール・データベースに
基づく16基準を統合し、人間キュレーター向けの下書き判定を生成します。残る7基準も
NOT_EVALUATEDとして保持し、1変異につき全28基準のVA-Spec EvidenceLineを返します。

設計方針・検証結果の詳細は [`ps3_bs3_ps4_implementation_v10.md`](ps3_bs3_ps4_implementation_v10.md)
を参照してください。

外部データの現在の取得経路、未対応・review待ちの情報、TogoVarへ置換する場合の
境界は [`doc/external_data_coverage_ja.md`](doc/external_data_coverage_ja.md) に分けて整理しています。

## ディレクトリ構成

```
acmg_pipeline/            判定パイプライン本体
  classification.py         ACMG/AMPカテゴリの最終分類ロジック
  common.py                 共通データ型(PaperContribution, FinalResult 等)
  export.py                 GA4GH VA-Spec EvidenceLine形式での出力
  gate.py                   ClinGen ERepoの既存キュレーション確認ゲート
  vcf_record.py             criterion共通入力のVariantRecord
  clinical_note.py          criterion共通入力のClinicalNoteExtraction
  inputs.py                 共通入力から遺伝子・HGVS等を読む補助関数
  pipeline.py                PubMed MCP + LLMを繋ぐメイン実行スクリプト
  criteria/                  基準ごとの判定ロジック(PS3/BS3, PS4, PP1/BS4)

acmg/                     移植した自動判定器(16基準)
config/                   自動判定の閾値・疾患文脈
tests/                    自動判定と統合インターフェースのpytest

democase/                 デモ用の臨床ノート・VCF・正解データ
doc/                       設計・参加者向け資料
examples/                 APIクライアントのスケルトンコード
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

```bash
pip install -r requirements.txt
```

- `ga4gh.va_spec`(`ga4gh.core`/`ga4gh.vrs`を自動導入) — `acmg_pipeline/export.py` がVA-Spec形式での出力に使用
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

## 統合インターフェース

他のスクリプトからは `acmg_pipeline.pipeline.evaluate_variant_evidence_lines()` を
呼び出します。入力型はmain側の共通クラスそのものです。

```python
lines = await evaluate_variant_evidence_lines(
    variant,                 # acmg_pipeline.vcf_record.VariantRecord
    clinical_note,           # acmg_pipeline.clinical_note.ClinicalNoteExtraction
    normalized_evidence=normalized_evidence,
    automated_config=config,
    mcp=pubmed_session,
    erepo_client=erepo_client,
)
```

`normalized_evidence` と `automated_config` は呼出し側が明示的に渡します。VCFのINFOを
暗黙の判定根拠へ変換しません。返却順はACMGの標準順で固定され、各EvidenceLineの
`bh26AssessmentDetails` に実状態、`referenceLink` に補助URLが入ります。補助URLは
エビデンスとして参照済みであることを意味しないため `reportedIn` には入れません。

### ネットワーク不要(ロジックのみ)

```bash
python3 test_classification.py
python3 test_ps3_bs3_judgment.py
python3 test_ps3_bs3_ps4_gate.py
python3 test_ps3_bs3_ps4_gate_full.py
python3 demo_ps3_bs3_judgment.py
pytest tests -q
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

## criterion入力インターフェース

実装済みの PS3、BS3、PS4、PP1、BS4 と、未実装criterionのstubを含む
全28 criterionの公開入力は、他チームとの結合用に次の2オブジェクトへ
統一しています。

```python
def judge(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
) -> ...:
    ...
```

文献プロンプトを作る `build_prompt()` も同じ2オブジェクトを受け取ります。
遺伝子名・HGVS・別表記は `VariantRecord.info` の `GENE`、`HGVSC`、
`HGVSP`、`EQUIVALENTS` から読みます。文献だけを評価し患者情報を使わない
criterionでも、呼び出し境界を揃えるため空の `ClinicalNoteExtraction()` を
渡してください。出力は従来どおり GA4GH VA-Spec `EvidenceLine` です。

## APIサーバーの起動

1.  ビルド

    ````bash
    docker build -t acmg-api:v1 .
    ```

2. 起動(ホストの8000番で公開)

    ````bash
    docker run -d --name acmg-api -p 8000:8000 acmg-api:v1
    ```

3. 動作確認(これまでと同じcurl手順がそのまま使える)

    ````bash
    curl -s http://localhost:8000/health
    ```

    ````bash
    JOB_ID=$(curl -s -X POST http://localhost:8000/v1/variant \
      -H "Content-Type: application/json" \
      -d @democase/case1_api_input_case1-noise2.json | jq -r .job_id)
    ````

    ````bash
    curl -s http://localhost:8000/v1/variant/$JOB_ID | jq .
    ```

4. 停止

    ````bash
    docker rm -f acmg-api
    ```

上記1〜2の代わりに [`compose.yaml`](compose.yaml) を使ってもよい
(`docker compose up -d --build` / 停止は `docker compose down`)。

## APIクライアントの例

APIクライアントの実装例は以下にあります。

* [`examples/client_example.py`](examples/client_example.py) 

APIサーバー起動後に以下のコマンドを実行してください。

```bash
python3 examples/client_example.py democase/case1_api_input_case1-noise2.json
```

## ライセンス・注意事項

- `.env` にはAPIキーが含まれるため、絶対にコミットしないでください
  (`.gitignore` で除外済み)。
- `logs/` 以下の実行ログ、`__pycache__/` はgit管理対象外です。
