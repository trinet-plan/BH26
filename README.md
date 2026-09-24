# BH26 — ACMG/AMP 統合EvidenceLine生成パイプライン

LLM(vLLM上のgemma-4)とPubMed MCPを組み合わせ、変異のACMG/AMP分類基準のうち
「文献読解が必要な5基準」(PS3, BS3, PS4, PP1, BS4)と、ルール・データベースに
基づく16基準を統合し、人間キュレーター向けの下書き判定を生成します。残る7基準も
NOT_EVALUATEDとして保持し、1変異につき全28基準のVA-Spec EvidenceLineを返します。

設計方針・検証結果の詳細は [`ps3_bs3_ps4_implementation_v10.md`](doc/ps3_bs3_ps4_implementation_v10.md)
を参照してください。

外部データの現在の取得経路、未対応・review待ちの情報、TogoVarへ置換する場合の
境界は [`doc/external_data_coverage_ja.md`](doc/external_data_coverage_ja.md) に分けて整理しています。
変異単体の集団頻度はTogoVar GRCh38 APIを取得口とします。利用するデータソースは
`config/demo-rules.json` の `population_sources` でプロバイダーとその内部ソースを選択します。
現在はTogoVar内の `gnomad` と `tommo`（gnomAD exomes/genomesとToMMo）だけを利用します。
公開済みClinGen/VCEP仕様から転記したBS1のレビュー用DRAFTは
[`config/bs1_thresholds_draft.json`](config/bs1_thresholds_draft.json) にあります。DRAFTは実行時の閾値ではありません。

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
config/                   自動判定の閾値・BA1例外・CLI互換用文脈
tests/                    自動判定と統合インターフェースのpytest
test_data/                28基準分のground truthデータ一式(判定ロジックは含まない)。
  collectors/               ClinGen ERepo/Ensembl MANEから生データを再取得するスクリプト + データ本体(full_criteria_ground_truth.py)
  fetched_data/             collectors/が生成する生JSONスナップショット。criterion_runner.pyがこのデータで任意の判定関数を採点する

democase/                 デモ用の臨床ノート・VCF・正解データ
doc/                       設計・参加者向け資料
examples/                 APIクライアントのスケルトンコード
va_spec_output/           パイプライン実行結果(VA-Spec JSON)
logs/                      実行ログ(git管理対象外)
ref_impl/                  参考実装アーカイブ

standalone_tests/         スタンドアロンのtest_*.py(pytest不要)
```

## セットアップ

### 前提

- Python 3.10+ 推奨(開発環境は Python 3.14.5)
- vLLMサーバー、およびPubMed MCPサーバーへのネットワーク到達性
  (`python -m acmg_pipeline.pipeline` / `scripts/check_mcp_llm_connection.py` の実行時のみ必要)

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
値が空の場合、`scripts/check_mcp_llm_connection.py` / `acmg_pipeline/pipeline.py` は起動時に
`RuntimeError` を送出します。

## 動作確認

`evaluation/`(ground truthとの一致率測定)と`scripts/`(実行・出力確認用)。リポジトリルートから`.venv/bin/python3`で実行。

```bash
.venv/bin/python3 evaluation/run_validation_64.py              # 文献3基準(PS3/BS3/PS4)をground truthと比較
.venv/bin/python3 evaluation/run_automated_validation_64.py    # 自動判定16基準をground truthと比較
.venv/bin/python3 evaluation/run_integrated_validation_64.py   # 統合28基準+classify()を64件で検証
.venv/bin/python3 evaluation/run_integrated_validation_demo.py # 同上、democase 4症例のみの高速版
.venv/bin/python3 evaluation/run_pvs1_validation.py            # PVS1を専門家パネルと比較(--contexts erepo|curated)
.venv/bin/python3 scripts/run_automated_api_va_spec.py         # APIパス(自動判定のみ)をin-processで実行しVA-Spec出力
.venv/bin/python3 scripts/run_all_tests.py                     # standalone_tests/配下のtest_*.pyを一括実行
```

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
python3 standalone_tests/test_classification.py
python3 standalone_tests/test_ps3_bs3_judgment.py
python3 standalone_tests/test_ps3_bs3_ps4_gate.py
python3 standalone_tests/test_ps3_bs3_ps4_gate_full.py
python3 scripts/demo_ps3_bs3_judgment.py
pytest tests -q
```

いずれもpytest不要のスタンドアロンスクリプトで、末尾に `N passed, M failed` の
ように結果が表示されます(`test_classification.py` は失敗があると exit code 1)。

### LLM/MCP経由の実行(`.env` 設定 + ネットワーク到達性が必要)

vLLM・PubMed MCP・TogoMCPへ実際に接続して動かす手順は
[`doc/llm_mcp_execution_ja.md`](doc/llm_mcp_execution_ja.md) を参照してください。

## 自動判定CLIの実行手順(democase)

`acmg_pipeline.automated_cli` が16本の自動化criterionを評価します。**4ステップで、
3番目を飛ばすとPP2/BP1がミスセンス全件で `unknown` になります。**

```bash
# 1. 原本VCFを監査(原本は変更しない)
python3 -m acmg_pipeline.automated_cli audit-demo   --input-dir democase --output-dir work/run/audit

# 2. identityを解決し、外部プロバイダから normalised evidence を集める
python3 -m acmg_pipeline.automated_cli prepare-demo-online   --input-dir democase   --cache-dir tests/fixtures/ensembl-cache   --evidence-cache-dir tests/fixtures/external-cache   --output-dir work/run/prepared --ensembl-release 116   --with-gnomad --gnomad-release 4.1.1   --with-clinvar --clinvar-release 2026-09-15   --with-pm1-hotspot --with-dbnsfp   --rules config/demo-rules.json --offline

# 3. 遺伝子--疾患の機序assessmentを展開して evidence に足す(PP2/BP1に必須)
python3 -m acmg_pipeline.automated_cli build-gene-disease-evidence   --input config/gene-disease-review-decisions.json   --base-evidence work/run/prepared/evidence.json   --output work/run/evidence.json

# 4. 評価する
python3 -m acmg_pipeline.automated_cli evaluate   --input work/run/prepared/variants.json   --evidence work/run/evidence.json   --config config/demo-rules.json   --context config/curated-context.json   --output-dir work/run/evaluated --offline
```

**ステップ3を省略できない理由。** ステップ2の `--with-clingen-dosage` /
`--with-gene2phenotype` が取得するのは遺伝子レベルのLoF機序で、PVS1のゲートには答えますが
「ミスセンスが疾患機序か」は記録していません。PP2とBP1はそれを読むため、transcriptスコープ
の機序・変異スペクトラムreviewである `config/gene-disease-review-decisions.json` の展開が
要ります。投入しない場合、PP2/BP1は推測せず `unknown` を返し、`missing_inputs` に不足
フィールド名を列挙します。

**APIサーバーでは不要です。** サーバー経路は同じ
`config/gene-disease-review-decisions.json` を `demo-rules.json` の
`gene_disease_assessments` から自動で読みます。CLIとサーバーが同じ変異に別の根拠で答えない
ようにするためで、review済みの判断が無い変異だけが統計的示唆
(`gene_disease_draft`)へフォールバックします。

reviewされた決定が無い遺伝子の変異には assessment が付きません(他遺伝子の判断を一般化
しないため)。democaseでは7グループが28件中15件をカバーします。

`--offline` はキャッシュ済みレスポンスのみを使い、ネットワークへ出ません。実行可能な形の
同じ手順が `tests/test_demo_pipeline.py` にあります。

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

## PVS1が判定に至らない場合の出力

PVS1は、決定木が判定に到達していても**疾患特異的なLoF機序が確認できていない**
場合は `met` を返しません。strength・direction・evidenceOutcome はいずれも
付かず、VA-Spec上は他の結論の出なかったcriterionと同じ
`not_met` + curator hint として報告されます。

**`status` だけで「適用されたか」が分かります。** 他のcriterionと同じ読み方で、
consumer 側で追加のフィールドを確認する必要はありません。

### 捨てていない情報

判定に至らなくても、variant側の評価結果は `bh26AssessmentDetails` に残ります。

| 場所 | 内容 |
|---|---|
| `provenance.preliminary_assessment` | 決定木がどのノードまで到達したか、疾患機序が確認できれば何になるか(`candidate_strength`) |
| `provenance.candidate_conditions` | 疾患未指定時の候補疾患(病名・機序・遺伝形式付き) |
| `missingInputs` | 何が足りないか |
| `curatorHints` | 人が何を判断すればよいか |

`candidate_strength` は criterion の `strength` には**なりません**。
「疾患機序が確認できればこうなる」という仮定の結論であって、確認された
という主張ではないためです。VA-Spec にも `PVS1 supports` としては出ません。

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
    JOB_ID=$(curl -s -X POST http://localhost:8000/v1/classify_criteria \
      -H "Content-Type: application/json" \
      -d @democase/case1_api_input_case1-noise2.json | jq -r .job_id)
    ````

    ````bash
    curl -s http://localhost:8000/v1/classify_criteria/$JOB_ID | jq .
    ```

    `/v1/get_evidence_line_by_target_criteria` は指定した基準だけを評価する。
    自動判定/未実装の基準のみなら同期で即座に結果が返る(`jq`でdemoケースの
    JSONに`criteria`を足して投げる)。


    ```bash
    curl -s -X POST http://localhost:8000/v1/get_evidence_line_by_target_criteria \
      -H "Content-Type: application/json" \
      -d "$(jq '. + {criteria: ["PM2", "BA1"]}' democase/case1_api_input_case1-noise2.json)" | jq .
    ```

    文献(LLM)判定基準(`PS3`/`BS3`/`PS4`)を1つでも含めると、
    `/v1/classify_criteria` と同様にjob_id + pollingになる。

    ```bash
    JOB_ID=$(curl -s -X POST http://localhost:8000/v1/get_evidence_line_by_target_criteria \
      -H "Content-Type: application/json" \
      -d "$(jq '. + {criteria: ["PS3", "BS3"]}' democase/case1_api_input_case1-noise2.json)" | jq -r .job_id)
    ```

    ```bash
    curl -s http://localhost:8000/v1/get_evidence_line_by_target_criteria/$JOB_ID | jq .
    ```

4. 停止

    ````bash
    docker rm -f acmg-api
    ```

上記1〜2の代わりに [`compose.yaml`](compose.yaml) を使ってもよい
(`docker compose up -d --build` / 停止は `docker compose down`)。

### 判定結果をファイルに残す

ジョブはプロセス内メモリにしか残らず(`api/job_store.py`)、同期応答に至っては
どこにも保存されない。`ACMG_API_OUTPUT_DIR` にディレクトリを指定すると、
完了した判定が1件1ファイルのJSONとして書き出される。未指定なら何も書かない。

```bash
docker run -d --name acmg-api -p 8000:8000 \
  -e ACMG_API_OUTPUT_DIR=/data/results \
  -v "$PWD/api_results:/data/results" \
  acmg-api:v1
```

ファイル名は `<完了時刻>-<gene>-<hgvsc>-<endpoint>-<id>.json`。中身は判定結果
(`result`)に加えて、どのリクエストが生んだものかを示す `variant` / `criteria` /
`endpoint` / `status` / `finished_at` を持つ。失敗したジョブも `status: "failed"`
として残る。

書き込みに失敗してもリクエスト自体は成功する(判定結果は返る)。失敗は
サーバーログに `Could not write API result` として出る。

## APIエンドポイント仕様

各エンドポイントのリクエスト/レスポンス形・エラー条件は
[`doc/api_spec_ja.md`](doc/api_spec_ja.md) を参照してください。

## APIクライアントの例

APIクライアントの実装例は以下にあります。

* [`examples/client_example.py`](examples/client_example.py) — `/v1/classify_criteria`
  (全28基準を評価しACMG分類まで出す。常に非同期)

  APIサーバー起動後に以下のコマンドを実行してください。

  ```bash
  python3 examples/client_example.py democase/case1_api_input_case1-noise2.json
  ```

* [`examples/target_criteria_client_example.py`](examples/target_criteria_client_example.py) —
  `/v1/get_evidence_line_by_target_criteria`
  (指定した基準だけを評価しVA-Spec EvidenceLineを返す。分類は行わない)

  自動判定/未実装の基準のみを指定した場合は同期で即座に結果が返ります。

  ```bash
  python3 examples/target_criteria_client_example.py democase/case1_api_input_case1-noise2.json PM2 BA1
  ```

  文献(LLM)判定基準(`PS3`/`BS3`/`PS4`)を1つでも含めると非同期(job_id + ポーリング)になります。

  ```bash
  python3 examples/target_criteria_client_example.py democase/case1_api_input_case1-noise2.json PS3 BS3
  ```

## ライセンス・注意事項

- `.env` にはAPIキーが含まれるため、絶対にコミットしないでください
  (`.gitignore` で除外済み)。
- `logs/` 以下の実行ログ、`__pycache__/` はgit管理対象外です。
