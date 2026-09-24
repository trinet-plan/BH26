# AGENTS.md

このファイルはAIコーディングエージェント向けの作業ガイドです。人間向けのセットアップ手順・アーキテクチャ概要は [`README.md`](README.md) を、設計判断の経緯は [`ps3_bs3_ps4_implementation_v10.md`](ps3_bs3_ps4_implementation_v10.md) を参照してください。

## プロジェクト概要

- LLM(vLLM上のgemma-4)とPubMed MCPを組み合わせ、ACMG/AMP 2015の28基準について、人間キュレーターの一次スクリーニングを高速化する下書き判定を生成するパイプライン
- fastapi + uvicorn をdockerで動かす形でAPIが実装されている
- 実装状況(mockではない): PS3/BS3/PS4はPubMed MCP + LLM(vLLM)による文献読解、PP1/BS4/PP4はClinGen 2024 Bayesian-pointsエンジン(診断イールド統計をPubMed+LLMで検索。見つからなければ正直にUNKNOWN)、残り16基準はルール・外部DB照合による自動判定(`acmg_pipeline/criteria/`, `acmg_pipeline/providers/`)。未実装分はNOT_EVALUATEDのstubとして返す。正確な内訳・件数は `acmg_pipeline/constants.py` の `LITERATURE_CODES` / `PHENOTYPE_SEGREGATION_CODES` / `AUTOMATED_CODES` / `STUB_CODES` を参照(READMEの「5基準」表記とは数が合わないので注意)

## test_data/ ディレクトリ

28基準分のground truthデータと検証ハーネス。`collectors/`にClinGen ERepo/Ensembl MANEから生データを再取得するスクリプトとデータ本体(`full_criteria_ground_truth.py`)、`fetched_data/`に収集済みの生JSONスナップショット、直下の`criterion_runner.py`に任意の判定関数をこのデータで自動採点するハーネスがある。判定ロジック自体はここには無い(`acmg_pipeline/criteria/`側)。

## 検証スクリプト

`evaluation/`(ground truth/専門家パネルとの一致率を測る)と`scripts/`(実行・出力確認用、正解比較なし)にある。すべてリポジトリルートから`.venv/bin/python3 evaluation/xxx.py`のように実行(内部で`config/` `test_data/` `democase/`等をルート基準の絶対パスで参照するため)。

- `evaluation/run_validation_64.py` — 文献3基準(PS3/BS3/PS4)をground truthと比較
- `evaluation/run_automated_validation_64.py` — 自動判定16基準をground truthと比較
- `evaluation/run_integrated_validation_64.py` — 統合28基準+`classify()`を64件で検証
- `evaluation/run_integrated_validation_demo.py` — 同上、democase 4症例のみの高速版
- `evaluation/run_pvs1_validation.py` — PVS1を専門家パネルと比較(`--contexts erepo|curated`)
- `scripts/run_automated_api_va_spec.py` — APIパス(自動判定のみ)をin-processで実行しVA-Spec出力
- `scripts/run_all_tests.py` — `standalone_tests/`配下の`test_*.py`を一括実行

## APIのテスト方法

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

4. 停止

    ````bash
    docker rm -f acmg-api
    ```


補足: `acmg_pipeline/pipeline_interface.py` の `run_pipeline()` と `run_selected_criteria()` はどちらも毎回 `extract_clinical_note()`(`clinical_note.py`が遅延importする`clinical_extraction.py`)と `hpo_mondo_extraction.resolve_diagnosis_mondo()` を呼ぶ。この2モジュールはimport時にVLLM_BASE_URL/VLLM_API_KEYが無いと`RuntimeError`を送出するため、`/v1/classify_criteria`・`/v1/get_evidence_line_by_target_criteria`のどちらも(自動判定基準のみを指定した場合を含め)環境変数なしでは動かない。`docker run -e VLLM_BASE_URL=... -e VLLM_API_KEY=...`または`.env`で必ず設定すること(doc 8章参照)。

※ `pipeline_interface.py` 冒頭のモジュールdocstring(「実際の判定ロジックが未接続で全28コードをNOT_EVALUATEDとして返す」)は古い記述のまま残っている。実装(`run_pipeline()`/`run_selected_criteria()`本体)はすでに`evaluate_variant_evidence_lines()`/`evaluate_selected_criteria()`に接続済みで、この点でもコード側のコメントとREADMEに乖離がある。


## 初期段階の構想

| 層 | 内容 | 対応するACMGコード | 実装 |
|---|---|---|---|
| 1. Automated Evidence | 構造化データから機械的に判定 | PVS1, PS1, PM1, PM2, PM4, PM5, PP2, PP3, PP5, BA1, BS1, BP1, BP3, BP4, BP6, BP7 | InterVar/AutoPVS1等既存の決定的ツール |
| 2. Semi-automated Evidence | 表現型を使う判定 | PP4 | PubCaseFinder(入力は自然文の臨床記述) |
| 3. Evidence Gap / Manual Review Required | 現行データでは判断できない点を明示 | PS2, PS3, PS4, PM3, PM6, PP1, BS2, BS3, BS4, BP2, BP5 | ルール+LLM(Evidence Gap Analysis) |
| 4. Documentation | Evidence・判断結果の自然文整理 | — | ローカルLLM |