# AGENTS.md

このファイルはAIコーディングエージェント向けの作業ガイドです。人間向けのセットアップ手順・アーキテクチャ概要は [`README.md`](README.md) を、設計判断の経緯は [`ps3_bs3_ps4_implementation_v10.md`](ps3_bs3_ps4_implementation_v10.md) を参照してください。

## プロジェクト概要

- LLM(vLLM上のgemma-4)とPubMed MCPを組み合わせ、ACMG/AMP 2015の28基準について、人間キュレーターの一次スクリーニングを高速化する下書き判定を生成するパイプライン
- fastapi + uvicorn をdockerで動かす形でAPIが実装されている
- 現状パイプラインの実装はmockである

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


補足: 今のrun_pipeline()はまだ暫定実装(vLLM/PubMed MCPを呼ばない)なので、VLLM_BASE_URL/VLLM_API_KEY等の環境変数は今回のコンテナには不要でした。実パイプラインが接続されたら、doc 8章にある通りdocker run -e VLLM_BASE_URL=... -e VLLM_API_KEY=...のように環境変数を渡す形に拡張することになります。


## 初期段階の構想

| 層 | 内容 | 対応するACMGコード | 実装 |
|---|---|---|---|
| 1. Automated Evidence | 構造化データから機械的に判定 | PVS1, PS1, PM1, PM2, PM4, PM5, PP2, PP3, PP5, BA1, BS1, BP1, BP3, BP4, BP6, BP7 | InterVar/AutoPVS1等既存の決定的ツール |
| 2. Semi-automated Evidence | 表現型を使う判定 | PP4 | PubCaseFinder(入力は自然文の臨床記述) |
| 3. Evidence Gap / Manual Review Required | 現行データでは判断できない点を明示 | PS2, PS3, PS4, PM3, PM6, PP1, BS2, BS3, BS4, BP2, BP5 | ルール+LLM(Evidence Gap Analysis) |
| 4. Documentation | Evidence・判断結果の自然文整理 | — | ローカルLLM |