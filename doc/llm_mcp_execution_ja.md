# LLM/MCP経由の実行

`scripts/` 配下のサンプル・デモスクリプトと `acmg_pipeline.pipeline` を、実際にvLLM・
PubMed MCP・TogoMCPへ接続して動かす手順。**`.env` の設定(`VLLM_BASE_URL` /
`VLLM_API_KEY`)とネットワーク到達性が必要**。ネットワーク不要のロジックのみの
確認方法は [`README.md`](../README.md) の該当節を参照。

```bash
# TogoMCP + PubMed MCPへの接続確認(smoke test)。固定の質問に対してLLMがツールを
# 呼び出しながら回答する一連の流れを実行。ログは logs/run_YYYYmmdd_HHMMSS.log に保存。
python3 scripts/check_mcp_llm_connection.py

# PS3/BS3・PS4・PP1/BS4の各判定エンジンで、コード内定義済みのテストケース
# (MYH7・PTEN等)を対象に、PubMed MCPから論文全文取得 → LLM判定 →
# ACMG/AMP分類までを一気通貫で実行。
# 結果は logs/ps3bs3_run_YYYYmmdd_HHMMSS.log と va_spec_output/ 以下のJSONに出力。
python3 -m acmg_pipeline.pipeline

# classify_variant_from_structured_input() を、APIが受け取るのと同じ
# dict形式("vcf"/"clinical_note"キー)から実際に叩くデモ。democaseの2症例
# (ERepoにPMIDがある変異/無い変異)で両方の分岐を確認する。
python3 scripts/demo_structured_input_run.py
```

CLI引数は用意されていません。対象の遺伝子/変異を変えたい場合は
`acmg_pipeline/pipeline.py` の `main()` 内 `test_cases`(505行目付近)を
直接編集してください。
