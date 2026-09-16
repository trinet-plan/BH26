# doc/docker_api_deployment_plan_v1_ja.md 8章。
# APIプロセス(api/main.py)のみを起動するイメージ。vLLM/PubMed MCP本体はこの
# コンテナの外にあり、含まれない(将来 run_pipeline() を実接続する際、
# VLLM_BASE_URL/VLLM_API_KEY等を環境変数として渡す想定)。
FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY acmg_pipeline/ acmg_pipeline/
COPY api/ api/

EXPOSE 8000
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
