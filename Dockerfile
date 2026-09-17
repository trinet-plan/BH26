# doc/docker_api_deployment_plan_v1_ja.md 8章。
FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY acmg_pipeline/ acmg_pipeline/
COPY api/ api/
COPY config/ config/

EXPOSE 8000
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
