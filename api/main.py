"""
main.py
doc/docker_api_deployment_plan_v1_ja.md の実装(2章エンドポイント設計、
7章エラーハンドリング)。ジョブ管理(job_store.py)とva_spec変換
(acmg_pipeline/va_spec_statement.py)のみを担い、判定ロジック自体は
acmg_pipeline/pipeline_interface.py の run_pipeline() に委譲する。

起動: uvicorn api.main:app --reload (9章のcurl手順を参照)
"""

from __future__ import annotations

from fastapi import BackgroundTasks, FastAPI, HTTPException
from pydantic import BaseModel

from acmg_pipeline.api_input import ApiCaseInput
from acmg_pipeline.pipeline_interface import run_pipeline
from acmg_pipeline.va_spec_statement import build_variant_statement
from api.job_store import Job, create_job, delete_job, get_job, mark_failed, mark_running, mark_succeeded

app = FastAPI(title="ACMG判定API")


class VariantRequest(BaseModel):
    vcf: str
    clinical_note: str = ""


def _execute_job(job_id: str, case: ApiCaseInput, gene: str, hgvsc: str, hgvsp: str) -> None:
    mark_running(job_id)
    try:
        output = run_pipeline(case)
        va_spec = build_variant_statement(output, gene=gene, hgvsc=hgvsc, hgvsp=hgvsp)
    except Exception as e:
        mark_failed(job_id, {"type": type(e).__name__, "message": str(e)})
        return
    mark_succeeded(job_id, va_spec)


@app.post("/v1/variant", status_code=202)
def submit_variant(request: VariantRequest, background_tasks: BackgroundTasks) -> dict:
    case = ApiCaseInput(vcf=request.vcf, clinical_note=request.clinical_note)
    try:
        record = case.parse_vcf().record
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    gene = record.info.get("GENE", "")
    hgvsc = record.info.get("HGVSC", "")
    hgvsp = record.info.get("HGVSP", "")
    job = create_job(variant={"gene": gene, "hgvsc": hgvsc, "hgvsp": hgvsp})
    background_tasks.add_task(_execute_job, job.job_id, case, gene, hgvsc, hgvsp)
    return {"job_id": job.job_id, "status": job.status, "poll_url": f"/v1/variant/{job.job_id}"}


@app.get("/v1/variant/{job_id}")
def get_variant(job_id: str) -> Job:
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found")
    return job


@app.delete("/v1/variant/{job_id}", status_code=204)
def delete_variant(job_id: str) -> None:
    if not delete_job(job_id):
        raise HTTPException(status_code=404, detail="job not found")


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}
