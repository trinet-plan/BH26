"""
main.py
doc/docker_api_deployment_plan_v1_ja.md の実装(2章エンドポイント設計、
7章エラーハンドリング)。ジョブ管理(job_store.py)とva_spec変換
(acmg_pipeline/va_spec_statement.py)のみを担い、判定ロジック自体は
acmg_pipeline/pipeline_interface.py の run_pipeline() に委譲する。

起動: uvicorn api.main:app --reload (9章のcurl手順を参照)
"""

from __future__ import annotations

from contextlib import AsyncExitStack

from fastapi import BackgroundTasks, FastAPI, HTTPException
from pydantic import BaseModel

from acmg_pipeline.constants import ALL_ACMG_CODES
from acmg_pipeline.gate import ERepoClient
from acmg_pipeline.pipeline import connect_pubmed
from acmg_pipeline.pipeline_interface import (
    needs_literature_workflow,
    parse_request_vcf,
    run_pipeline,
    run_selected_criteria,
)
from acmg_pipeline.va_spec_statement import build_variant_statement
from api.job_store import Job, create_job, delete_job, get_job, mark_failed, mark_running, mark_succeeded

app = FastAPI(title="ACMG判定API")


class ClassifyCriteriaRequest(BaseModel):
    vcf: str
    clinical_note: str = ""


class EvidenceLineByTargetCriteriaRequest(BaseModel):
    vcf: str
    clinical_note: str = ""
    criteria: list[str]


async def _execute_classify_job(job_id: str, record, clinical_note: str, gene: str, hgvsc: str, hgvsp: str) -> None:
    mark_running(job_id)
    try:
        async with AsyncExitStack() as stack:
            mcp = await connect_pubmed(stack)
            output = await run_pipeline(
                record, clinical_note,
                mcp=mcp, erepo_client=ERepoClient(),
            )
        va_spec = build_variant_statement(output, gene=gene, hgvsc=hgvsc, hgvsp=hgvsp)
    except Exception as e:
        mark_failed(job_id, {"type": type(e).__name__, "message": str(e)})
        return
    mark_succeeded(job_id, va_spec)


@app.post("/v1/classify_criteria", status_code=202)
def submit_classify_criteria(request: ClassifyCriteriaRequest, background_tasks: BackgroundTasks) -> dict:
    """全28基準を評価し、ACMG分類(P/LP/VUS/LB/B)まで出す。常に非同期(job_id + polling)。"""
    try:
        record = parse_request_vcf(request.vcf)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    gene = record.info.get("GENE", "")
    hgvsc = record.info.get("HGVSC", "")
    hgvsp = record.info.get("HGVSP", "")
    job = create_job(variant={"gene": gene, "hgvsc": hgvsc, "hgvsp": hgvsp})
    background_tasks.add_task(
        _execute_classify_job, job.job_id, record, request.clinical_note,
        gene, hgvsc, hgvsp,
    )
    return {"job_id": job.job_id, "status": job.status, "poll_url": f"/v1/classify_criteria/{job.job_id}"}


@app.get("/v1/classify_criteria/{job_id}")
def get_classify_criteria(job_id: str) -> Job:
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found")
    return job


@app.delete("/v1/classify_criteria/{job_id}", status_code=204)
def delete_classify_criteria(job_id: str) -> None:
    if not delete_job(job_id):
        raise HTTPException(status_code=404, detail="job not found")


def _canonicalize_criteria(criteria: list[str]) -> tuple[str, ...]:
    if not criteria:
        raise HTTPException(status_code=400, detail="criteria must not be empty")
    unknown = sorted(set(criteria) - set(ALL_ACMG_CODES))
    if unknown:
        raise HTTPException(status_code=400, detail=f"Unrecognized ACMG code(s): {unknown}")
    requested = set(criteria)
    return tuple(code for code in ALL_ACMG_CODES if code in requested)


async def _execute_target_criteria_job(job_id: str, record, clinical_note: str, criteria: tuple[str, ...]) -> None:
    mark_running(job_id)
    try:
        async with AsyncExitStack() as stack:
            mcp = await connect_pubmed(stack)
            evidence_lines = await run_selected_criteria(
                record, clinical_note, criteria,
                mcp=mcp, erepo_client=ERepoClient(),
            )
    except Exception as e:
        mark_failed(job_id, {"type": type(e).__name__, "message": str(e)})
        return
    mark_succeeded(job_id, evidence_lines)


@app.post("/v1/get_evidence_line_by_target_criteria")
async def submit_evidence_line_by_target_criteria(
    request: EvidenceLineByTargetCriteriaRequest, background_tasks: BackgroundTasks
) -> dict:
    """指定した基準だけを評価し、VA-Spec EvidenceLineを返す(分類は行わない)。

    自動判定基準/未実装基準のみなら同期で200即返し、文献(LLM)判定基準を
    1つでも含む場合は非同期(job_id + polling)になる。
    """
    try:
        record = parse_request_vcf(request.vcf)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    criteria = _canonicalize_criteria(request.criteria)

    if not needs_literature_workflow(criteria):
        evidence_lines = await run_selected_criteria(record, request.clinical_note, criteria)
        return {"status": "succeeded", "criteria": list(criteria), "evidence_lines": evidence_lines}

    gene = record.info.get("GENE", "")
    hgvsc = record.info.get("HGVSC", "")
    hgvsp = record.info.get("HGVSP", "")
    job = create_job(variant={"gene": gene, "hgvsc": hgvsc, "hgvsp": hgvsp, "criteria": list(criteria)})
    background_tasks.add_task(_execute_target_criteria_job, job.job_id, record, request.clinical_note, criteria)
    return {
        "job_id": job.job_id,
        "status": job.status,
        "poll_url": f"/v1/get_evidence_line_by_target_criteria/{job.job_id}",
    }


@app.get("/v1/get_evidence_line_by_target_criteria/{job_id}")
def get_evidence_line_by_target_criteria(job_id: str) -> Job:
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found")
    return job


@app.delete("/v1/get_evidence_line_by_target_criteria/{job_id}", status_code=204)
def delete_evidence_line_by_target_criteria(job_id: str) -> None:
    if not delete_job(job_id):
        raise HTTPException(status_code=404, detail="job not found")


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}
