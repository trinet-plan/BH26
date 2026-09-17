"""
job_store.py
doc/docker_api_deployment_plan_v1_ja.md 3章のジョブモデル・状態遷移
(queued -> running -> succeeded/failed)をプロセス内メモリで管理する。

6章の通りv1(BackgroundTasks構成)向けの最小実装。バックグラウンドタスクは
Starletteのスレッドプール上で動くため、GETによる読み取りと同時に走りうる
書き込みに備えて素朴なロックだけ入れてある。永続化やマルチワーカー対応は
将来タスクキューへ移行する際に差し替える(6章「未決事項」参照)。
"""

from __future__ import annotations

import threading
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class Job:
    job_id: str
    status: str  # queued / running / succeeded / failed
    created_at: str
    variant: dict
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    va_spec: Optional[dict] = None
    error: Optional[dict] = None


_jobs: dict[str, Job] = {}
_lock = threading.Lock()


def create_job(variant: dict) -> Job:
    job = Job(job_id=str(uuid.uuid4()), status="queued", created_at=_now(), variant=variant)
    with _lock:
        _jobs[job.job_id] = job
    return job


def get_job(job_id: str) -> Optional[Job]:
    with _lock:
        return _jobs.get(job_id)


def delete_job(job_id: str) -> bool:
    with _lock:
        return _jobs.pop(job_id, None) is not None


def mark_running(job_id: str) -> None:
    with _lock:
        _jobs[job_id].status = "running"
        _jobs[job_id].started_at = _now()


def mark_succeeded(job_id: str, va_spec: dict) -> None:
    with _lock:
        job = _jobs[job_id]
        job.status = "succeeded"
        job.va_spec = va_spec
        job.completed_at = _now()


def mark_failed(job_id: str, error: dict) -> None:
    with _lock:
        job = _jobs[job_id]
        job.status = "failed"
        job.error = error
        job.completed_at = _now()
