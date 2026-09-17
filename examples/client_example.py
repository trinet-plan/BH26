"""
client_example.py
api/main.py が公開するACMG判定API(POST/GET /v1/variant)を呼び出す
最小構成のクライアントスケルトン。ジョブ登録 → ポーリング → 結果取得の
一連の流れのみを示す。認証やリトライ等は含まれていないため、実運用では
そのまま使わず、必要な処理を足して利用すること。

使い方(あらかじめ api/main.py を起動しておく。README.md の
「APIのテスト方法」参照):

    python3 examples/client_example.py democase/case1_api_input_case1-noise2.json
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import requests

DEFAULT_BASE_URL = "http://localhost:8000"


def submit_variant(base_url: str, vcf: str, clinical_note: str = "") -> str:
    """POST /v1/variant でジョブを登録し、job_idを返す。"""
    resp = requests.post(f"{base_url}/v1/variant", json={"vcf": vcf, "clinical_note": clinical_note})
    resp.raise_for_status()
    return resp.json()["job_id"]


def wait_for_result(base_url: str, job_id: str, interval_sec: float = 2.0, timeout_sec: float = 300.0) -> dict:
    """GET /v1/variant/{job_id} を status が succeeded/failed になるまでポーリングする。"""
    deadline = time.monotonic() + timeout_sec
    while True:
        resp = requests.get(f"{base_url}/v1/variant/{job_id}")
        resp.raise_for_status()
        job = resp.json()
        if job["status"] in ("succeeded", "failed"):
            return job
        if time.monotonic() > deadline:
            raise TimeoutError(f"job {job_id} did not finish within {timeout_sec}s (status={job['status']})")
        time.sleep(interval_sec)


def main() -> None:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <api_input.json>", file=sys.stderr)
        sys.exit(1)

    case = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    job_id = submit_variant(DEFAULT_BASE_URL, vcf=case["vcf"], clinical_note=case.get("clinical_note", ""))
    print(f"submitted job_id={job_id}")

    job = wait_for_result(DEFAULT_BASE_URL, job_id)
    if job["status"] == "failed":
        print(f"job failed: {job['error']}", file=sys.stderr)
        sys.exit(1)

    print(json.dumps(job["va_spec"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
