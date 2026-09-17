"""
target_criteria_client_example.py
api/main.py が公開するACMG判定API(POST/GET /v1/get_evidence_line_by_target_criteria)を
呼び出す最小構成のクライアントスケルトン。指定したcriteriaが自動判定/未実装のみなら
サーバーは200を同期即返しし、文献(LLM)判定criteriaを1つでも含む場合は202 + job_idを
返して非同期ジョブになる - この2パターンをどちらも扱う流れを示す。認証やリトライ等は
含まれていないため、実運用ではそのまま使わず、必要な処理を足して利用すること。

使い方(あらかじめ api/main.py を起動しておく。README.md の
「APIのテスト方法」参照):

    # 自動判定criteriaのみ -> 同期で即返る
    python3 examples/target_criteria_client_example.py democase/case1_api_input_case1-noise2.json PM2 BA1

    # 文献(LLM)判定criteriaを含む -> job_id + ポーリングになる
    python3 examples/target_criteria_client_example.py democase/case1_api_input_case1-noise2.json PS3 BS3
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import requests

DEFAULT_BASE_URL = "http://localhost:8000"


def submit_target_criteria(base_url: str, vcf: str, criteria: list[str], clinical_note: str = "") -> dict:
    """POST /v1/get_evidence_line_by_target_criteria を呼び出す。

    自動判定/未実装のcriteriaのみなら {"status": "succeeded", "evidence_lines": ...}
    がそのまま返る(ジョブ登録なし)。文献判定criteriaを1つでも含む場合は
    {"job_id": ..., "status": "queued", ...} が返るので、呼び出し側が
    wait_for_result() でポーリングする必要がある。
    """
    resp = requests.post(
        f"{base_url}/v1/get_evidence_line_by_target_criteria",
        json={"vcf": vcf, "clinical_note": clinical_note, "criteria": criteria},
    )
    resp.raise_for_status()
    return resp.json()


def wait_for_result(base_url: str, job_id: str, interval_sec: float = 2.0, timeout_sec: float = 300.0) -> dict:
    """GET /v1/get_evidence_line_by_target_criteria/{job_id} を
    status が succeeded/failed になるまでポーリングする。"""
    deadline = time.monotonic() + timeout_sec
    while True:
        resp = requests.get(f"{base_url}/v1/get_evidence_line_by_target_criteria/{job_id}")
        resp.raise_for_status()
        job = resp.json()
        if job["status"] in ("succeeded", "failed"):
            return job
        if time.monotonic() > deadline:
            raise TimeoutError(f"job {job_id} did not finish within {timeout_sec}s (status={job['status']})")
        time.sleep(interval_sec)


def main() -> None:
    if len(sys.argv) < 3:
        print(f"usage: {sys.argv[0]} <api_input.json> <CRITERION> [<CRITERION> ...]", file=sys.stderr)
        sys.exit(1)

    case = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    criteria = sys.argv[2:]

    result = submit_target_criteria(
        DEFAULT_BASE_URL, vcf=case["vcf"], criteria=criteria, clinical_note=case.get("clinical_note", ""),
    )

    if "job_id" in result:
        print(f"submitted job_id={result['job_id']} (criteria require the literature/LLM workflow)")
        result = wait_for_result(DEFAULT_BASE_URL, result["job_id"])
        if result["status"] == "failed":
            print(f"job failed: {result['error']}", file=sys.stderr)
            sys.exit(1)
        evidence_lines = result["va_spec"]
    else:
        print("answered synchronously (automated/stub criteria only)")
        evidence_lines = result["evidence_lines"]

    print(json.dumps(evidence_lines, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
