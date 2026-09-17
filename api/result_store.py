"""
result_store.py
Write each completed API result to a directory, so a run leaves an artifact
behind instead of only a response body.

[Why this is opt-in]
  job_store.py keeps results in process memory: a restart, a DELETE, or a
  second worker loses them, and the synchronous
  /v1/get_evidence_line_by_target_criteria response is never stored at all.
  Persisting is a deployment decision, not a pipeline one - the container
  needs a mounted volume to write into - so nothing is written unless
  ACMG_API_OUTPUT_DIR names a directory. Unset, the API behaves exactly as
  before.

[Why a failure here does not fail the request]
  The API's obligation is to answer the caller. A full disk or a read-only
  mount must not turn a completed evaluation into an error response, so a
  write failure is logged and the result is still returned. It is logged
  rather than swallowed: a deployment that asked for artifacts and is not
  getting them should be able to find out why.

[What is written]
  One JSON file per completed evaluation, holding the payload plus the
  request context needed to identify it later - the variant, the criteria
  asked for, the endpoint, and when it finished. The payload alone does not
  say which request produced it.
"""

from __future__ import annotations

import json
import logging
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

LOGGER = logging.getLogger(__name__)
OUTPUT_DIR_ENV = "ACMG_API_OUTPUT_DIR"

# A filename has to survive both a Linux container and a Windows bind mount,
# and HGVS is full of characters neither agrees on (">", ":", "*"). Everything
# outside this set becomes "_", so c.3382G>A reads as c_3382G_A rather than
# being dropped or quietly renaming one variant onto another.
_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")


def output_dir() -> Optional[Path]:
    """The configured artifact directory, or None when persistence is off."""
    value = os.environ.get(OUTPUT_DIR_ENV, "").strip()
    return Path(value) if value else None


def _slug(value: str, fallback: str) -> str:
    cleaned = _UNSAFE.sub("_", (value or "").strip()).strip("._-")
    return cleaned[:60] or fallback


def filename(*, endpoint: str, variant: dict, identifier: str, finished_at: str) -> str:
    """A name that sorts by time and still says what it holds.

    The job id (or a request id for a synchronous call) is kept last so two
    evaluations of the same variant in the same second cannot collide.
    """
    stamp = _UNSAFE.sub("", finished_at)[:15]
    gene = _slug(variant.get("gene", ""), "nogene")
    hgvsc = _slug(variant.get("hgvsc", ""), "novariant")
    return f"{stamp}-{gene}-{hgvsc}-{_slug(endpoint, 'result')}-{identifier[:8]}.json"


def save(*, endpoint: str, variant: dict, criteria, payload, status: str,
         identifier: str) -> Optional[Path]:
    """Write one completed result. Returns the path written, or None.

    None means either that persistence is off or that the write failed; the
    caller carries on in both cases, because neither is the caller's problem
    to report to the person who asked for a judgment.
    """
    directory = output_dir()
    if directory is None:
        return None
    finished_at = datetime.now(timezone.utc).isoformat()
    document = {
        "schema_version": "1.0",
        "endpoint": endpoint,
        "status": status,
        "finished_at": finished_at,
        "identifier": identifier,
        "variant": variant,
        "criteria": list(criteria) if criteria is not None else None,
        "result": payload,
    }
    path = directory / filename(endpoint=endpoint, variant=variant,
                                identifier=identifier, finished_at=finished_at)
    try:
        directory.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(document, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
    except (OSError, TypeError, ValueError) as error:
        LOGGER.warning("Could not write API result to %s: %s", path, error)
        return None
    return path
