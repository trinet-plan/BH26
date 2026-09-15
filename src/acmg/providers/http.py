"""Replayable, content-checked response cache with bounded HTTP retries."""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


class FetchError(ValueError):
    """Data unavailable; callers must not interpret this as a negative observation."""


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      allow_nan=False)


class CachedHttpClient:
    def __init__(self, cache_dir, *, offline=False, timeout=20, attempts=3, delay=0.4,
                 opener=urlopen, sleeper=time.sleep):
        if not 1 <= attempts <= 5 or timeout <= 0 or delay < 0:
            raise ValueError("Invalid HTTP retry policy")
        self.cache_dir = Path(cache_dir)
        self.offline = offline
        self.timeout = timeout
        self.attempts = attempts
        self.delay = delay
        self.opener = opener
        self.sleeper = sleeper
        self.used = {}
        self.network_used = False

    def fetch(self, url, *, data=None, response_format="json", dataset_version=None,
              allow_application_errors=False):
        if not url.startswith("https://") or response_format not in {"json", "text"}:
            raise ValueError("HTTPS and JSON/text response formats are required")
        request_key = {"url": url, "data": data, "format": response_format,
                       "dataset_version": dataset_version}
        # Preserve existing cache identities for ordinary requests.
        if allow_application_errors:
            request_key["allow_application_errors"] = True
        key = hashlib.sha256(canonical_json(request_key).encode()).hexdigest()
        path = self.cache_dir / f"{key}.json"
        if path.exists():
            document = json.loads(path.read_text(encoding="utf-8"))
            body = document.get("body")
            digest = hashlib.sha256(canonical_json(body).encode()).hexdigest()
            if document.get("request") != request_key or digest != document.get("body_sha256"):
                raise FetchError("CACHE_INTEGRITY_FAILURE")
            self.used[key] = {"path": str(path), "sha256": digest, "retrieved_at": document["retrieved_at"]}
            return document
        if self.offline:
            raise FetchError(f"OFFLINE_CACHE_MISS:{key}")
        payload = canonical_json(data).encode() if data is not None else None
        request = Request(url, data=payload, headers={"Accept": "application/json" if response_format == "json" else "text/plain",
                          "Content-Type": "application/json", "User-Agent": "BH26-ACMG/0.1"})
        for attempt in range(self.attempts):
            self.sleeper(self.delay * (2 ** attempt))
            try:
                with self.opener(request, timeout=self.timeout) as response:
                    raw = response.read().decode("utf-8")
                self.network_used = True
                body = json.loads(raw) if response_format == "json" else raw
                if (not allow_application_errors and isinstance(body, dict)
                        and (body.get("errors") or body.get("error"))):
                    raise FetchError("REMOTE_APPLICATION_ERROR")
                break
            except HTTPError as exc:
                if exc.code not in {429, 500, 502, 503, 504} or attempt == self.attempts - 1:
                    raise FetchError(f"HTTP_ERROR:{exc.code}") from exc
            except (URLError, TimeoutError) as exc:
                if attempt == self.attempts - 1:
                    raise FetchError(f"NETWORK_UNAVAILABLE:{type(exc).__name__}") from exc
            except json.JSONDecodeError as exc:
                raise FetchError("INVALID_REMOTE_JSON") from exc
        document = {"request": request_key, "body": body,
                    "body_sha256": hashlib.sha256(canonical_json(body).encode()).hexdigest(),
                    "retrieved_at": datetime.now(timezone.utc).isoformat()}
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        # One client per CLI run; no background writers mutate the cache.
        path.write_text(canonical_json(document) + "\n", encoding="utf-8")
        self.used[key] = {"path": str(path), "sha256": document["body_sha256"],
                          "retrieved_at": document["retrieved_at"]}
        return document
