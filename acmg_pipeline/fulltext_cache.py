"""
acmg_pipeline/fulltext_cache.py

A disk-backed, dict-like cache for pipeline.fetch_full_text()'s
(full_text, note) results, keyed by PMID. Added 2026-09-16 in response to
two limits of the original per-variant in-memory cache (test_case_ground_
truth.md's 227-pair validation run):

  1. A fresh dict per variant means the same PMID cited by two DIFFERENT
     variants (not uncommon - one functional-assay paper often covers
     several missense variants) still gets re-fetched from PubMed MCP once
     per variant.
  2. An in-memory dict doesn't survive past one run of the script, so
     re-running it (e.g. after adding more ground-truth variants, or after
     a crash) re-fetches everything from scratch.

A published paper's full text never changes, and "this PMID has no PMCID"
is an equally permanent fact - so entries here have no expiry. One JSON
file per PMID (cache/pubmed_fulltext/<pmid>.json by convention) rather
than one big file, so a single write only touches one small file and a
crash mid-run can't corrupt entries already on disk.

fetch_full_text()'s `cache` parameter only ever uses the dict protocol
(`in`, `[...]`, `[...] = ...`), so this class is a drop-in replacement for
a plain dict - pass ONE instance for the whole run (not one per variant)
to get cross-variant AND cross-run sharing.

[Size cap]
  Real scale here is tiny (order of 100 PMIDs), so no eviction would ever
  be needed for correctness or disk space. `max_files`/`evict_batch` exist
  purely as a safety net against unbounded growth if this cache is reused
  for a much larger variant set later: once the cache exceeds `max_files`
  entries, the oldest `evict_batch` (by file mtime - a fetched paper is
  only ever written once, so mtime already reflects fetch order) are
  deleted. Eviction is checked only on write, never on read, so lookups
  stay cheap.
"""

from __future__ import annotations

import json
from pathlib import Path


class DiskBackedFullTextCache:
    def __init__(self, directory: Path | str, max_files: int = 5000, evict_batch: int = 100):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        if max_files <= 0:
            raise ValueError("max_files must be positive")
        if evict_batch <= 0:
            raise ValueError("evict_batch must be positive")
        self.max_files = max_files
        self.evict_batch = evict_batch

    def _path(self, pmid: str) -> Path:
        return self.directory / f"{pmid}.json"

    def __contains__(self, pmid: str) -> bool:
        return self._path(pmid).exists()

    def __getitem__(self, pmid: str) -> tuple[str | None, str]:
        data = json.loads(self._path(pmid).read_text(encoding="utf-8"))
        return data["full_text"], data["note"]

    def __setitem__(self, pmid: str, value: tuple[str | None, str]) -> None:
        full_text, note = value
        self._path(pmid).write_text(
            json.dumps({"full_text": full_text, "note": note}, ensure_ascii=False),
            encoding="utf-8",
        )
        self._evict_if_over_capacity()

    def _evict_if_over_capacity(self) -> None:
        files = list(self.directory.glob("*.json"))
        if len(files) <= self.max_files:
            return
        files.sort(key=lambda p: p.stat().st_mtime)  # oldest first
        for stale in files[: self.evict_batch]:
            stale.unlink(missing_ok=True)

    def __len__(self) -> int:
        return len(list(self.directory.glob("*.json")))
