"""
llm_cache.py

Disk-backed cache for LLM judgment calls (acmg_pipeline.pipeline.
call_llm_for_judgment()), modeled on fulltext_cache.DiskBackedFullTextCache -
same dict-protocol shape (__contains__/__getitem__/__setitem__), same
one-file-per-key persistence, same size-cap eviction.

[Why this exists, 2026-09-17]
  DiskBackedFullTextCache already avoids re-fetching a paper's full text
  across runs and variants, but the LLM inference call itself (the actual
  slow, costly step - 8-35s per paper observed in practice, dominated by
  the model, not the network) was never cached: two runs asking the exact
  same question (same prompt) always paid for a fresh completion. Added
  while planning a 64-variant integrated validation run
  (run_integrated_validation_64.py) where many papers and variants repeat
  across runs (re-running after a code fix, or a --limit smoke test
  followed by the full run) - without this, every one of those repeats
  pays full LLM latency again for an answer already computed.

[Why keying on (model, prompt) is safe here specifically]
  The LLM endpoint's temperature is left at its default (not 0 - see
  pipeline.call_llm_for_judgment()), so two calls with an identical prompt
  are not guaranteed to return byte-identical completions in general. This
  cache deliberately freezes ONE sampled answer per (model, prompt) pair,
  trading a small amount of sampling variance for reproducibility and
  speed across repeated validation runs - the same tradeoff this project's
  own design already makes elsewhere (e.g. run_validation_64.py's
  full-text cache), and appropriate here because this cache is opt-in
  (callers must construct and pass one explicitly - see cache=None default
  on call_llm_for_judgment()) and is meant for validation/development
  reruns, not for treating a single cached sample as more authoritative
  than a fresh call would be in production curator-facing use.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


class DiskBackedLLMCache:
    def __init__(self, directory: str | Path, max_files: int = 5000, evict_batch: int = 100):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.max_files = max_files
        self.evict_batch = evict_batch

    @staticmethod
    def _key(model: str, prompt: str) -> str:
        return hashlib.sha256(f"{model}\x00{prompt}".encode("utf-8")).hexdigest()

    def _path(self, model: str, prompt: str) -> Path:
        return self.directory / f"{self._key(model, prompt)}.json"

    def __contains__(self, item: tuple[str, str]) -> bool:
        model, prompt = item
        return self._path(model, prompt).is_file()

    def __getitem__(self, item: tuple[str, str]) -> tuple[dict, dict]:
        model, prompt = item
        path = self._path(model, prompt)
        path.touch()  # refresh mtime for the LRU-by-mtime eviction below
        data = json.loads(path.read_text(encoding="utf-8"))
        return data["response"], data["stats"]

    def __setitem__(self, item: tuple[str, str], value: tuple[dict, dict]) -> None:
        model, prompt = item
        response, stats = value
        path = self._path(model, prompt)
        path.write_text(
            json.dumps({"response": response, "stats": stats}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        self._evict_if_over_capacity()

    def __len__(self) -> int:
        return sum(1 for _ in self.directory.glob("*.json"))

    def _evict_if_over_capacity(self) -> None:
        files = list(self.directory.glob("*.json"))
        if len(files) <= self.max_files:
            return
        files.sort(key=lambda p: p.stat().st_mtime)
        for path in files[: self.evict_batch]:
            path.unlink(missing_ok=True)
