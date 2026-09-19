"""
llm_client.py

One place that decides which LLM backend every module in this project talks
to, so switching providers means changing environment variables, not code.

[Why an OpenAI-compatible client works for both]
  Every call site in this project (acmg_pipeline.pipeline, .clinical_
  extraction, .hpo_mondo_extraction, and everything downstream of them -
  ps3_bs3.py, pp4_literature_search.py, pp1_segregation_search.py,
  condition_from_literature_search.py) only ever does the same simple thing:
  `client.chat.completions.create(model=..., messages=[{"role": "system",
  ...}, {"role": "user", ...}])` and reads `response.choices[0].message.
  content`. Anthropic's Claude API offers an OpenAI SDK compatibility
  endpoint that accepts exactly this shape, so switching from the
  self-hosted vLLM server to Claude needs no new SDK and no call-site
  changes - only which base_url/api_key/model this one factory hands out.

[Why this project would want to switch at all]
  Empirically (2026-09-19), every LLM call in this project's logs takes
  9-16 seconds even for a short prompt (~1200 input tokens, ~250 output
  tokens) - a floor that does not track prompt size, pointing at the
  self-hosted vLLM server itself (not PubMed fetching, which is either an
  instant cache hit or, when live, has never been the slow step in any log
  this project has produced) as the bottleneck. Claude's API is a
  candidate for a faster, more reliable backend for the exact same prompts.

[What still needs to be true for LLM_PROVIDER=claude to work]
  A real ANTHROPIC_API_KEY the user controls, added to .env by the user
  themselves (never pasted into a chat/agent session) - this module raises
  a clear error rather than silently falling back when it is missing so a
  misconfigured run fails loudly instead of quietly hitting the wrong
  backend. Not yet verified end-to-end against a real ANTHROPIC_API_KEY in
  this environment - the vLLM path (the default, unchanged) has been
  exercised all hackathon; test with a real key before relying on this for
  anything that matters.
"""

from __future__ import annotations

import os
from pathlib import Path

from openai import OpenAI

ROOT_DIR = Path(__file__).resolve().parent.parent


def _load_env_file(path: Path = ROOT_DIR / ".env") -> None:
    if not path.is_file():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip())


_load_env_file()

DEFAULT_VLLM_MODEL = "google/gemma-4-26B-A4B-it"
DEFAULT_CLAUDE_MODEL = "claude-sonnet-5"
DEFAULT_CLAUDE_BASE_URL = "https://api.anthropic.com/v1/"


def get_provider() -> str:
    """'vllm' (default, this project's original self-hosted backend) or
    'claude' - set via LLM_PROVIDER in the environment or .env."""
    return os.environ.get("LLM_PROVIDER", "vllm").strip().lower()


def make_client() -> tuple[OpenAI, str]:
    """Returns (client, model) for whichever provider is configured.

    Every caller already does its own `client.chat.completions.create(
    model=model, ...)` - this only decides which server that goes to.
    """
    provider = get_provider()
    if provider in ("claude", "anthropic"):
        api_key = os.environ.get("ANTHROPIC_API_KEY", "")
        if not api_key:
            raise RuntimeError(
                "LLM_PROVIDER=claude/anthropic requires ANTHROPIC_API_KEY to be set "
                "(add it to .env yourself - never paste an API key into a chat/agent session)."
            )
        base_url = os.environ.get("ANTHROPIC_BASE_URL_OPENAI_COMPAT", DEFAULT_CLAUDE_BASE_URL)
        model = os.environ.get("ANTHROPIC_MODEL", DEFAULT_CLAUDE_MODEL)
        return OpenAI(base_url=base_url, api_key=api_key), model

    base_url = os.environ.get("VLLM_BASE_URL", "")
    api_key = os.environ.get("VLLM_API_KEY", "")
    if not base_url or not api_key:
        raise RuntimeError(
            "VLLM_BASE_URL / VLLM_API_KEY が設定されていません。"
            "リポジトリ直下に .env を作成してください(.env.example を参照)。"
        )
    model = os.environ.get("VLLM_MODEL", DEFAULT_VLLM_MODEL)
    return OpenAI(base_url=base_url, api_key=api_key), model
