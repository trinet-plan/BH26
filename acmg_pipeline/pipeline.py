"""
ps3_bs3_llm_pipeline.py
Integration script that connects the PS3/BS3 judgment pipeline to a real LLM
(gemma-4 on vLLM).

Reuses the connection pattern from mcp_sample_multi.py (PubMed MCP + an
OpenAI-compatible client) as-is, and wires together the following:

  Given a PMID
    -> fetch full text via PubMed MCP (convert_article_ids -> get_full_text_article)
    -> build the prompt via ps3_bs3_judgment.build_prompt()
    -> send it to gemma-4 (vLLM) and get structured JSON output back
    -> parse it with PS3BS3Judgment.from_json()
    -> run ps3_bs3_judgment.finalize() to attach the approved-assay
       cross-check and human-curator hints
    -> display and log the result, and compare against ground truth if known

[Important note]
  This script assumes reachability to the vLLM server (113.43.212.98:8000)
  and the PubMed MCP server (pubmed.mcp.claude.com). Neither is reachable
  from the development sandbox, so actual runtime verification has to happen
  in the user's own environment (see the design doc appendix on known
  network-allowlist constraints).

[Requirements]
  pip install mcp openai
"""

import asyncio
import json
import os
import re
import sys
import time
from contextlib import AsyncExitStack
from datetime import datetime
from functools import lru_cache
from pathlib import Path

from mcp import ClientSession
try:
    from mcp.client.streamable_http import streamable_http_client
except ImportError:
    # Depending on the mcp version, the function may be named
    # streamablehttp_client (no underscore). An ImportError was actually
    # observed in the user's environment (mcp.client.streamable_http), so
    # both spellings are supported here.
    from mcp.client.streamable_http import streamablehttp_client as streamable_http_client
from mcp import types as mcp_types

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import types

from acmg_pipeline.criteria import ps3_bs3 as ps3bs3_mod
from acmg_pipeline.criteria import ps4 as ps4_mod
from acmg_pipeline.criteria import segregation as seg_mod
from acmg_pipeline.criteria import stubs as stubs_mod
from acmg_pipeline.common import (
    PaperContribution, AggregatedJudgment, FinalResult, CuratorHint,
    MatchStatus, VariantMatchingResult,
)
from acmg_pipeline.gate import ERepoClient, _protein_equivalents
from acmg_pipeline.export import (
    build_automated_evidence_line,
    build_evidence_line,
    build_stub_evidence_line,
    build_workflow_evidence_line,
)
from acmg_pipeline.classification import (
    ALL_ACMG_CODES,
    AUTOMATED_CODES,
    IMPLEMENTED_CODES,
    LITERATURE_CODES,
    PHENOTYPE_SEGREGATION_CODES,
    ClassificationResult,
    classify,
    from_aggregated_judgment,
)
from acmg_pipeline.criteria import stubs
from acmg_pipeline.criteria import pp1_bs4_pp4_engine
from acmg_pipeline.api_input import ApiCaseInput
from acmg_pipeline.clinical_note import ClinicalNoteExtraction, extract_clinical_note
from acmg_pipeline.inputs import empty_clinical_note
from acmg_pipeline.vcf_record import VariantRecord
from acmg_pipeline.automated_core.identity import reconcile
from acmg_pipeline.automated_core.models import CRITERIA as AUTOMATED_CRITERIA
from acmg_pipeline.automated_core.models import Variant as AutomatedVariant
from acmg_pipeline.gene_disease import build_assessment_document
from acmg_pipeline.services.resolve import PROVIDER_ERRORS, ProviderEvidenceResolver
from acmg_pipeline.automated_engine import evaluate_record as evaluate_automated_record, make_services

VA_SPEC_OUTPUT_DIR = Path("va_spec_output")
VA_SPEC_OUTPUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Judgment engines: one per criterion family (PS3/BS3, PS4, PP1/BS4)
# ---------------------------------------------------------------------------
#
# Added 2026-09-15 alongside ps4_judgment.py / segregation_judgment.py. Each
# criterion family has its own extraction schema, prompt, and finalize()
# safety nets, but judge_single_paper()/judge_variant() below don't need to
# know any of that - they just need, for whichever family is running: how
# to build a prompt, how to parse the LLM's JSON, how to finalize one
# paper's judgment, and how to aggregate multiple papers. A JudgmentEngine
# bundles exactly those four functions so the same judge_single_paper()/
# judge_variant() work for all three families unmodified.
#
# finalize() is wrapped in a same-shaped lambda for all three engines, even
# though ps4/segregation's own finalize() only takes (judgment, pmid=...)
# and ignores gene/vcep_name/criterion - PS3/BS3's finalize() needs all
# three (for the approved-assay cross-check), so the call site
# (judge_single_paper) always passes the full set and lets each engine's
# adapter drop what it doesn't need.

class JudgmentEngine:
    def __init__(self, name, build_prompt, from_json, finalize, aggregate, not_clear):
        self.name = name
        self.build_prompt = build_prompt
        self.from_json = from_json
        self.finalize = finalize
        self.aggregate = aggregate
        # This engine's own "no confident direction" Enum member (e.g.
        # OverallDirection.NOT_CLEAR) - needed by judge_single_paper() to
        # build a placeholder contribution when full text can't be fetched
        # at all (see the full_text-unavailable branch below), so that
        # PMID's "why it was never evaluated" reason still shows up in the
        # aggregated EvidenceLine's curator hints instead of the PMID
        # silently vanishing.
        self.not_clear = not_clear


PS3BS3_ENGINE = JudgmentEngine(
    name="PS3/BS3",
    build_prompt=ps3bs3_mod.build_prompt,
    from_json=ps3bs3_mod.PS3BS3Judgment.from_json,
    finalize=lambda judgment, gene, vcep_name, criterion, pmid: ps3bs3_mod.finalize(
        judgment, gene=gene, vcep_name=vcep_name, criterion=criterion, pmid=pmid,
    ),
    aggregate=ps3bs3_mod.aggregate_multi_paper_results,
    not_clear=ps3bs3_mod.OverallDirection.NOT_CLEAR,
)

PS4_ENGINE = JudgmentEngine(
    name="PS4",
    build_prompt=ps4_mod.build_prompt,
    from_json=ps4_mod.PS4Judgment.from_json,
    finalize=lambda judgment, gene, vcep_name, criterion, pmid: ps4_mod.finalize(judgment, pmid=pmid),
    aggregate=ps4_mod.aggregate_multi_paper_results,
    not_clear=ps4_mod.CaseControlDirection.NOT_CLEAR,
)

SEGREGATION_ENGINE = JudgmentEngine(
    name="PP1/BS4",
    build_prompt=seg_mod.build_prompt,
    from_json=seg_mod.SegregationJudgment.from_json,
    finalize=lambda judgment, gene, vcep_name, criterion, pmid: seg_mod.finalize(judgment, pmid=pmid),
    aggregate=seg_mod.aggregate_multi_paper_results,
    not_clear=seg_mod.SegregationDirection.NOT_CLEAR,
)

ENGINE_BY_CRITERION = {
    "PS3": PS3BS3_ENGINE, "BS3": PS3BS3_ENGINE,
    "PS4": PS4_ENGINE,
    "PP1": SEGREGATION_ENGINE, "BS4": SEGREGATION_ENGINE,
}

# ---------------------------------------------------------------------------
# Configuration (reused from mcp_sample_multi.py)
# ---------------------------------------------------------------------------
ROOT_DIR = Path(__file__).resolve().parent.parent


def _load_env_file(path: Path = ROOT_DIR / ".env") -> None:
    """.env があれば os.environ に読み込む(python-dotenv 不使用の簡易実装)。"""
    if not path.exists():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip())


_load_env_file()

from acmg_pipeline.llm_client import make_client as _make_llm_client

MCP_SERVERS = {
    "pubmed": "https://pubmed.mcp.claude.com/mcp",
}

# Which LLM backend every module downstream of this one talks to (PS3/BS3/
# PS4, PP4/PP1 literature search, MONDO resolution, ...) is decided in one
# place now - acmg_pipeline.llm_client - via LLM_PROVIDER in .env. Default
# stays this project's original self-hosted vLLM server; set
# LLM_PROVIDER=claude (plus ANTHROPIC_API_KEY) to switch to Claude instead,
# see that module's own docstring.
client, MODEL = _make_llm_client()

LOG_DIR = Path("logs")
LOG_DIR.mkdir(exist_ok=True)
LOG_PATH = LOG_DIR / f"ps3bs3_run_{datetime.now():%Y%m%d_%H%M%S}.log"


def log(text: str = ""):
    with LOG_PATH.open("a", encoding="utf-8") as f:
        f.write(text + "\n")


def show(text: str):
    print(text)
    log(text)


# ---------------------------------------------------------------------------
# A safe wrapper around tool calls
# ---------------------------------------------------------------------------
#
# In mcp==2.2.0, ClientSession.call_tool() internally calls send_request()
# and then re-validates the structured output via validate_tool_result(),
# which raises a RuntimeError if the tool declares an outputSchema but
# result.structured_content is None (confirmed directly against the source
# of mcp/client/session.py). This happens when the PubMed MCP server does
# not return structuredContent.
#
# Calling send_request() directly bypasses this re-validation
# (validate_tool_result), since it is never invoked on that path. This
# construction (how to build CallToolRequest / CallToolRequestParams) was
# confirmed directly against the mcp==2.2.0 source. It has not, however,
# been tested against the real PubMed MCP server (unreachable from this
# sandbox).

async def call_tool_safe(mcp: ClientSession, name: str, arguments: dict):
    try:
        return await mcp.call_tool(name, arguments)
    except RuntimeError as e:
        if "did not return structured content" not in str(e):
            raise
        # Confirmed in mcp==2.2.0: call_tool() calls send_request()
        # internally and then re-validates the structured output in
        # validate_tool_result(), which is where this RuntimeError is
        # raised. Calling send_request() directly skips that re-validation
        # (no need to wrap in ClientRequest; the method field is set
        # automatically inside CallToolRequest).
        return await mcp.send_request(
            mcp_types.CallToolRequest(
                params=mcp_types.CallToolRequestParams(name=name, arguments=arguments),
            ),
            mcp_types.CallToolResult,
        )


# ---------------------------------------------------------------------------
# Fetching full text via PubMed MCP
# ---------------------------------------------------------------------------

async def fetch_full_text(
    mcp: ClientSession, pmid: str, cache: dict[str, tuple[str | None, str]] | None = None,
) -> tuple[str | None, str]:
    """
    Fetches the full text for a PMID. Returns None with an explanatory note
    if the article is not in PMC. As noted in design doc section 2-3, only
    about 20% of articles are in PMC, so this branch is mandatory.

    `cache`, when passed, is checked first and populated on the way out -
    added 2026-09-16 so that testing a single variant against all 5
    implemented criteria (test_case_ground_truth.md's 227-pair validation
    run) doesn't re-fetch the same PMID's full text from PubMed MCP once
    per criterion. The caller owns the cache's lifetime: pass a fresh dict
    per variant (not a single dict for the whole run) so a variant's PMIDs
    don't linger in memory or get reused for an unrelated variant that
    happens to cite the same paper.
    """
    if cache is not None and pmid in cache:
        cached_text, cached_note = cache[pmid]
        return cached_text, f"{cached_note} [full_text_cache hit - no PubMed MCP call made]"

    result = await _fetch_full_text_uncached(mcp, pmid)
    if cache is not None:
        cache[pmid] = result
    return result


async def _fetch_full_text_uncached(mcp: ClientSession, pmid: str) -> tuple[str | None, str]:
    unavailable_reason = None

    convert_result = await call_tool_safe(mcp, "convert_article_ids", {"ids": [pmid], "id_type": "pmid"})
    convert_text = "\n".join(getattr(b, "text", str(b)) for b in convert_result.content)
    try:
        convert_data = json.loads(convert_text)
        pmcid = convert_data["records"][0].get("pmcid")
    except (json.JSONDecodeError, KeyError, IndexError):
        pmcid = None

    if not pmcid:
        unavailable_reason = f"PMID:{pmid} is not in PMC (no PMCID)."
    else:
        ft_result = await call_tool_safe(mcp, "get_full_text_article", {"pmc_ids": [pmcid]})
        ft_text = "\n".join(getattr(b, "text", str(b)) for b in ft_result.content)
        try:
            ft_data = json.loads(ft_text)
            article = ft_data["articles"][0]
            full_text = article.get("full_text", "")
            doi = article.get("doi", "")
            if full_text:
                return full_text, f"Fetched from PubMed. PMCID={pmcid}, DOI={doi}"
            # Real bug found 2026-09-15 running the democase variants: the
            # PubMed MCP server can return a 200-ish "articles" record for a
            # PMCID with an empty/missing full_text field (e.g. PMID:8282798,
            # PMID:19645038 - both have a PMCID but no body text came back).
            # The caller only ever checked `full_text is None`, so an empty
            # string silently passed as "fetched successfully" and got
            # substituted into the prompt as a blank "Full text of the
            # paper:" section - the LLM then correctly complained it had no
            # text to work with ("Please provide the full text..."), which
            # broke JSON parsing and looked like a flaky LLM failure rather
            # than the real cause (no text was ever sent). Treat this the
            # same as "not in PMC" so it falls through to the abstract
            # fallback below instead of being prompted with blank text.
            unavailable_reason = f"PMID:{pmid} has a PMCID ({pmcid}) but no full_text came back from PubMed MCP (empty article body)."
        except (json.JSONDecodeError, KeyError, IndexError) as e:
            unavailable_reason = f"Failed to parse the full-text fetch result for PMID:{pmid}: {e}"

    # Abstract fallback, added 2026-09-16: PMC full text is only available
    # for ~20% of PMIDs (design doc section 2-3) - before this fallback
    # existed, the other ~80% contributed literally nothing, even when
    # get_article_metadata's abstract directly discusses the target variant
    # (confirmed against real data: PMID:17351073, cited for MYH7 c.1594T>C's
    # PP1 evidence, was skipped as "not in PMC" in every run to date, but its
    # abstract explicitly reports the S532P mutant's force-generation data -
    # exactly the kind of PS3-relevant finding this pipeline was missing).
    # The `[ABSTRACT ONLY ...]` marker is prepended to the text itself
    # (rather than added as a separate return value) so every build_prompt()
    # caller sees the caveat with no signature change anywhere downstream -
    # a real fingerprint of what evidence quality this judgment rests on.
    abstract, abstract_note = await _fetch_abstract(mcp, pmid)
    if abstract:
        marked_text = (
            "[ABSTRACT ONLY - full text was not available for this paper; the excerpt "
            "below is the PubMed abstract, not the full article. Numeric details, "
            "specific experiment counts, and per-family/per-patient data that would "
            "normally only appear in the full text may be absent. Judge accordingly - "
            "prefer not_clear over inferring specifics the abstract doesn't state.]\n\n"
            + abstract
        )
        return marked_text, f"{unavailable_reason} {abstract_note}"

    return None, f"{unavailable_reason} Abstract fallback also unavailable ({abstract_note}). Judgment marked insufficient_data."


async def _fetch_abstract(mcp: ClientSession, pmid: str) -> tuple[str | None, str]:
    """
    Fetches just the title/abstract via get_article_metadata() - used by
    _fetch_full_text_uncached() when PMC full text isn't available. Returns
    (abstract, note); abstract is None (with an explanatory note) if
    get_article_metadata has no abstract for this PMID either (e.g. a very
    old record, or a non-journal source).
    """
    try:
        result = await call_tool_safe(mcp, "get_article_metadata", {"pmids": [pmid]})
        text = "\n".join(getattr(b, "text", str(b)) for b in result.content)
        data = json.loads(text)
        article = data["articles"][0]
        abstract = article.get("abstract")
    except (json.JSONDecodeError, KeyError, IndexError) as e:
        return None, f"get_article_metadata failed to parse for PMID:{pmid}: {e}"

    if not abstract:
        return None, f"get_article_metadata returned no abstract for PMID:{pmid}."
    return abstract, "Using the PubMed abstract instead (via get_article_metadata)."


# ---------------------------------------------------------------------------
# LLM call + JSON extraction
# ---------------------------------------------------------------------------

def extract_json(raw_text: str) -> dict:
    """
    Extracts a JSON object from the LLM's response, handling both the case
    where it is wrapped in a code fence (```json ... ```) and the case where
    it has explanatory text before/after it.
    """
    fence_match = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", raw_text, re.DOTALL)
    if fence_match:
        candidate = fence_match.group(1)
    else:
        # naively take everything from the first '{' to the matching last '}'
        start = raw_text.find("{")
        end = raw_text.rfind("}")
        if start == -1 or end == -1 or end <= start:
            raise ValueError(f"Could not extract JSON from the response: {raw_text[:300]!r}")
        candidate = raw_text[start:end + 1]
    return json.loads(candidate)


def call_llm_for_judgment(prompt: str, cache=None) -> tuple[dict, dict]:
    """
    Sends the prompt to gemma-4 and gets structured JSON output back.
    Returns (parsed JSON, raw API response stats).

    `cache`, when passed (an acmg_pipeline.llm_cache.DiskBackedLLMCache or
    anything with the same dict protocol), is checked first keyed on
    (MODEL, prompt) - see that module's own docstring for why freezing one
    sampled answer per prompt is an acceptable tradeoff here (opt-in,
    validation/development reruns, not production curator-facing use). A
    cache hit returns instantly with elapsed_sec=0.0 and a cache_hit=True
    marker in stats so callers/logs can tell it apart from a fresh call.
    """
    if cache is not None and (MODEL, prompt) in cache:
        parsed, stats = cache[(MODEL, prompt)]
        return parsed, {**stats, "elapsed_sec": 0.0, "cache_hit": True}
    t0 = time.perf_counter()
    res = client.chat.completions.create(
        model=MODEL,
        messages=[
            {
                "role": "system",
                "content": (
                    "You are an assistant supporting a clinical genetics "
                    "curator. Output only the JSON as instructed, with no "
                    "other explanatory text."
                ),
            },
            {"role": "user", "content": prompt},
        ],
    )
    dt = time.perf_counter() - t0
    msg = res.choices[0].message
    stats = {
        "elapsed_sec": dt,
        "prompt_tokens": res.usage.prompt_tokens,
        "completion_tokens": res.usage.completion_tokens,
        "cache_hit": False,
    }
    parsed = extract_json(msg.content or "")
    if cache is not None:
        cache[(MODEL, prompt)] = (parsed, stats)
    return parsed, stats


# ---------------------------------------------------------------------------
# Main function tying the whole flow together
# ---------------------------------------------------------------------------
#
# Restructured 2026-09-15 to support multiple papers per variant, in
# response to a design gap the user pointed out: the original version only
# ever looked at one PMID per variant, whereas a real ClinGen classification
# typically cites several papers for a single criterion (e.g., 4 PMIDs for
# RUNX1 c.601C>T, seen during design doc section 9). judge_single_paper()
# now handles exactly one paper and returns a PaperContribution (or None if
# that paper could not be used at all); judge_variant() calls it once per
# PMID and then runs aggregate_multi_paper_results() over the results.

async def judge_single_paper(
    engine: JudgmentEngine,
    mcp: ClientSession,
    pmid: str,
    gene: str,
    hgvsc: str,
    hgvsp: str,
    equivalents: list[str],
    vcep_name: str | None,
    criterion: str,
    full_text_cache: dict[str, tuple[str | None, str]] | None = None,
    llm_cache=None,
) -> PaperContribution | None:
    show(f"\n--- PMID:{pmid} ---")

    full_text, note = await fetch_full_text(mcp, pmid, cache=full_text_cache)
    show(f"[Full text] {note}")
    if not full_text:
        show(f"  -> Full text unavailable; PMID:{pmid} contributes no evidence (skipping the LLM call)")
        # Still return a PaperContribution (not None) so this PMID and the
        # reason it could never be evaluated survive into the aggregation
        # and the VA-Spec output's curatorHints, instead of silently
        # vanishing. match_status=UNSUCCESSFUL keeps it correctly excluded
        # from the "relevant" (usable-evidence) set in aggregate_multi_
        # paper_results, so it never counts toward direction agreement or
        # the paper-count strength heuristic - it only adds a record of
        # "this paper was cited but could not be read."
        placeholder_judgment = types.SimpleNamespace(
            variant_matching=VariantMatchingResult(match_status=MatchStatus.UNSUCCESSFUL, notes=note),
            overall_evidence=types.SimpleNamespace(direction=engine.not_clear, rationale=note),
        )
        result = FinalResult(
            judgment=placeholder_judgment,
            curator_hints=[CuratorHint("info", f"Full text could not be obtained for PMID:{pmid}: {note}")],
            effective_direction=engine.not_clear,
        )
        return PaperContribution(pmid=pmid, result=result)

    # engine.build_prompt() takes (variant, clinical_note, full_text) as of
    # h.muroda's 2026-09-16 "unify criterion input interfaces" refactor
    # (ps3_bs3.py/ps4.py/segregation.py all match this now) - this function
    # only has the already-flattened gene/hgvsc/hgvsp/equivalents (its own
    # signature predates that refactor and several callers, e.g.
    # run_validation_64.py, only ever have those flat values, not a real
    # VariantRecord/ClinicalNoteExtraction), so a minimal VariantRecord is
    # reconstructed here rather than threading a real one through every
    # caller. EQUIVALENTS is stored in INFO because variant_identity()
    # (which every build_prompt() calls) reads it from there - without it,
    # protein-equivalent notation matching within a paper would be lost.
    variant = VariantRecord(
        chrom="", pos=1, id="", ref="", alt="", qual="", filter="",
        info={"GENE": gene, "HGVSC": hgvsc, "HGVSP": hgvsp, "EQUIVALENTS": equivalents},
    )
    prompt = engine.build_prompt(variant, empty_clinical_note(), full_text)
    log(f"--- Prompt (first 1000 chars) ---\n{prompt[:1000]}")

    try:
        raw_json, stats = call_llm_for_judgment(prompt, cache=llm_cache)
    except Exception as e:
        show(f"  -> Error during the LLM call / JSON parsing: {e}")
        return None

    cache_note = " [llm_cache hit - no LLM call made]" if stats.get("cache_hit") else ""
    show(f"[LLM response] {stats['elapsed_sec']:.1f}s, "
         f"prompt={stats['prompt_tokens']}tok, completion={stats['completion_tokens']}tok{cache_note}")
    log(f"--- Raw LLM JSON ---\n{json.dumps(raw_json, ensure_ascii=False, indent=2)}")

    try:
        judgment = engine.from_json(raw_json)
    except Exception as e:
        show(f"  -> Failed to parse the structured output: {e} (the LLM may not have followed the schema)")
        show(f"     Raw JSON (first 500 chars): {json.dumps(raw_json, ensure_ascii=False)[:500]}")
        return None

    result = engine.finalize(judgment, gene=gene, vcep_name=vcep_name, criterion=criterion, pmid=pmid)

    show(f"LLM's raw judgment: direction = {judgment.overall_evidence.direction.value}")
    show(f"Rationale: {judgment.overall_evidence.rationale}")
    if result.effective_direction != judgment.overall_evidence.direction:
        show(f"-> Forced override due to conflict detection: PMID:{pmid}'s judgment = {result.effective_direction.value}")
    for hint in result.curator_hints:
        show(f"  [{hint.severity:7s}] {hint.message}")

    return PaperContribution(pmid=pmid, result=result)


async def judge_variant(
    engine: JudgmentEngine,
    mcp: ClientSession,
    pmids: list[str],
    gene: str,
    hgvsc: str,
    hgvsp: str,
    equivalents: list[str],
    vcep_name: str | None,
    criterion: str,
    ground_truth: str | None = None,
    full_text_cache: dict[str, tuple[str | None, str]] | None = None,
    llm_cache=None,
) -> AggregatedJudgment:
    """
    Judges a variant using ALL of the given PMIDs (not just the first one),
    then aggregates the per-paper results into a single variant-level
    direction via `engine.aggregate`. Pass a single-element list for the
    common case of only having one known citing paper.

    `engine` selects which criterion family is being run (PS3/BS3, PS4, or
    PP1/BS4 - see ENGINE_BY_CRITERION) - everything else about this
    function is identical across all three.

    `full_text_cache`, when passed, is forwarded to fetch_full_text() via
    judge_single_paper() - added 2026-09-16 for the 227-pair validation run
    (test_case_ground_truth.md) where the same variant is judged once per
    implemented criterion it has ground truth for (up to 5x), each call
    citing the same PMIDs. Any object satisfying the dict protocol (`in`,
    `[...]`, `[...] = ...`) works here - a plain dict for simple ad hoc
    reuse across a handful of judge_variant() calls, or (preferred for a
    real batch run) acmg_pipeline.fulltext_cache.DiskBackedFullTextCache
    passed as ONE instance for the WHOLE run: since a paper's full text
    never changes, sharing it across variants (not just within one) and
    across separate runs of the script is strictly better, not a risk -
    see that module's docstring for why a per-variant-only cache still
    left cross-variant and cross-run reuse on the table.

    Returns the AggregatedJudgment itself (not just its direction), so the
    caller can both score it against ground_truth (see score_direction())
    and build a VA-Spec EvidenceLine from it (see va_spec_export.py) -
    the aggregate function already handles an empty contributions list
    correctly (falls through to not_clear with an explanatory hint), so
    there is no separate early-return branch here for "no paper yielded
    usable output".
    """
    show(f"\n{'='*70}\n[{engine.name}] {gene} {hgvsc} ({hgvsp}) - {len(pmids)} paper(s): {', '.join(pmids)}\n{'='*70}")

    contributions = []
    for pmid in pmids:
        contribution = await judge_single_paper(
            engine, mcp, pmid, gene, hgvsc, hgvsp, equivalents, vcep_name, criterion,
            full_text_cache=full_text_cache, llm_cache=llm_cache,
        )
        if contribution is not None:
            contributions.append(contribution)

    aggregated = engine.aggregate(contributions)
    show(f"\n[Aggregated across {len(pmids)} paper(s), {len(contributions)} usable]: "
         f"direction = {aggregated.aggregated_direction.value}")
    for hint in aggregated.aggregation_hints:
        show(f"  [{hint.severity:7s}] {hint.message}")

    if ground_truth is not None:
        show(f"\n[Scoring] Ground truth: {ground_truth} / Adopted judgment: {aggregated.aggregated_direction.value}")

    return aggregated


# ---------------------------------------------------------------------------
# Driving the literature criteria from structured API-shaped input
# ---------------------------------------------------------------------------
#
# Added 2026-09-16: the first wiring between acmg_pipeline.api_input.
# ApiCaseInput (this project's real input contract, built from a plain dict
# - e.g. an HTTP request body already run through json.loads(), never a
# file - see api_input.py's own from_dict()/from_json_file() split) and the
# existing literature JudgmentEngines. Takes VariantRecord directly (via
# case_input.parse_vcf().record) rather than a separate wrapper class, per
# the user's direction (2026-09-16): a criterion's input is one of the
# pulled-in entity classes (VariantRecord / ClinicalNoteExtraction), not a
# new parallel dataclass, so a later new VCF INFO key or clinical_note
# field needs no call-site signature change here.
#
# Only drives PS3/BS3/PS4/PP1/BS4 (ENGINE_BY_CRITERION) - the literature
# path. It does NOT touch case_input.clinical_note at all: none of these 5
# codes' literature path needs it (see test_case_ground_truth.md's
# criterion-input-group design). The clinical-note-only paths (PP1/BS4's
# alternative "patient's own family" route via acmg_pipeline.criteria.
# segregation.from_clinical_note_family(), and PS2/PM6's de-novo rule) are
# separate, not-yet-wired call sites - this function covers literature
# only.
#
# PMIDs come from ERepoClient first; when ERepo has zero citations for this
# variant (e.g. it's genuinely novel/unregistered - the gap raised by the
# user 2026-09-16, "PMIDありきの実装を変更しないといけない"), fall back to
# search_candidate_pmids() below, a live PubMed search. See that function's
# docstring for the query strategy and its real-data validation.

async def search_candidate_pmids(
    mcp: ClientSession, gene: str, hgvsp: str | None = None, disease: str | None = None,
    max_results: int = 10,
) -> list[str]:
    """
    Live PubMed search fallback for when ERepo has no curated citations at
    all for a variant. Tries several query strategies, in order, stopping
    once max_results candidates are collected - validated 2026-09-16
    against real cases in this project's own demo data:

      - "<gene> AND <protein change>" (e.g. "MYH7 AND Arg719Trp") works
        well for missense variants - 5/5 real hits for a known MYH7
        variant, including PMIDs already in this project's own ERepo-
        sourced ground truth.
      - Exact HGVS c./p. notation (e.g. "MYBPC3 AND c.278delA",
        "MYBPC3 AND Lys93Argfs") returns ZERO results even for a variant
        whose own source paper is indexed in PubMed - confirmed
        empirically, not assumed. PubMed's indexing does not reliably
        match this notation as literal text.
      - "<gene>[Title] AND <disease, keywords>" works well as a fallback,
        especially for indel/frameshift variants where the protein-change
        query above doesn't apply well - found the exact real source paper
        for MYBPC3 c.278delA (democase Case 1) via "MYBPC3[Title] AND
        hypertrophic cardiomyopathy AND apical aneurysm" -> PMID 42428486,
        the correct paper.
      - A last-resort "<gene> AND novel variant" query, if the above
        didn't find enough.

    Results are CANDIDATES, not confirmed citations - unlike ERepo's
    evidenceLinks (a VCEP already verified these papers discuss the
    variant), a search hit is not guaranteed relevant until judge_single_
    paper()'s own variant-matching/not_clear handling checks it against
    the actual text. That's an acceptable tradeoff: the alternative is
    finding nothing at all for a variant ERepo has never curated.
    """
    queries: list[str] = []
    if hgvsp and hgvsp != "N/A":
        aa_change = hgvsp.strip().removeprefix("p.").strip("()")
        if aa_change:
            queries.append(f"{gene} AND {aa_change}")
    if disease:
        queries.append(f"{gene}[Title] AND {disease.replace('_', ' ')}")
    queries.append(f"{gene} AND novel variant")

    seen: list[str] = []
    for query in queries:
        if len(seen) >= max_results:
            break
        try:
            result = await call_tool_safe(mcp, "search_articles", {"query": query, "max_results": max_results})
            text = "\n".join(getattr(b, "text", str(b)) for b in result.content)
            data = json.loads(text)
        except Exception as e:
            show(f"  [PubMed search] query {query!r} failed: {e}")
            continue
        for pmid in data.get("pmids", []):
            if pmid not in seen:
                seen.append(pmid)
    return seen[:max_results]


async def resolve_pmids_for_variant(
    mcp: ClientSession, erepo_client: ERepoClient, gene: str, hgvsc: str,
    hgvsp: str | None = None, disease: str | None = None,
) -> tuple[list[str], str]:
    """
    The single, shared "how do we get PMIDs for this variant" resolver:
    ERepo's curated evidenceLinks first, a live PubMed search (search_
    candidate_pmids() above) only when ERepo has nothing. Returns
    (pmids, source) where source is "erepo" or "pubmed_search" - an empty
    pmids list means neither source found anything at all.

    Extracted 2026-09-16 so every caller shares one resolution path rather
    than each hand-rolling its own ERepo-only lookup: run_validation_64.py
    used to call erepo_client.lookup() directly and never got the search
    fallback judge_variant_from_structured_input() (below) already had,
    silently skipping the exact variants (e.g. MYBPC3 c.1000G>A, DSG2
    c.1592T>G - no ERepo record at all) that most needed it. Per the
    user's explicit direction (2026-09-16): no caller should reference
    ERepo directly anymore - always go through this function instead, so
    a future change to the resolution strategy (e.g. a better search query,
    or a third PMID source) only needs to happen in one place.
    """
    erepo_result = erepo_client.lookup(gene, hgvsc)
    pmids = erepo_result.evidence_pmids
    show(f"\n[ERepo] {gene} {hgvsc}: found_in_erepo={erepo_result.found_in_erepo}, evidence_pmids={pmids}")
    if pmids:
        return pmids, "erepo"

    show("  -> No evidence PMIDs from ERepo; falling back to a live PubMed search "
         "(this variant has no curated citations, e.g. it may be novel/unregistered)")
    pmids = await search_candidate_pmids(mcp, gene, hgvsp, disease)
    show(f"[PubMed search] candidate PMIDs: {pmids}")
    if pmids:
        show("  [caution] These PMIDs came from a live PubMed search, not a VCEP-curated "
             "citation list - unlike ERepo's evidenceLinks, a search hit is not confirmed to "
             "actually discuss this variant until the per-paper judgment checks it.")
    else:
        show("  -> PubMed search also found nothing; this variant cannot be evaluated by the literature path at all")
    return pmids, "pubmed_search"


async def judge_variant_from_structured_input(
    case_input: ApiCaseInput,
    mcp: ClientSession,
    erepo_client: ERepoClient,
    criteria: tuple[str, ...] = ("PS3", "BS3", "PS4"),
    vcep_name: str | None = None,
    full_text_cache: dict[str, tuple[str | None, str]] | None = None,
    llm_cache=None,
) -> dict[str, AggregatedJudgment]:
    """
    Parses `case_input`'s embedded VCF (exactly 1 variant, per ApiCaseInput's
    own contract) and runs judge_variant() once per requested literature
    criterion, using resolve_pmids_for_variant() above (ERepo first, live
    PubMed search fallback) to find citing PMIDs. Returns {criterion:
    AggregatedJudgment}; a criterion is omitted from the result (not given a
    not_clear placeholder) only when neither PMID source found anything for
    this variant at all, in which case the whole result dict is empty and
    the caller should treat this variant as "nothing this pipeline could
    evaluate".

    clinical_note.condition_id is resolved here (EBI OLS4, real candidates
    only) rather than trusted from extract_clinical_note() alone, which
    never sets it - see clinical_extraction._force_condition_id_null()'s
    own docstring for why that path was removed (2026-09-18). The import is
    local because
    hpo_mondo_extraction requires VLLM_BASE_URL/VLLM_API_KEY at import
    time, the same reason acmg_pipeline.criteria.pp1_bs4_pp4_engine.
    evaluate() imports it lazily too.
    """
    variant = case_input.parse_vcf().record
    clinical_note = extract_clinical_note(case_input.clinical_note)
    from acmg_pipeline import hpo_mondo_extraction
    clinical_note = await hpo_mondo_extraction.resolve_diagnosis_mondo(clinical_note, case_input.clinical_note)
    return await judge_variant_from_shared_input(
        variant,
        clinical_note,
        mcp,
        erepo_client,
        criteria=criteria,
        vcep_name=vcep_name,
        full_text_cache=full_text_cache,
        llm_cache=llm_cache,
    )


async def judge_variant_from_shared_input(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    mcp: ClientSession,
    erepo_client: ERepoClient,
    criteria: tuple[str, ...] = ("PS3", "BS3", "PS4"),
    vcep_name: str | None = None,
    full_text_cache: dict[str, tuple[str | None, str]] | None = None,
    llm_cache=None,
) -> dict[str, AggregatedJudgment]:
    """Run the five literature criteria directly from main's shared input class."""
    if not isinstance(variant, VariantRecord):
        raise TypeError("variant must be acmg_pipeline.vcf_record.VariantRecord")
    if not isinstance(clinical_note, ClinicalNoteExtraction):
        raise TypeError("clinical_note must be acmg_pipeline.clinical_note.ClinicalNoteExtraction")
    unknown = set(criteria) - LITERATURE_CODES
    if unknown:
        raise ValueError(f"Unsupported literature criteria: {sorted(unknown)}")
    gene = variant.info.get("GENE", "")
    hgvsc = variant.info.get("HGVSC", "")
    hgvsp = variant.info.get("HGVSP", "N/A")
    equivalents = list(_protein_equivalents(hgvsp)) if hgvsp and hgvsp != "N/A" else [hgvsc]

    disease = clinical_note.diagnosis
    pmids, _pmid_source = await resolve_pmids_for_variant(mcp, erepo_client, gene, hgvsc, hgvsp, disease)
    if not pmids:
        return {}

    results: dict[str, AggregatedJudgment] = {}
    for criterion in criteria:
        engine = ENGINE_BY_CRITERION[criterion]
        results[criterion] = await judge_variant(
            engine, mcp, pmids=pmids, gene=gene, hgvsc=hgvsc, hgvsp=hgvsp,
            equivalents=equivalents, vcep_name=vcep_name, criterion=criterion,
            full_text_cache=full_text_cache, llm_cache=llm_cache,
        )
    return results


_IDENTITY_INFO_KEYS = ("GENE", "TRANSCRIPT", "HGVSC", "HGVSP", "CLNVARIATIONID")


def _identity_from_info(variant: VariantRecord) -> dict:
    """The INFO subset the providers need to look this variant up.

    Deliberately an allowlist, not the whole INFO dict: CLNSIG and
    ACMG_CODES sit in the same column and are conclusions, not lookup
    keys - nothing downstream should be able to reach them by accident.
    Disease context is deliberately absent: it is supplied separately by the clinical-note
    parser, so a VCF CONDITION or DISEASE_ASSOCIATION value cannot select the disease.
    """
    wanted = {key.casefold(): key for key in _IDENTITY_INFO_KEYS}
    identity = {}
    for name, value in variant.info.items():
        key = wanted.get(name.casefold())
        if key is not None and value not in (None, ""):
            identity[key] = str(value)
    return identity


def _provider_identity(
    variant: VariantRecord, clinical_note: ClinicalNoteExtraction,
) -> dict:
    """Provider lookup identity plus the parser-owned MONDO disease context."""
    identity = _identity_from_info(variant)
    if clinical_note.condition_id:
        identity["CONDITION"] = clinical_note.condition_id
    return identity


_VCF_DISEASE_INFO_KEYS = {
    "condition", "condition_label", "condition_mapping", "condition_ancestors",
    "disease_association",
}


def _variant_without_vcf_disease_context(variant: VariantRecord) -> VariantRecord:
    """Copy a request variant while dropping every VCF-supplied disease field."""
    return VariantRecord(
        chrom=variant.chrom,
        pos=variant.pos,
        id=variant.id,
        ref=variant.ref,
        alt=variant.alt,
        qual=variant.qual,
        filter=variant.filter,
        info={name: value for name, value in variant.info.items()
              if name.casefold() not in _VCF_DISEASE_INFO_KEYS},
    )


_curated_context_cache: dict[str, dict] = {}


def _apply_curated_context(variant: VariantRecord, automated_config: dict) -> None:
    """Merge non-disease curated context such as the BA1 exception check.

    Opt-in via automated_config["curated_context_path"] (unset by default, so existing
    callers are unaffected). Before this, InitiationProvider/UpstreamPathogenicProvider's
    fix pattern repeated: acmg_pipeline.automated_core.context's load_context()/
    apply_context() already resolve the BA1 exception list precisely for every variant
    (not just the ones it lists - "complete": true means absence from it is itself
    resolved as is_exception=False), but only automated_cli.py's separate batch path ever
    called them - a live BA1 evaluation through this module always saw
    ba1_exception_assessment as absent, regardless of what curated-context.json already
    had recorded. See acmg_pipeline.automated_core.context's own docstring for why this
    lives in a separate versioned document rather than being read off VCF INFO directly.
    """
    path = automated_config.get("curated_context_path")
    if not path:
        return
    if path not in _curated_context_cache:
        from acmg_pipeline.automated_core.context import load_context
        document = json.loads(Path(path).read_text(encoding="utf-8"))
        _curated_context_cache[path] = load_context(document)
    from acmg_pipeline.automated_core.context import apply_context
    record = {
        "variant": {"assembly": "GRCh38", "chrom": variant.chrom, "pos": variant.pos,
                    "ref": variant.ref, "alt": variant.alt},
        "record_id": variant.id,
    }
    updated = apply_context(record, _curated_context_cache[path])
    # condition/condition_label entries in legacy documents are intentionally ignored. The
    # clinical-note parser is the only owner of the disease selected for this case.
    for key in ("inheritance", "disease_frequency_threshold", "ba1_exception_assessment"):
        if key in updated:
            variant.info[key] = updated[key]


# What a request carries about which variant it means, in the shape reconcile() audits.
_IDENTITY_KEYS = ("TRANSCRIPT", "HGVSC", "CLNVARIATIONID")


def _resolve_identity(variant: VariantRecord, resolver) -> list[str]:
    """Settle which variant a request means, and say so, before anything is evaluated.

    A request names a variant twice: as coordinates, and as a transcript HGVS or a ClinVar
    accession. Those can disagree, and one of them can be missing - the project's own demo
    request bodies leave ALT as "." and put the identity in TRANSCRIPT/HGVSC, which nothing
    on this path was resolving, so evaluation stopped before it began.

    acmg_pipeline.automated_core.identity.reconcile() is what settles it, and it is an audit
    rather than a lookup: a candidate has to match the identifiers this request supplied,
    version included, and carry its own provenance, and a coordinate the request got wrong is
    corrected only on corroborated evidence. It already backs the batch command; this is the
    same function, on the same Ensembl provider the evidence comes from.

    Evaluation continues either way. An unresolved identity is reported, not raised: the
    criteria will see whatever coordinates the request gave and answer from those, and a
    caller that cannot tell a verified variant from an unverified one is worse off than one
    holding an answer it has been told to check. The issues are returned for the caller to
    surface, and the status travels on the variant.
    """
    record = {
        "record_id": variant.id,
        "raw_variant": {"assembly": "GRCh38", "chrom": variant.chrom, "pos": variant.pos,
                        "ref": variant.ref, "alt": variant.alt},
        "parsed_variant": None,
        "identity": {key: value for key, value in variant.info.items()
                     if key.upper() in _IDENTITY_KEYS},
        "issues": [],
    }
    try:
        record["parsed_variant"] = AutomatedVariant(
            assembly="GRCh38", chrom=variant.chrom, pos=int(variant.pos),
            ref=variant.ref, alt=variant.alt).to_dict()
    except (ValueError, TypeError):
        # ALT "." and the like: nothing to compare a candidate against, which is exactly the
        # case reconcile() reports as CORRECTED rather than VERIFIED.
        pass
    provider = getattr(resolver, "identity_provider", None)
    if provider is None:
        # A caller that injected its own resolver supplied the evidence itself and has said,
        # by doing so, which variant it is about. There is nothing independent to audit
        # against, so the request's own coordinates stand.
        return
    try:
        candidate, _annotation, _predictions = provider.map_record_with_evidence(record)
        candidates = [candidate]
    except PROVIDER_ERRORS as exc:
        record["issues"].append(f"IDENTITY_PROVIDER_ERROR: {exc}")
        candidates = []
    outcome = reconcile(record, candidates, provider.reference)
    variant.info["identity_status"] = outcome["identity_status"]
    resolution = outcome.get("resolution")
    if resolution:
        settled = resolution["variant"]
        variant.chrom, variant.pos = settled["chrom"], settled["pos"]
        variant.ref, variant.alt = settled["ref"], settled["alt"]
        variant.info["identity_provenance"] = resolution["evidence"]
    if outcome["issues"]:
        variant.info["identity_issues"] = list(outcome["issues"])


def _automated_variant(variant: VariantRecord) -> AutomatedVariant:
    """The evidence-cli Variant for provider lookups (GRCh38 unless INFO says otherwise)."""
    assembly = next(
        (str(value) for name, value in variant.info.items() if name.casefold() == "assembly"),
        "GRCh38",
    )
    return AutomatedVariant(
        assembly=assembly,
        chrom=variant.chrom,
        pos=variant.pos,
        ref=variant.ref,
        alt=variant.alt,
    )


@lru_cache(maxsize=4)
def _reviewed_gene_disease_document(path: str) -> tuple[dict, ...]:
    """The reviewed gene-disease assessments, expanded once per server process.

    Cached on the path: this is server-owned configuration, read on every request and
    identical between them, so re-expanding 7 groups per criterion evaluation buys nothing.
    A change to the file takes effect on restart, like the rest of load_automated_config().
    """
    document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    return tuple(build_assessment_document(document)["evidence"])


def _reviewed_gene_disease_records(variant: AutomatedVariant, automated_config: dict) -> list[dict]:
    """The reviewed mechanism assessments for exactly this variant, or [].

    PP2 and BP1 read a transcript-scoped mechanism and variant-spectrum review. No provider
    produces one: ClinGen Dosage and Gene2Phenotype state whether loss of function is a
    mechanism, which answers PVS1's gate and nothing else. Until this was wired in, the
    server had no way to reach a reviewed record at all and always fell through to
    mechanism.py's gene_disease_draft suggestion - so the same variant could be answered
    from a curator's review through the CLI and from gnomAD constraint through the API.

    Matching is on the exact variant key, never on the gene: a reviewed decision about one
    variant in a gene is not a decision about another one. A variant with no reviewed record
    still reaches the draft suggestion, which is the behavior this replaces only where a
    review exists.
    """
    path = automated_config.get("gene_disease_assessments")
    if not path:
        return []
    key = variant.key
    return [dict(item) for item in _reviewed_gene_disease_document(str(path))
            if item.get("variant_key") == key]


async def evaluate_variant_evidence_lines(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    *,
    automated_config: dict,
    mcp: ClientSession,
    erepo_client: ERepoClient,
    evidence_resolver=None,
    vcep_name: str | None = None,
    full_text_cache: dict[str, tuple[str | None, str]] | None = None,
    llm_cache=None,
) -> list[dict]:
    """Return exactly one VA-Spec EvidenceLine for each of the 28 ACMG codes.

    Evidence for the automated criteria is retrieved and normalized by the
    server-side ProviderEvidenceResolver. The VCF INFO column supplies variant
    identity only (GENE, TRANSCRIPT, HGVS and ClinVar accession). Disease context
    comes exclusively from ClinicalNoteExtraction. INFO is never read as evidence:
    the demo VCFs' CLNSIG/ACMG_CODES are already-reached conclusions, the
    AM_*/AG_* scores have no calibration entry, and no population
    frequency is present at all. See acmg/services/resolve.py.
    """
    if not isinstance(variant, VariantRecord):
        raise TypeError("variant must be acmg_pipeline.vcf_record.VariantRecord")
    if not isinstance(clinical_note, ClinicalNoteExtraction):
        raise TypeError(
            "clinical_note must be acmg_pipeline.clinical_note.ClinicalNoteExtraction"
        )
    if not isinstance(automated_config, dict):
        raise TypeError("automated_config must be a dictionary")

    resolver = evidence_resolver or ProviderEvidenceResolver(
        automated_config.get("evidence_cache_dir", "cache/evidence"),
        offline=bool(automated_config.get("offline")),
        ensembl_release=automated_config.get("ensembl_release"),
        population_sources=automated_config.get("population_sources"),
        hotspot_policy=automated_config.get("PM1", {}).get("hotspot"),
        with_clingen_dosage=bool(automated_config.get("with_clingen_dosage")),
        with_pvs1_transcript_gates=bool(automated_config.get("with_pvs1_transcript_gates")),
        with_gene_disease_draft=bool(automated_config.get("gene_disease_draft")),
        gene_disease_draft_policy=automated_config.get("gene_disease_draft"),
        with_clinvar_spectrum=bool(automated_config.get("with_clinvar_spectrum")),
        with_splice_default=bool(automated_config.get("PVS1", {}).get("splice_default_policy_version")),
        splice_default_policy_version=automated_config.get("PVS1", {}).get("splice_default_policy_version"),
        with_initiation_assessment=bool(automated_config.get("PVS1", {}).get("with_initiation_assessment")),
        with_gene2phenotype=bool(automated_config.get("with_gene2phenotype")),
        cspec_applicability_path=automated_config.get("cspec_applicability"),
        with_disease_matching=bool(automated_config.get("with_disease_matching")),
        with_gene_disease_associations=bool(
            automated_config.get("with_gene_disease_associations")),
    )
    # Identity first: the curated context is looked up by the variant's own key, and a
    # request that named its variant only as a transcript HGVS has no key until this runs.
    _resolve_identity(variant, resolver)
    _apply_curated_context(variant, automated_config)
    automated_variant = _automated_variant(variant)
    resolved = resolver.resolve(
        _provider_identity(variant, clinical_note), automated_variant)
    services = make_services(
        [*resolved.records,
         *_reviewed_gene_disease_records(automated_variant, automated_config)],
        automated_config.get("population_providers"),
        failures=resolved.failures,
    )
    criterion_variant = _variant_without_vcf_disease_context(variant)
    automated_results = evaluate_automated_record(
        criterion_variant,
        clinical_note,
        services,
        automated_config,
        criteria=AUTOMATED_CRITERIA,
    )
    automated_by_code = {result.criterion: result for result in automated_results}

    literature_results = await judge_variant_from_shared_input(
        variant,
        clinical_note,
        mcp,
        erepo_client,
        vcep_name=vcep_name,
        full_text_cache=full_text_cache,
        llm_cache=llm_cache,
    )
    gene = str(variant.info.get("GENE", ""))
    hgvsc = str(variant.info.get("HGVSC", ""))

    by_code: dict[str, dict] = {}
    for code in AUTOMATED_CODES:
        result = automated_by_code.get(code)
        if result is None:
            by_code[code] = build_workflow_evidence_line(
                code,
                variant,
                status="unknown",
                description=f"{code} automated evaluation returned no result.",
                details={"missingInputs": ["automated criterion result"]},
            )
        else:
            by_code[code] = build_automated_evidence_line(result, variant)

    for code in LITERATURE_CODES:
        aggregated = literature_results.get(code)
        if aggregated is None:
            by_code[code] = build_workflow_evidence_line(
                code,
                variant,
                status="unknown",
                description=f"{code} was not evaluated because no literature was resolved.",
                details={"missingInputs": ["resolvable literature PMID"]},
            )
        else:
            by_code[code] = build_evidence_line(
                aggregated,
                gene,
                hgvsc,
                code,
                vcep_name=vcep_name,
                variant=variant,
            )

    phenotype_segregation_results = await pp1_bs4_pp4_engine.evaluate(variant, clinical_note, automated_config)
    for code in PHENOTYPE_SEGREGATION_CODES:
        by_code[code] = pp1_bs4_pp4_engine.build_evidence_line(
            code, phenotype_segregation_results[code], variant,
        )

    for code in ALL_ACMG_CODES:
        if code not in IMPLEMENTED_CODES:
            by_code[code] = build_stub_evidence_line(code, variant)

    lines = [by_code[code] for code in ALL_ACMG_CODES]
    ids = [line["id"] for line in lines]
    method_codes = [line["specifiedBy"]["methodType"] for line in lines]
    if len(lines) != 28 or len(set(ids)) != 28 or tuple(method_codes) != ALL_ACMG_CODES:
        raise RuntimeError("Integrated ACMG output must contain 28 ordered, unique criteria")
    return lines


async def evaluate_selected_criteria(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    criteria: tuple[str, ...],
    *,
    automated_config: dict,
    mcp: ClientSession | None = None,
    erepo_client: ERepoClient | None = None,
    evidence_resolver=None,
    vcep_name: str | None = None,
    full_text_cache: dict[str, tuple[str | None, str]] | None = None,
    llm_cache=None,
) -> dict[str, dict]:
    """Return one VA-Spec EvidenceLine per requested code only.

    Unlike evaluate_variant_evidence_lines() (always all 28), this only runs
    the automated and/or literature machinery actually needed for `criteria`
    - a stub-only or automated-only request never touches PubMed/the LLM,
    and mcp/erepo_client may be left None in that case. Raises ValueError up
    front if `criteria` needs LITERATURE_CODES but mcp/erepo_client weren't
    given, so a caller's sync/async routing mistake fails loudly instead of
    silently returning UNKNOWN literature lines.
    """
    if not isinstance(variant, VariantRecord):
        raise TypeError("variant must be acmg_pipeline.vcf_record.VariantRecord")
    if not isinstance(clinical_note, ClinicalNoteExtraction):
        raise TypeError("clinical_note must be acmg_pipeline.clinical_note.ClinicalNoteExtraction")
    if not isinstance(automated_config, dict):
        raise TypeError("automated_config must be a dictionary")

    unknown = [code for code in criteria if code not in ALL_ACMG_CODES]
    if unknown:
        raise ValueError(f"Unrecognized ACMG code(s): {unknown}")

    requested = set(criteria)
    automated_subset = tuple(code for code in ALL_ACMG_CODES if code in requested and code in AUTOMATED_CODES)
    literature_subset = tuple(code for code in ALL_ACMG_CODES if code in requested and code in LITERATURE_CODES)
    phenotype_segregation_subset = tuple(
        code for code in ALL_ACMG_CODES if code in requested and code in PHENOTYPE_SEGREGATION_CODES
    )
    stub_subset = [code for code in criteria if code not in IMPLEMENTED_CODES]

    if literature_subset and (mcp is None or erepo_client is None):
        raise ValueError(
            f"criteria {literature_subset} require the literature workflow "
            "(mcp + erepo_client), but none were provided"
        )

    by_code: dict[str, dict] = {}

    if automated_subset:
        resolver = evidence_resolver or ProviderEvidenceResolver(
            automated_config.get("evidence_cache_dir", "cache/evidence"),
            offline=bool(automated_config.get("offline")),
            ensembl_release=automated_config.get("ensembl_release"),
            population_sources=automated_config.get("population_sources"),
            hotspot_policy=automated_config.get("PM1", {}).get("hotspot"),
            with_clingen_dosage=bool(automated_config.get("with_clingen_dosage")),
            with_pvs1_transcript_gates=bool(automated_config.get("with_pvs1_transcript_gates")),
            with_gene_disease_draft=bool(automated_config.get("gene_disease_draft")),
            gene_disease_draft_policy=automated_config.get("gene_disease_draft"),
            with_clinvar_spectrum=bool(automated_config.get("with_clinvar_spectrum")),
            with_splice_default=bool(automated_config.get("PVS1", {}).get("splice_default_policy_version")),
            splice_default_policy_version=automated_config.get("PVS1", {}).get("splice_default_policy_version"),
            with_initiation_assessment=bool(automated_config.get("PVS1", {}).get("with_initiation_assessment")),
            with_gene2phenotype=bool(automated_config.get("with_gene2phenotype")),
            cspec_applicability_path=automated_config.get("cspec_applicability"),
            with_disease_matching=bool(automated_config.get("with_disease_matching")),
            with_gene_disease_associations=bool(
                automated_config.get("with_gene_disease_associations")),
        )
        # Identity first - see evaluate_variant_evidence_lines() for why.
        _resolve_identity(variant, resolver)
        _apply_curated_context(variant, automated_config)
        automated_variant = _automated_variant(variant)
        resolved = resolver.resolve(
            _provider_identity(variant, clinical_note), automated_variant)
        # The same reviewed assessments the all-28 path gets: a caller asking only for PP2
        # must not be answered from a weaker source than one asking for everything.
        services = make_services(
            [*resolved.records,
             *_reviewed_gene_disease_records(automated_variant, automated_config)],
            automated_config.get("population_providers"),
            failures=resolved.failures)
        criterion_variant = _variant_without_vcf_disease_context(variant)
        automated_results = evaluate_automated_record(
            criterion_variant, clinical_note, services, automated_config, criteria=automated_subset,
        )
        automated_by_code = {result.criterion: result for result in automated_results}
        for code in automated_subset:
            result = automated_by_code.get(code)
            if result is None:
                by_code[code] = build_workflow_evidence_line(
                    code, variant, status="unknown",
                    description=f"{code} automated evaluation returned no result.",
                    details={"missingInputs": ["automated criterion result"]},
                )
            else:
                by_code[code] = build_automated_evidence_line(result, variant)

    if literature_subset:
        literature_results = await judge_variant_from_shared_input(
            variant, clinical_note, mcp, erepo_client, criteria=literature_subset,
            vcep_name=vcep_name, full_text_cache=full_text_cache, llm_cache=llm_cache,
        )
        gene = str(variant.info.get("GENE", ""))
        hgvsc = str(variant.info.get("HGVSC", ""))
        for code in literature_subset:
            aggregated = literature_results.get(code)
            if aggregated is None:
                by_code[code] = build_workflow_evidence_line(
                    code, variant, status="unknown",
                    description=f"{code} was not evaluated because no literature was resolved.",
                    details={"missingInputs": ["resolvable literature PMID"]},
                )
            else:
                by_code[code] = build_evidence_line(
                    aggregated, gene, hgvsc, code, vcep_name=vcep_name, variant=variant,
                )

    if phenotype_segregation_subset:
        phenotype_segregation_results = await pp1_bs4_pp4_engine.evaluate(variant, clinical_note, automated_config)
        for code in phenotype_segregation_subset:
            by_code[code] = pp1_bs4_pp4_engine.build_evidence_line(
                code, phenotype_segregation_results[code], variant,
            )

    for code in stub_subset:
        by_code[code] = build_stub_evidence_line(code, variant)

    return {code: by_code[code] for code in criteria}


# ---------------------------------------------------------------------------
# Full 28-code classification from structured input
# ---------------------------------------------------------------------------
#
# Added 2026-09-16, per the user's decision that criteria this project
# This older classification-only entry point runs the five literature
# criteria. Use evaluate_variant_evidence_lines() above for the integrated
# 28-code VA-Spec output (five literature, sixteen automated, seven stubs).
# classify() exactly the same way test_full_criteria_ground_truth.py
# already does for the ground-truth dataset.
#
# Execution order (cheapest-first, matching classification.classify()'s own
# BA1-short-circuit philosophy): stub lookups are instant, so in practice
# they're computed after the (comparatively expensive) literature judgment
# here - there's nothing today for a stub to short-circuit, since this
# project doesn't compute the Layer-1 values (gnomAD AF, ClinVar assertions,
# etc.) that would make e.g. a BA1 short-circuit meaningful. Once another
# team's real Layer-1/PP4 modules exist, swap their real CriterionEvidence
# in place of stubs.stub_evidence(code) below - registry.get_criterion_
# evidence() already supports exactly that swap without any other change
# here.

async def classify_variant_from_structured_input(
    case_input: ApiCaseInput,
    mcp: ClientSession,
    erepo_client: ERepoClient,
    vcep_name: str | None = None,
    full_text_cache: dict[str, tuple[str | None, str]] | None = None,
    llm_cache=None,
) -> ClassificationResult:
    literature_results = await judge_variant_from_structured_input(
        case_input, mcp, erepo_client, vcep_name=vcep_name,
        full_text_cache=full_text_cache, llm_cache=llm_cache,
    )

    evidence = [
        from_aggregated_judgment(aggregated, code)
        for code, aggregated in literature_results.items()
    ]
    evaluated_codes = set(literature_results)
    for code in ALL_ACMG_CODES:
        if code in evaluated_codes:
            continue
        if code in IMPLEMENTED_CODES:
            # Implemented by this project, but not evaluated for this
            # particular variant (e.g. ERepo had zero PMIDs at all, so
            # judge_variant_from_structured_input() returned {} - see its
            # own docstring). Distinct from "not our team's job" (the stub
            # branch below): this is "should have been run, wasn't", so
            # classify()'s not_evaluated_codes correctly flags it rather
            # than silently treating it as a stub.
            continue
        evidence.append(stubs.stub_evidence(code))

    result = classify(evidence)
    show(f"\n{'='*70}\n[Classification]\n{'='*70}")
    show(f"category = {result.category.value} (score={result.score}"
         f"{', BA1 override' if result.ba1_override else ''})")
    show(f"  met: {', '.join(f'{e.code}({e.strength.value})' for e in result.met) or '(none)'}")
    show(f"  not_met: {', '.join(e.code for e in result.not_met) or '(none)'}")
    if result.not_evaluated_codes:
        show(f"  not_evaluated (implemented but not run for this variant): "
             f"{', '.join(c for c in result.not_evaluated_codes if c in IMPLEMENTED_CODES) or '(none)'}")
    return result


def parse_ground_truth(ground_truth: str) -> tuple[str, bool]:
    """Parses a 'PS3-Moderate = met' style string into (code, is_met)."""
    code, status = (s.strip() for s in ground_truth.split("="))
    return code, status == "met"


def score_direction(direction, criterion: str, gt_is_met: bool) -> str:
    """
    Buckets one variant's outcome as 'match' / 'mismatch' / 'reserved',
    following the same 3-way scheme used for the manual CGBench pilot in
    the design doc (section 7): a not_clear direction is always 'reserved'
    (the pipeline declined to commit, which is the safe failure mode, not
    a scored error). Otherwise, 'match' iff the committed direction
    matching `criterion` agrees with whether ground truth says the
    criterion was met.

    Checks `direction.value == "not_clear"` rather than comparing against a
    specific Enum's NOT_CLEAR member, so this works for any of the three
    engines' direction Enums (OverallDirection / CaseControlDirection /
    SegregationDirection) - same convention as va_spec_export.py.
    """
    if direction.value == "not_clear":
        return "reserved"
    implied_met = (direction.value == criterion)
    return "match" if implied_met == gt_is_met else "mismatch"


async def connect_pubmed(stack: AsyncExitStack) -> ClientSession:
    """
    Connects to the PubMed MCP server. The number of values returned by
    streamablehttp_client differs by mcp version (confirmed with mcp==1.9.4:
    3 values - read_stream, write_stream, get_session_id_callback. Directly
    unpacking as `read, write = ...` assuming 2 values raises
    "ValueError: too many values to unpack", which actually occurred). Using
    only the first two elements makes this work whether the installed mcp
    version returns 2 or 3 values.
    """
    ctx = await stack.enter_async_context(streamable_http_client(MCP_SERVERS["pubmed"]))
    read, write = ctx[0], ctx[1]
    mcp = await stack.enter_async_context(ClientSession(read, write))
    await mcp.initialize()
    return mcp


async def main():
    # One timestamp for the whole run, prefixed onto every output filename so
    # files from different runs sort together and never silently clobber an
    # earlier run's output for the same variant.
    run_ts = datetime.now().strftime("%Y%m%d%H%M%S")
    async with AsyncExitStack() as stack:
        mcp = await connect_pubmed(stack)
        show(f"[MCP] Connected to PubMed")

        erepo_client = ERepoClient()

        # Test cases reproducing, with a real LLM (gemma-4), the "genuinely
        # blind" verification done manually by Claude earlier in the design
        # conversation (3 cases), so the results can be compared directly.
        #
        # "pmids" is no longer hard-coded here: it is looked up live from
        # ClinGen ERepo (ERepoLookupResult.evidence_pmids) for each variant
        # just before judge_variant() runs, in the loop below. This can
        # return more (or fewer/different) PMIDs than the single
        # PS3/BS3-evidence-code-row citation used previously, since
        # evidenceLinks is the full citation list for the variant's overall
        # classification, not filtered to one criterion.
        test_cases = [
            dict(
                gene="MYH7", hgvsc="c.1594T>C", hgvsp="p.Ser532Pro",
                equivalents=["S532P", "p.(Ser532Pro)"],
                vcep_name="Cardiomyopathy VCEP", criterion="PS3",
                ground_truth="PS3 = met",
            ),
            dict(
                gene="PTEN", hgvsc="c.112C>T", hgvsp="p.Pro38Ser",
                equivalents=["P38S", "p.(Pro38Ser)"],
                vcep_name="PTEN VCEP", criterion="PS3",
                # Corrected 2026-09-15: originally "PS3-Moderate = not_met",
                # carried over from pre-ERepo-integration design-doc work and
                # never cross-checked against live ERepo. A user-prompted
                # audit of all main() ground_truth values against live
                # ERepoClient found this (and 3 other cases below) to be the
                # same CGBench met_status labeling bug documented in design
                # doc sections 8-9, recurring here independently of the
                # audited 36-case list. Live ERepo: PS3_Moderate = met.
                ground_truth="PS3_Moderate = met",
            ),
            dict(
                gene="TP53", hgvsc="c.875A>G", hgvsp="p.Lys292Arg",
                equivalents=["K292R", "p.(Lys292Arg)"],
                vcep_name="TP53 VCEP", criterion="BS3",
                ground_truth="BS3 = met",
            ),
            # Added 2026-09-15: 4 fresh variants from the previously-untouched
            # CGBench pool (clingen_vci_pubmed_var_na_filtered.csv), selected
            # to expand the sample size beyond the original 3 (see design doc
            # section 14, "next step"). Chosen from candidates confirmed to
            # have PMC full text.
            dict(
                gene="PTEN", hgvsc="c.170T>G", hgvsp="p.Leu57Trp",
                equivalents=["L57W", "p.(Leu57Trp)"],
                vcep_name="PTEN VCEP", criterion="PS3",
                # Corrected 2026-09-15: live ERepo has no bare "PS3" key for
                # this variant, only "PS3_Moderate" (met) - the strength
                # suffix was simply missing before, not a met/not_met error.
                ground_truth="PS3_Moderate = met",
            ),
            dict(
                gene="PTEN", hgvsc="c.278A>G", hgvsp="p.His93Arg",
                equivalents=["H93R", "p.(His93Arg)"],
                vcep_name="PTEN VCEP", criterion="PS3",
                # Corrected 2026-09-15 (see PTEN c.112C>T above for context):
                # live ERepo: PS3_Supporting = met, not not_met. This was the
                # sole "mismatch" reported in the earlier 19-case LLM run -
                # the LLM's PS3 call was actually correct; the recorded
                # ground truth was the error.
                ground_truth="PS3_Supporting = met",
            ),
            dict(
                gene="MTOR", hgvsc="c.4447T>C", hgvsp="p.Cys1483Arg",
                equivalents=["C1483R", "p.(Cys1483Arg)"],
                vcep_name="Brain Malformations VCEP", criterion="PS3",
                # Corrected 2026-09-15: live ERepo: PS3_Supporting = met.
                ground_truth="PS3_Supporting = met",
            ),
            dict(
                gene="MYOC", hgvsc="c.1187_1188insCCCAGA", hgvsp="N/A",
                equivalents=["c.1187_1188insCCCAGA"],
                vcep_name="Glaucoma VCEP", criterion="PS3",
                # Corrected 2026-09-15: live ERepo: PS3_Moderate = met.
                ground_truth="PS3_Moderate = met",
            ),
            # Added 2026-09-15 (2nd expansion): 12 more variants, mined from
            # the same CGBench pool (clingen_vci_pubmed_fulltext_dedup_pmid_
            # CORRECTED.csv) for PS3/BS3 rows with a clear single-criterion
            # met status, across genes/VCEPs not already covered above. Each
            # one was checked against the now-working live ERepoClient before
            # being added here (found_in_erepo=True with non-empty
            # evidence_pmids for all 12).
            dict(
                gene="BRCA1", hgvsc="c.191G>A", hgvsp="p.Cys64Tyr",
                equivalents=["C64Y", "p.(Cys64Tyr)"],
                vcep_name="ENIGMA BRCA1 and BRCA2 VCEP", criterion="PS3",
                ground_truth="PS3 = met",
            ),
            dict(
                gene="BRCA2", hgvsc="c.831T>G", hgvsp="p.Asn277Lys",
                equivalents=["N277K", "p.(Asn277Lys)"],
                vcep_name="ENIGMA BRCA1 and BRCA2 VCEP", criterion="BS3",
                ground_truth="BS3 = met",
            ),
            dict(
                gene="CDH1", hgvsc="c.1057G>A", hgvsp="p.Glu353Lys",
                equivalents=["E353K", "p.(Glu353Lys)"],
                vcep_name="CDH1 VCEP", criterion="PS3",
                ground_truth="PS3 = met",
            ),
            dict(
                gene="CDH1", hgvsc="c.387+5G>A", hgvsp="N/A",
                equivalents=["c.387+5G>A"],
                vcep_name="CDH1 VCEP", criterion="BS3",
                ground_truth="BS3 = met",
            ),
            dict(
                gene="DYSF", hgvsc="c.1717C>T", hgvsp="p.Arg573Trp",
                equivalents=["R573W", "p.(Arg573Trp)"],
                vcep_name="Limb Girdle Muscular Dystrophy VCEP", criterion="PS3",
                ground_truth="PS3-Moderate = met",
            ),
            dict(
                gene="GUCY2D", hgvsc="c.1762C>T", hgvsp="p.Arg588Trp",
                equivalents=["R588W", "p.(Arg588Trp)"],
                vcep_name="Leber Congenital Amaurosis/early onset Retinal Dystrophy VCEP",
                criterion="PS3",
                ground_truth="PS3-Supporting = met",
            ),
            dict(
                gene="LDLR", hgvsc="c.2575G>A", hgvsp="p.Val859Met",
                equivalents=["V859M", "p.(Val859Met)"],
                vcep_name="Familial Hypercholesterolemia VCEP", criterion="BS3",
                ground_truth="BS3 = met",
            ),
            dict(
                gene="MSH2", hgvsc="c.1012G>A", hgvsp="p.Gly338Arg",
                equivalents=["G338R", "p.(Gly338Arg)"],
                vcep_name="InSiGHT Hereditary Colorectal Cancer/Polyposis VCEP",
                criterion="PS3",
                ground_truth="PS3 = met",
            ),
            dict(
                gene="PTPN11", hgvsc="c.184T>G", hgvsp="p.Tyr62Asp",
                equivalents=["Y62D", "p.(Tyr62Asp)"],
                vcep_name="RASopathy VCEP", criterion="PS3",
                ground_truth="PS3 = met",
            ),
            dict(
                gene="RIT1", hgvsc="c.268A>G", hgvsp="p.Met90Val",
                equivalents=["M90V", "p.(Met90Val)"],
                vcep_name="RASopathy VCEP", criterion="PS3",
                ground_truth="PS3-Supporting = met",
            ),
            dict(
                gene="SCN2A", hgvsc="c.5317G>A", hgvsp="p.Ala1773Thr",
                equivalents=["A1773T", "p.(Ala1773Thr)"],
                vcep_name="Epilepsy Sodium Channel VCEP", criterion="PS3",
                ground_truth="PS3 = met",
            ),
            dict(
                gene="VHL", hgvsc="c.273C>A", hgvsp="p.Phe91Leu",
                equivalents=["F91L", "p.(Phe91Leu)"],
                vcep_name="VHL VCEP", criterion="PS3",
                ground_truth="PS3-Supporting = met",
            ),
            # Added 2026-09-15 (3rd expansion): PS4 (case-control) and PP1/BS4
            # (family segregation) - the two other literature-backed criteria
            # in scope per the user's decision (see doc/BH26_participant_
            # briefing_v3_en.md line 66 for the full "manual" Layer-3 list;
            # PS2/PM3/PM6/BS2/BP2/BP5 depend on the patient's own clinical/
            # genetic-testing records rather than papers, so stay out of
            # scope for this LLM pipeline). Each variant/criterion pair below
            # was checked against the live ERepoClient before being added
            # (found_in_erepo=True, non-empty evidence_pmids, and the named
            # evidence code present in evidence_code_status).
            dict(
                gene="RUNX1", hgvsc="c.601C>T", hgvsp="p.Arg201Ter",
                equivalents=["R201*", "p.(Arg201Ter)"],
                vcep_name="Myeloid Malignancy VCEP", criterion="PS4",
                ground_truth="PS4 = met",
            ),
            dict(
                gene="GJB2", hgvsc="c.109G>A", hgvsp="p.Val37Ile",
                equivalents=["V37I", "p.(Val37Ile)"],
                vcep_name="Hearing Loss VCEP", criterion="PS4",
                ground_truth="PS4 = met",
            ),
            dict(
                gene="RUNX1", hgvsc="c.601C>T", hgvsp="p.Arg201Ter",
                equivalents=["R201*", "p.(Arg201Ter)"],
                vcep_name="Myeloid Malignancy VCEP", criterion="PP1",
                ground_truth="PP1_Strong = met",
            ),
            dict(
                gene="USH2A", hgvsc="c.8559-2A>G", hgvsp="N/A",
                equivalents=["c.8559-2A>G"],
                vcep_name="Hearing Loss VCEP", criterion="PP1",
                ground_truth="PP1_Strong = met",
            ),
            dict(
                gene="RUNX1", hgvsc="c.601C>T", hgvsp="p.Arg201Ter",
                equivalents=["R201*", "p.(Arg201Ter)"],
                vcep_name="Myeloid Malignancy VCEP", criterion="BS4",
                ground_truth="BS4 = not_met",
            ),
            dict(
                gene="GJB2", hgvsc="c.109G>A", hgvsp="p.Val37Ile",
                equivalents=["V37I", "p.(Val37Ile)"],
                vcep_name="Hearing Loss VCEP", criterion="BS4",
                ground_truth="BS4 = not_met",
            ),
        ]

        tally = {"match": 0, "mismatch": 0, "reserved": 0}
        mismatches = []
        # Accumulates this run's real CriterionEvidence per variant (gene,
        # hgvsc), keyed by which of the 5 implemented codes were actually
        # evaluated for it - test_cases above often covers only one or two
        # of the 5 per variant (e.g. MYH7 c.1594T>C is only ever run for
        # PS3 here), not all 5, so this can't reuse registry.
        # get_criterion_evidence()'s all-28-codes contract (it requires real
        # evidence for every implemented code, by design - see registry.py).
        # Instead the final classification pass below combines whatever was
        # actually run with the 23 stub codes directly, and lets classify()
        # itself report the implemented-but-not-run-here codes as
        # not_evaluated_codes - an honest reflection of this demo's partial
        # coverage, not a claim that this project doesn't implement them.
        variant_evidence: dict[tuple[str, str], dict[str, object]] = {}

        for case in test_cases:
            gene, hgvsc = case["gene"], case["hgvsc"]
            engine = ENGINE_BY_CRITERION[case["criterion"]]
            erepo_result = erepo_client.lookup(gene, hgvsc)
            pmids = erepo_result.evidence_pmids
            show(f"\n[ERepo] {gene} {hgvsc}: found_in_erepo={erepo_result.found_in_erepo}, "
                 f"evidence_pmids={pmids}")
            if not pmids:
                show("  -> No evidence PMIDs from ERepo; skipping this variant")
                continue

            aggregated = await judge_variant(engine, mcp, pmids=pmids, **case)
            direction = aggregated.aggregated_direction

            criterion_evidence = from_aggregated_judgment(aggregated, case["criterion"])
            variant_evidence.setdefault((gene, hgvsc), {})[case["criterion"]] = criterion_evidence

            evidence_line = build_evidence_line(
                aggregated, gene=gene, hgvsc=hgvsc,
                criterion=case["criterion"], vcep_name=case.get("vcep_name"),
            )
            out_path = VA_SPEC_OUTPUT_DIR / f"{run_ts}_{evidence_line['id'].replace('evline:', '')}.json"
            with out_path.open("w", encoding="utf-8") as f:
                json.dump(evidence_line, f, ensure_ascii=False, indent=2)
            show(f"[VA-Spec] Wrote EvidenceLine -> {out_path}")

            ground_truth = case.get("ground_truth")
            if ground_truth is not None:
                _, gt_is_met = parse_ground_truth(ground_truth)
                bucket = score_direction(direction, case["criterion"], gt_is_met)
                tally[bucket] += 1
                if bucket == "mismatch":
                    mismatches.append(
                        f"{gene} {hgvsc} ({case['criterion']}): "
                        f"ground_truth={ground_truth}, LLM direction={direction.value}"
                    )

        show(f"\n{'='*70}\n[Final tally over {sum(tally.values())} scored variant(s)]\n{'='*70}")
        show(f"match={tally['match']}  mismatch={tally['mismatch']}  reserved(not_clear)={tally['reserved']}")
        committed = tally["match"] + tally["mismatch"]
        if committed:
            show(f"Accuracy among committed judgments: {tally['match']}/{committed} = "
                 f"{tally['match'] / committed:.1%}")
        else:
            show("No committed judgments to score (every case resolved to not_clear).")
        if mismatches:
            show("\nMismatches:")
            for m in mismatches:
                show(f"  - {m}")

        # Final ACMG/AMP classification per variant (classification.classify(),
        # see acmg_pipeline/classification.py). Combines whichever of the 5
        # implemented codes were actually run above for that variant with
        # NOT_EVALUATED placeholders for the 23 codes this project doesn't
        # implement (acmg_pipeline/criteria/stubs.py) - this is a pipeline
        # demo, not a real curation, so the category shown here is only as
        # complete as the evidence gathered above.
        show(f"\n{'='*70}\n[Final ACMG/AMP classification per variant]\n{'='*70}")
        stub_evidence = stubs_mod.all_stub_evidence()
        for (v_gene, v_hgvsc), real_by_code in variant_evidence.items():
            full_evidence = list(real_by_code.values()) + stub_evidence
            result = classify(full_evidence)
            show(f"\n{v_gene} {v_hgvsc}:")
            met_str = ", ".join(f"{e.code}({e.strength.value})" for e in result.met) or "(none)"
            not_met_str = ", ".join(e.code for e in result.not_met) or "(none)"
            show(f"  category = {result.category.value}  (score={result.score}"
                 f"{', BA1 override' if result.ba1_override else ''})")
            show(f"  met: {met_str}")
            show(f"  not_met: {not_met_str}")
            implemented_not_run = [c for c in result.not_evaluated_codes if c not in stubs_mod.STUB_CODES]
            if implemented_not_run:
                show(f"  not run in this demo (implemented, but no test case above covered it): {', '.join(implemented_not_run)}")


if __name__ == "__main__":
    asyncio.run(main())
