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
from openai import OpenAI

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
from acmg_pipeline.gate import ERepoClient
from acmg_pipeline.export import build_evidence_line
from acmg_pipeline.classification import classify, from_aggregated_judgment

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

VLLM_BASE_URL = os.environ.get("VLLM_BASE_URL", "")
VLLM_API_KEY = os.environ.get("VLLM_API_KEY", "")
if not VLLM_BASE_URL or not VLLM_API_KEY:
    raise RuntimeError(
        "VLLM_BASE_URL / VLLM_API_KEY が設定されていません。"
        "リポジトリ直下に .env を作成してください(.env.example を参照)。"
    )
MODEL = "google/gemma-4-26B-A4B-it"

MCP_SERVERS = {
    "pubmed": "https://pubmed.mcp.claude.com/mcp",
}

client = OpenAI(base_url=VLLM_BASE_URL, api_key=VLLM_API_KEY)

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
    convert_result = await call_tool_safe(mcp, "convert_article_ids", {"ids": [pmid], "id_type": "pmid"})
    convert_text = "\n".join(getattr(b, "text", str(b)) for b in convert_result.content)
    try:
        convert_data = json.loads(convert_text)
        pmcid = convert_data["records"][0].get("pmcid")
    except (json.JSONDecodeError, KeyError, IndexError):
        pmcid = None

    if not pmcid:
        return None, f"PMID:{pmid} is not in PMC (no PMCID). Judgment should fall back to abstract-only or be marked insufficient_data."

    ft_result = await call_tool_safe(mcp, "get_full_text_article", {"pmc_ids": [pmcid]})
    ft_text = "\n".join(getattr(b, "text", str(b)) for b in ft_result.content)
    try:
        ft_data = json.loads(ft_text)
        article = ft_data["articles"][0]
        full_text = article.get("full_text", "")
        doi = article.get("doi", "")
        if not full_text:
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
            # same as "not in PMC" so it's skipped instead of prompted.
            return None, f"PMID:{pmid} has a PMCID ({pmcid}) but no full_text came back from PubMed MCP (empty article body). Judgment should fall back to abstract-only or be marked insufficient_data."
        return full_text, f"Fetched from PubMed. PMCID={pmcid}, DOI={doi}"
    except (json.JSONDecodeError, KeyError, IndexError) as e:
        return None, f"Failed to parse the full-text fetch result: {e}"


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


def call_llm_for_judgment(prompt: str) -> tuple[dict, dict]:
    """
    Sends the prompt to gemma-4 and gets structured JSON output back.
    Returns (parsed JSON, raw API response stats).
    """
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
    }
    parsed = extract_json(msg.content or "")
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

    prompt = engine.build_prompt(gene, hgvsc, hgvsp, equivalents, full_text)
    log(f"--- Prompt (first 1000 chars) ---\n{prompt[:1000]}")

    try:
        raw_json, stats = call_llm_for_judgment(prompt)
    except Exception as e:
        show(f"  -> Error during the LLM call / JSON parsing: {e}")
        return None

    show(f"[LLM response] {stats['elapsed_sec']:.1f}s, "
         f"prompt={stats['prompt_tokens']}tok, completion={stats['completion_tokens']}tok")
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
            full_text_cache=full_text_cache,
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
            out_path = VA_SPEC_OUTPUT_DIR / f"{evidence_line['id'].replace('evline:', '')}.json"
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
