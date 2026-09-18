"""
hpo_mondo_extraction.py

Normalize ClinicalFeature.label values in ClinicalNoteExtraction to official HPO
IDs (via TogoMCP), and ClinicalNoteExtraction.diagnosis to a MONDO disease class
(via EBI OLS4). Renamed from hpo_extraction.py (2026-09-18) when MONDO
resolution was added alongside the original HPO resolution.

Primary API:
    await normalize_hpo(extraction) -> ClinicalNoteExtraction
    await resolve_diagnosis_mondo(extraction) -> ClinicalNoteExtraction

The source-grounded ClinicalFeature.label/diagnosis text is preserved in both
cases. Only hpo_id / condition_id is filled. The current shared ClinicalFeature
class has no official-HPO-label field, so the official label is used only
during candidate validation and is not stored.

[Why MONDO resolution uses OLS4, not TogoMCP - switched 2026-09-18]
  The original MONDO implementation used the same TogoMCP SPARQL
  search-then-choose pattern as HPO below, but an LLM has to author a fresh
  SPARQL query on every call, and that hit two real, reproducible query-
  construction bugs in production testing: a stray `PREFIX bif:`
  declaration (Virtuoso reserves that name and rejects the whole query,
  HTTP 400, even when bif:contains is never actually called) and
  `ORDER BY STRLEN(?label) ASC` (SQL syntax, not valid SPARQL - the correct
  form is `ORDER BY ASC(STRLEN(?label))`). EBI's OLS4 (Ontology Lookup
  Service) exposes its own MCP server (https://www.ebi.ac.uk/ols4/api/mcp)
  with a purpose-built `searchClasses` tool - a real, relevance-ranked
  ontology search endpoint, not a query language the LLM has to write
  correctly from scratch each call - so neither failure mode can occur
  here. Verified directly: `searchClasses(query=..., ontologyId="mondo")`
  returns the correct MONDO term as the first result for both
  "hypertrophic cardiomyopathy" (MONDO:0005045) and "arrhythmogenic right
  ventricular cardiomyopathy" (MONDO:0016587), without needing the
  shortest-label-first workaround the TogoMCP version needed.
"""

from __future__ import annotations

import argparse
import asyncio
import copy
import csv
import io
import json
import os
import re
from contextlib import AsyncExitStack
from dataclasses import asdict
from pathlib import Path
from typing import Any, Optional

from mcp import ClientSession
try:
    from mcp.client.streamable_http import streamable_http_client
except ImportError:
    # Same mcp-version naming difference already worked around in
    # pipeline.py's connect_pubmed() (2026-09-17 fix, applied here too).
    from mcp.client.streamable_http import streamablehttp_client as streamable_http_client
from openai import OpenAI

from acmg_pipeline.clinical_note import ClinicalFeature, ClinicalNoteExtraction


# parents[1]: this file lives at <repo root>/acmg_pipeline/hpo_mondo_extraction.py,
# so one parent up is <repo root> (where .env lives) - same fix as
# clinical_extraction.py's ROOT_DIR (2026-09-17), for the same reason: the
# pp4_pp1_bs4 branch's parents[3] only resolves correctly under a
# different, deeper checkout layout than this repo's own.
ROOT_DIR = Path(__file__).resolve().parents[1]


def _load_env_file(path: Path = ROOT_DIR / ".env") -> None:
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
MODEL = os.environ.get("VLLM_MODEL", "google/gemma-4-26B-A4B-it")
TOGOMCP_URL = os.environ.get("TOGOMCP_URL", "https://togomcp.rdfportal.org/mcp")
OLS4_MCP_URL = os.environ.get("OLS4_MCP_URL", "https://www.ebi.ac.uk/ols4/api/mcp")

if not VLLM_BASE_URL or not VLLM_API_KEY:
    raise RuntimeError(
        "VLLM_BASE_URL / VLLM_API_KEY are not configured in the repository-root .env"
    )

client = OpenAI(base_url=VLLM_BASE_URL, api_key=VLLM_API_KEY)


def to_openai_tools(mcp_tools: list[Any]) -> list[dict[str, Any]]:
    return [
        {
            "type": "function",
            "function": {
                "name": tool.name,
                "description": tool.description or "",
                # mcp.types.Tool's JSON-schema field name differs by installed
                # mcp version - camelCase inputSchema in some releases,
                # snake_case input_schema in others (confirmed: mcp==2.2.0
                # exposes input_schema only). Same class of version drift
                # already worked around for streamable_http_client's naming
                # above and in pipeline.py's connect_pubmed().
                "parameters": getattr(tool, "inputSchema", None) or getattr(tool, "input_schema", None),
            },
        }
        for tool in mcp_tools
    ]


def extract_text(result: Any) -> str:
    return "\n".join(getattr(block, "text", str(block)) for block in result.content)


def parse_sparql_csv(text: str) -> list[dict[str, str]]:
    """Parse a run_sparql CSV result with ?term/?notation/?label columns.

    HPO-only now: MONDO resolution moved to EBI OLS4's own JSON search API
    (2026-09-18, see this module's docstring), which needs no SPARQL CSV
    parsing at all.
    """
    rows: list[dict[str, str]] = []
    lines = [
        line
        for line in text.splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    if not lines:
        return rows

    try:
        reader = csv.DictReader(io.StringIO("\n".join(lines)))
        for row in reader:
            notation = (row.get("notation") or "").strip()
            label = (row.get("label") or "").strip()
            term = (row.get("term") or "").strip()
            if notation and label:
                rows.append(
                    {
                        "term": term,
                        "hpo_id": notation,
                        "hpo_label": label,
                    }
                )
    except Exception:
        pass

    return rows


async def get_ontology_mie(mcp: ClientSession) -> str:
    result = await mcp.call_tool("get_MIE_file", {"database": "ontology"})
    return extract_text(result)


async def search_hpo_candidates(
    mcp: ClientSession,
    run_sparql_tool: dict[str, Any],
    mie: str,
    expression: str,
) -> list[dict[str, str]]:
    # Keep the prompt that previously succeeded in the BH26 demo.
    system_prompt = f"""
You are resolving a clinical phenotype expression against
The Human Phenotype Ontology (HPO) using TogoMCP.

Before writing SPARQL, follow the ontology MIE below.

----- BEGIN ONTOLOGY MIE -----
{mie}
----- END ONTOLOGY MIE -----

You have one available tool: run_sparql.

Your job:
- Search the HPO graph for HPO terms relevant to the supplied clinical expression.
- Use database="ontology".
- Pin the HPO graph:
  GRAPH <http://rdfportal.org/ontology/hp>
- Retrieve:
  ?term
  ?notation
  ?label
- ?notation must be the HPO CURIE such as HP:0001250.
- Exclude obsolete terms.
- Use the label-handling rules from the MIE.
- Keep the query bounded.
- Do not search any other ontology.
- Do not invent an HPO ID from memory.
- Actually call run_sparql.

Use the verified HPO pattern from the MIE, adapted to the supplied expression.
The result should provide a small candidate set, not every matching HPO term.
"""

    question = f"""
Find HPO term candidates for this clinical expression:

{expression}
"""

    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": question},
        ],
        tools=[run_sparql_tool],
        tool_choice={"type": "function", "function": {"name": "run_sparql"}},
    )

    message = response.choices[0].message
    if not message.tool_calls:
        raise RuntimeError(f"run_sparql was not called for: {expression}")

    args = json.loads(message.tool_calls[0].function.arguments)
    args["database"] = "ontology"

    result = await mcp.call_tool("run_sparql", args)
    return parse_sparql_csv(extract_text(result))


def choose_best_hpo(
    expression: str,
    candidates: list[dict[str, str]],
) -> Optional[dict[str, str]]:
    if not candidates:
        return None

    candidate_text = "\n".join(
        f"{index}. {item['hpo_id']}|||{item['hpo_label']}"
        for index, item in enumerate(candidates, start=1)
    )

    system_prompt = """
Select the best Human Phenotype Ontology match for a clinical expression.

You MUST choose only from the supplied candidates.

Rules:
- Do not create a new HPO label.
- Do not create or change an HPO ID.
- If none of the candidates adequately represents the clinical expression,
  return exactly: NOT_FOUND
- Otherwise return exactly one candidate in this format:
  HPO_ID|||HPO_LABEL
- Do not add explanations.
"""

    question = f"""
Clinical expression:

{expression}

HPO candidates returned by TogoMCP:

{candidate_text}
"""

    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": question},
        ],
    )

    answer = (response.choices[0].message.content or "").strip()
    if answer == "NOT_FOUND":
        return None

    for item in candidates:
        if answer == f"{item['hpo_id']}|||{item['hpo_label']}":
            return item
    return None


async def search_mondo_candidates(
    mcp: ClientSession,
    expression: str,
    *,
    page_size: int = 8,
) -> list[dict[str, str]]:
    """Search EBI OLS4's MONDO index for disease class candidates matching a
    free-text diagnosis, via OLS4's own searchClasses tool (real network
    call, no SPARQL an LLM has to author - see this module's own docstring
    for why that replaced the earlier TogoMCP approach).
    """
    result = await mcp.call_tool(
        "searchClasses", {"query": expression, "ontologyId": "mondo", "pageSize": page_size},
    )
    try:
        data = json.loads(extract_text(result))
    except (json.JSONDecodeError, ValueError):
        return []

    candidates: list[dict[str, str]] = []
    seen_ids: set[str] = set()
    for item in data.get("items") or []:
        curie = item.get("curie")
        labels = item.get("label") or []
        if not curie or not curie.startswith("MONDO:") or not labels or curie in seen_ids:
            continue
        seen_ids.add(curie)
        candidates.append({"mondo_id": curie, "mondo_label": labels[0]})
    return candidates


# Splits a compound diagnosis string ("hypertrophic cardiomyopathy (HCM)
# complicated by a left ventricular apical aneurysm") at its first
# qualifier/complication clause or parenthetical, for _primary_diagnosis_clause()
# below.
_DIAGNOSIS_QUALIFIER_SPLIT = re.compile(
    r"\s*\(|\s+(?:complicated by|with|due to|secondary to|associated with)\s+",
    re.IGNORECASE,
)


def _primary_diagnosis_clause(diagnosis: str) -> Optional[str]:
    """The leading clause of a compound diagnosis string, or None if it is
    already a single clause.

    OLS4's searchClasses does phrase-relevance matching, not the free
    substring match TogoMCP's SPARQL CONTAINS() did - confirmed empirically:
    the full string "hypertrophic cardiomyopathy (HCM) complicated by a
    left ventricular apical aneurysm" returns zero candidates, even though
    "hypertrophic cardiomyopathy" alone finds the right one immediately.
    This is used only to build a broader FOLLOW-UP search query when the
    first search returns nothing - never assigned to diagnosis/condition_id
    directly, and choose_best_mondo() is still shown the full original
    diagnosis text for its final pick, so a bad split can only mean fewer
    candidates are found, never a wrong disease getting chosen.
    """
    head = _DIAGNOSIS_QUALIFIER_SPLIT.split(diagnosis, maxsplit=1)[0].strip()
    return head if head and head != diagnosis else None


def choose_best_mondo(
    expression: str,
    candidates: list[dict[str, str]],
) -> Optional[dict[str, str]]:
    if not candidates:
        return None

    candidate_text = "\n".join(
        f"{index}. {item['mondo_id']}|||{item['mondo_label']}"
        for index, item in enumerate(candidates, start=1)
    )

    system_prompt = """
Select the best Mondo Disease Ontology (MONDO) match for a clinical
diagnosis expression.

You MUST choose only from the supplied candidates.

Rules:
- Do not create a new MONDO label.
- Do not create or change a MONDO ID.
- If none of the candidates adequately represents the diagnosis,
  return exactly: NOT_FOUND
- Otherwise return exactly one candidate in this format:
  MONDO_ID|||MONDO_LABEL
- Do not add explanations.
"""

    question = f"""
Clinical diagnosis:

{expression}

MONDO candidates returned by EBI OLS4:

{candidate_text}
"""

    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": question},
        ],
    )

    answer = (response.choices[0].message.content or "").strip()
    if answer == "NOT_FOUND":
        return None

    for item in candidates:
        if answer == f"{item['mondo_id']}|||{item['mondo_label']}":
            return item
    return None


async def resolve_diagnosis_mondo(
    extraction: ClinicalNoteExtraction,
) -> ClinicalNoteExtraction:
    """
    Return a deep-copied ClinicalNoteExtraction with condition_id filled,
    resolving extraction.diagnosis (free text) to a MONDO disease class via
    EBI OLS4 - the same search-then-choose shape normalize_hpo() uses for
    clinical features (real candidates only, a final pick validated to
    exactly match one of them - see choose_best_mondo()), applied to the
    diagnosis field instead.

    This is the ONLY producer of ClinicalNoteExtraction.condition_id
    (2026-09-18): clinical_extraction.py's LLM used to also guess this
    field directly from the note text with no grounding at all beyond a
    regex format check, which produced a confirmed, reproducible wrong
    answer in real testing (a demo case's diagnosis "arrhythmogenic right
    ventricular cardiomyopathy" came back condition_id=MONDO:0005045, which
    is hypertrophic cardiomyopathy - a different disease). That prompt
    instruction was removed; this function is now the only path that can
    set condition_id, and it can only return a real, existing MONDO term
    that OLS4 itself returned as a candidate, never a fabricated one.

    extraction.diagnosis is never replaced or overwritten; only condition_id
    is added.
    """
    result = copy.deepcopy(extraction)
    diagnosis = (result.diagnosis or "").strip()
    if not diagnosis:
        return result

    best: Optional[dict[str, str]] = None

    async with AsyncExitStack() as stack:
        ctx = await stack.enter_async_context(streamable_http_client(OLS4_MCP_URL))
        read, write = ctx[0], ctx[1]
        mcp = await stack.enter_async_context(ClientSession(read, write))
        await mcp.initialize()
        tools = (await mcp.list_tools()).tools
        tool_names = {tool.name for tool in tools}
        if "searchClasses" not in tool_names:
            raise RuntimeError("Required OLS4 tool 'searchClasses' is missing")

        print(f"[MONDO] {diagnosis}")
        candidates = await search_mondo_candidates(mcp, diagnosis)
        if not candidates:
            # OLS4 does phrase-relevance matching, not a substring search -
            # a compound diagnosis ("X complicated by Y") can return zero
            # candidates for the full text even though its primary clause
            # alone finds the right term. Broaden the QUERY only; the final
            # choose_best_mondo() call below still sees the full, original
            # diagnosis text (see _primary_diagnosis_clause()'s own
            # docstring for why this cannot introduce a wrong disease).
            fallback_query = _primary_diagnosis_clause(diagnosis)
            if fallback_query:
                print(f"        no candidates for full text, retrying with: {fallback_query}")
                candidates = await search_mondo_candidates(mcp, fallback_query)
        # choose_best_mondo() is an LLM call, not a deterministic lookup -
        # confirmed empirically (2026-09-18) that a borderline diagnosis
        # (a specific sub-phenotype whose exact MONDO term isn't among the
        # candidates, e.g. "hypertrophic cardiomyopathy with apical
        # ventricular aneurysm" vs. plain "hypertrophic cardiomyopathy")
        # returns NOT_FOUND on roughly 1 in 5 calls, purely from sampling
        # variance, even though the same candidates are judged sufficient
        # on the other 4. A missing condition_id here cascades into PVS1's
        # (and other criteria's) mechanism gate going unevaluated, so a
        # few independent retries meaningfully reduces the odds a real,
        # available candidate gets lost to a single unlucky call - this
        # does not change which candidates are offered, just how many
        # chances the same judgment gets to land on one of them.
        best = None
        for attempt in range(3):
            best = choose_best_mondo(diagnosis, candidates)
            if best is not None:
                break
        if best is None:
            print("        selected: NOT_FOUND")
        else:
            print(f"        selected: {best['mondo_id']} {best['mondo_label']}")

    result.condition_id = best["mondo_id"] if best else None
    return result


def _all_features(extraction: ClinicalNoteExtraction) -> list[ClinicalFeature]:
    features: list[ClinicalFeature] = []
    features.extend(extraction.proband.phenotype.clinical_features)
    for relative in extraction.family.relatives:
        features.extend(relative.clinical_features)
    return features


async def resolve_hpo_labels(
    labels: list[str],
) -> dict[str, Optional[dict[str, str]]]:
    unique_labels: list[str] = []
    seen: set[str] = set()

    for label in labels:
        cleaned = (label or "").strip()
        if not cleaned:
            continue
        key = cleaned.casefold()
        if key in seen:
            continue
        seen.add(key)
        unique_labels.append(cleaned)

    if not unique_labels:
        return {}

    async with AsyncExitStack() as stack:
        # The number of values streamable_http_client() yields differs by
        # mcp version (2 vs 3 - see acmg_pipeline.pipeline.connect_pubmed()'s
        # own docstring for the same issue, hit there first); take only the
        # first two (read, write) regardless, rather than unpacking directly.
        ctx = await stack.enter_async_context(streamable_http_client(TOGOMCP_URL))
        read, write = ctx[0], ctx[1]
        mcp = await stack.enter_async_context(ClientSession(read, write))
        await mcp.initialize()
        tools = (await mcp.list_tools()).tools
        tool_names = {tool.name for tool in tools}
        required = {"get_MIE_file", "run_sparql"}
        missing = required - tool_names
        if missing:
            raise RuntimeError(
                "Required TogoMCP tools are missing: " + ", ".join(sorted(missing))
            )

        run_sparql_mcp_tool = next(tool for tool in tools if tool.name == "run_sparql")
        run_sparql_tool = to_openai_tools([run_sparql_mcp_tool])[0]
        mie = await get_ontology_mie(mcp)

        resolved: dict[str, Optional[dict[str, str]]] = {}
        for index, label in enumerate(unique_labels, start=1):
            print(f"[HPO] [{index}/{len(unique_labels)}] {label}")
            candidates = await search_hpo_candidates(
                mcp=mcp,
                run_sparql_tool=run_sparql_tool,
                mie=mie,
                expression=label,
            )
            best = choose_best_hpo(label, candidates)
            if best is None:
                print("      selected: NOT_FOUND")
            else:
                print(f"      selected: {best['hpo_id']} {best['hpo_label']}")
            resolved[label.casefold()] = best

        return resolved


async def normalize_hpo(
    extraction: ClinicalNoteExtraction,
) -> ClinicalNoteExtraction:
    """
    Return a deep-copied ClinicalNoteExtraction with ClinicalFeature.hpo_id filled.

    The original source-grounded ClinicalFeature.label is never replaced.
    """
    result = copy.deepcopy(extraction)
    features = _all_features(result)
    resolved = await resolve_hpo_labels([feature.label for feature in features])

    for feature in features:
        best = resolved.get(feature.label.casefold())
        feature.hpo_id = best["hpo_id"] if best else None

    return result


def proband_hpo_ids(extraction: ClinicalNoteExtraction) -> list[str]:
    """Convenience helper for the later PP4 phenotype-matching step."""
    return [
        feature.hpo_id
        for feature in extraction.proband.phenotype.clinical_features
        if feature.hpo_id is not None
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize ClinicalNoteExtraction clinical features to official HPO IDs."
    )
    parser.add_argument(
        "clinical_extraction_json",
        type=Path,
        help="Debug JSON representation of ClinicalNoteExtraction.",
    )
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw = json.loads(args.clinical_extraction_json.read_text(encoding="utf-8"))
    extraction = ClinicalNoteExtraction.from_json(raw)
    normalized = asyncio.run(normalize_hpo(extraction))

    payload = json.dumps(asdict(normalized), ensure_ascii=False, indent=2)
    print(payload)

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
        print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
