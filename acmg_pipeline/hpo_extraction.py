"""
hpo_extraction.py

Normalize ClinicalFeature.label values in ClinicalNoteExtraction to official HPO
IDs using TogoMCP.

Primary API:
    await normalize_hpo(extraction) -> ClinicalNoteExtraction

The source-grounded ClinicalFeature.label is preserved. Only hpo_id is filled.
The current shared ClinicalFeature class has no official-HPO-label field, so the
official label is used only during candidate validation and is not stored.
"""

from __future__ import annotations

import argparse
import asyncio
import copy
import csv
import io
import json
import os
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


# parents[1]: this file lives at <repo root>/acmg_pipeline/hpo_extraction.py,
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
                "parameters": tool.inputSchema,
            },
        }
        for tool in mcp_tools
    ]


def extract_text(result: Any) -> str:
    return "\n".join(getattr(block, "text", str(block)) for block in result.content)


def parse_sparql_csv(text: str) -> list[dict[str, str]]:
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
