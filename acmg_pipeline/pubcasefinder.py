"""
pubcasefinder.py

Automatic PP4 phenotype-specificity signal via PubCaseFinder (through TogoMCP),
per Layer 2 ("表現型データが必要") of the BH26 ACMG criteria definition v2
(2026-09-14): "臨床記述からのHPO抽出→疾患との一致判定が必要。PubCaseFinder等を
想定". This replaces the need to hand-curate a required-HPO-term list per gene
(PP4ReferenceRecord.phenotype_hpo) for every config/pp4_reference_records.json
entry - the patient's own normalized HPO profile (acmg_pipeline.hpo_extraction)
is enough to ask PubCaseFinder directly whether this gene is the best
phenotype match.

Primary API:
    await rank_genes_by_phenotype(hpo_ids) -> dict   (raw PubCaseFinder JSON)
    phenotype_match_from_gene_ranking(ranking, gene_symbol) -> PhenotypeMatchResult

What this does and does not decide
-----------------------------------
PubCaseFinder returns an information-content-weighted phenotype-SIMILARITY
RANKING over genes, not a calibrated diagnostic-yield/posterior probability.
It is therefore used only to answer PP4's "is the patient's phenotype specific
to this locus" gate (phenotype_match) - never to invent a diagnostic_yield /
ClinGen 2024 Bayesian-points number. That number stays a literature-curated
fact in config/pp4_reference_records.json (see
acmg_pipeline.criteria.evaluator.DIAGNOSTIC_YIELD_POINT_TABLE); PubCaseFinder's
own tool description is explicit that it is "decision support over
annotations, NOT a diagnosis".

[Rank-1-only, ties included, is a provisional design choice]
  matched=True requires the candidate gene to be *among* PubCaseFinder's
  top-ranked genes for the patient's HPO profile (rank == 1; ties with other
  genes are allowed, since real disease loci are frequently heterogeneous and
  that heterogeneity is already captured by reference.locus_model /
  diagnostic_yield, not by this gate). Any lower rank counts as matched=False.
  This interpretation has not been reviewed by whoever curates
  pp4_reference_records.json going forward - flagged here the same way
  acmg_pipeline.criteria.pp1_pp4_strength_table flags its own tie-break
  rule for Table 4.
"""

from __future__ import annotations

import json
import os
from contextlib import AsyncExitStack
from pathlib import Path
from typing import Optional, Sequence

from mcp import ClientSession
try:
    from mcp.client.streamable_http import streamable_http_client
except ImportError:
    # Same mcp-version naming difference already worked around in
    # pipeline.py's connect_pubmed() and hpo_extraction.py's resolve_hpo_labels().
    from mcp.client.streamable_http import streamablehttp_client as streamable_http_client

from acmg_pipeline.criteria.pp4_pp1_bs4 import PatientHpoTerm, PhenotypeMatchResult

# parents[1]: this file lives at <repo root>/acmg_pipeline/pubcasefinder.py,
# so one parent up is <repo root> (where .env lives) - same convention as
# hpo_extraction.py's ROOT_DIR.
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

TOGOMCP_URL = os.environ.get("TOGOMCP_URL", "https://togomcp.rdfportal.org/mcp")
PUBCASEFINDER_RANK_TOOL = "pubcasefinder_rank_by_phenotypes"
_MAX_LIMIT = 100


def _extract_text(result) -> str:
    return "\n".join(getattr(block, "text", str(block)) for block in result.content)


async def rank_genes_by_phenotype(hpo_ids: Sequence[str], *, limit: int = _MAX_LIMIT) -> dict:
    """Ask PubCaseFinder (via TogoMCP) to rank candidate genes for these HPO IDs.

    Returns the tool's parsed JSON object (rank/score/gene_symbol/matched_hpo_ids/...
    per row) as-is. Propagates whatever the tool itself raises unchanged:
    ValueError on a malformed HPO ID, when none of the supplied IDs is
    recognized, when the shared PubCaseFinder rate budget is exhausted, or on
    an upstream HTTP error - callers must treat those as "not evaluable", not
    as a negative phenotype-match observation.
    """
    hpo_ids = list(hpo_ids)
    if not hpo_ids:
        raise ValueError("hpo_ids must not be empty")

    async with AsyncExitStack() as stack:
        # AsyncExitStack + take-first-two, not direct unpacking: the number of
        # values streamable_http_client() yields differs by installed mcp
        # version (same issue already hit in hpo_extraction.resolve_hpo_labels()).
        ctx = await stack.enter_async_context(streamable_http_client(TOGOMCP_URL))
        read, write = ctx[0], ctx[1]
        mcp = await stack.enter_async_context(ClientSession(read, write))
        await mcp.initialize()

        tools = (await mcp.list_tools()).tools
        if not any(tool.name == PUBCASEFINDER_RANK_TOOL for tool in tools):
            raise RuntimeError(f"TogoMCP does not expose {PUBCASEFINDER_RANK_TOOL!r}")

        result = await mcp.call_tool(
            PUBCASEFINDER_RANK_TOOL,
            {"hpo_ids": hpo_ids, "target": "gene", "limit": limit},
        )
        return json.loads(_extract_text(result))


def phenotype_match_from_gene_ranking(
    ranking: dict,
    gene_symbol: str,
    *,
    patient_terms: Sequence[PatientHpoTerm] = (),
) -> PhenotypeMatchResult:
    """Decide PP4's phenotype_match gate from a PubCaseFinder gene ranking.

    matched=True  - gene_symbol is one of the rank-1 (best-matching) genes for
                     the patient's HPO profile.
    matched=False - gene_symbol was ranked, but not at rank 1.
    matched=None  - gene_symbol never appears among PubCaseFinder's ranked
                     candidates for this phenotype set (e.g. it has no HPO
                     annotation in PubCaseFinder's data, or every candidate
                     scored fell outside ``limit``) - a real "cannot evaluate",
                     not a fabricated non-match.
    """
    results = ranking.get("results") or []
    gene_row = next((row for row in results if row.get("gene_symbol") == gene_symbol), None)

    if gene_row is None:
        return PhenotypeMatchResult(
            matched=None,
            patient_terms=list(patient_terms),
            reason="gene_not_present_in_pubcasefinder_ranking",
        )

    matched = gene_row.get("rank") == 1
    return PhenotypeMatchResult(
        matched=matched,
        patient_terms=list(patient_terms),
        reason=(
            "pubcasefinder_gene_top_ranked_for_phenotype"
            if matched
            else "pubcasefinder_gene_not_top_ranked_for_phenotype"
        ),
    )
