"""
pp4_literature_search.py

Curator-support literature search for PP4's diagnostic_yield input, per
the user's explicit direction (2026-09-17): search PubMed (the same MCP
server acmg_pipeline.pipeline already uses for PS3/BS3/PS4) and let the
project's own LLM (vLLM gemma-4, same client pipeline.py already
configures) extract a candidate diagnostic-yield statistic, rather than
requiring a human to have already typed the number into
config/pp4_reference_records.json before PP4 can produce anything at all.

[Why this doesn't need a separate approval gate before it can be used]
  Confirmed with the user (2026-09-17): this project has no "AI drafts,
  a human approves before use" principle. PS3/BS3/PS4 already show the
  actual pattern - a live LLM judgment against real literature becomes
  the CriterionEvidence (MET/NOT_MET/strength) immediately, carrying a
  disclosure ("heuristic ... requires human curator confirmation") rather
  than being withheld pending separate approval. A literature-search-
  derived PP4 entry follows the same pattern - see
  acmg_pipeline.criteria.pp1_bs4_pp4_engine.evaluate() and
  acmg_pipeline.export._strength_blocks()'s own disclosure message (folded
  into curatorHints, 2026-09-18) for the precedent.

[Why the denominator is checked explicitly, not just "a percentage"]
  Two real experiments (case3/MYH7 against both a GeneReviews entry and a
  PubMed search hit) found the exact same trap: a gene's percentage is
  very often reported as "share among patients where SOME gene was
  already found positive," not "share of ALL patients with this
  phenotype." Feeding the former into evaluator.evaluate_pp4() as
  diagnostic_yield would badly overstate the evidence (confirmed:
  GeneReviews' MYH7=33%-of-solved-cases understates the true ~10%
  overall yield for sporadic presentations by more than 3x). The
  extraction prompt below requires the model to explicitly confirm a
  percentage is against ALL tested/reviewed patients with the phenotype
  (`is_overall_yield`); this module discards any percentage not
  explicitly confirmed as such, rather than guessing which denominator
  applies.

[Why this connects to the PubMed MCP server, not TogoMCP]
  acmg_pipeline.hpo_mondo_extraction opens its own TogoMCP connection (used
  for HPO/MONDO normalization); this module instead reuses acmg_pipeline.pipeline's
  existing connect_pubmed() / search_candidate_pmids() / fetch_full_text()
  / extract_json() / client / MODEL - the
  same building blocks the PS3/BS3/PS4 literature workflow already uses -
  rather than duplicating a second PubMed client. The import is lazy
  (inside search_diagnostic_yield(), not at module load) because
  acmg_pipeline.pipeline requires VLLM_BASE_URL/VLLM_API_KEY at import
  time, the same reason acmg_pipeline.clinical_note lazily imports
  clinical_extraction.
"""

from __future__ import annotations

from contextlib import AsyncExitStack
from dataclasses import dataclass
from typing import Optional

EXTRACTION_SYSTEM_PROMPT = """\
You are a clinical genetics literature-review assistant helping a human
curator find candidate PP4 (ACMG/AMP) diagnostic-yield evidence.

Given the text of one paper, decide whether it reports a DIAGNOSTIC YIELD
statistic for the described gene-phenotype pair: i.e. "of X patients
tested/reviewed with this phenotype, Y were found to have a pathogenic
variant in this gene" (a fraction/percentage with a clear denominator).

Rules:
- Only extract a number that is EXPLICITLY stated in the text. Never
  compute, estimate, or infer a percentage that is not written down.
- Decide whether the percentage is computed against ALL patients tested/
  reviewed for this phenotype (is_overall_yield=true), or only against a
  narrower group such as patients already known to carry SOME pathogenic
  variant (is_overall_yield=false) - e.g. "MYH7 causes 33% of genetically
  solved HCM cases" is NOT an overall yield (it excludes patients with no
  variant found at all); "36.4% of the 225 patients tested had a
  pathogenic MYH7 variant" IS an overall yield. When genuinely unsure,
  set is_overall_yield=false rather than guessing true.
- If the paper does not report this kind of statistic (e.g. it is a
  single case report with no cohort denominator), say so explicitly - do
  not invent one.
- Return ONLY valid JSON, no markdown fences, in this exact shape:
{
  "yield_percent_stated": <number 0-100, or null>,
  "is_overall_yield": <true/false>,
  "denominator_description": "<what population this percentage is over, or null>",
  "sample_size": <integer, or null>,
  "quote": "<the exact sentence(s) the number/statement came from, or null>",
  "no_yield_statistic_found": <true/false>
}
"""


@dataclass
class LiteratureYieldResult:
    found: bool
    yield_fraction: Optional[float] = None
    sample_size: Optional[int] = None
    denominator_description: Optional[str] = None
    quote: Optional[str] = None
    pmid: Optional[str] = None
    reason: str = ""


async def _judge_paper_for_yield(pipeline, gene: str, phenotype_description: str, text: str, pmid: str) -> dict:
    response = pipeline.client.chat.completions.create(
        model=pipeline.MODEL,
        messages=[
            {"role": "system", "content": EXTRACTION_SYSTEM_PROMPT},
            {"role": "user", "content": (
                f"Gene: {gene}\nPhenotype of interest: {phenotype_description}\n\n"
                f"----- BEGIN PAPER TEXT (PMID:{pmid}) -----\n{text}\n----- END PAPER TEXT -----"
            )},
        ],
    )
    raw = response.choices[0].message.content or ""
    try:
        return pipeline.extract_json(raw)
    except Exception:
        return {"no_yield_statistic_found": True}


async def search_diagnostic_yield(
    gene: str, phenotype_description: str, *, max_candidates: int = 5,
) -> LiteratureYieldResult:
    """Search PubMed for gene/phenotype and ask the project's own LLM to
    extract a confirmed OVERALL diagnostic-yield statistic.

    Returns the first candidate paper whose extraction confirms
    is_overall_yield=true with a numeric percentage; LiteratureYieldResult.
    found=False (with `reason`) if nothing usable turned up, including
    when `phenotype_description` is empty (nothing to search for).
    """
    if not phenotype_description.strip():
        return LiteratureYieldResult(found=False, reason="no_phenotype_description_available")

    from acmg_pipeline import pipeline

    async with AsyncExitStack() as stack:
        mcp = await pipeline.connect_pubmed(stack)
        pmids = await pipeline.search_candidate_pmids(
            mcp, gene, disease=f"{phenotype_description} diagnostic yield genetic testing",
            max_results=max_candidates,
        )
        if len(pmids) < 3:
            more = await pipeline.search_candidate_pmids(mcp, gene, disease=phenotype_description, max_results=max_candidates)
            pmids = list(dict.fromkeys(pmids + more))[:max_candidates]

        for pmid in pmids:
            text, _note = await pipeline.fetch_full_text(mcp, pmid)
            if not text:
                continue
            data = await _judge_paper_for_yield(pipeline, gene, phenotype_description, text, pmid)
            if data.get("no_yield_statistic_found") or not data.get("is_overall_yield"):
                continue
            percent = data.get("yield_percent_stated")
            if not isinstance(percent, (int, float)) or not 0 <= percent <= 100:
                continue
            return LiteratureYieldResult(
                found=True,
                yield_fraction=float(percent) / 100.0,
                sample_size=data.get("sample_size"),
                denominator_description=data.get("denominator_description"),
                quote=data.get("quote"),
                pmid=pmid,
                reason="confirmed_overall_yield_statistic_found",
            )

    return LiteratureYieldResult(
        found=False,
        reason="no_confirmed_overall_yield_statistic_found_in_pubmed_search",
    )
