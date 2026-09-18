"""
condition_from_literature_search.py

Curator-support literature search for a variant's disease/condition context,
for callers that have no ERepo record and no free-text clinical note to draw
a diagnosis from at all (e.g. run_integrated_validation_64.py's demo-case-
sourced variants, such as MYBPC3 c.278delA/c.2905+1G>A/c.836del - real,
well-established HCM variants that simply are not curated in ERepo). Mirrors
acmg_pipeline.pp4_literature_search's and .pp1_segregation_search's design:
live PubMed search + this project's own LLM, no persistence.

[Why this doesn't just reuse pp4_literature_search's query]
  PP4's search asks "what fraction of ALL patients with this phenotype carry
  a pathogenic variant in this gene" - it already assumes a phenotype to
  search FOR. This module answers a different, prior question: "what disease
  does the literature report THIS variant causing," so its query is built
  from gene+variant alone, with no phenotype input required.

[Why the found condition still goes through resolve_diagnosis_mondo()]
  This module only extracts a free-text disease NAME from a paper (e.g.
  "hypertrophic cardiomyopathy") - never a MONDO ID. Grounding that name to a
  real MONDO term (so PVS1's mechanism gate and PP1/BS4/PP4 can use it) is
  acmg_pipeline.hpo_mondo_extraction.resolve_diagnosis_mondo()'s job, exactly
  as it already is for a diagnosis pulled from a real clinical note - kept
  as a separate step here rather than duplicated, so there is only one place
  that can ever produce a condition_id.
"""

from __future__ import annotations

from contextlib import AsyncExitStack
from dataclasses import dataclass
from typing import Optional

EXTRACTION_SYSTEM_PROMPT = """\
You are a clinical genetics literature-review assistant helping a human
curator identify what disease/condition a specific genetic variant is
reported to cause.

Given the text of one paper and a target variant, decide whether the paper
identifies this SPECIFIC variant and states what disease/condition it is
associated with or reported to cause.

Rules:
- Only extract a condition name EXPLICITLY stated as associated with THIS
  variant - never a disease merely associated with the gene in general if
  the variant itself is not discussed, and never inferred from what would
  be "typical" for the gene.
- Report the condition using its full disease name, not an abbreviation or
  acronym (e.g. "hypertrophic cardiomyopathy", not "HCM") - this name is
  looked up against a disease ontology afterwards, which indexes full names,
  not acronyms. Expand only a standard, unambiguous abbreviation the paper
  itself defines or that is common medical usage; never guess an expansion
  you are not confident in.
- If the paper does not mention this variant at all, or mentions it without
  stating an associated disease, say so explicitly.
- Return ONLY valid JSON, no markdown fences, in this exact shape:
{
  "variant_condition_described": <true/false>,
  "condition_name": "<string, or null>",
  "quote": "<the sentence(s) the condition name came from, or null>",
  "no_condition_found": <true/false>
}
"""


@dataclass
class LiteratureConditionResult:
    found: bool
    condition_name: Optional[str] = None
    pmid: Optional[str] = None
    quote: Optional[str] = None
    reason: str = ""


async def _judge_paper_for_condition(pipeline, gene: str, hgvsc: str, text: str, pmid: str) -> dict:
    response = pipeline.client.chat.completions.create(
        model=pipeline.MODEL,
        messages=[
            {"role": "system", "content": EXTRACTION_SYSTEM_PROMPT},
            {"role": "user", "content": (
                f"Target variant: {gene} {hgvsc}\n\n"
                f"----- BEGIN PAPER TEXT (PMID:{pmid}) -----\n{text}\n----- END PAPER TEXT -----"
            )},
        ],
    )
    raw = response.choices[0].message.content or ""
    try:
        return pipeline.extract_json(raw)
    except Exception:
        return {"no_condition_found": True}


async def search_condition_from_literature(
    gene: str, hgvsc: str, *, hgvsp: str | None = None, max_candidates: int = 8,
) -> LiteratureConditionResult:
    """Search PubMed for this variant and ask the project's own LLM to
    extract the disease/condition it reports the variant as causing.

    Returns the first candidate paper whose extraction confirms
    variant_condition_described=true with a non-empty condition_name;
    LiteratureConditionResult.found=False (with `reason`) if nothing usable
    turned up.

    [Query strategy - a known, real limitation, not an oversight]
      acmg_pipeline.pipeline.search_candidate_pmids() itself documents that
      exact HGVS c./p. notation (e.g. "MYBPC3 AND c.278delA") returns ZERO
      PubMed hits even for a variant whose real source paper IS indexed -
      its own best fallback for a variant with no known phenotype yet is a
      gene-only query ("<gene> AND novel variant"), which is broad and
      relies on the per-paper LLM judge below to reject every paper that
      does not actually discuss this exact variant. hgvsp, when the caller
      has it (a missense change reads far better, e.g. "MYH7 AND
      Arg719Trp"), is passed through and improves the query; for a
      frameshift/indel with no hgvsp this can legitimately find nothing,
      which is reported honestly as not found rather than guessed.
    """
    if not gene or not hgvsc:
        return LiteratureConditionResult(found=False, reason="no_gene_or_hgvsc_available")

    from acmg_pipeline import pipeline

    async with AsyncExitStack() as stack:
        mcp = await pipeline.connect_pubmed(stack)
        pmids = await pipeline.search_candidate_pmids(mcp, gene, hgvsp=hgvsp, max_results=max_candidates)

        for pmid in pmids:
            text, _note = await pipeline.fetch_full_text(mcp, pmid)
            if not text:
                continue
            data = await _judge_paper_for_condition(pipeline, gene, hgvsc, text, pmid)
            if data.get("no_condition_found") or not data.get("variant_condition_described"):
                continue
            condition_name = data.get("condition_name")
            if not condition_name or not str(condition_name).strip():
                continue
            return LiteratureConditionResult(
                found=True,
                condition_name=str(condition_name).strip(),
                pmid=pmid,
                quote=data.get("quote"),
                reason="confirmed_variant_condition_found",
            )

    return LiteratureConditionResult(
        found=False,
        reason="no_confirmed_variant_condition_found_in_pubmed_search",
    )
