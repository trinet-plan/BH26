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
- Decide how the tested cohort was ASCERTAINED (cohort_ascertainment):
  "gene_specific" when patients were selected because of a phenotype
  characteristic OF THIS GENE'S OWN DISEASE specifically (e.g. "suspected
  Lynch syndrome", "clinically diagnosed X-linked retinitis pigmentosa",
  "hereditary breast and ovarian cancer referral criteria") - the
  denominator is a population where this gene is already a strong
  candidate cause; "broad_panel" when patients were selected via a large
  multigene/NGS panel covering many DIFFERENT, genetically heterogeneous
  conditions (e.g. "266-gene inherited retinal dystrophy panel", "hereditary
  cancer multigene panel", "all patients referred for genetic testing") and
  this gene is only one of many possible findings - this yield is diluted
  by every other gene the panel could have found instead; "unclear" if the
  text does not describe the ascertainment clearly enough to tell.
- If the paper does not report this kind of statistic (e.g. it is a
  single case report with no cohort denominator), say so explicitly - do
  not invent one.
- Return ONLY valid JSON, no markdown fences, in this exact shape:
{
  "yield_percent_stated": <number 0-100, or null>,
  "is_overall_yield": <true/false>,
  "cohort_ascertainment": "<gene_specific, broad_panel, or unclear>",
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
    cohort_ascertainment: Optional[str] = None
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
    gene: str, phenotype_description: str, *,
    preferred_pmids: tuple[str, ...] = (), max_candidates: int = 5,
) -> LiteratureYieldResult:
    """Search PubMed for gene/phenotype and ask the project's own LLM to
    extract a confirmed OVERALL diagnostic-yield statistic.

    Returns the first candidate paper whose extraction confirms
    is_overall_yield=true with a numeric percentage; LiteratureYieldResult.
    found=False (with `reason`) if nothing usable turned up, including
    when `phenotype_description` is empty (nothing to search for).

    `preferred_pmids`, when given (this variant's own ERepo evidenceLinks),
    are tried BEFORE any live search query - the same "curated citations
    first, live search only as fallback" pattern already used by
    resolve_pmids_for_variant() (PS3/BS3/PS4/BP5) and
    acmg_pipeline.pp1_segregation_search.search_family_segregation(). Even
    when a preferred paper turns out not to carry a diagnostic-yield
    statistic (ERepo's citations are usually the variant's own
    classification evidence, not a cohort-yield paper), trying it first is
    free - it just falls through to the live search below.

    [Broad-panel dilution - confirmed empirically (2026-09-25)]
      Real ground-truth PP4-MET misses (BRCA1 c.191G>A, MSH2 c.1012G>A,
      RPGR c.492G>T) all landed below PP4_MIN_POSTERIOR not because the
      search failed, but because it correctly found a genuine, correctly-
      extracted gene-specific overall-yield statistic - just diluted by a
      broad multigene/NGS panel denominator (e.g. RPGR's yield among ALL
      5201 patients on a 266-gene retinal-dystrophy panel: 4.5%, versus
      the much higher yield expected among patients whose phenotype is
      specifically suggestive of RPGR-related (X-linked) retinopathy). Two
      changes address this: (1) the query itself now tries a
      panel-excluding phrasing first ('NOT "multigene panel" NOT "gene
      panel" NOT "NGS panel"'), before falling back to the original,
      panel-inclusive phrasing; (2) the extraction prompt now also
      classifies each candidate's `cohort_ascertainment` as gene_specific
      vs broad_panel, and the loop below prefers a gene_specific hit,
      continuing to scan further candidates rather than stopping at the
      first broad_panel one - only falling back to a broad_panel result
      (disclosed as such via `cohort_ascertainment`) if no gene_specific
      candidate exists among all `max_candidates` tried.
    """
    if not phenotype_description.strip():
        return LiteratureYieldResult(found=False, reason="no_phenotype_description_available")

    from acmg_pipeline import pipeline

    async with AsyncExitStack() as stack:
        mcp = await pipeline.connect_pubmed(stack)
        pmids = list(dict.fromkeys(preferred_pmids))
        narrow = await pipeline.search_candidate_pmids(
            mcp, gene,
            disease=(
                f'{phenotype_description} diagnostic yield '
                'NOT "multigene panel" NOT "gene panel" NOT "NGS panel"'
            ),
            max_results=max_candidates,
        )
        pmids = list(dict.fromkeys([*pmids, *narrow]))
        # Always run the panel-inclusive query too, even when the narrow query
        # above already filled max_candidates on its own - confirmed (2026-09-25)
        # that skipping it whenever the narrow query alone reaches max_candidates
        # can silently crowd out a real, previously-found hit (BRCA1 c.191G>A:
        # the narrow query's own 5 candidates didn't include PMID:42738331,
        # which only the panel-inclusive query below ever found).
        searched = await pipeline.search_candidate_pmids(
            mcp, gene, disease=f"{phenotype_description} diagnostic yield genetic testing",
            max_results=max_candidates,
        )
        pmids = list(dict.fromkeys([*pmids, *searched]))
        if len(pmids) < 3:
            more = await pipeline.search_candidate_pmids(mcp, gene, disease=phenotype_description, max_results=max_candidates)
            pmids = list(dict.fromkeys([*pmids, *more]))
        # Judge more candidates than a single query's max_results, now that two
        # independent queries (narrow-first, then panel-inclusive) are merged -
        # otherwise the panel-inclusive query's own hits would never get a turn.
        pmids = pmids[:max(2 * max_candidates, len(preferred_pmids))]

        broad_panel_fallback: Optional[LiteratureYieldResult] = None
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
            ascertainment = data.get("cohort_ascertainment")
            result = LiteratureYieldResult(
                found=True,
                yield_fraction=float(percent) / 100.0,
                sample_size=data.get("sample_size"),
                denominator_description=data.get("denominator_description"),
                quote=data.get("quote"),
                pmid=pmid,
                cohort_ascertainment=ascertainment,
                reason="confirmed_overall_yield_statistic_found",
            )
            if ascertainment == "gene_specific":
                return result
            if broad_panel_fallback is None:
                broad_panel_fallback = result

        if broad_panel_fallback is not None:
            return broad_panel_fallback

    return LiteratureYieldResult(
        found=False,
        reason="no_confirmed_overall_yield_statistic_found_in_pubmed_search",
    )
