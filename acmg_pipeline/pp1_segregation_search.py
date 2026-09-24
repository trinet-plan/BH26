"""
pp1_segregation_search.py

Curator-support literature search for PP1/BS4's family-segregation input,
mirroring acmg_pipeline.pp4_literature_search's design (2026-09-18): search
PubMed (the same MCP server acmg_pipeline.pipeline already uses for
PS3/BS3/PS4) and let the project's own LLM extract a candidate family
pedigree/segregation description, rather than requiring
clinical_note.family.relatives to already be populated before PP1/BS4 can
produce anything at all.

[Why this is a separate search from PP4's]
  PP4's search asks "what fraction of ALL patients with this phenotype
  carry a pathogenic variant in this gene" - a cohort statistic. PP1/BS4
  need something structurally different: a paper that reports THIS
  variant segregating (or not) through a specific family - each relative's
  affected status and whether they carry the variant. Checked directly
  against the three papers already cached for today's PP4-MET subset
  (PMID 35028538, 36274938, 32849802): none contain pedigree/segregation
  language at all - they are the same case-control/functional papers
  ERepo cites for the variant's overall classification, not family
  studies. A PP4-oriented search query would not reliably surface a
  family study, so this module runs its own PubMed search biased toward
  "family"/"pedigree"/"segregation" terms instead of reusing PP4's query.

[No curated registry, same reasoning as PP4]
  Same precedent as acmg_pipeline.pp4_literature_search: PS3/BS3/PS4/PP4
  already show this project's real pattern (a live LLM judgment against
  real literature becomes the CriterionEvidence immediately, carrying a
  disclosure, not withheld pending human approval). A family study is
  rare enough (most papers ERepo already cites don't have one) that
  repeating the search per call is an acceptable cost.

[Only what acmg_pipeline.criteria.pp4_pp1_bs4._score_family_segregation()
 actually reads]
  The extraction only asks for what that scorer consumes: per relative,
  relationship / affected_status / variant_status / zygosity, plus the
  family's overall inheritance_pattern and (for autosomal recessive)
  ar_case_mode. fully_penetrant/low_phenocopy are deliberately NOT
  extracted here - those are conservative judgment calls this project
  already treats as "leave as None unless a human curator establishes
  them" (see evaluate_locus_evidence()'s own docstring), not something to
  infer from one paper's wording.
"""

from __future__ import annotations

from contextlib import AsyncExitStack
from dataclasses import dataclass, field
from typing import Optional

EXTRACTION_SYSTEM_PROMPT = """\
You are a clinical genetics literature-review assistant helping a human
curator find candidate PP1/BS4 (ACMG/AMP family co-segregation) evidence.

Given the text of one paper and a target variant, decide whether it
describes a SPECIFIC FAMILY (pedigree) in which the target variant was
tested in multiple relatives - i.e. for at least one relative besides the
proband, both (a) whether they are affected/unaffected and (b) whether
they carry the target variant are stated.

Rules:
- Only extract facts EXPLICITLY stated in the text. Never infer a
  relative's status from population statistics or from what would be
  "typical" for the disease.
- If the paper is a case report/cohort study that does not test relatives
  for the variant (e.g. a solo proband, or a cohort paper with no
  per-family segregation data), say so explicitly - do not invent one.
- inheritance_pattern should be one of "autosomal dominant",
  "autosomal recessive", "x-linked recessive", or null if not stated/
  determinable for this gene-family.
- ar_case_mode applies only when inheritance_pattern is "autosomal
  recessive": "homozygous" or "compound_heterozygous" for the proband, or
  null if not applicable/not stated.
- Each relative's relationship should be relative to the proband (e.g.
  "father", "mother", "sibling", "son", "daughter", "aunt"), affected_status
  and variant_status as true/false/null (null if not stated for that
  relative), zygosity as "heterozygous"/"homozygous"/null.
- Return ONLY valid JSON, no markdown fences, in this exact shape:
{
  "family_segregation_described": <true/false>,
  "inheritance_pattern": "<string, or null>",
  "ar_case_mode": "<string, or null>",
  "relatives": [
    {"relationship": "<string>", "affected_status": <bool/null>,
     "variant_status": <bool/null>, "zygosity": "<string, or null>"}
  ],
  "quote": "<the sentence(s) the relatives' data came from, or null>",
  "no_segregation_data_found": <true/false>
}
"""


@dataclass
class SegregationRelative:
    relationship: str
    affected_status: Optional[bool] = None
    variant_status: Optional[bool] = None
    zygosity: Optional[str] = None


@dataclass
class LiteratureSegregationResult:
    found: bool
    inheritance_pattern: Optional[str] = None
    ar_case_mode: Optional[str] = None
    relatives: list[SegregationRelative] = field(default_factory=list)
    pmid: Optional[str] = None
    quote: Optional[str] = None
    reason: str = ""


async def _judge_paper_for_segregation(pipeline, gene: str, hgvsc: str, text: str, pmid: str) -> dict:
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
        return {"no_segregation_data_found": True}


def _parse_relatives(raw_relatives: list) -> list[SegregationRelative]:
    relatives = []
    for r in raw_relatives or []:
        if not isinstance(r, dict) or not r.get("relationship"):
            continue
        relatives.append(SegregationRelative(
            relationship=str(r["relationship"]),
            affected_status=r.get("affected_status"),
            variant_status=r.get("variant_status"),
            zygosity=r.get("zygosity"),
        ))
    return relatives


async def search_family_segregation(
    gene: str, hgvsc: str, hgvsp: Optional[str] = None, *,
    preferred_pmids: tuple[str, ...] = (), max_candidates: int = 5,
) -> LiteratureSegregationResult:
    """Search PubMed for a family/pedigree study of this variant and ask
    the project's own LLM to extract per-relative segregation data.

    Returns the first candidate paper whose extraction reports at least
    one relative with both affected_status and variant_status stated;
    LiteratureSegregationResult.found=False (with `reason`) if nothing
    usable turned up.

    Query construction deliberately never embeds the raw HGVS c. notation
    (e.g. "c.1594T>C") into a PubMed query string - confirmed empirically
    in acmg_pipeline.pipeline.search_candidate_pmids's own docstring
    (2026-09-16) that exact c./p. notation returns ZERO results even for a
    variant whose own source paper is indexed in PubMed, which is exactly
    what the previous version of this query did. `hgvsp`, when available,
    lets search_candidate_pmids build its own protein-notation query
    ("{gene} AND {aa_change}", validated 5/5 on real cases) as the first
    candidate query instead.

    `preferred_pmids`, when given (this variant's own ERepo evidenceLinks -
    real VCEP-curated citations for its overall classification, not search
    results), are tried BEFORE any live search query. Confirmed empirically
    (2026-09-24) that a live search here, regardless of query wording,
    returns recent/generic gene papers with zero overlap against real
    ERepo citations for the same variant (checked directly: MYH7
    c.1594T>C, MSH2 c.1012G>A, RUNX1 c.601C>T all had 0/N overlap) - a
    family/pedigree paper is rare enough that when a VCEP already cited
    something for this variant, checking those first is far more likely
    to surface the real family study than a fresh keyword search, the same
    "curated citations first, live search only as fallback" pattern
    resolve_pmids_for_variant() already uses for PS3/BS3/PS4/BP5.
    """
    if not gene or not hgvsc:
        return LiteratureSegregationResult(found=False, reason="no_gene_or_hgvsc_available")

    from acmg_pipeline import pipeline

    async with AsyncExitStack() as stack:
        mcp = await pipeline.connect_pubmed(stack)
        pmids = list(dict.fromkeys(preferred_pmids))
        if len(pmids) < max_candidates:
            searched = await pipeline.search_candidate_pmids(
                mcp, gene, hgvsp=hgvsp, disease="family segregation pedigree", max_results=max_candidates,
            )
            pmids = list(dict.fromkeys([*pmids, *searched]))
        if len(pmids) < 3:
            more = await pipeline.search_candidate_pmids(
                mcp, gene, disease="family pedigree", max_results=max_candidates,
            )
            pmids = list(dict.fromkeys(pmids + more))[:max(max_candidates, len(preferred_pmids))]

        for pmid in pmids:
            text, _note = await pipeline.fetch_full_text(mcp, pmid)
            if not text:
                continue
            data = await _judge_paper_for_segregation(pipeline, gene, hgvsc, text, pmid)
            if data.get("no_segregation_data_found") or not data.get("family_segregation_described"):
                continue
            relatives = _parse_relatives(data.get("relatives"))
            informative = [r for r in relatives if r.affected_status is not None and r.variant_status is not None]
            if not informative:
                continue
            return LiteratureSegregationResult(
                found=True,
                inheritance_pattern=data.get("inheritance_pattern"),
                ar_case_mode=data.get("ar_case_mode"),
                relatives=relatives,
                pmid=pmid,
                quote=data.get("quote"),
                reason="confirmed_family_segregation_data_found",
            )

    return LiteratureSegregationResult(
        found=False,
        reason="no_confirmed_family_segregation_data_found_in_pubmed_search",
    )
