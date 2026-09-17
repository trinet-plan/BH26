"""
phenotype_matchers.py

Pluggable strategies for PP4's phenotype_match gate ("does the patient's
phenotype match the cohort the curated diagnostic_yield was derived
from?") - the input evaluate_locus_evidence()'s phenotype_match_override
seam (acmg_pipeline/criteria/pp4_pp1_bs4.py) accepts.

Why this exists as its own module
----------------------------------
Until now, acmg_pipeline.criteria.pp1_bs4_pp4_engine.evaluate() called
acmg_pipeline.pubcasefinder directly and unconditionally. Per the user's
explicit direction (2026-09-17): PubCaseFinder's differential-diagnosis
RANKING is a reasonable automatic proxy, but it is not what Biesecker et
al., 2024 (the ClinGen PP1/BS4/PP4 guidance this project's algorithm
otherwise follows exactly - see pp4_pp1_bs4.py/evaluator.py) actually
specifies. That paper's own Table 2 footnote is explicit that
diagnostic_yield "is wholly dependent upon the phenotype criteria
matching" - i.e. the match is against the SAME curated phenotype
definition (e.g. "typical cystinosis," Ghent criteria for Marfan
syndrome) the cited diagnostic-yield study itself used to enroll its
cohort, not a general "which gene best explains these symptoms" query.
That is exactly what PP4ReferenceRecord.phenotype_hpo +
match_phenotype_constellation() (pp4_pp1_bs4.py, left unmodified) already
implements - curated_hpo_list_matcher below is a thin wrapper around it.

PubCaseFinder is kept here as a second, selectable option rather than
removed - a real gene may end up curated with no phenotype_hpo definition
recorded, or a future case may want an automatic cross-check - but it is
no longer the default, since it is not the paper's own method.

Adding a new matcher
---------------------
Any callable matching the PhenotypeMatcher signature below, added to
PHENOTYPE_MATCHERS under a new name, becomes selectable the same way -
see pp1_bs4_pp4_engine.evaluate()'s docstring for how a config key /
per-gene `gates.phenotype_matcher` entry picks one.
"""

from __future__ import annotations

from typing import Awaitable, Callable

from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.criteria.pp4_pp1_bs4 import (
    PP4ReferenceRecord,
    PhenotypeMatchResult,
    match_phenotype_constellation,
    patient_hpo_terms,
)

PhenotypeMatcher = Callable[
    [ClinicalNoteExtraction, PP4ReferenceRecord, str], Awaitable[PhenotypeMatchResult]
]


async def curated_hpo_list_matcher(
    clinical_note: ClinicalNoteExtraction,
    reference: PP4ReferenceRecord,
    gene: str,
) -> PhenotypeMatchResult:
    """Biesecker et al., 2024's own method: exact match against the curated
    phenotype definition (PP4ReferenceRecord.phenotype_hpo) the cited
    diagnostic-yield study itself used. `gene` is unused - kept only so
    every PhenotypeMatcher shares one call signature.

    A thin, synchronous-underneath wrapper: pp4_pp1_bs4.match_phenotype_
    constellation() is pure and already returns a PhenotypeMatchResult, so
    there is nothing to await here besides the coroutine wrapper itself.
    """
    return match_phenotype_constellation(clinical_note, reference)


async def pubcasefinder_matcher(
    clinical_note: ClinicalNoteExtraction,
    reference: PP4ReferenceRecord,
    gene: str,
) -> PhenotypeMatchResult:
    """Automatic proxy via PubCaseFinder's HPO-based gene ranking (TogoMCP).

    `reference` is unused (no curated phenotype_hpo definition is needed) -
    kept only for the shared PhenotypeMatcher signature. See
    acmg_pipeline.pubcasefinder's own module docstring for what this does
    and does not decide, and why matched=None (not a guessed non-match) is
    returned when PubCaseFinder itself is unavailable or the patient has no
    normalized HPO terms.
    """
    from acmg_pipeline import pubcasefinder

    patient_terms = patient_hpo_terms(clinical_note)
    hpo_ids = [term.hpo_id for term in patient_terms]
    if not hpo_ids:
        return PhenotypeMatchResult(
            matched=None, patient_terms=patient_terms, reason="no_normalized_proband_hpo",
        )

    try:
        ranking = await pubcasefinder.rank_genes_by_phenotype(hpo_ids)
    except (ValueError, RuntimeError) as exc:
        return PhenotypeMatchResult(
            matched=None, patient_terms=patient_terms,
            reason=f"pubcasefinder_unavailable:{exc}",
        )
    return pubcasefinder.phenotype_match_from_gene_ranking(ranking, gene, patient_terms=patient_terms)


PHENOTYPE_MATCHERS: dict[str, PhenotypeMatcher] = {
    "curated_hpo_list": curated_hpo_list_matcher,
    "pubcasefinder": pubcasefinder_matcher,
}

# Biesecker et al., 2024's own method - see this module's docstring for why
# PubCaseFinder is available but not the default.
DEFAULT_PHENOTYPE_MATCHER = "curated_hpo_list"
