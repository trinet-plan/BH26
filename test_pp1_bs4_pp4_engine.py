"""
test_pp1_bs4_pp4_engine.py

Tests acmg_pipeline.criteria.pp1_bs4_pp4_engine - the glue connecting
pp4_pp1_bs4.py's judgment logic (pulled in from r-kobayashi's branch,
unmodified, already covered by test_pp4_pp1_bs4.py) to this project's
classify()/VA-Spec pipeline. PP4's diagnostic-yield input always comes
from a live acmg_pipeline.pp4_literature_search call now (2026-09-17, no
curated registry - see pp1_bs4_pp4_engine.py's own docstring for why).
This file fakes that call to stay offline/deterministic, the same way it
fakes hpo_mondo_extraction.normalize_hpo().
"""

import asyncio
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import acmg_pipeline.hpo_mondo_extraction as hpo_mondo_extraction
import acmg_pipeline.pp1_segregation_search as pp1_segregation_search
import acmg_pipeline.pp4_literature_search as pp4_literature_search
from acmg_pipeline.classification import CriterionStatus, Strength
from acmg_pipeline.clinical_note import (
    ClinicalFeature, ClinicalNoteExtraction, Family, Proband, ProbandPhenotype, Relative,
)
from acmg_pipeline.criteria import pp1_bs4_pp4_engine as engine
from acmg_pipeline.vcf_record import VariantRecord
from test_harness import Harness

h = Harness()
check = h.check


async def _identity_normalize_hpo(extraction):
    # engine.evaluate() calls hpo_mondo_extraction.normalize_hpo() (a real
    # TogoMCP + LLM round trip) before evaluate_locus_evidence(). Patching
    # it to a pass-through keeps this file fast, deterministic, and
    # offline - the phenotype_match itself no longer depends on HPO terms
    # at all (see [3] below), so there is nothing left to fake there.
    return extraction


hpo_mondo_extraction.normalize_hpo = _identity_normalize_hpo


def _fake_search(*, found, yield_fraction=None, sample_size=None,
                  denominator="all patients tested with this phenotype",
                  pmid="99999999", quote=None):
    """Builds a fake acmg_pipeline.pp4_literature_search.search_diagnostic_yield()."""
    async def _fake(gene, phenotype_description):
        if not found:
            return pp4_literature_search.LiteratureYieldResult(found=False, reason="test_not_found")
        return pp4_literature_search.LiteratureYieldResult(
            found=True, yield_fraction=yield_fraction, sample_size=sample_size,
            denominator_description=denominator, pmid=pmid, quote=quote,
        )
    return _fake


def _counting(fake):
    """Wraps a fake search coroutine to also record how many times it ran."""
    async def _wrapped(gene, phenotype_description):
        search_call_count["n"] += 1
        return await fake(gene, phenotype_description)
    return _wrapped


search_call_count = {"n": 0}
pp4_literature_search.search_diagnostic_yield = _fake_search(found=True, yield_fraction=0.70, sample_size=100)


def _fake_segregation_search(*, found, inheritance_pattern="autosomal dominant",
                              ar_case_mode=None, relatives=None, pmid="88888888"):
    """Builds a fake acmg_pipeline.pp1_segregation_search.search_family_segregation()
    (added 2026-09-18 alongside the live search itself) - same offline/
    deterministic-by-default convention as _fake_search() above, so every
    existing "no family data" fixture in this file keeps meaning exactly
    that, rather than silently making a real PubMed/LLM call."""
    async def _fake(gene, hgvsc):
        if not found:
            return pp1_segregation_search.LiteratureSegregationResult(found=False, reason="test_not_found")
        return pp1_segregation_search.LiteratureSegregationResult(
            found=True, inheritance_pattern=inheritance_pattern, ar_case_mode=ar_case_mode,
            relatives=relatives or [], pmid=pmid,
        )
    return _fake


# Default: no segregation study found, so a fixture with relatives=[] means
# exactly that (not "a live search hasn't been mocked yet").
pp1_segregation_search.search_family_segregation = _fake_segregation_search(found=False)


def run_evaluate(variant, note, config=None):
    # engine.evaluate() is async (hpo_mondo_extraction/pp4_literature_search are
    # both real TogoMCP/PubMed round trips in production) - this test file
    # stays plain top-to-bottom script style (test_harness.py convention),
    # so each call site just drives its own event loop.
    return asyncio.run(engine.evaluate(variant, note, config or {}))


def _variant(gene: str) -> VariantRecord:
    return VariantRecord(chrom="1", pos=1, id="", ref="A", alt="G", qual="", filter="",
                         info={"GENE": gene, "HGVSC": "c.1A>G"})


def _note(relatives=None, inheritance="autosomal dominant", diagnosis="some well-studied disease") -> ClinicalNoteExtraction:
    return ClinicalNoteExtraction(
        proband=Proband(phenotype=ProbandPhenotype(affected_status=True)),
        family=Family(inheritance_pattern=inheritance, relatives=relatives or []),
        diagnosis=diagnosis,
    )


# ============================================================================
# [1] No diagnosis at all -> UNKNOWN, search never attempted
# ============================================================================
print("[1] No diagnosis -> PP4 UNKNOWN, search never attempted")
pp4_literature_search.search_diagnostic_yield = _counting(
    _fake_search(found=True, yield_fraction=0.70, sample_size=100)
)
no_diagnosis_results = run_evaluate(_variant("NO_DIAGNOSIS_GENE"), _note(diagnosis=None))
check("PP4 is UNKNOWN with no diagnosis", no_diagnosis_results["PP4"].status == CriterionStatus.UNKNOWN)
check("PP1 is UNKNOWN too (no family data in this fixture)", no_diagnosis_results["PP1"].status == CriterionStatus.UNKNOWN)
check("BS4 is UNKNOWN too (no family data in this fixture)", no_diagnosis_results["BS4"].status == CriterionStatus.UNKNOWN)
check("search was never attempted", search_call_count["n"] == 0)

# The important case: PP1/BS4 must NOT be dragged down to UNKNOWN just
# because PP4's literature search has nothing to work with - family
# co-segregation (Table 3) is independent of PP4's diagnostic yield.
no_diagnosis_with_family_results = run_evaluate(
    _variant("NO_DIAGNOSIS_WITH_FAMILY_GENE"),
    _note(diagnosis=None, relatives=[Relative("sister", affected_status=True, variant_status=True)]),
)
check("PP4 UNKNOWN with no diagnosis (still true)",
      no_diagnosis_with_family_results["PP4"].status == CriterionStatus.UNKNOWN)
check("PP1 is MET from family data alone, despite PP4 having no diagnosis to search for",
      no_diagnosis_with_family_results["PP1"].status == CriterionStatus.MET
      and no_diagnosis_with_family_results["PP1"].strength == Strength.SUPPORTING)


# ============================================================================
# [2] Diagnosis present, search finds nothing usable -> PP4 UNKNOWN, PP1/BS4 unaffected
# ============================================================================
print("\n[2] Diagnosis present, nothing found -> PP4 UNKNOWN, PP1/BS4 still scored")
pp4_literature_search.search_diagnostic_yield = _fake_search(found=False)
not_found_results = run_evaluate(_variant("NOTHING_FOUND_GENE"), _note())
check("PP4 UNKNOWN when nothing usable is found", not_found_results["PP4"].status == CriterionStatus.UNKNOWN)
check("PP1 UNKNOWN too (no family data in this fixture)", not_found_results["PP1"].status == CriterionStatus.UNKNOWN)

not_found_with_family_results = run_evaluate(
    _variant("NOTHING_FOUND_WITH_FAMILY_GENE"),
    _note(relatives=[Relative("sister", affected_status=True, variant_status=True)]),
)
check("PP4 UNKNOWN when nothing usable is found (with family data present)",
      not_found_with_family_results["PP4"].status == CriterionStatus.UNKNOWN)
check("PP1 is MET from family data alone, despite PP4's search finding nothing",
      not_found_with_family_results["PP1"].status == CriterionStatus.MET
      and not_found_with_family_results["PP1"].strength == Strength.SUPPORTING)


# ============================================================================
# [3] Diagnosis present, search finds a confirmed overall yield -> real evidence
# ============================================================================
print("\n[3] Confirmed yield found -> real MET/NOT_MET evidence")
pp4_literature_search.search_diagnostic_yield = _fake_search(
    found=True, yield_fraction=0.70, sample_size=100, pmid="99999999",
    quote="70 of 100 patients tested had a pathogenic variant.",
)

# 70% -> 4.0 points (PP4) + one co-segregating sibling, autosomal dominant
# -> +1.0 (PP1), combined and capped at +5 - reproduces test_pp4_pp1_bs4.py's
# own case (2), through the full engine.evaluate() -> to_criterion_evidence()
# path. phenotype_match is assumed true (the search query WAS "this gene +
# this diagnosis"). yield=70% is below the 90% locus_model="homogeneous"
# approximation threshold (_build_reference_from_literature(), 2026-09-18),
# so locus_model="heterogeneous" here and PP1 is not suppressed.
note_with_relative = _note(relatives=[Relative("sister", affected_status=True, variant_status=True)])
sibling_results = run_evaluate(_variant("GENE1"), note_with_relative)
check("PP4 MET (4.0 points -> STRONG)",
      sibling_results["PP4"].status == CriterionStatus.MET
      and sibling_results["PP4"].strength == Strength.STRONG)
check("PP1 MET (1.0 point of family segregation -> SUPPORTING)",
      sibling_results["PP1"].status == CriterionStatus.MET
      and sibling_results["PP1"].strength == Strength.SUPPORTING)
check("BS4 NOT_MET (no non-segregation observed)",
      sibling_results["BS4"].status == CriterionStatus.NOT_MET)
check("disclosure is present in PP4's source",
      "auto-extracted via literature search" in sibling_results["PP4"].source)
check("PMID is present in PP4's source", "99999999" in sibling_results["PP4"].source)

# High yield (>90%) -> locus_model="homogeneous" approximation
# (_build_reference_from_literature(), 2026-09-18) suppresses PP1 even
# though family co-segregation data is present, matching the paper's own
# CTNS/cystinosis argument that PP1 must not be added on top of a
# locus-homogeneous, high-yield PP4 (double counting the same evidence).
pp4_literature_search.search_diagnostic_yield = _fake_search(
    found=True, yield_fraction=0.958, sample_size=100, pmid="22222222",
)
high_yield_results = run_evaluate(
    _variant("HIGH_YIELD_GENE"),
    _note(relatives=[Relative("sister", affected_status=True, variant_status=True)]),
)
check("PP4 MET at the +5.0 cap (95.8% yield -> STRONG)",
      high_yield_results["PP4"].status == CriterionStatus.MET
      and high_yield_results["PP4"].strength == Strength.STRONG)
check("PP1 NOT_MET - suppressed as redundant with a locus-homogeneous, high-yield PP4",
      high_yield_results["PP1"].status == CriterionStatus.NOT_MET)

# Reset back to the 70%-yield fake before the remaining sections.
pp4_literature_search.search_diagnostic_yield = _fake_search(
    found=True, yield_fraction=0.70, sample_size=100, pmid="99999999",
    quote="70 of 100 patients tested had a pathogenic variant.",
)

# No family data at all -> PP4 alone (segregation not evaluable).
note_no_family = _note(relatives=[])
solo_results = run_evaluate(_variant("GENE1"), note_no_family)
check("PP4 still MET without family data", solo_results["PP4"].status == CriterionStatus.MET)
check("PP1 UNKNOWN without family data (not evaluable, not a fabricated NOT_MET)",
      solo_results["PP1"].status == CriterionStatus.UNKNOWN)
check("BS4 UNKNOWN without family data", solo_results["BS4"].status == CriterionStatus.UNKNOWN)


# ============================================================================
# [4] Small sample size gets a caution note, without being blocked
# ============================================================================
print("\n[4] Small sample size -> caution note, not blocked")
pp4_literature_search.search_diagnostic_yield = _fake_search(
    found=True, yield_fraction=0.70, sample_size=3, pmid="11111111",
)
small_n_results = run_evaluate(_variant("SMALL_N_GENE"), _note(relatives=[]))
check("small sample size still produces a result (not blocked)",
      small_n_results["PP4"].status == CriterionStatus.MET)
check("small sample size gets a caution note in the source",
      "CAUTION" in small_n_results["PP4"].source and "n=3" in small_n_results["PP4"].source)


# ============================================================================
# [5] Every call re-searches - no persistent cache (deliberate, 2026-09-17)
# ============================================================================
print("\n[5] Every call re-searches; nothing is persisted to disk")
# An earlier design persisted found entries back into a config/ registry so
# a later call for the same gene would skip re-searching. The user's
# explicit direction (2026-09-17, after weighing PS3/BS3/PS4's own
# per-call, non-persistent pattern against the maintenance cost of any
# registry) was to drop persistence entirely: every call searches live,
# same as PS3/BS3/PS4 already do.
search_call_count["n"] = 0
pp4_literature_search.search_diagnostic_yield = _counting(
    _fake_search(found=True, yield_fraction=0.70, sample_size=100)
)
run_evaluate(_variant("GENE1"), note_no_family)
run_evaluate(_variant("GENE1"), note_no_family)
check("the same gene is searched again on a second call (no caching)",
      search_call_count["n"] == 2)


# ============================================================================
# [6] build_evidence_line() produces a valid, real VA-Spec EvidenceLine
# ============================================================================
print("\n[6] EvidenceLine construction")
variant = _variant("GENE1")
pp4_line = engine.build_evidence_line("PP4", sibling_results["PP4"], variant)
check("PP4 line has the right criterion id", pp4_line["specifiedBy"]["methodType"] == "PP4")
check("PP4 line direction is 'supports' (MET, pathogenic code)",
      pp4_line["directionOfEvidenceProvided"] == "supports")
check("PP4 line carries a strengthOfEvidenceProvided",
      pp4_line.get("strengthOfEvidenceProvided", {}).get("primaryCoding", {}).get("code") == "strong")
check("PP4 line's status extension round-trips status=met",
      any(e["name"] == "status" and e["value"] == "met"
          for e in pp4_line["extensions"]))

bs4_line = engine.build_evidence_line("BS4", sibling_results["BS4"], variant)
check("BS4 (NOT_MET) line direction is 'neutral'", bs4_line["directionOfEvidenceProvided"] == "neutral")
check("BS4 (NOT_MET) line has no strengthOfEvidenceProvided", "strengthOfEvidenceProvided" not in bs4_line)

unknown_line = engine.build_evidence_line("PP1", no_diagnosis_results["PP1"], _variant("NO_DIAGNOSIS_GENE"))
check("UNKNOWN line direction is 'neutral'", unknown_line["directionOfEvidenceProvided"] == "neutral")
check("UNKNOWN line has no evidenceOutcome", "evidenceOutcome" not in unknown_line)


# ============================================================================
# [7] Full pipeline wiring: evaluate_variant_evidence_lines()/evaluate_
#     selected_criteria() actually call this engine for PP1/BS4/PP4
# ============================================================================
print("\n[7] Wired into acmg_pipeline.constants.IMPLEMENTED_CODES")
from acmg_pipeline.constants import IMPLEMENTED_CODES, PHENOTYPE_SEGREGATION_CODES, STUB_CODES

check("PP1/BS4/PP4 are in PHENOTYPE_SEGREGATION_CODES",
      PHENOTYPE_SEGREGATION_CODES == {"PP1", "BS4", "PP4"})
check("PP1/BS4/PP4 are now IMPLEMENTED, not STUB",
      PHENOTYPE_SEGREGATION_CODES <= IMPLEMENTED_CODES
      and PHENOTYPE_SEGREGATION_CODES.isdisjoint(STUB_CODES))


# ============================================================================
# [8] PP1/BS4's own live literature search (added 2026-09-18): when the
#     caller supplies no family.relatives at all, engine.evaluate() tries
#     acmg_pipeline.pp1_segregation_search before giving up as UNKNOWN -
#     the same "search live rather than stay UNKNOWN forever" pattern PP4
#     already uses for diagnostic yield.
# ============================================================================
print("\n[8] PP1/BS4 fall back to a live family-segregation search")
from acmg_pipeline.pp1_segregation_search import SegregationRelative

pp4_literature_search.search_diagnostic_yield = _fake_search(found=False)

pp1_segregation_search.search_family_segregation = _fake_segregation_search(found=False)
no_segregation_results = run_evaluate(_variant("NO_SEGREGATION_PAPER_GENE"), _note(relatives=[]))
check("PP1 stays UNKNOWN when the segregation search also finds nothing",
      no_segregation_results["PP1"].status == CriterionStatus.UNKNOWN)

pp1_segregation_search.search_family_segregation = _fake_segregation_search(
    found=True, inheritance_pattern="autosomal dominant",
    relatives=[SegregationRelative("sister", affected_status=True, variant_status=True)],
)
found_segregation_results = run_evaluate(_variant("SEGREGATION_PAPER_GENE"), _note(relatives=[]))
check("PP1 is MET once a literature search finds a segregating relative",
      found_segregation_results["PP1"].status == CriterionStatus.MET
      and found_segregation_results["PP1"].strength == Strength.SUPPORTING)

# Caller-supplied family data always wins - the search only fills a gap,
# never overrides real input the caller already has.
pp1_segregation_search.search_family_segregation = _counting(_fake_segregation_search(found=True))
search_call_count["n"] = 0
run_evaluate(_variant("GENE1"), note_with_relative)
check("the search is skipped when the caller already supplied family.relatives",
      search_call_count["n"] == 0)


h.report_and_exit()
