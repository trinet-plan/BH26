"""
test_pp1_bs4_pp4_engine.py

Tests acmg_pipeline.criteria.pp1_bs4_pp4_engine - the glue connecting
pp4_pp1_bs4.py's judgment logic (pulled in from r-kobayashi's branch,
unmodified, already covered by test_pp4_pp1_bs4.py) to this project's
classify()/VA-Spec pipeline. This file tests the ENGINE (curated-reference
lookup, points-to-Strength mapping, EvidenceLine construction), not the
judgment logic itself - see test_pp4_pp1_bs4.py for that.
"""

import asyncio
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import acmg_pipeline.hpo_extraction as hpo_extraction
import acmg_pipeline.pubcasefinder as pubcasefinder
from acmg_pipeline.classification import CriterionStatus, Strength
from acmg_pipeline.clinical_note import (
    ClinicalFeature, ClinicalNoteExtraction, Family, Proband, ProbandPhenotype, Relative,
)
from acmg_pipeline.criteria import pp1_bs4_pp4_engine as engine
from acmg_pipeline.criteria.pp4_pp1_bs4 import evaluate_locus_evidence
from acmg_pipeline.vcf_record import VariantRecord
from test_harness import Harness

h = Harness()
check = h.check


async def _identity_normalize_hpo(extraction):
    # engine.evaluate() now calls hpo_extraction.normalize_hpo() (a real
    # TogoMCP + LLM round trip) before evaluate_locus_evidence(). This
    # test's fixtures already hand-set ClinicalFeature.hpo_id to synthetic
    # values (HP:0000001/HP:0000002) matching the synthetic reference
    # records below, purely to test the ENGINE's own logic - not TogoMCP
    # resolution of real phenotype labels, which is hpo_extraction.py's own
    # concern (exercised directly, live, in the earlier manual check of
    # resolve_hpo_labels()). Patching normalize_hpo() to a pass-through
    # keeps this file fast, deterministic, and offline, as it was before
    # engine.evaluate() became async.
    return extraction


hpo_extraction.normalize_hpo = _identity_normalize_hpo


def _fake_rank_genes_by_phenotype(*, top_gene: str = "GENE1"):
    # engine.evaluate() now also calls pubcasefinder.rank_genes_by_phenotype()
    # (a real TogoMCP round trip to PubCaseFinder) for PP4's phenotype_match
    # gate instead of matching reference.phenotype_hpo directly. Faking the
    # ranking here keeps this file offline/deterministic, the same way
    # _identity_normalize_hpo() stands in for the real TogoMCP HPO lookup.
    # `top_gene` is whichever gene this fixture puts at rank 1; the other
    # gene of the pair is always ranked 2, to test both PP4's matched=True
    # (rank 1) and matched=False (ranked, but not rank 1) paths.
    async def _fake(hpo_ids):
        other = "OTHER_GENE" if top_gene != "OTHER_GENE" else "GENE1"
        return {"results": [
            {"rank": 1, "score": 1.0, "gene_symbol": top_gene, "matched_hpo_ids": hpo_ids},
            {"rank": 2, "score": 0.5, "gene_symbol": other, "matched_hpo_ids": []},
        ]}

    return _fake


pubcasefinder.rank_genes_by_phenotype = _fake_rank_genes_by_phenotype()


def run_evaluate(variant, note, config):
    # engine.evaluate() became async on 2026-09-17 once it started calling
    # acmg_pipeline.hpo_extraction.normalize_hpo() (a TogoMCP + LLM round
    # trip) for genes with a curated reference record - this test file
    # stays plain top-to-bottom script style (test_harness.py convention),
    # so each call site just drives its own event loop instead of the
    # whole file becoming async.
    return asyncio.run(engine.evaluate(variant, note, config))


def _variant(gene: str) -> VariantRecord:
    return VariantRecord(chrom="1", pos=1, id="", ref="A", alt="G", qual="", filter="",
                         info={"GENE": gene, "HGVSC": "c.1A>G"})


def _note(relatives=None, inheritance="autosomal dominant") -> ClinicalNoteExtraction:
    # HPO IDs must match the synthetic reference's phenotype_hpo (see
    # section [2] below) for match_phenotype_constellation() to return
    # matched=True - without any patient HPO terms at all it returns
    # matched=None ("not evaluable"), same as test_pp4_pp1_bs4.py's own
    # case (5) fixture.
    return ClinicalNoteExtraction(
        proband=Proband(phenotype=ProbandPhenotype(
            affected_status=True,
            clinical_features=[
                ClinicalFeature("feature A", "HP:0000001"),
                ClinicalFeature("feature B", "HP:0000002"),
            ],
        )),
        family=Family(inheritance_pattern=inheritance, relatives=relatives or []),
    )


# ============================================================================
# [1] No curated reference record for this gene -> all three UNKNOWN
# ============================================================================
print("[1] No curated reference -> UNKNOWN")
empty_registry = Path(tempfile.mkdtemp()) / "empty.json"
empty_registry.write_text(json.dumps({"schema_version": "1.0", "entries": []}), encoding="utf-8")

results = run_evaluate(
    _variant("UNCURATED_GENE"), _note(),
    {"pp4_reference_records_path": str(empty_registry)},
)
check("all 3 codes returned", set(results) == {"PP1", "BS4", "PP4"})
check("PP4 is UNKNOWN", results["PP4"].status == CriterionStatus.UNKNOWN)
check("PP1 is UNKNOWN", results["PP1"].status == CriterionStatus.UNKNOWN)
check("BS4 is UNKNOWN", results["BS4"].status == CriterionStatus.UNKNOWN)

# A DRAFT (not APPROVED) entry must not load either.
draft_registry = Path(tempfile.mkdtemp()) / "draft.json"
draft_registry.write_text(json.dumps({
    "schema_version": "1.0",
    "entries": [{
        "status": "DRAFT", "id": "draft-1", "gene": "GENE1",
        "phenotype_label": "x", "phenotype_hpo": [], "locus_model": "heterogeneous",
        "diagnostic_yield": 0.7, "testing_method": "sequencing", "source_citation": "x",
    }],
}), encoding="utf-8")
draft_results = run_evaluate(_variant("GENE1"), _note(), {"pp4_reference_records_path": str(draft_registry)})
check("a DRAFT (unapproved) entry does not load", draft_results["PP4"].status == CriterionStatus.UNKNOWN)


# ============================================================================
# [2] Curated reference present -> real evidence, via a synthetic APPROVED entry
# ============================================================================
print("\n[2] Curated (synthetic) reference -> real MET/NOT_MET evidence")
approved_registry = Path(tempfile.mkdtemp()) / "approved.json"
approved_registry.write_text(json.dumps({
    "schema_version": "1.0",
    "entries": [{
        "status": "APPROVED", "id": "test-gene1-v1", "gene": "GENE1",
        "phenotype_label": "Demo phenotype", "phenotype_hpo": ["HP:0000001", "HP:0000002"],
        "locus_model": "heterogeneous", "diagnostic_yield": 0.70,
        "testing_method": "sequencing", "source_citation": "demo source",
        "gates": {"method_comparable": True},
    }],
}), encoding="utf-8")

config = {"pp4_reference_records_path": str(approved_registry)}

# Reproduces test_pp4_pp1_bs4.py case (2): heterogeneous disease, PP4 + PP1, capped at +5 -
# but here going through the FULL engine.evaluate() -> to_criterion_evidence() path, not
# evaluate_locus_evidence() directly.
note_with_relative = _note(relatives=[Relative("sister", affected_status=True, variant_status=True)])
sibling_results = run_evaluate(_variant("GENE1"), note_with_relative, config)
check("PP4 MET (matches case (1)'s 4.0 points -> Strength tier below 8, at/above 4 -> STRONG)",
      sibling_results["PP4"].status == CriterionStatus.MET
      and sibling_results["PP4"].strength == Strength.STRONG)
check("PP1 MET (1.0 point of family segregation, at/above 1 -> SUPPORTING)",
      sibling_results["PP1"].status == CriterionStatus.MET
      and sibling_results["PP1"].strength == Strength.SUPPORTING)
check("BS4 NOT_MET (no non-segregation observed)",
      sibling_results["BS4"].status == CriterionStatus.NOT_MET)

# No family data at all -> PP4 alone (segregation not evaluable).
note_no_family = _note(relatives=[])
solo_results = run_evaluate(_variant("GENE1"), note_no_family, config)
check("PP4 still MET without family data", solo_results["PP4"].status == CriterionStatus.MET)
check("PP1 UNKNOWN without family data (not evaluable, not a fabricated NOT_MET)",
      solo_results["PP1"].status == CriterionStatus.UNKNOWN)
check("BS4 UNKNOWN without family data", solo_results["BS4"].status == CriterionStatus.UNKNOWN)


# ============================================================================
# [2b] PubCaseFinder supplies phenotype_match, not reference.phenotype_hpo
# ============================================================================
print("\n[2b] PubCaseFinder phenotype-specificity gate")

# GENE1 ranks below OTHER_GENE for this phenotype profile -> not a phenotype
# match, even though reference.phenotype_hpo (still present in the fixture
# above) would have matched under the old required-HPO-list check.
pubcasefinder.rank_genes_by_phenotype = _fake_rank_genes_by_phenotype(top_gene="OTHER_GENE")
not_top_ranked_results = run_evaluate(_variant("GENE1"), note_no_family, config)
check("PP4 NOT_MET when GENE1 is not PubCaseFinder's top-ranked gene",
      not_top_ranked_results["PP4"].status == CriterionStatus.NOT_MET)
check("reason cites PubCaseFinder, not the curated HPO list",
      "pubcasefinder" in not_top_ranked_results["PP4"].source)

# PubCaseFinder itself unavailable (rate limit / no recognized HPO IDs / HTTP
# error) -> "not evaluable", not a fabricated non-match.
async def _raising_rank(hpo_ids):
    raise ValueError("PUBCASEFINDER_RATE_LIMIT_EXCEEDED")


pubcasefinder.rank_genes_by_phenotype = _raising_rank
unavailable_results = run_evaluate(_variant("GENE1"), note_no_family, config)
check("PP4 UNKNOWN when PubCaseFinder raises",
      unavailable_results["PP4"].status == CriterionStatus.UNKNOWN)
check("reason cites pubcasefinder_unavailable",
      "pubcasefinder_unavailable" in unavailable_results["PP4"].source)

pubcasefinder.rank_genes_by_phenotype = _fake_rank_genes_by_phenotype()


# ============================================================================
# [3] build_evidence_line() produces a valid, real VA-Spec EvidenceLine
# ============================================================================
print("\n[3] EvidenceLine construction")
variant = _variant("GENE1")
pp4_line = engine.build_evidence_line("PP4", sibling_results["PP4"], variant)
check("PP4 line has the right criterion id", pp4_line["specifiedBy"]["methodType"] == "PP4")
check("PP4 line direction is 'supports' (MET, pathogenic code)",
      pp4_line["directionOfEvidenceProvided"] == "supports")
check("PP4 line carries a strengthOfEvidenceProvided",
      pp4_line.get("strengthOfEvidenceProvided", {}).get("primaryCoding", {}).get("code") == "strong")
check("PP4 line's bh26AssessmentDetails round-trips status=met",
      any(e["name"] == "bh26AssessmentDetails" and e["value"]["status"] == "met"
          for e in pp4_line["extensions"]))

bs4_line = engine.build_evidence_line("BS4", sibling_results["BS4"], variant)
check("BS4 (NOT_MET) line direction is 'neutral'", bs4_line["directionOfEvidenceProvided"] == "neutral")
check("BS4 (NOT_MET) line has no strengthOfEvidenceProvided", "strengthOfEvidenceProvided" not in bs4_line)

unknown_line = engine.build_evidence_line("PP1", results["PP1"], _variant("UNCURATED_GENE"))
check("UNKNOWN line direction is 'neutral'", unknown_line["directionOfEvidenceProvided"] == "neutral")
check("UNKNOWN line has no evidenceOutcome", "evidenceOutcome" not in unknown_line)


# ============================================================================
# [4] Full pipeline wiring: evaluate_variant_evidence_lines()/evaluate_
#     selected_criteria() actually call this engine for PP1/BS4/PP4
# ============================================================================
print("\n[4] Wired into acmg_pipeline.constants.IMPLEMENTED_CODES")
from acmg_pipeline.constants import IMPLEMENTED_CODES, PHENOTYPE_SEGREGATION_CODES, STUB_CODES

check("PP1/BS4/PP4 are in PHENOTYPE_SEGREGATION_CODES",
      PHENOTYPE_SEGREGATION_CODES == {"PP1", "BS4", "PP4"})
check("PP1/BS4/PP4 are now IMPLEMENTED, not STUB",
      PHENOTYPE_SEGREGATION_CODES <= IMPLEMENTED_CODES
      and PHENOTYPE_SEGREGATION_CODES.isdisjoint(STUB_CODES))


h.report_and_exit()
