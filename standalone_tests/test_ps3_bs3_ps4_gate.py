import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from acmg_pipeline.gate import (
    match_variant_in_text, MatchStatus,
    OfflineERepoClient, CriterionStatus,
    check_approved_assay, AssayApplicability,
    run_gate, GateMode,
)
from test_harness import Harness

h = Harness()
check = h.check


# --- 1. Variant matching (the MT-ND6 off-by-one-coordinate case, section 6-2) ---
print("[1] Variant matching (coordinate tolerance)")
text = "single fiber PCR analysis shows m.14512_14513del segregation..."
result = match_variant_in_text("m.14513_14514del", text)
check("MT-ND6 off-by-one is caught as heuristic", result.status == MatchStatus.HEURISTIC)
check("the drift is recorded as 1", result.tolerance_used == 1)

# a genuine mismatch (a different CDH1 variant, +1 vs +5)
result2 = match_variant_in_text("c.387+5G>A", "the CDH1 c.387+1G>A variant was...")
check("a genuinely different variant (+1 vs +5) does not match", result2.status == MatchStatus.UNSUCCESSFUL)

# --- 2. ERepo lookup (real RIT1 data, idx 8336) ---
print("\n[2] ERepo lookup (real RIT1 data)")
rit1_response = {
    "variantInterpretations": [{
        "gene": {"label": "RIT1"},
        "caid": "CAR:CA342802812",
        "hgvs": ["NM_006912.6(RIT1):c.268A>G (p.Met90Val)"],
        "guidelines": [{
            "outcome": {"label": "Pathogenic"},
            "agents": [{"evidenceCodes": [
                {"label": "PS3_Supporting", "status": "Met"},
                {"label": "PS2_Very Strong", "status": "Met"},
            ]}],
        }],
    }]
}
client = OfflineERepoClient([rit1_response])
lookup = client.lookup("RIT1", "c.268A>G")
check("the RIT1 variant is found in ERepo", lookup.found_in_erepo)
check("PS3_Supporting reads as MET", lookup.evidence_code_status.get("PS3_Supporting") == CriterionStatus.MET)

# --- 3. The gate itself: in production mode, a known variant does not start the pipeline ---
print("\n[3] The gate itself (production mode)")
decision = run_gate("RIT1", "c.268A>G", "PS3_Supporting", client, mode=GateMode.PRODUCTION)
check("pipeline is not started in production", decision.should_run_pipeline is False)
check("resolved_status is MET", decision.resolved_status == CriterionStatus.MET)

# --- 4. The gate itself: in validation mode, the pipeline is forced to run with the truth hidden ---
print("\n[4] The gate itself (validation mode)")
decision_v = run_gate("RIT1", "c.268A>G", "PS3_Supporting", client, mode=GateMode.VALIDATION)
check("pipeline is forced to run in validation", decision_v.should_run_pipeline is True)
check("resolved_status is not shown downstream (None)", decision_v.resolved_status is None)
check("the ground truth is held only in _hidden_ground_truth", decision_v._hidden_ground_truth == CriterionStatus.MET)

# --- 5. An uncurated variant starts the pipeline in either mode ---
print("\n[5] An uncurated variant")
decision_unknown = run_gate("BRCA2", "c.9999X>Y", "PS3", client, mode=GateMode.PRODUCTION)
check("pipeline starts even in production when uncurated", decision_unknown.should_run_pipeline is True)

# --- 6. Approved-assay cross-check ---
print("\n[6] Approved-assay cross-check")
r1 = check_approved_assay("RASopathy VCEP", "RIT1", "ERK1/2 phosphorylation assay in HEK293T cells", "PS3")
check("RIT1's ERK assay is APPROVED", r1.applicability == AssayApplicability.APPROVED)

r2 = check_approved_assay("RASopathy VCEP", "PIK3CA", "PI3K activity and S6 phosphorylation assay", "PS3")
check("PIK3CA is not on the RASopathy approved list -> ASSAY_NOT_APPROVED", r2.applicability == AssayApplicability.ASSAY_NOT_APPROVED)

r3 = check_approved_assay("RASopathy VCEP", "RIT1", "any assay", "BS3")
check("BS3 is BS3_NOT_APPLICABLE for RASopathy", r3.applicability == AssayApplicability.BS3_NOT_APPLICABLE)

r4 = check_approved_assay("Unknown VCEP", "GENE1", "some assay", "PS3")
check("an unregistered VCEP yields NO_SPEC_AVAILABLE", r4.applicability == AssayApplicability.NO_SPEC_AVAILABLE)

r5 = check_approved_assay("RASopathy VCEP", "LZTR1", "ERK1/2 phosphorylation assay", "PS3")
check("ERK assay + LZTR1 (within the applicable-gene list) is APPROVED", r5.applicability == AssayApplicability.APPROVED)

# Regression test: discovered when translating descriptions to English -
# "ERK Activation Assay (... after serum stimulation)" was incorrectly
# matched against "MEK Activation Assay" because both descriptions contain
# the generic word "after" (MEK is checked first in the list).
r6 = check_approved_assay(
    "RASopathy VCEP", "RIT1",
    "ERK Activation Assay (phospho-ERK1/2 western blot after serum stimulation)",
    "PS3",
)
check("the RIT1 ERK assay is correctly matched to ERK, not MEK (stopword-filter regression test)",
      r6.applicability == AssayApplicability.APPROVED and r6.matched_assay.name == "ERK Activation Assay")

# --- 7. Multi-paper evidence extraction from ERepo's evidenceLinks (real RUNX1 data, 4 PMIDs) ---
print("\n[7] Multi-paper evidence extraction (real RUNX1 data, idx 3859 area)")
runx1_multi_pmid_response = {
    "variantInterpretations": [{
        "gene": {"label": "RUNX1"},
        "caid": "CAR:CA16602487",
        "hgvs": ["NM_001754.4(RUNX1):c.601C>T (p.Arg201Ter)"],
        "guidelines": [{
            "outcome": {"label": "Pathogenic"},
            "agents": [{"evidenceCodes": [
                {"label": "PVS1", "status": "Met"},
                {"label": "PS4", "status": "Met"},
            ]}],
        }],
        "evidenceLinks": [
            {"@id": "https://www.ncbi.nlm.nih.gov/pubmed/20549580"},
            {"@id": "https://www.ncbi.nlm.nih.gov/pubmed/10508512"},
            {"@id": "https://www.ncbi.nlm.nih.gov/pubmed/28513614"},
            {"@id": "https://www.ncbi.nlm.nih.gov/pubmed/19387465"},
        ],
    }]
}
multi_client = OfflineERepoClient([runx1_multi_pmid_response])
lookup_multi = multi_client.lookup("RUNX1", "c.601C>T")
check("all 4 real PMIDs are extracted from evidenceLinks",
      lookup_multi.evidence_pmids == ["20549580", "10508512", "28513614", "19387465"])

h.report_and_exit()
