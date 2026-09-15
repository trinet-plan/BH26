"""
test_ps3_bs3_ps4_gate_full.py
A regression test covering all the real data gathered so far.

Coverage:
  A. Variant-matching gate: matches all 102 PS3/BS3-family rows in the
     corrected CGBench CSV against the actual paper full text (the ft
     column), and confirms this reproduces exactly the results from design
     doc section 6 (matched 57 / heuristic 1 / unsuccessful 44).
  B. ERepo lookup + the gate itself: for the same 102 rows, runs the gate in
     production mode via CsvBackedERepoClient and confirms every result is
     consistent with the CSV's own met_status. Also individually verifies
     that the 26 corrections (not_met -> met) actually confirmed via the
     live ClinGen API in section 9 are correctly reflected.
  C. Approved-assay cross-check: checks against real-world examples actually
     encountered (RIT1's ERK assay, the PIK3CA/PIK3R2 PI3K-pathway assay).
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.gate import (
    match_variant_in_text, MatchStatus,
    CsvBackedERepoClient, OfflineERepoClient, CriterionStatus,
    check_approved_assay, AssayApplicability,
    run_gate, GateMode,
)

# References the CSV placed alongside this test file (no absolute-path dependency)
CORRECTED_CSV = str(Path(__file__).resolve().parent / "clingen_vci_pubmed_fulltext_dedup_pmid_CORRECTED.csv")

passed = 0
failed = 0
failures: list[str] = []


def check(label, cond):
    global passed, failed
    if cond:
        passed += 1
    else:
        failed += 1
        failures.append(label)
        print(f"  FAIL {label}")


def base_code(code: str) -> str:
    return code.split("-")[0]


# ============================================================================
# A. Full regression test of the variant-matching gate (102 rows)
# ============================================================================
print("[A] Variant-matching gate: full regression test over all 102 PS3/BS3 rows")

with open(CORRECTED_CSV, newline="", encoding="utf-8") as f:
    all_rows = list(csv.DictReader(f))

ps3_bs3_rows = [r for r in all_rows if base_code(r["evidence_code"]) in ("PS3", "BS3")]
check("the target row count is 102", len(ps3_bs3_rows) == 102)

status_counts = {"matched": 0, "heuristic": 0, "unsuccessful": 0}
mismatches = []
for row in ps3_bs3_rows:
    result = match_variant_in_text(row["variant"], row["ft"])
    status_counts[result.status.value] += 1

print(f"  Results: matched={status_counts['matched']}, heuristic={status_counts['heuristic']}, "
      f"unsuccessful={status_counts['unsuccessful']}")

check("matched count matches section 6's result (57)", status_counts["matched"] == 57)
check("heuristic count matches section 6's result (1)", status_counts["heuristic"] == 1)
check("unsuccessful count matches section 6's result (44)", status_counts["unsuccessful"] == 44)


# ============================================================================
# B. ERepo lookup + the gate itself: all 102 rows + individual checks on the
#    26 corrections confirmed in section 9
# ============================================================================
print("\n[B] ERepo lookup + the gate itself: all rows + individual checks on the 26 corrected cases")

erepo = CsvBackedERepoClient(CORRECTED_CSV)

# B-1: for all 102 rows, confirm the gate's resolved_status is never
# inconsistent with the CSV's own met_status
n_checked = 0
for row in ps3_bs3_rows:
    decision = run_gate(row["hgnc_gene"], row["variant"], row["evidence_code"], erepo, mode=GateMode.PRODUCTION)
    expected = CriterionStatus.MET if row["met_status"] == "met" else CriterionStatus.NOT_MET
    if decision.resolved_status is not None:
        n_checked += 1
        check(
            f"{row['hgnc_gene']} {row['variant']} ({row['evidence_code']}) judgment matches",
            decision.resolved_status == expected,
        )
check("the gate resolved more than 0 of the 102 rows", n_checked > 0)
print(f"  The gate resolved {n_checked} of the 102 rows (pipeline not started, in production mode)")

# B-2: individually confirm that the 26 "corrections" confirmed via the live
# ClinGen API in section 9 are correctly reflected
# (entry_index, evidence_code, gene, variant) -> confirmed ground truth (all MET after correction)
VERIFIED_CORRECTIONS = [
    ("3786", "MT-ND1", "m.3890G>A", "PS3-Supporting"),
    ("3964", "DICER1", "c.2642T>C", "PS3-Supporting"),
    ("4046", "LDLR", "c.-136C>G", "PS3-Moderate"),
    ("5072", "LDLR", "c.1003G>T", "PS3-Moderate"),
    ("8146", "LDLR", "c.1444G>C", "PS3-Moderate"),
    ("8149", "LDLR", "c.1739C>T", "PS3-Supporting"),
    ("4263", "MT-TL1", "m.3291T>C", "PS3-Supporting"),
    ("4870", "MT-TN", "m.5690A>G", "PS3-Supporting"),
    ("5109", "MT-TS1", "m.7497G>A", "PS3-Supporting"),
    ("5545", "HNF1A", "c.66C>G", "PS3-Supporting"),
    ("6600", "MT-ND6", "m.14513_14514del", "PS3-Supporting"),
    ("4852", "MT-ND6", "m.14487T>C", "PS3-Moderate"),
    ("6612", "GCK", "c.447C>A", "PS3-Moderate"),
    ("7035", "RAG1", "c.2690G>A", "PS3-Moderate"),
    ("7223", "VHL", "c.273C>A", "PS3-Supporting"),
    ("7456", "MT-TK", "m.8340G>A", "PS3-Supporting"),
    ("7655", "MT-ND6", "m.14597A>G", "PS3-Supporting"),
    ("8336", "RIT1", "c.268A>G", "PS3-Supporting"),
    ("9142", "DYSF", "c.953T>A", "PS3-Moderate"),
    ("9144", "DYSF", "c.1906G>A", "PS3-Moderate"),
    ("9155", "DYSF", "c.1717C>T", "PS3-Moderate"),
    ("9200", "SGCB", "c.452C>G", "PS3-Moderate"),
    ("9523", "GUCY2D", "c.1762C>T", "PS3-Supporting"),
    ("9551", "GUCY2D", "c.1724C>T", "BS3-Supporting"),
    ("1587", "TP53", "c.329G>A", "BS3-Supporting"),
    ("3208", "MYOC", "c.898G>A", "PS3"),
]
check("the corrections list has 26 entries", len(VERIFIED_CORRECTIONS) == 26)

for entry_index, gene, variant, code in VERIFIED_CORRECTIONS:
    decision = run_gate(gene, variant, code, erepo, mode=GateMode.PRODUCTION)
    check(
        f"[corrected] {gene} {variant} ({code}, idx {entry_index}) resolves to MET",
        decision.resolved_status == CriterionStatus.MET,
    )

# B-2b: verify behavior for PS4 (since this gate handles the criterion
# string generically, it works for PS4 as-is, not just PS3/BS3; confirmed
# with real data: RUNX1 c.601C>T, PS4:Met)
runx1_response = {
    "variantInterpretations": [{
        "gene": {"label": "RUNX1"},
        "caid": "CAR:CA16602487",
        "hgvs": ["NM_001754.4(RUNX1):c.601C>T (p.Arg201Ter)"],
        "guidelines": [{
            "outcome": {"label": "Pathogenic"},
            "agents": [{"evidenceCodes": [
                {"label": "PVS1", "status": "Met"},
                {"label": "PS4", "status": "Met"},
                {"label": "PM2_Supporting", "status": "Met"},
            ]}],
        }],
    }]
}
runx1_client = OfflineERepoClient([runx1_response])
decision_ps4 = run_gate("RUNX1", "c.601C>T", "PS4", runx1_client, mode=GateMode.PRODUCTION)
check("the gate works for PS4 unmodified (real RUNX1 data)", decision_ps4.resolved_status == CriterionStatus.MET)
check("pipeline is not started for PS4 in production mode either", decision_ps4.should_run_pipeline is False)

# B-3: also confirm that the 10 cases confirmed as "already correct" in
# section 9 still resolve to NOT_MET in production mode
VERIFIED_CORRECT_AS_IS = [
    ("985", "USH2A", "c.12295-3T>A", "PS3"),
    ("3162", "MYOC", "c.731G>T", "BS3"),
    ("3184", "MYOC", "c.1412A>G", "BS3"),
    ("3859", "RUNX1", "c.496C>T", "PS3"),
    ("7053", "HNF1A", "c.347C>T", "PS3"),
    ("9547", "GUCY2D", "c.3271C>T", "BS3"),
    ("4479", "MYOC", "c.1037G>C", "PS3"),
]
check("the already-correct list has 7 entries (the 3 PI3K-pathway cases are checked individually)", len(VERIFIED_CORRECT_AS_IS) == 7)
for entry_index, gene, variant, code in VERIFIED_CORRECT_AS_IS:
    decision = run_gate(gene, variant, code, erepo, mode=GateMode.PRODUCTION)
    check(
        f"[verified] {gene} {variant} ({code}, idx {entry_index}) resolves to NOT_MET as-is",
        decision.resolved_status == CriterionStatus.NOT_MET,
    )

# B-4: in validation mode, even a known variant must force the pipeline to
# run without ever exposing the ground truth downstream
decision_val = run_gate("RIT1", "c.268A>G", "PS3-Supporting", erepo, mode=GateMode.VALIDATION)
check("validation mode: pipeline forced to run even for a known variant", decision_val.should_run_pipeline is True)
check("validation mode: resolved_status is not exposed downstream (None)", decision_val.resolved_status is None)
check("validation mode: the ground truth is held as MET in hidden_ground_truth", decision_val._hidden_ground_truth == CriterionStatus.MET)


# ============================================================================
# C. Approved-assay cross-check: real-world examples actually encountered
# ============================================================================
print("\n[C] Approved-assay cross-check: based on real-world examples")

real_cases = [
    ("RASopathy VCEP", "RIT1", "elevated and prolonged ERK1/2 phosphorylation in HEK293 cells",
     "PS3", AssayApplicability.APPROVED),
    ("RASopathy VCEP", "PIK3CA", "increased PI3K activity and S6 phosphorylation",
     "PS3", AssayApplicability.ASSAY_NOT_APPROVED),
    ("RASopathy VCEP", "PIK3R2", "increased PI3K activity and S6 phosphorylation",
     "PS3", AssayApplicability.ASSAY_NOT_APPROVED),
    ("RASopathy VCEP", "RIT1", "n/a", "BS3", AssayApplicability.BS3_NOT_APPLICABLE),
    ("Brain Malformations VCEP", "PIK3CA", "increased PI3K activity and S6 phosphorylation",
     "PS3", AssayApplicability.NO_SPEC_AVAILABLE),  # the actually-correct VCEP, but its table is not yet registered
]
for vcep, gene, desc, criterion, expected in real_cases:
    result = check_approved_assay(vcep, gene, desc, criterion)
    check(f"{vcep}/{gene}/{criterion}: {desc[:30]}... -> {expected.value}", result.applicability == expected)


# ============================================================================
print(f"\n{'='*50}\n{passed} passed, {failed} failed (total {passed + failed})\n{'='*50}")
if failures:
    print("Failed cases:")
    for f_ in failures:
        print(f"  - {f_}")
sys.exit(1 if failed else 0)
