import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.criteria.ps3_bs3 import (
    PS3BS3Judgment, finalize, MatchStatus, OverallDirection,
    detect_experiment_conflict, detect_narrative_categorical_mismatch,
    detect_empty_experiments_with_definitive_direction,
    detect_no_quantitative_evidence,
    PaperContribution, aggregate_multi_paper_results,
)

passed = 0
failed = 0


def check(label, cond):
    global passed, failed
    if cond:
        passed += 1
        print(f"  OK   {label}")
    else:
        failed += 1
        print(f"  FAIL {label}")


minimal_json = {
    "variant_matching": {"match_status": "matched", "confidence": "high", "notes": ""},
    "experiments": [{
        "assay_type": "ERK Activation Assay",
        "experimental_system": "HEK293T",
        "readout": "pERK/ERK",
        "comparator": "WT",
        "result_direction": "functionally_abnormal",
    }],
    "multiple_assay_types_in_paper": False,
    "single_study_only": True,
    "overall_evidence": {"direction": "PS3", "rationale": "test"},
}

print("[1] JSON parsing")
j = PS3BS3Judgment.from_json(minimal_json)
check("match_status parses correctly", j.variant_matching.match_status == MatchStatus.MATCHED)
check("overall_evidence.direction parses correctly", j.overall_evidence.direction == OverallDirection.PS3)
check("1 experiment parses correctly", len(j.experiments) == 1)

print("\n[2] Hint generation from the approved-assay cross-check")
result = finalize(j, gene="RIT1", vcep_name="RASopathy VCEP", criterion="PS3")
check("an APPROVED info hint is generated", any("matches an approved assay" in h.message for h in result.curator_hints))
check("a single_study_only caution hint is generated", any("single study only" in h.message for h in result.curator_hints))

print("\n[3] Hint generation for an unapproved assay")
j2 = PS3BS3Judgment.from_json({**minimal_json, "experiments": [{
    "assay_type": "totally unknown assay type",
    "experimental_system": "X", "readout": "Y", "comparator": "Z",
    "result_direction": "functionally_abnormal",
}]})
result2 = finalize(j2, gene="RIT1", vcep_name="RASopathy VCEP", criterion="PS3")
check("an ASSAY_NOT_APPROVED warning hint is generated", any("did not match any" in h.message for h in result2.curator_hints))

print("\n[4] No hints when unsuccessful")
j3 = PS3BS3Judgment.from_json({
    "variant_matching": {"match_status": "unsuccessful"},
    "experiments": [],
    "overall_evidence": {"direction": "not_clear"},
})
result3 = finalize(j3, gene="RIT1", vcep_name="RASopathy VCEP", criterion="PS3")
check("0 hints when unsuccessful", len(result3.curator_hints) == 0)

print("\n[5] BS3-not-applicable hint")
j4 = PS3BS3Judgment.from_json({**minimal_json, "overall_evidence": {"direction": "BS3", "rationale": "t"}})
result4 = finalize(j4, gene="RIT1", vcep_name="RASopathy VCEP", criterion="BS3")
check("a BS3_NOT_APPLICABLE warning hint is generated", any("BS3 itself is not applicable" in h.message for h in result4.curator_hints))

print("\n[6] Defensive parsing against the 'PS3(damaging)' pattern seen in real gemma-4 output (regression test)")
gemma_raw_json = {
    "variant_matching": {"match_status": "matched", "confidence": "high", "notes": ""},
    "experiments": [{
        "assay_type": "test",
        "experimental_system": "test",
        "readout": "test",
        "comparator": "test",
        "result_direction": "functionally_abnormal(gain-of-function)",  # simulating the same kind of habit
    }],
    "multiple_assay_types_in_paper": False,
    "single_study_only": True,
    "overall_evidence": {"direction": "PS3(damaging)", "rationale": "test"},  # the actual format gemma-4 produced
}
j_gemma = PS3BS3Judgment.from_json(gemma_raw_json)
check("'PS3(damaging)' parses as OverallDirection.PS3", j_gemma.overall_evidence.direction == OverallDirection.PS3)
check("'functionally_abnormal(gain-of-function)' also parses as a ResultDirection",
      j_gemma.experiments[0].result_direction.value == "functionally_abnormal")

print("\n[7] Automatic detection of experiment conflict -> forced override to not_clear (reproducing the real PTEN misjudgment)")
pten_gemma_like = PS3BS3Judgment.from_json({
    "variant_matching": {"match_status": "matched", "confidence": "high", "notes": ""},
    "experiments": [
        {
            "assay_type": "VAMP-seq", "experimental_system": "HEK293T",
            "readout": "abundance score", "comparator": "WT",
            "result_direction": "functionally_normal",
        },
        {
            "assay_type": "Akt phosphorylation assay", "experimental_system": "HEK293T",
            "readout": "pAkt", "comparator": "WT",
            "result_direction": "functionally_abnormal",
        },
    ],
    "multiple_assay_types_in_paper": True,
    "single_study_only": True,
    # gemma-4 actually concluded PS3 here and was wrong (ground truth: not_met)
    "overall_evidence": {"direction": "PS3", "rationale": "dominant-negative activity observed"},
})
result7 = finalize(pten_gemma_like, gene="PTEN", vcep_name=None, criterion="PS3")
check("the LLM's raw judgment stays PS3 (kept for audit)", pten_gemma_like.overall_evidence.direction == OverallDirection.PS3)
check("the adopted judgment (effective_direction) is forced to not_clear", result7.effective_direction == OverallDirection.NOT_CLEAR)
check("a forced-override warning hint is generated", any("forced to not_clear" in h.message for h in result7.curator_hints))

print("\n[8] Automatic detection of a narrative/categorical-label mismatch (reproducing the second real PTEN misjudgment)")
# The actual raw gemma-4 JSON obtained when the user ran this on 2026-09-15
pten_gemma_run2 = PS3BS3Judgment.from_json({
    "variant_matching": {"match_status": "matched"},
    "experiments": [
        {
            "assay_type": "VAMP-seq (Steady-state protein abundance assay via fluorescent reporters)",
            "experimental_system": "Engineered HEK 293T cells",
            "readout": "EGFP:mCherry ratio (abundance score)",
            "comparator": "Wild-type PTEN",
            "result_direction": "functionally_abnormal",
        },
        {
            "assay_type": "Akt phosphorylation assay (Functional signaling assay)",
            "experimental_system": "Cell lines (presence of endogenous wild type PTEN)",
            "readout": "Increased Akt phosphorylation (dominant negative effect)",
            "comparator": "Wild-type PTEN",
            "result_direction": "functionally_abnormal",
        },
    ],
    "multiple_assay_types_in_paper": True,
    "single_study_only": True,
    "overall_evidence": {
        "direction": "PS3",
        "rationale": (
            "The variant p.Pro38Ser was explicitly tested using VAMP-seq, showing a "
            "WT-like abundance score (1.14), and a functional assay demonstrating a "
            "dominant negative effect through increased Akt phosphorylation."
        ),
    },
})
check("the simple experiment-conflict check does not catch this (as intended)", not detect_experiment_conflict(pten_gemma_run2))
check("it is caught as a narrative/categorical-label mismatch", detect_narrative_categorical_mismatch(pten_gemma_run2))
result8 = finalize(pten_gemma_run2, gene="PTEN", vcep_name=None, criterion="PS3")
check("the adopted judgment is forced to not_clear", result8.effective_direction == OverallDirection.NOT_CLEAR)

print("\n[9] Automatic detection of a definitive direction with zero extracted experiments (reproducing the real TP53 misjudgment)")
# Reconstructed from the actual gemma-4 run (2026-09-15, English prompt): the
# rationale explicitly admitted the variant-specific data was only in
# supplementary material, with no experiments extracted, yet the LLM still
# confidently concluded PS3. Ground truth was BS3 = met.
tp53_gemma_run = PS3BS3Judgment.from_json({
    "variant_matching": {"match_status": "heuristic"},
    "experiments": [],
    "multiple_assay_types_in_paper": False,
    "single_study_only": True,
    "overall_evidence": {
        "direction": "PS3",
        "rationale": (
            "The target variant (p.Lys292Arg) is part of a comprehensive "
            "MITE library... Although the specific Z-score for K292R is in "
            "the supplementary data rather than the main text, the paper's "
            "methodology... support a damaging effect."
        ),
    },
})
check("detected as a definitive direction with no supporting experiments",
      detect_empty_experiments_with_definitive_direction(tp53_gemma_run))
result9 = finalize(tp53_gemma_run, gene="TP53", vcep_name=None, criterion="PS3")
check("the adopted judgment is forced to not_clear", result9.effective_direction == OverallDirection.NOT_CLEAR)
check("an explanatory warning hint is generated", any("no concrete supporting experiment" in h.message for h in result9.curator_hints))

print("\n[10] Defensive fallback when variant_matching is flattened to a bare string (reproducing the real MYH7/PTEN crash)")
# Reconstructed from the actual gemma-4 run (2026-09-15, English prompt):
# variant_matching was returned as "matched" (a bare string) instead of
# {"match_status": "matched", ...}, causing "string indices must be
# integers, not 'str'" when the old code tried vm["match_status"].
flattened_vm_json = {
    "variant_matching": "matched",
    "experiments": [{
        "assay_type": "in vitro motility and biochemical assays (isolated myosin)",
        "experimental_system": "mouse",
        "readout": "actin filament velocities and actin-activated myosin ATPases",
        "comparator": "wild type controls",
        "result_direction": "functionally_abnormal",
    }],
    "multiple_assay_types_in_paper": False,
    "single_study_only": True,
    "overall_evidence": {"direction": "PS3", "rationale": "test"},
}
j10 = PS3BS3Judgment.from_json(flattened_vm_json)
check("a bare-string variant_matching no longer crashes and parses as MATCHED",
      j10.variant_matching.match_status == MatchStatus.MATCHED)

print("\n[11] Negation-aware narrative check: 'lack WT-like activity' must NOT count as a normal finding (real TP53 data)")
# The exact key_findings text from the real gemma-4 run (2026-09-15) that
# caused a false-positive narrative-mismatch trigger for the wrong reason.
tp53_real_run = PS3BS3Judgment.from_json({
    "variant_matching": {
        "match_status": "heuristic",
        "notes": "The target variant p.Lys292Arg (K292R) is contained within the MITE library...",
    },
    "experiments": [
        {
            "assay_type": "Pooled positive-selection screen (Nutlin-3)",
            "experimental_system": "p53-null A549 human lung carcinoma cell line",
            "readout": "Allele enrichment/depletion (Z-score) via massively-parallel sequencing",
            "comparator": "Wild-type p53 and silent mutations",
            "result_direction": "functionally_abnormal",
            "key_findings": [
                "The screen identifies alleles that interfere with p53 function (dominant-negative) or lack function (loss-of-function) by measuring resistance to Nutlin-3 in a p53-null background.",
            ],
        },
        {
            "assay_type": "Pooled positive-selection screen (Etoposide)",
            "experimental_system": "p53-null A549 human lung carcinoma cell line",
            "readout": "Allele enrichment/depletion (Z-score) via massively-parallel sequencing",
            "comparator": "Wild-type p53 and silent mutations",
            "result_direction": "functionally_abnormal",
            "key_findings": [
                "The screen identifies alleles that lack WT-like activity (loss-of-function) by measuring depletion under DNA damage-induced stress.",
            ],
        },
    ],
    "multiple_assay_types_in_paper": True,
    "single_study_only": True,
    "overall_evidence": {
        "direction": "PS3",
        "rationale": (
            "The paper utilizes a large-scale, high-throughput MITE library "
            "screen in A549 cells to quantify the functional impact of "
            "missense mutations on p53 activity (Nutlin-3 and Etoposide "
            "response). The target variant K292R is part of this "
            "comprehensive functional dataset showing loss of "
            "function/dominant-negative activity."
        ),
    },
})
check("'lack WT-like activity' is correctly NOT counted as a normal finding (negation-aware)",
      not detect_narrative_categorical_mismatch(tp53_real_run))
check("this case IS caught instead by the no-quantitative-evidence check",
      detect_no_quantitative_evidence(tp53_real_run))
result11 = finalize(tp53_real_run, gene="TP53", vcep_name="TP53 VCEP", criterion="PS3")
check("the adopted judgment is forced to not_clear", result11.effective_direction == OverallDirection.NOT_CLEAR)

print("\n[12] The no-quantitative-evidence check must NOT fire on genuine variant-specific data (real MYH7 data)")
myh7_real_run = PS3BS3Judgment.from_json({
    "variant_matching": {"match_status": "matched"},
    "experiments": [
        {
            "assay_type": "In vivo whole heart mechanics (isovolumic contraction)",
            "experimental_system": "transgenic mouse (S532P/+)",
            "readout": "LV developed pressure",
            "comparator": "wild-type (+/+) mice",
            "result_direction": "functionally_abnormal",
            "key_findings": [
                "S532P/+ mice demonstrated depressed contractile function with lower developed pressure at 31 and 35 \u03bcL balloon volumes compared to controls.",
            ],
        },
    ],
    "multiple_assay_types_in_paper": True,
    "single_study_only": True,
    "overall_evidence": {"direction": "PS3", "rationale": "test"},
})
check("real MYH7 data (contains '31 and 35') is NOT flagged as lacking quantitative evidence",
      not detect_no_quantitative_evidence(myh7_real_run))

print("\n[13] The no-quantitative-evidence check must NOT fire on genuine variant-specific data (real PTEN data)")
pten_real_run = PS3BS3Judgment.from_json({
    "variant_matching": {"match_status": "matched"},
    "experiments": [
        {
            "assay_type": "VAMP-seq (Variant Abundance by Massively Parallel Sequencing)",
            "experimental_system": "engineered HEK 293T cell line",
            "readout": "steady-state protein abundance (EGFP:mCherry ratio)",
            "comparator": "wild type PTEN",
            "result_direction": "functionally_normal",
            "key_findings": [
                "p.Pro38Ser had a slightly higher abundance score than WT (1.14)",
                "classified as WT-like abundance",
            ],
        },
        {
            "assay_type": "Akt phosphorylation assay",
            "experimental_system": "cell culture (presence of endogenous wild type PTEN)",
            "readout": "Akt phosphorylation levels",
            "comparator": "wild type PTEN",
            "result_direction": "functionally_abnormal",
            "key_findings": [
                "p.Pro38Ser drove increased Akt phosphorylation, suggesting a dominant negative effect",
            ],
        },
    ],
    "multiple_assay_types_in_paper": True,
    "single_study_only": True,
    "overall_evidence": {"direction": "PS3", "rationale": "test"},
})
check("real PTEN data (contains '1.14') is NOT flagged as lacking quantitative evidence",
      not detect_no_quantitative_evidence(pten_real_run))
check("real PTEN data is still caught by the experiment-conflict check (unaffected by this change)",
      detect_experiment_conflict(pten_real_run))

print("\n[14] Multi-paper aggregation: 0 papers yield usable evidence")
empty_agg = aggregate_multi_paper_results([
    PaperContribution("11111111", finalize(
        PS3BS3Judgment.from_json({"variant_matching": {"match_status": "unsuccessful"}, "experiments": [], "overall_evidence": {"direction": "not_clear"}}),
        gene="X", vcep_name=None, criterion="PS3",
    )),
])
check("0 usable papers -> not_clear", empty_agg.aggregated_direction == OverallDirection.NOT_CLEAR)
check("an explanatory info hint is generated", any("None of the" in h.message for h in empty_agg.aggregation_hints))

print("\n[15] Multi-paper aggregation: a single usable paper (real MYH7 result)")
myh7_result = finalize(
    PS3BS3Judgment.from_json({
        "variant_matching": {"match_status": "matched"},
        "experiments": [{
            "assay_type": "In vivo whole heart mechanics", "experimental_system": "mouse",
            "readout": "LV developed pressure", "comparator": "wild-type",
            "result_direction": "functionally_abnormal",
            "key_findings": ["lower developed pressure at 31 and 35 uL balloon volumes"],
        }],
        "multiple_assay_types_in_paper": False,
        "single_study_only": True,
        "overall_evidence": {"direction": "PS3", "rationale": "test"},
    }),
    gene="MYH7", vcep_name="Cardiomyopathy VCEP", criterion="PS3",
)
single_agg = aggregate_multi_paper_results([PaperContribution("23313350", myh7_result)])
check("1 usable paper -> that paper's direction is adopted", single_agg.aggregated_direction == OverallDirection.PS3)
check("a caution hint about the single-paper count is generated",
      any("Only 1 of 1 paper" in h.message for h in single_agg.aggregation_hints))

print("\n[16] Multi-paper aggregation: 2 papers agree -> stronger support")
agree_agg = aggregate_multi_paper_results([
    PaperContribution("23313350", myh7_result),
    PaperContribution("99999999", myh7_result),  # simulating a second independent paper with the same conclusion
])
check("2 agreeing papers -> the agreed direction is adopted", agree_agg.aggregated_direction == OverallDirection.PS3)
check("an info hint about agreement across papers is generated",
      any("2 independent papers agree" in h.message for h in agree_agg.aggregation_hints))

print("\n[17] Multi-paper aggregation: papers disagree -> forced to not_clear")
bs3_result = finalize(
    PS3BS3Judgment.from_json({
        "variant_matching": {"match_status": "matched"},
        "experiments": [{
            "assay_type": "some assay", "experimental_system": "cells",
            "readout": "activity", "comparator": "WT",
            "result_direction": "functionally_normal",
            "key_findings": ["activity was 98% of wild-type (p=0.6)"],
        }],
        "multiple_assay_types_in_paper": False,
        "single_study_only": True,
        "overall_evidence": {"direction": "BS3", "rationale": "test"},
    }),
    gene="MYH7", vcep_name="Cardiomyopathy VCEP", criterion="BS3",
)
conflict_agg = aggregate_multi_paper_results([
    PaperContribution("23313350", myh7_result),   # PS3
    PaperContribution("88888888", bs3_result),     # BS3 - conflicts
])
check("conflicting papers -> forced to not_clear", conflict_agg.aggregated_direction == OverallDirection.NOT_CLEAR)
check("a warning hint naming the conflicting PMIDs is generated",
      any("23313350=PS3" in h.message and "88888888=BS3" in h.message for h in conflict_agg.aggregation_hints))

print(f"\n{'='*40}\n{passed} passed, {failed} failed\n{'='*40}")
sys.exit(1 if failed else 0)
