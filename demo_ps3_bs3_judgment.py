"""
demo_ps3_bs3_judgment.py
An end-to-end demo using two pieces of real data.

[Important note] At this stage, the actual LLM API call is not yet wired up
(ps3_bs3_judgment.build_prompt() only returns the prompt string). Because of
this, the demo feeds in, directly as JSON, "the structured output the LLM
should return" - matching exactly the manual judgment made in design doc
section 6 by actually reading the paper text. In other words, what is being
verified here is "does the post-processing (approved-assay cross-check,
curator-hint generation) work correctly", not "can the LLM judge correctly".
Verifying the latter requires an actual API connection (the next step).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.criteria.ps3_bs3 import PS3BS3Judgment, finalize, build_prompt

# ----------------------------------------------------------------------------
# Example 1: RIT1 c.268A>G (p.Met90Val), idx 8336
# JSON representation of what was found by actually reading the paper text
# (pilot/rit1_fulltext.txt).
# Confirmed ground truth: PS3-Supporting = Met (live ClinGen ERepo lookup,
# design doc section 9)
# ----------------------------------------------------------------------------
rit1_judgment_json = {
    "variant_matching": {
        "match_status": "matched",
        "match_type": "protein_notation",
        "confidence": "high",
        "notes": "Explicitly tested as one of 6 NS-associated RIT1 variants",
    },
    "experiments": [
        {
            "assay_type": "ERK Activation Assay (phospho-ERK1/2 western blot after serum stimulation)",
            "experimental_system": "HEK293T cells, transient transfection",
            "readout": "phospho-ERK1/2 / total ERK1/2 ratio (5/15/30 min after serum stimulation)",
            "comparator": "RIT1 wild-type",
            "result_direction": "functionally_abnormal",
            "key_findings": [
                "M90V showed a statistically significant increase in pERK1/2 vs WT at 5, 15, and 30 min after serum stimulation",
            ],
        },
        {
            "assay_type": "PAK1 binding / co-immunoprecipitation (PAK[CRIB] pulldown)",
            "experimental_system": "HEK293T cells",
            "readout": "RIT1 co-precipitation via GST-PAK[CRIB] pulldown and co-IP with endogenous PAK1",
            "comparator": "RIT1 wild-type",
            "result_direction": "functionally_abnormal",
            "key_findings": [
                "M90V showed an increased trend in PAK1 co-precipitation across all experiments, with a statistically significant increase in co-IP under serum-starved conditions",
            ],
        },
    ],
    "multiple_assay_types_in_paper": True,
    "single_study_only": True,
    "overall_evidence": {
        "direction": "PS3",
        "strength_hint": "not_clear",
        "rationale": (
            "M90V significantly increases ERK1/2 phosphorylation after serum "
            "stimulation relative to wild-type, showing a gain-of-function "
            "phenotype similar to the other 5 NS-associated RIT1 variants. "
            "The PAK1-binding assay corroborates this. Supports PS3 (damaging)."
        ),
    },
}

# ----------------------------------------------------------------------------
# Example 2: MT-TN m.5690A>G, idx 4870
# From the paper text (pilot/mtTN_fulltext.txt). Single-fiber segregation
# analysis.
# Confirmed ground truth: PS3-Supporting = Met (live ClinGen ERepo lookup,
# design doc section 9). No approved-assay table is registered yet for the
# Mitochondrial Diseases VCEP.
# ----------------------------------------------------------------------------
mt_tn_judgment_json = {
    "variant_matching": {
        "match_status": "matched",
        "match_type": "hgvsc",
        "confidence": "high",
        "notes": "The target variant is directly examined via single-fiber analysis of a patient-derived muscle biopsy",
    },
    "experiments": [
        {
            "assay_type": "Single fiber PCR segregation analysis",
            "experimental_system": "Patient-derived skeletal-muscle biopsy, single muscle fibers",
            "readout": "Heteroplasmy rate in COX-deficient fibers vs COX-normal fibers",
            "comparator": "COX-positive fibers (within the same patient)",
            "result_direction": "functionally_abnormal",
            "key_findings": [
                "Heteroplasmy rate was 86.8% in COX-deficient fibers vs 27.9% in COX-positive fibers (p<0.0001)",
                "The paper's own conclusion classifies this variant as 'definitely pathogenic'",
            ],
        },
    ],
    "multiple_assay_types_in_paper": False,
    "single_study_only": True,
    "overall_evidence": {
        "direction": "PS3",
        "strength_hint": "not_clear",
        "rationale": (
            "There is a significant difference in heteroplasmy rate between "
            "COX-deficient and COX-positive fibers (p<0.0001), showing clear "
            "segregation with the biochemical defect. Supports PS3."
        ),
    },
}


def run_example(title: str, judgment_json: dict, gene: str, vcep_name, criterion: str):
    print(f"\n{'='*70}\n{title}\n{'='*70}")
    judgment = PS3BS3Judgment.from_json(judgment_json)
    result = finalize(judgment, gene=gene, vcep_name=vcep_name, criterion=criterion)

    print(f"Judgment: direction={judgment.overall_evidence.direction.value}, "
          f"multiple_assay_types={judgment.multiple_assay_types_in_paper}, "
          f"single_study_only={judgment.single_study_only}")
    print(f"\nCurator hints ({len(result.curator_hints)}):")
    for hint in result.curator_hints:
        print(f"  [{hint.severity:7s}] {hint.message}")


if __name__ == "__main__":
    run_example(
        "Example 1: RIT1 c.268A>G (p.Met90Val) - RASopathy VCEP",
        rit1_judgment_json,
        gene="RIT1",
        vcep_name="RASopathy VCEP",
        criterion="PS3",
    )

    run_example(
        "Example 2: MT-TN m.5690A>G - Mitochondrial Diseases VCEP (not yet registered)",
        mt_tn_judgment_json,
        gene="MT-TN",
        vcep_name="Mitochondrial Diseases VCEP",
        criterion="PS3",
    )

    # Also check what the prompt itself looks like (just the beginning)
    print(f"\n{'='*70}\nSample output of build_prompt() (first 500 chars)\n{'='*70}")
    prompt = build_prompt(
        gene="RIT1", hgvsc="c.268A>G", hgvsp="p.Met90Val",
        equivalents=["M90V", "p.(Met90Val)"],
        full_text="(the full paper text would go here)",
    )
    print(prompt[:500] + "...")
