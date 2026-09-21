"""Regression checks for the shared criterion input interface."""

from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.common import AggregatedJudgment
from acmg_pipeline.constants import STUB_CODES
from acmg_pipeline.criteria import ps3_bs3, ps4, segregation, stubs
from acmg_pipeline.export import build_evidence_line
from acmg_pipeline.vcf_record import VariantRecord


variant = VariantRecord(
    chrom="14",
    pos=23425971,
    id="case3-var1",
    ref="G",
    alt="A",
    qual="99",
    filter="PASS",
    info={
        "GENE": "MYH7",
        "HGVSC": "c.2155C>T",
        "HGVSP": "p.(Arg719Trp)",
        "EQUIVALENTS": ["R719W"],
    },
)
clinical_note = ClinicalNoteExtraction.from_json({
    "proband": {
        "phenotype": {
            "affected_status": True,
            "clinical_features": [{"label": "Hypertrophic cardiomyopathy", "hpo_id": "HP:0001639"}],
        },
        "genotype": {"variant_status": True, "zygosity": "heterozygous"},
    },
})


criterion_modules = {
    "PS3": ps3_bs3,
    "BS3": ps3_bs3,
    "PS4": ps4,
    "PP1": segregation,
    "BS4": segregation,
}
for code, module in criterion_modules.items():
    prompt = module.build_prompt(variant, clinical_note, "paper body")
    assert "MYH7" in prompt
    assert "c.2155C>T" in prompt
    assert "R719W" in prompt
    assert code in {"PS3", "BS3", "PS4", "PP1", "BS4"}

stub_results = stubs.all_stub_evidence()
assert len(stub_results) == len(STUB_CODES)

assert clinical_note.proband.phenotype.clinical_features[0].hpo_id == "HP:0001639"

aggregated = AggregatedJudgment(
    contributions=[],
    aggregated_direction=ps3_bs3.OverallDirection.NOT_CLEAR,
    aggregation_hints=[],
)
evidence_line = build_evidence_line(
    aggregated, gene="MYH7", hgvsc="c.2155C>T", criterion="PS3",
    vcep_name="Cardiomyopathy VCEP", variant=variant,
)
assert evidence_line["id"] == "evline:MYH7_c_2155C_T_PS3"
assert evidence_line["specifiedBy"]["methodType"] == "PS3"

print("shared input interface: passed")
