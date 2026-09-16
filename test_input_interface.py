"""Regression checks for the shared criterion input interface."""

from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.common import AggregatedJudgment
from acmg_pipeline.criteria import ps3_bs3, ps4, segregation
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


for module in (ps3_bs3, ps4, segregation):
    prompt = module.build_prompt(variant, clinical_note, "paper body")
    assert "MYH7" in prompt
    assert "c.2155C>T" in prompt
    assert "R719W" in prompt

assert clinical_note.proband.phenotype.clinical_features[0].hpo_id == "HP:0001639"

aggregated = AggregatedJudgment(
    contributions=[],
    aggregated_direction=ps3_bs3.OverallDirection.NOT_CLEAR,
    aggregation_hints=[],
)
evidence_line = build_evidence_line(aggregated, variant, criterion="PS3", vcep_name="Cardiomyopathy VCEP")
assert evidence_line["id"] == "evline:MYH7_c_2155C_T_PS3"
assert evidence_line["specifiedBy"]["methodType"] == "PS3"

print("shared input interface: passed")
