from acmg_pipeline.clinical_note import (
    ClinicalFeature,
    ClinicalNoteExtraction,
    Family,
    Proband,
    ProbandPhenotype,
    Relative,
)
from acmg_pipeline.criteria.pp4_pp1_bs4 import (
    PP4ReferenceRecord,
    evaluate_locus_evidence,
)


def ref(*, locus_model="heterogeneous", yield_=0.70):
    return PP4ReferenceRecord(
        reference_id="demo",
        gene="GENE1",
        phenotype_label="Demo phenotype",
        phenotype_hpo=("HP:0000001", "HP:0000002"),
        locus_model=locus_model,
        diagnostic_yield=yield_,
        testing_method="sequencing",
        source_citation="demo source",
    )


def case(relatives=None, inheritance="autosomal dominant"):
    return ClinicalNoteExtraction(
        proband=Proband(
            phenotype=ProbandPhenotype(
                affected_status=True,
                clinical_features=[
                    ClinicalFeature("feature A", "HP:0000001"),
                    ClinicalFeature("feature B", "HP:0000002"),
                ],
            )
        ),
        family=Family(
            inheritance_pattern=inheritance,
            relatives=relatives or [],
        ),
    )


# 1) PP4 works without family data.
r = evaluate_locus_evidence(case(), ref(), method_comparable=True)
assert r.pp4 is not None and r.pp4.points == 4.0
assert r.segregation.pp1_points_raw == 0.0
assert r.combined_positive_locus_points == 4.0

# 2) Heterogeneous disease: PP4 + PP1, cap at +5.
r = evaluate_locus_evidence(
    case([Relative("sister", affected_status=True, variant_status=True)]),
    ref(),
    method_comparable=True,
)
assert r.pp4 is not None and r.pp4.points == 4.0
assert r.segregation.pp1_points_raw == 1.0
assert r.combined_positive_locus_points == 5.0

# 3) High-yield homogeneous disease: PP1 is not added to PP4.
r = evaluate_locus_evidence(
    case([Relative("sister", affected_status=True, variant_status=True)]),
    ref(locus_model="homogeneous", yield_=0.958),
    method_comparable=True,
)
assert r.pp4 is not None and r.pp4.points == 5.0
assert r.segregation.pp1_points_raw == 1.0
assert r.segregation.pp1_points_used == 0.0
assert r.combined_positive_locus_points == 5.0

# 4) Clear AD non-segregation under explicit assumptions -> BS4.
r = evaluate_locus_evidence(
    case([Relative("sister", affected_status=True, variant_status=False)]),
    ref(),
    method_comparable=True,
    low_phenocopy=True,
)
assert r.segregation.bs4_met is True
assert r.bs4_points == -4.0
assert r.combined_positive_locus_points == 0.0
assert r.net_locus_points == -4.0

# 5) Missing curated reference HPO definition -> PP4 is not auto-scored,
#    but segregation can still be evaluated.
reference_no_hpo = PP4ReferenceRecord(
    reference_id="no_hpo",
    gene="GENE1",
    phenotype_label="Demo phenotype",
    phenotype_hpo=(),
    locus_model="heterogeneous",
    diagnostic_yield=0.70,
    testing_method="sequencing",
    source_citation="demo source",
)
r = evaluate_locus_evidence(
    case([Relative("sister", affected_status=True, variant_status=True)]),
    reference_no_hpo,
    method_comparable=True,
)
assert r.pp4 is None
assert r.segregation.pp1_points_raw == 1.0
assert r.combined_positive_locus_points == 1.0

print("5 passed, 0 failed")
