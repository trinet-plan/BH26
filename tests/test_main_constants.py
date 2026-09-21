from acmg_pipeline.constants import (
    ALL_ACMG_CODES, AUTOMATED_CODES, DE_NOVO_CODES, IMPLEMENTED_CODES,
    LITERATURE_CODES, PHENOTYPE_SEGREGATION_CODES, STUB_CODES, CriterionStatus,
)


def test_main_constants_cover_acmg_codes_once():
    assert len(ALL_ACMG_CODES) == 28
    # +3 (PP1/BS4/PP4) on 2026-09-17 once acmg_pipeline.criteria.
    # pp1_bs4_pp4_engine connected the ClinGen 2024 PP4+PP1/BS4 evaluator,
    # +2 (PS2/PM6) on 2026-09-19 once acmg_pipeline.criteria.ps2_pm6
    # connected a rule-based de-novo evaluator - see each module's own
    # docstring.
    assert len(IMPLEMENTED_CODES) == 24
    assert len(STUB_CODES) == 4
    assert IMPLEMENTED_CODES == (
        LITERATURE_CODES | AUTOMATED_CODES | PHENOTYPE_SEGREGATION_CODES | DE_NOVO_CODES
    )
    assert IMPLEMENTED_CODES.isdisjoint(STUB_CODES)
    assert IMPLEMENTED_CODES | STUB_CODES == set(ALL_ACMG_CODES)


def test_only_three_criterion_statuses_exist():
    assert set(CriterionStatus) == {
        CriterionStatus.MET, CriterionStatus.NOT_MET, CriterionStatus.UNKNOWN,
    }
