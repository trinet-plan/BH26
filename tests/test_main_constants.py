from acmg_pipeline.constants import (
    ALL_ACMG_CODES, AUTOMATED_CODES, IMPLEMENTED_CODES, LITERATURE_CODES,
    STUB_CODES, CriterionStatus,
)


def test_main_constants_cover_acmg_codes_once():
    assert len(ALL_ACMG_CODES) == 28
    assert len(IMPLEMENTED_CODES) == 19
    assert len(STUB_CODES) == 9
    assert IMPLEMENTED_CODES == LITERATURE_CODES | AUTOMATED_CODES
    assert IMPLEMENTED_CODES.isdisjoint(STUB_CODES)
    assert IMPLEMENTED_CODES | STUB_CODES == set(ALL_ACMG_CODES)


def test_only_three_criterion_statuses_exist():
    assert set(CriterionStatus) == {
        CriterionStatus.MET, CriterionStatus.NOT_MET, CriterionStatus.UNKNOWN,
    }
