"""acmg_pipeline.export._direction_of_evidence() - directionOfEvidenceProvided
means "for/against pathogenicity", not "for/against the criterion under test".

Found 2026-09-18: a real BS3 line (LDLR c.2575G>A, live run_integrated_validation_64
run) whose literature evidence genuinely established BS3 (benign) raised
_check_acmg_semantics()'s ValueError, because this function unconditionally
returned SUPPORTS whenever the established direction matched the criterion
being tested - correct for a pathogenic criterion (PS3 evidence supporting a
PS3 test) but backwards for a benign one (BS3 evidence disputes pathogenicity,
it does not support it). acmg_pipeline.criteria.pp1_bs4_pp4_engine's own
build_evidence_line() already got this right via PATHOGENIC_CODES; export.py's
literature-engine counterpart did not.
"""

import unittest

from ga4gh.va_spec.base.core import Direction

from acmg_pipeline.criteria.ps3_bs3 import OverallDirection
from acmg_pipeline.export import _direction_of_evidence


class DirectionOfEvidenceTests(unittest.TestCase):
    def test_ps3_evidence_for_a_ps3_test_supports_pathogenicity(self):
        self.assertEqual(
            _direction_of_evidence(OverallDirection.PS3, "PS3"), Direction.SUPPORTS)

    def test_bs3_evidence_for_a_bs3_test_disputes_pathogenicity(self):
        self.assertEqual(
            _direction_of_evidence(OverallDirection.BS3, "BS3"), Direction.DISPUTES)

    def test_bs3_evidence_found_while_testing_ps3_disputes_the_test(self):
        self.assertEqual(
            _direction_of_evidence(OverallDirection.BS3, "PS3"), Direction.DISPUTES)

    def test_ps3_evidence_found_while_testing_bs3_disputes_the_test(self):
        self.assertEqual(
            _direction_of_evidence(OverallDirection.PS3, "BS3"), Direction.DISPUTES)

    def test_not_clear_is_neutral_regardless_of_criterion(self):
        for criterion in ("PS3", "BS3"):
            with self.subTest(criterion=criterion):
                self.assertEqual(
                    _direction_of_evidence(OverallDirection.NOT_CLEAR, criterion),
                    Direction.NEUTRAL)


if __name__ == "__main__":
    unittest.main()
