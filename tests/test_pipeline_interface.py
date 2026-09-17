"""_evidence_from_line()'s GA4GH ACMG strength code -> Strength enum conversion.

Found 2026-09-17: a real PVS1 "very strong" VA-Spec line (GA4GH's own
space-separated ACMG strength coding, acmg_pipeline.automated_va_spec.STRENGTHS)
was passed directly to Strength("very strong") and raised ValueError, since
this project's own Strength.VERY_STRONG enum value is "very_strong" - this
had never been exercised before because no criterion had reached MET at
very_strong/stand_alone strength through the integrated pipeline until PVS1's
gates were wired up that same day.
"""

import unittest

from acmg_pipeline.classification import CriterionStatus, Strength
from acmg_pipeline.pipeline_interface import _evidence_from_line


def _line(status, ga4gh_strength_code=None):
    line = {"extensions": [{"name": "bh26AssessmentDetails", "value": {"status": status}}]}
    if ga4gh_strength_code is not None:
        line["strengthOfEvidenceProvided"] = {
            "primaryCoding": {"system": "ACMG Guidelines, 2015", "code": ga4gh_strength_code}
        }
    return line


class EvidenceFromLineTests(unittest.TestCase):
    def test_ga4gh_very_strong_maps_to_the_underscored_enum_value(self):
        evidence = _evidence_from_line("PVS1", _line("met", "very strong"))
        self.assertEqual(evidence.strength, Strength.VERY_STRONG)

    def test_ga4gh_standalone_maps_to_stand_alone(self):
        evidence = _evidence_from_line("BA1", _line("met", "standalone"))
        self.assertEqual(evidence.strength, Strength.STAND_ALONE)

    def test_the_remaining_strengths_round_trip_unchanged(self):
        for code, expected in (("strong", Strength.STRONG), ("moderate", Strength.MODERATE),
                              ("supporting", Strength.SUPPORTING)):
            with self.subTest(code=code):
                self.assertEqual(_evidence_from_line("PS1", _line("met", code)).strength, expected)

    def test_a_not_met_line_carries_no_strength_and_never_looks_at_the_code(self):
        evidence = _evidence_from_line("PS1", _line("not_met"))
        self.assertEqual(evidence.status, CriterionStatus.NOT_MET)
        self.assertIsNone(evidence.strength)

    def test_an_unrecognized_code_is_a_loud_error_not_a_silent_guess(self):
        with self.assertRaises(ValueError):
            _evidence_from_line("PVS1", _line("met", "extremely strong"))


if __name__ == "__main__":
    unittest.main()
