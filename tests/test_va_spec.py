import unittest

from ga4gh.va_spec.acmg_2015 import VariantPathogenicityEvidenceLine

from acmg.core.models import CriterionResult, Status
from acmg.va_spec.mapper import export_document, to_evidence_line


VARIANT = {"assembly": "GRCh38", "chrom": "1", "pos": 2, "ref": "C", "alt": "T"}
EVIDENCE = [{"evidence_id": "test:frequency", "category": "population"}]


class VaSpecTests(unittest.TestCase):
    def test_met_is_validated_by_reference_model(self):
        result = CriterionResult("PM2", Status.MET, VARIANT, "rare", "supporting", "supports",
                                 "PM2_supporting", evidence=EVIDENCE)
        line = to_evidence_line(result)
        validated = VariantPathogenicityEvidenceLine.model_validate(line)
        self.assertEqual(validated.directionOfEvidenceProvided, "supports")
        self.assertEqual(validated.evidenceOutcome.primaryCoding.code.root, "PM2_supporting")

    def test_not_met_maps_to_machine_readable_neutral(self):
        result = CriterionResult("PM2", Status.NOT_MET, VARIANT, "not rare", None, "none",
                                 "PM2_not_met", evidence=EVIDENCE)
        line = to_evidence_line(result)
        self.assertEqual(line["directionOfEvidenceProvided"], "neutral")
        self.assertNotIn("strengthOfEvidenceProvided", line)
        VariantPathogenicityEvidenceLine.model_validate(line)

    def test_workflow_statuses_are_not_exported(self):
        for code, status in (("PVS1", Status.MANUAL_REVIEW), ("PP5", Status.DEPRECATED),
                             ("BP6", Status.DEPRECATED), ("PM1", Status.NOT_EVALUATED)):
            result = CriterionResult(code, status, VARIANT, "workflow state")
            self.assertIsNone(to_evidence_line(result))

    def test_all_supported_criteria_map(self):
        from acmg.core.models import CRITERIA
        from acmg.criteria.common import DEFAULT_STRENGTH

        for code in CRITERIA:
            if code in {"PVS1", "PP5", "BP6"}:
                continue
            strength = DEFAULT_STRENGTH[code]
            outcome = code
            result = CriterionResult(code, Status.MET, VARIANT, "met", strength,
                                     "disputes" if code.startswith("B") else "supports", outcome,
                                     evidence=EVIDENCE)
            with self.subTest(code=code):
                VariantPathogenicityEvidenceLine.model_validate(to_evidence_line(result))

    def test_referenced_evidence_id_conflict_is_rejected(self):
        first = CriterionResult("PM2", Status.NOT_MET, VARIANT, "a", None, "none", "PM2_not_met",
                                evidence=[{"evidence_id": "test:x", "value": 1}]).to_dict()
        second = CriterionResult("BA1", Status.NOT_MET, VARIANT, "b", None, "none", "BA1_not_met",
                                 evidence=[{"evidence_id": "test:x", "value": 2}]).to_dict()
        with self.assertRaisesRegex(ValueError, "Conflicting evidence"):
            export_document([{"record_id": "test:1", "variant": VARIANT,
                              "results": [first, second]}])
