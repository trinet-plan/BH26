import unittest
from copy import deepcopy

from acmg.core.models import CriterionResult, Status
from acmg.va_spec.mapper import SCHEMA_ID, export_document, to_evidence_line, validate_1_0_1


VARIANT = {"assembly": "GRCh38", "chrom": "1", "pos": 2, "ref": "C", "alt": "T"}
EVIDENCE = [{
    "evidence_id": "https://example.org/evidence/frequency",
    "category": "population",
    "source": "gnomAD",
    "source_version": "4.1.1",
    "population": "exome:global",
    "AC": 1,
    "AN": 100000,
    "AF": "0.00001",
    "quality_status": "PASS",
    "callable": True,
    "filters": [],
    "variant_flags": [],
    "retrieved_at": "2026-09-15T00:00:00Z",
}]


class VaSpecTests(unittest.TestCase):
    def test_met_is_validated_by_reference_model(self):
        result = CriterionResult("PM2", Status.MET, VARIANT, "rare", "supporting", "supports",
                                 "PM2_supporting", evidence=EVIDENCE)
        line = to_evidence_line(result)
        validate_1_0_1(line, "PM2")
        self.assertEqual(line["directionOfEvidenceProvided"], "supports")
        self.assertEqual(line["specifiedBy"]["methodType"], "PM2")
        self.assertEqual(line["evidenceOutcome"]["primaryCoding"]["code"], "PM2_supporting")
        self.assertNotIn("type", line["evidenceOutcome"])
        self.assertNotIn("type", line["strengthOfEvidenceProvided"])
        study = line["hasEvidenceItems"][0]
        self.assertEqual(study["type"], "CohortAlleleFrequencyStudyResult")
        self.assertEqual(study["focusAlleleCount"], 1)
        self.assertEqual(study["locusAlleleCount"], 100000)
        self.assertEqual(study["focusAlleleFrequency"], 0.00001)
        self.assertEqual(study["sourceDataSet"]["type"], "DataSet")
        self.assertEqual(study["cohort"]["type"], "StudyGroup")
        self.assertEqual(study["specifiedBy"]["type"], "Method")

    def test_not_met_maps_to_machine_readable_neutral(self):
        result = CriterionResult("PM2", Status.NOT_MET, VARIANT, "not rare", None, "none",
                                 "PM2_not_met", evidence=EVIDENCE)
        line = to_evidence_line(result)
        self.assertEqual(line["directionOfEvidenceProvided"], "neutral")
        self.assertNotIn("strengthOfEvidenceProvided", line)
        validate_1_0_1(line, "PM2")

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
                validate_1_0_1(to_evidence_line(result), code)

    def test_document_records_pinned_official_schema(self):
        document = export_document([])
        self.assertEqual(document["validated_by"]["schema_version"], "1.0.1")
        self.assertEqual(document["validated_by"]["schema_id"], SCHEMA_ID)

    def test_1_0_1_contract_rejects_newer_gks_discriminator_and_mismatched_method(self):
        result = CriterionResult("PM2", Status.NOT_MET, VARIANT, "not rare", None, "none",
                                 "PM2_not_met", evidence=EVIDENCE)
        line = to_evidence_line(result)
        invalid = deepcopy(line)
        invalid["evidenceOutcome"]["type"] = "MappableConcept"
        with self.assertRaisesRegex(ValueError, "schema validation failed"):
            validate_1_0_1(invalid, "PM2")
        invalid = deepcopy(line)
        invalid["specifiedBy"]["methodType"] = "BA1"
        with self.assertRaisesRegex(ValueError, "criterion mapping mismatch"):
            validate_1_0_1(invalid, "PM2")

    def test_referenced_evidence_id_conflict_is_rejected(self):
        first = CriterionResult("PM2", Status.NOT_MET, VARIANT, "a", None, "none", "PM2_not_met",
                                evidence=[{"evidence_id": "test:x", "value": 1}]).to_dict()
        second = CriterionResult("BA1", Status.NOT_MET, VARIANT, "b", None, "none", "BA1_not_met",
                                 evidence=[{"evidence_id": "test:x", "value": 2}]).to_dict()
        with self.assertRaisesRegex(ValueError, "Conflicting evidence"):
            export_document([{"record_id": "test:1", "variant": VARIANT,
                              "results": [first, second]}])
