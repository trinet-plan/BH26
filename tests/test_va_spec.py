from acmg_pipeline.constants import CriterionStatus
import unittest
from copy import deepcopy

from acmg_pipeline.automated_core.models import CriterionResult
from acmg_pipeline.automated_va_spec import (
    SCHEMA_ID,
    export_document,
    to_evidence_line,
    validate_1_0_1,
    validate_envelope,
)


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
        result = CriterionResult("PM2", CriterionStatus.MET, VARIANT, "rare", "supporting", "supports",
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

    def test_togovar_population_provenance_is_exported(self):
        evidence = deepcopy(EVIDENCE)
        evidence[0].update({
            "evidence_id": "https://grch38.togovar.org/variant/tgv123#frequency:tommo",
            "source": "TogoVar",
            "source_version": "API 0.9.1",
            "population": "tommo:global",
        })
        result = CriterionResult(
            "PM2", CriterionStatus.MET, VARIANT, "rare", "supporting", "supports",
            "PM2_supporting", evidence=evidence,
        )
        study = to_evidence_line(result)["hasEvidenceItems"][0]
        self.assertEqual(study["sourceDataSet"]["id"], "https://grch38.togovar.org/")
        self.assertEqual(study["sourceDataSet"]["version"], "API 0.9.1")
        self.assertEqual(study["cohort"]["name"], "tommo:global")
        self.assertEqual(study["specifiedBy"]["name"],
                         "TogoVar API allele frequency aggregation")

    def test_not_met_maps_to_machine_readable_neutral(self):
        result = CriterionResult("PM2", CriterionStatus.NOT_MET, VARIANT, "not rare", None, "none",
                                 "PM2_not_met", evidence=EVIDENCE)
        line = to_evidence_line(result)
        self.assertEqual(line["directionOfEvidenceProvided"], "neutral")
        self.assertNotIn("strengthOfEvidenceProvided", line)
        validate_1_0_1(line, "PM2")

    def test_workflow_statuses_are_not_exported(self):
        for code, status in (("PVS1", CriterionStatus.UNKNOWN), ("PP5", CriterionStatus.UNKNOWN),
                             ("BP6", CriterionStatus.UNKNOWN), ("PM1", CriterionStatus.UNKNOWN)):
            result = CriterionResult(code, status, VARIANT, "workflow state")
            self.assertIsNone(to_evidence_line(result))

    def test_all_supported_criteria_map(self):
        from acmg_pipeline.automated_core.models import CRITERIA
        from acmg_pipeline.criteria.common import DEFAULT_STRENGTH

        for code in CRITERIA:
            if code in {"PVS1", "PP5", "BP6"}:
                continue
            strength = DEFAULT_STRENGTH[code]
            outcome = code
            result = CriterionResult(code, CriterionStatus.MET, VARIANT, "met", strength,
                                     "disputes" if code.startswith("B") else "supports", outcome,
                                     evidence=EVIDENCE)
            with self.subTest(code=code):
                validate_1_0_1(to_evidence_line(result), code)

    def test_all_pvs1_strengths_map(self):
        expected = {
            "very_strong": "PVS1",
            "strong": "PVS1_strong",
            "moderate": "PVS1_moderate",
            "supporting": "PVS1_supporting",
        }
        for strength, outcome in expected.items():
            value = CriterionResult("PVS1", CriterionStatus.MET, VARIANT, "met", strength,
                                    "supports", outcome, evidence=EVIDENCE)
            with self.subTest(strength=strength):
                line = to_evidence_line(value)
                validate_1_0_1(line, "PVS1")
                self.assertEqual(line["evidenceOutcome"]["primaryCoding"]["code"], outcome)

    def test_document_records_pinned_official_schema(self):
        document = export_document([])
        self.assertEqual(document["validated_by"]["schema_version"], "1.0.1")
        self.assertEqual(document["validated_by"]["schema_id"], SCHEMA_ID)
        self.assertEqual(document["validated_by"]["audit_envelope_schema_version"], "1.1")
        self.assertTrue(document["validated_by"]["audit_envelope_schema_sha256"])

    def test_1_0_1_contract_rejects_newer_gks_discriminator_and_mismatched_method(self):
        result = CriterionResult("PM2", CriterionStatus.NOT_MET, VARIANT, "not rare", None, "none",
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
        first = CriterionResult("PM2", CriterionStatus.NOT_MET, VARIANT, "a", None, "none", "PM2_not_met",
                                evidence=[{"evidence_id": "test:x", "value": 1}]).to_dict()
        second = CriterionResult("BA1", CriterionStatus.NOT_MET, VARIANT, "b", None, "none", "BA1_not_met",
                                 evidence=[{"evidence_id": "test:x", "value": 2}]).to_dict()
        with self.assertRaisesRegex(ValueError, "Conflicting evidence"):
            export_document([{"record_id": "test:1", "variant": VARIANT,
                              "results": [first, second]}])

    def test_envelope_resolves_iri_to_structured_study_result(self):
        result = CriterionResult("PM2", CriterionStatus.NOT_MET, VARIANT, "not rare", None, "none",
                                 "PM2_not_met", evidence=EVIDENCE).to_dict()
        document = export_document([{
            "record_id": "test:study-result", "variant": VARIANT, "results": [result],
        }])
        self.assertEqual(document["envelope_schema_version"], "1.1")
        record = document["records"][0]
        wrapped = record["evidence_lines"][0]
        identifier = EVIDENCE[0]["evidence_id"]
        self.assertEqual(identifier, wrapped["evidence_line"]["hasEvidenceItems"][0]["id"])
        item = record["referenced_evidence"][identifier]
        self.assertEqual(item["type"], "StudyResult")
        self.assertEqual(item["sourceDataSet"]["version"], "4.1.1")
        observations = next(extension["value"] for extension in item["extensions"]
                            if extension["name"] == "observations")
        evidence_hash = next(extension["value"] for extension in item["extensions"]
                             if extension["name"] == "normalizedEvidenceSha256")
        self.assertEqual(observations["AC"], 1)
        self.assertEqual(observations["AN"], 100000)
        self.assertEqual(len(evidence_hash), 64)
        self.assertEqual(wrapped["assessment_details"]["status"], "not_met")
        # evidenceItemIds lives only in the audit-index entry now, not in the
        # embedded VA-Spec content (see automated_va_spec.assessment_details()'s
        # docstring, 2026-09-18).
        self.assertEqual(
            record["criterion_assessments"][0]["evidenceItemIds"], [identifier])
        # status is its own top-level extension now, not grouped under one
        # bh26AssessmentDetails object (2026-09-18, per the user's direction).
        line_extensions = {item["name"]: item["value"]
                           for item in wrapped["evidence_line"]["extensions"]}
        self.assertEqual(line_extensions["status"], "not_met")
        self.assertNotIn("evidenceItemIds", line_extensions)

    def test_pvs1_assessment_details_keep_decision_trace(self):
        result = CriterionResult(
            "PVS1", CriterionStatus.MET, VARIANT, "NMD expected", "very_strong", "supports", "PVS1",
            evidence=EVIDENCE,
            evaluation_context={"gene": "TEST", "condition_status": "NOT_PROVIDED"},
            decision_trace=[{"node_id": "NF02", "result": "PASS", "value": True}],
            rules_used=[{"source": "ClinGen PVS1 2018"}],
            warnings=["Condition was not provided."],
        ).to_dict()
        wrapped = export_document([{
            "record_id": "test:pvs1", "variant": VARIANT, "results": [result],
        }])["records"][0]["evidence_lines"][0]
        details = wrapped["assessment_details"]
        self.assertEqual(details["decisionTrace"][0]["node_id"], "NF02")
        self.assertEqual(details["evaluationContext"]["gene"], "TEST")
        self.assertEqual(details["rulesUsed"][0]["source"], "ClinGen PVS1 2018")
        self.assertEqual(wrapped["evidence_line"]["extensions"][0]["value"], details)
        self.assertNotIn("warnings", details)
        hints = wrapped["evidence_line"]["extensions"][1]
        self.assertEqual(hints, {
            "name": "curatorHints",
            "value": [{
                "severity": "warning",
                "category": "warning",
                "message": "Condition was not provided.",
            }],
        })

    def test_automated_review_messages_use_shared_curator_hints_extension(self):
        result = CriterionResult(
            "PM2", CriterionStatus.MET, VARIANT, "rare", "supporting", "supports",
            "PM2_supporting", evidence=EVIDENCE,
            conflict_flags=["Conflicting automated assessments"],
            review_points=["Confirm disease-specific frequency threshold"],
            warnings=["Population coverage is limited"],
        )
        line = to_evidence_line(result)
        extensions = {item["name"]: item["value"] for item in line["extensions"]}
        self.assertNotIn("reviewPoints", extensions["bh26AssessmentDetails"])
        self.assertNotIn("conflictFlags", extensions["bh26AssessmentDetails"])
        self.assertNotIn("warnings", extensions["bh26AssessmentDetails"])
        self.assertEqual(extensions["curatorHints"], [
            {"severity": "warning", "category": "conflict",
             "message": "Conflicting automated assessments"},
            {"severity": "caution", "category": "review",
             "message": "Confirm disease-specific frequency threshold"},
            {"severity": "warning", "category": "warning",
             "message": "Population coverage is limited"},
        ])

    def test_workflow_only_assessment_and_its_evidence_remain_auditable(self):
        evidence = [{
            "evidence_id": "https://example.org/evidence/annotation",
            "category": "annotation", "source": "Ensembl VEP", "source_version": "116",
            "retrieved_at": "2026-09-16T00:00:00Z", "quality_status": "PASS",
            "variant_key": "GRCh38:1:2:C:T", "transcript": "NM_TEST.1",
        }]
        result = CriterionResult(
            "PM1", CriterionStatus.UNKNOWN, VARIANT, "Reviewed region unavailable",
            evidence=evidence, missing_inputs=["region"],
            review_points=["Curate a disease-relevant functional region"],
        ).to_dict()
        record = export_document([{
            "record_id": "test:workflow", "variant": VARIANT, "results": [result],
        }])["records"][0]
        self.assertEqual(record["evidence_lines"], [])
        assessment = record["criterion_assessments"][0]
        self.assertEqual(assessment["status"], "unknown")
        self.assertEqual(assessment["missingInputs"], ["region"])
        self.assertNotIn("reviewPoints", assessment)
        self.assertNotIn("reviewPoints", assessment)
        self.assertIn(evidence[0]["evidence_id"], record["referenced_evidence"])

    def test_audit_envelope_rejects_an_unresolved_evidence_reference(self):
        result = CriterionResult("PM2", CriterionStatus.NOT_MET, VARIANT, "not rare", None, "none",
                                 "PM2_not_met", evidence=EVIDENCE).to_dict()
        document = export_document([{
            "record_id": "test:broken-reference", "variant": VARIANT, "results": [result],
        }])
        del document["records"][0]["referenced_evidence"][EVIDENCE[0]["evidence_id"]]
        with self.assertRaisesRegex(ValueError, "Unresolved Evidence Item references"):
            validate_envelope(document)
