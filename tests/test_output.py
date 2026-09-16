from pathlib import Path
import unittest
from unittest.mock import patch
import uuid

from acmg.core.models import CRITERIA
from acmg.output import load_object, run_internal


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


class OutputTests(unittest.TestCase):
    def output_dir(self):
        return ROOT / ".work" / "test-output" / uuid.uuid4().hex

    def run_fixture(self, path):
        return run_internal(FIXTURES / "synthetic-prepared.json", FIXTURES / "synthetic-evidence.json",
                            FIXTURES / "synthetic-rules.json", path)

    def test_offline_replay_is_identical(self):
        one, two = self.output_dir(), self.output_dir()
        first, second = self.run_fixture(one), self.run_fixture(two)
        self.assertEqual(first, second)
        self.assertEqual((one / "results.json").read_bytes(), (two / "results.json").read_bytes())
        first_manifest = load_object(one / "run-manifest.json")
        second_manifest = load_object(two / "run-manifest.json")
        self.assertEqual(first_manifest["result_sha256"], second_manifest["result_sha256"])
        self.assertFalse(first_manifest["network_used"])
        self.assertEqual(first_manifest["va_spec_export"], "NOT_PERFORMED")

    def test_all_criteria_and_summary_rows(self):
        output = self.output_dir()
        payload = self.run_fixture(output)
        self.assertEqual(payload["schema_version"], "1.1")
        results = payload["records"][0]["results"]
        self.assertEqual([item["criterion"] for item in results], list(CRITERIA))
        self.assertEqual(len((output / "summary.tsv").read_text().splitlines()), 17)
        pm2 = next(item for item in results if item["criterion"] == "PM2")
        self.assertEqual(pm2["evidence_outcome"], "PM2_supporting")
        for code in ("PP5", "BP6"):
            self.assertEqual(next(r for r in results if r["criterion"] == code)["status"], "DEPRECATED")

    def test_old_run_cannot_be_overwritten(self):
        output = self.output_dir()
        self.run_fixture(output)
        with self.assertRaises(FileExistsError):
            self.run_fixture(output)

    def test_va_spec_export_is_validated_and_replayable(self):
        one, two = self.output_dir(), self.output_dir()
        run_internal(FIXTURES / "synthetic-prepared.json", FIXTURES / "synthetic-evidence.json",
                     FIXTURES / "synthetic-rules.json", one, va_spec=True)
        run_internal(FIXTURES / "synthetic-prepared.json", FIXTURES / "synthetic-evidence.json",
                     FIXTURES / "synthetic-rules.json", two, va_spec=True)
        self.assertEqual((one / "evidence-lines.json").read_bytes(),
                         (two / "evidence-lines.json").read_bytes())
        manifest = load_object(one / "run-manifest.json")
        self.assertEqual(manifest["va_spec_export"], "VALIDATED")
        self.assertEqual(manifest["va_spec"]["schema_version"], "1.0.1")
        self.assertGreater(manifest["va_spec"]["instance_count"], 0)
        self.assertEqual(len(list((one / "va-spec-1.0.1").glob("*.json"))),
                         manifest["va_spec"]["instance_count"])

    def test_population_met_fixtures_emit_official_example_shaped_evidence_lines(self):
        output = self.output_dir()
        payload = run_internal(
            FIXTURES / "population-met-prepared.json",
            FIXTURES / "population-met-evidence.json",
            FIXTURES / "population-met-rules.json",
            output,
            criteria=("PM2", "BA1"),
            va_spec=True,
        )
        by_record = {
            record["record_id"]: {result["criterion"]: result for result in record["results"]}
            for record in payload["records"]
        }
        self.assertEqual(by_record["fixture:pm2-met"]["PM2"]["status"], "MET")
        self.assertEqual(by_record["fixture:ba1-met"]["BA1"]["status"], "MET")

        pm2 = load_object(output / "va-spec-1.0.1" / "fixture_pm2-met--PM2.json")
        self.assertEqual(pm2["directionOfEvidenceProvided"], "supports")
        self.assertEqual(pm2["strengthOfEvidenceProvided"]["primaryCoding"]["code"],
                         "supporting")
        self.assertEqual(pm2["evidenceOutcome"]["primaryCoding"]["code"], "PM2_supporting")
        self.assertEqual(pm2["hasEvidenceItems"][0]["type"],
                         "CohortAlleleFrequencyStudyResult")
        self.assertEqual(pm2["hasEvidenceItems"][0]["focusAlleleFrequency"], 0)

        ba1 = load_object(output / "va-spec-1.0.1" / "fixture_ba1-met--BA1.json")
        self.assertEqual(ba1["directionOfEvidenceProvided"], "disputes")
        self.assertEqual(ba1["strengthOfEvidenceProvided"]["primaryCoding"]["code"],
                         "standalone")
        self.assertEqual(ba1["evidenceOutcome"]["primaryCoding"]["code"], "BA1")
        study = ba1["hasEvidenceItems"][0]
        self.assertEqual(study["type"], "CohortAlleleFrequencyStudyResult")
        self.assertEqual(study["focusAlleleFrequency"], 0.0501)
        self.assertEqual(study["focusAlleleCount"], 501)
        self.assertEqual(study["locusAlleleCount"], 10000)
        self.assertIn("https://example.org/evidence/ba1-exception-check-negative",
                      ba1["hasEvidenceItems"])

    def test_va_spec_failure_leaves_no_partial_run(self):
        output = self.output_dir()
        with patch("acmg.va_spec.mapper.export_document", side_effect=ValueError("invalid")):
            with self.assertRaisesRegex(ValueError, "invalid"):
                run_internal(FIXTURES / "synthetic-prepared.json",
                             FIXTURES / "synthetic-evidence.json",
                             FIXTURES / "synthetic-rules.json", output, va_spec=True)
        self.assertFalse(output.exists())
