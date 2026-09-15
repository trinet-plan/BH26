import json
import unittest
import uuid
from pathlib import Path

from acmg.cli import main
from acmg.core.models import CRITERIA


ROOT = Path(__file__).resolve().parents[1]


class DemoPipelineTests(unittest.TestCase):
    def test_all_demo_records_evaluate_from_committed_cache(self):
        base = ROOT / ".work" / "demo-pipeline" / uuid.uuid4().hex
        prepared = base / "prepared"
        evaluated = base / "evaluated"
        code = main([
            "prepare-demo-online", "--input-dir", str(ROOT / "demo-data"),
            "--cache-dir", str(ROOT / "tests" / "fixtures" / "ensembl-cache"),
            "--evidence-cache-dir", str(ROOT / "tests" / "fixtures" / "external-cache"),
            "--output-dir", str(prepared), "--ensembl-release", "116",
            "--with-gnomad", "--gnomad-release", "4.1.1",
            "--with-clinvar", "--clinvar-release", "2026-09-15", "--offline",
        ])
        self.assertEqual(code, 0)
        audit = json.loads((prepared / "audit.json").read_text(encoding="utf-8"))
        self.assertEqual(len(audit["records"]), 28)
        self.assertFalse(any(row["identity_status"] == "PENDING" for row in audit["records"]))
        manifest = json.loads((prepared / "identity-manifest.json").read_text(encoding="utf-8"))
        self.assertFalse(manifest["network_used"])
        self.assertEqual(manifest["annotation_evidence"], 17)
        self.assertEqual(manifest["external_providers"][0]["queried_variants"], 17)
        self.assertEqual(manifest["external_providers"][0]["observed_variants"], 15)
        self.assertEqual(manifest["external_providers"][1]["matched_records"], 13)

        code = main([
            "evaluate", "--input", str(prepared / "variants.json"),
            "--evidence", str(prepared / "evidence.json"), "--offline",
            "--config", str(ROOT / "config" / "demo-rules.json"),
            "--output-dir", str(evaluated),
        ])
        self.assertEqual(code, 0)
        results = json.loads((evaluated / "results.json").read_text(encoding="utf-8"))
        self.assertEqual(len(results["records"]), 28)
        self.assertTrue(all(len(record["results"]) == len(CRITERIA)
                            for record in results["records"]))
        self.assertEqual(results["va_spec_export"], "VALIDATED")
        self.assertTrue((evaluated / "evidence-lines.json").is_file())
        pm2 = [result for record in results["records"] for result in record["results"]
               if result["criterion"] == "PM2"]
        self.assertTrue(any(result["status"] == "NOT_MET" for result in pm2))


if __name__ == "__main__":
    unittest.main()
