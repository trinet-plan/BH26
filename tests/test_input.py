from acmg_pipeline.constants import CriterionStatus
import unittest
from pathlib import Path

from acmg_pipeline.automated_core.input import audit_demo, read_xlsx_rows
from acmg_pipeline.automated_core.models import CRITERIA, CriterionResult, Variant


ROOT = Path(__file__).resolve().parents[1]


class DemoAuditTests(unittest.TestCase):
    def setUp(self):
        self.records = audit_demo(ROOT / "demo-data")

    def test_every_row_is_retained(self):
        expected = sum(
            len(line.split("\t")[4].split(","))
            for path in (ROOT / "demo-data").glob("*.vcf")
            for line in path.read_text().splitlines() if line and not line.startswith("#")
        )
        self.assertEqual(len(self.records), expected)
        self.assertEqual(len({r["record_id"] for r in self.records}), expected)
        self.assertEqual({r["source"]["case_id"] for r in self.records},
                         {"case1", "case2", "case3", "case4"})

    def test_missing_alt_is_never_assumed_deletion(self):
        for identifier in ("case1-var1", "case2-var2"):
            row = next(r for r in self.records if r["source"]["variant_id"] == identifier)
            self.assertIsNone(row["parsed_variant"])
            self.assertEqual(row["identity_status"], "PENDING")

    def test_labels_do_not_enter_identity(self):
        for row in self.records:
            self.assertNotIn("CLNSIG", row["identity"])
            self.assertNotIn("ACMG_CODES", row["identity"])
            self.assertNotIn("NOTE", row["identity"])

    def test_hgvs_conflicting_coordinates_flagged(self):
        rows = [r for r in self.records if r["identity"].get("HGVSC") == "c.861T>C"]
        self.assertGreaterEqual(len(rows), 2)
        for row in rows:
            self.assertIn("CONFLICTING_GENOMIC_REPRESENTATIONS_FOR_HGVS", row["issues"])

    def test_no_unverified_row_claims_validity(self):
        self.assertTrue(all(r["identity_status"] == "PENDING" for r in self.records))

    def test_spreadsheet_rows_are_joined_losslessly(self):
        path = ROOT / "demo-data" / "annotation_alphamissense_alphagenome_v1.xlsx"
        self.assertEqual(len(read_xlsx_rows(path)), 28)
        self.assertTrue(all("spreadsheet_source" in row for row in self.records))
        self.assertFalse(any("VCF_SPREADSHEET_IDENTITY_MISMATCH" in row["issues"]
                             for row in self.records))
        case1_var2 = next(row for row in self.records
                          if row["source"]["variant_id"] == "case1-var2")
        self.assertEqual(case1_var2["spreadsheet_annotations"]["am_pathogenicity"], 0.9485)
        self.assertNotIn("clnsig", case1_var2["spreadsheet_annotations"])


class DomainTests(unittest.TestCase):
    def test_variant_validation(self):
        self.assertEqual(Variant("GRCh38", "chr1", 1, "A", "G").chrom, "1")
        for alt in (".", "<DEL>", "*", "N", "A"):
            with self.assertRaises(ValueError):
                Variant("GRCh38", "1", 1, "A", alt)

    def test_deprecated_cannot_fire_and_pvs1_requires_strength(self):
        self.assertEqual(len(CRITERIA), 16)
        for code in ("PP5", "BP6"):
            with self.assertRaises(ValueError):
                CriterionResult(code, CriterionStatus.MET, {}, "invalid")
        with self.assertRaises(ValueError):
            CriterionResult("PVS1", CriterionStatus.MET, {}, "invalid")
        value = CriterionResult("PVS1", CriterionStatus.MET, {}, "valid", "very_strong")
        self.assertEqual(value.strength, "very_strong")


if __name__ == "__main__":
    unittest.main()
