import unittest
from pathlib import Path

from acmg.core.vcf_adapter import inputs_from_vcf, parse_vcf
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import ParsedVcf, VariantRecord


ROOT = Path(__file__).resolve().parents[1]


class VcfAdapterTests(unittest.TestCase):
    def test_current_demo_vcfs_convert_to_shared_inputs(self):
        converted = []
        for path in sorted((ROOT / "demo-data").glob("*.vcf")):
            converted.extend(inputs_from_vcf(path))

        self.assertEqual(len(converted), 28)
        self.assertTrue(all(isinstance(variant, VariantRecord) for variant, _ in converted))
        self.assertTrue(all(isinstance(note, ClinicalNoteExtraction) for _, note in converted))
        self.assertTrue(all(variant.info["ASSEMBLY"] == "GRCh38" for variant, _ in converted))

        case1_var2, note = next(
            pair for pair in converted if pair[0].id == "case1-var2"
        )
        self.assertEqual(case1_var2.info["AM_PATHOGENICITY"], 0.9485)
        self.assertEqual(note.proband.genotype.zygosity, "heterozygous")
        self.assertTrue(note.proband.genotype.variant_status)

    def test_explicit_clinical_note_is_preserved(self):
        clinical_note = ClinicalNoteExtraction.from_json({
            "proband": {"phenotype": {"affected_status": True}},
            "family": {"inheritance_pattern": "autosomal_dominant"},
        })
        converted = inputs_from_vcf(
            ROOT / "demo-data" / "case1_variants_v2.vcf",
            clinical_note=clinical_note,
        )
        self.assertTrue(all(note is clinical_note for _, note in converted))

    def test_multiallelic_number_a_info_is_selected_per_alt(self):
        parsed = parse_vcf(ROOT / "tests" / "fixtures" / "multi-alt.vcf")

        self.assertEqual(len(parsed), 2)
        self.assertTrue(all(isinstance(item, ParsedVcf) for item in parsed))
        self.assertEqual([item.record.alt for item in parsed], ["C", "G"])
        self.assertEqual([item.record.info["AF"] for item in parsed], [0.1, 0.2])

    def test_reference_conflict_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "conflicts"):
            parse_vcf(
                ROOT / "demo-data" / "case1_variants_v2.vcf",
                assembly="GRCh37",
            )


if __name__ == "__main__":
    unittest.main()
