"""Regression tests for the API-side VCF parser (acmg_pipeline.vcf_record).

The CLI's own parser (acmg_pipeline.automated_core.vcf_adapter._convert_scalar)
already maps VCF's "." missing value to None. This module's parser - the one
parse_request_vcf() calls for an API request - did not, so a declared
Number=1,Type=Float field written as "." raised ValueError out of float(".")
and the variant could not be parsed at all. Three real demo rows carry one
(SG10K_AF=. in case3-var1, case4-var2, case4-var3).
"""
import unittest
from pathlib import Path

from acmg_pipeline.vcf_record import parse_vcf

ROOT = Path(__file__).resolve().parents[1]

HEADER = "\n".join([
    "##fileformat=VCFv4.2",
    '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">',
    '##INFO=<ID=SG10K_AF,Number=1,Type=Float,Description="SG10K allele frequency">',
    '##INFO=<ID=DEPTHS,Number=.,Type=Integer,Description="Per-sample depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
])


def _parse(info: str):
    line = f"14\t23425971\tv1\tG\tA\t99\tPASS\t{info}"
    return parse_vcf(HEADER + "\n" + line + "\n").record


class MissingInfoValueTests(unittest.TestCase):
    def test_missing_float_becomes_none_instead_of_raising(self):
        record = _parse("GENE=MYH7;SG10K_AF=.")
        self.assertIsNone(record.info["SG10K_AF"])

    def test_present_float_is_still_cast(self):
        record = _parse("GENE=MYH7;SG10K_AF=0.0027")
        self.assertEqual(record.info["SG10K_AF"], 0.0027)

    def test_missing_value_in_a_list_field_becomes_none(self):
        record = _parse("GENE=MYH7;DEPTHS=30,.,12")
        self.assertEqual(record.info["DEPTHS"], [30, None, 12])

    def test_missing_value_of_an_undeclared_key_becomes_none(self):
        record = _parse("GENE=MYH7;GNOMAD_AF=.")
        self.assertIsNone(record.info["GNOMAD_AF"])


class DemoVcfTests(unittest.TestCase):
    def test_every_demo_variant_parses(self):
        """case3-var1 (the MYH7 c.2155C>T this was found on) and case4-var2/var3
        carry SG10K_AF=. and used to raise here."""
        for path in sorted((ROOT / "demo-data").glob("case*_variants_v2.vcf")):
            text = path.read_text(encoding="utf-8")
            header = [l for l in text.splitlines() if l.startswith("#")]
            for body in [l for l in text.splitlines() if l and not l.startswith("#")]:
                with self.subTest(vcf=path.name, variant=body.split("\t")[2]):
                    record = parse_vcf("\n".join([*header, body]) + "\n").record
                    self.assertTrue(record.info.get("GENE"))


if __name__ == "__main__":
    unittest.main()
