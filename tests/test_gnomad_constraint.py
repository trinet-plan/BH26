"""gnomAD per-gene constraint (mis_z / pLI / LOEUF) parsing and error handling."""

import unittest

from acmg_pipeline.providers.gnomad_constraint import GnomadConstraintProvider


class FakeClient:
    def __init__(self, body):
        self.body = body
        self.calls = []

    def fetch(self, url, *, data=None, **kwargs):
        self.calls.append(data)
        return {"body": self.body, "retrieved_at": "2026-09-17T00:00:00Z"}


def gene_response(mis_z=1.3, p_li=0.9, oe_lof_upper=0.4):
    constraint = {"mis_z": mis_z, "pLI": p_li, "oe_lof_upper": oe_lof_upper}
    return {"data": {"gene": {"gene_id": "ENSG00000000001", "gnomad_constraint": constraint}}}


class GnomadConstraintTests(unittest.TestCase):
    def test_a_found_gene_returns_the_three_constraint_metrics(self):
        provider = GnomadConstraintProvider(FakeClient(gene_response(mis_z=7.4, p_li=1e-20, oe_lof_upper=0.66)))
        record = provider.get_constraint("MYH7")
        self.assertEqual(record["mis_z"], 7.4)
        self.assertEqual(record["p_li"], 1e-20)
        self.assertEqual(record["loeuf"], 0.66)
        for field in ("source", "source_version", "retrieved_at", "evidence_id"):
            self.assertTrue(record[field], field)

    def test_a_gene_not_found_returns_none_not_zero(self):
        provider = GnomadConstraintProvider(FakeClient(
            {"data": {"gene": None}, "errors": [{"message": "Gene 'NOT_A_GENE' not found"}]}
        ))
        self.assertIsNone(provider.get_constraint("NOT_A_GENE"))

    def test_a_real_error_is_not_read_as_absence(self):
        provider = GnomadConstraintProvider(FakeClient(
            {"data": {"gene": None}, "errors": [{"message": "Internal server error"}]}
        ))
        with self.assertRaises(ValueError):
            provider.get_constraint("MYH7")

    def test_a_gene_with_no_constraint_block_returns_none(self):
        provider = GnomadConstraintProvider(FakeClient(
            {"data": {"gene": {"gene_id": "ENSG00000000001", "gnomad_constraint": None}}}
        ))
        self.assertIsNone(provider.get_constraint("MYH7"))

    def test_a_non_numeric_constraint_value_is_rejected(self):
        provider = GnomadConstraintProvider(FakeClient(
            {"data": {"gene": {"gene_id": "x", "gnomad_constraint": {"mis_z": "not-a-number",
                                                                     "pLI": 0.9, "oe_lof_upper": 0.4}}}}
        ))
        with self.assertRaises(ValueError):
            provider.get_constraint("MYH7")

    def test_release_must_be_recorded(self):
        with self.assertRaises(ValueError):
            GnomadConstraintProvider(FakeClient({}), release="")


if __name__ == "__main__":
    unittest.main()
