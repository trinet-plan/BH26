import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.providers.gnomad import GnomadProvider


class Client:
    def __init__(self, body):
        self.body = body

    def fetch(self, url, **kwargs):
        self.url = url
        self.kwargs = kwargs
        return {"body": self.body, "retrieved_at": "2026-09-15T00:00:00Z"}


class GnomadProviderTests(unittest.TestCase):
    def test_frequency_and_population_observations(self):
        body = {"data": {"variant": {
            "variant_id": "1-2-C-T", "reference_genome": "GRCh38",
            "chrom": "1", "pos": 2, "ref": "C", "alt": "T", "flags": [],
            "exome": {"ac": 2, "an": 10000, "af": 0.0002, "filters": [],
                      "populations": [{"id": "eas", "ac": 1, "an": 1000}]},
            "genome": None,
        }}}
        client = Client(body)
        result = GnomadProvider(client).get_frequency(Variant("GRCh38", "1", 2, "C", "T"))
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]["AF"], "0.0002")
        self.assertEqual(result[1]["population"], "exome:eas")
        self.assertEqual(result[1]["AF"], "0.001")
        self.assertTrue(client.kwargs["allow_application_errors"])

    def test_not_found_is_not_zero_frequency(self):
        body = {"errors": [{"message": "Variant not found"}], "data": {"variant": None}}
        result = GnomadProvider(Client(body)).get_frequency(
            Variant("GRCh38", "1", 2, "C", "T"))
        self.assertIsNone(result)

    def test_filtered_variant_is_not_pass(self):
        body = {"data": {"variant": {
            "variant_id": "1-2-C-T", "reference_genome": "GRCh38",
            "chrom": "1", "pos": 2, "ref": "C", "alt": "T", "flags": ["lcr"],
            "exome": {"ac": 1, "an": 100, "af": 0.01, "filters": [],
                      "populations": []}, "genome": None,
        }}}
        result = GnomadProvider(Client(body)).get_frequency(
            Variant("GRCh38", "1", 2, "C", "T"))
        self.assertEqual(result[0]["quality_status"], "FILTERED")

    def test_batch_keeps_hits_when_not_found_errors_are_unscoped(self):
        present = {
            "variant_id": "1-2-C-T", "reference_genome": "GRCh38",
            "chrom": "1", "pos": 2, "ref": "C", "alt": "T", "flags": [],
            "exome": {"ac": 0, "an": 10000, "af": 0, "filters": [],
                      "populations": []}, "genome": None,
        }
        body = {"errors": [{"message": "Variant not found"}],
                "data": {"v0": present, "v1": None}}
        variants = [Variant("GRCh38", "1", 2, "C", "T"),
                    Variant("GRCh38", "1", 3, "G", "A")]
        result = GnomadProvider(Client(body)).get_frequencies(variants)
        self.assertEqual(len(result[variants[0].key]), 1)
        self.assertIsNone(result[variants[1].key])


if __name__ == "__main__":
    unittest.main()
