import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.providers.dbnsfp import DbnsfpProvider
from acmg_pipeline.providers.http import FetchError


METADATA = {"build_version": "20260901", "src": {"dbnsfp": {"version": "4.8a"}}}
VARIANT_BODY = {
    "_id": "chr14:g.23425971G>A",
    "dbnsfp": {
        "revel": {"score": 0.814, "rankscore": 0.94008},
        # dbNSFP reports one entry per transcript.
        "alphamissense": {"score": [0.8811, 0.8654], "pred": ["P", "P"], "rankscore": 0.81537},
        "sift": {"score": 0.0, "pred": "D"},
    },
}


class Client:
    def __init__(self, body=VARIANT_BODY, metadata=METADATA, error=None):
        self.body = body
        self.metadata = metadata
        self.error = error
        self.urls = []

    def fetch(self, url, **kwargs):
        self.urls.append(url)
        if url.endswith("/metadata"):
            return {"body": self.metadata, "retrieved_at": "2026-09-15T00:00:00Z"}
        if self.error:
            raise self.error
        return {"body": self.body, "retrieved_at": "2026-09-15T00:00:01Z"}


class DbnsfpProviderTests(unittest.TestCase):
    variant = Variant("GRCh38", "14", 23425971, "G", "A")

    def test_scores_carry_the_dbnsfp_release(self):
        client = Client()
        records = DbnsfpProvider(client).get_predictions(self.variant, "NM_000257.4")
        self.assertEqual([item["predictor"] for item in records], ["REVEL", "AlphaMissense"])
        for item in records:
            self.assertEqual(item["predictor_version"], "dbNSFP-4.8a")
            self.assertTrue(item["calibration_eligible"])
            self.assertEqual(item["mechanism"], "protein")
            self.assertEqual(item["transcript"], "NM_000257.4")
            self.assertEqual(item["quality_status"], "PASS")
        self.assertEqual(records[0]["score"], "0.814")
        # A per-transcript list collapses to the first reported value, never silently to zero.
        self.assertEqual(records[1]["score"], "0.8811")
        self.assertEqual(records[1]["classification"], "P")

    def test_uncalibrated_predictors_are_not_carried_forward(self):
        records = DbnsfpProvider(Client()).get_predictions(self.variant)
        self.assertNotIn("SIFT", [item["predictor"] for item in records])

    def test_missing_variant_is_not_an_error_and_not_a_score(self):
        client = Client(error=FetchError("HTTP_ERROR:404"))
        self.assertEqual(DbnsfpProvider(client).get_predictions(self.variant), [])

    def test_other_transport_failures_are_raised(self):
        client = Client(error=FetchError("NETWORK_UNAVAILABLE:URLError"))
        with self.assertRaises(FetchError):
            DbnsfpProvider(client).get_predictions(self.variant)

    def test_absent_scores_are_skipped(self):
        client = Client(body={"dbnsfp": {"revel": {"rankscore": 0.5}}})
        self.assertEqual(DbnsfpProvider(client).get_predictions(self.variant), [])

    def test_version_must_be_reported(self):
        client = Client(metadata={"src": {"dbnsfp": {}}})
        with self.assertRaisesRegex(ValueError, "dbNSFP version"):
            DbnsfpProvider(client).get_predictions(self.variant)

    def test_metadata_is_fetched_once(self):
        client = Client()
        provider = DbnsfpProvider(client)
        provider.get_predictions(self.variant)
        provider.get_predictions(Variant("GRCh38", "1", 2, "C", "T"))
        self.assertEqual(sum(url.endswith("/metadata") for url in client.urls), 1)

    def test_indels_are_skipped_without_a_request(self):
        client = Client()
        records = DbnsfpProvider(client).get_predictions(Variant("GRCh38", "11", 47351252, "CT", "C"))
        self.assertEqual(records, [])
        self.assertEqual(client.urls, [])

    def test_only_grch38_is_supported(self):
        with self.assertRaises(ValueError):
            DbnsfpProvider(Client(), assembly="hg19")


if __name__ == "__main__":
    unittest.main()
