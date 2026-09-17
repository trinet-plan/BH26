import unittest
from decimal import Decimal

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.reference_links import togovar_variant_url
from acmg_pipeline.providers.togovar import API_URL, API_VERSION, TogoVarProvider
from acmg_pipeline.vcf_record import VariantRecord


class Client:
    def __init__(self, body):
        self.body = body

    def fetch(self, url, **kwargs):
        self.url = url
        self.kwargs = kwargs
        return {
            "body": self.body,
            "retrieved_at": "2026-09-17T00:00:00Z",
        }


class TogoVarProviderTests(unittest.TestCase):
    def test_curator_reference_uses_togovar_coordinate_url(self):
        variant = VariantRecord(
            chrom="12", pos=111803962, id="rs671", ref="G", alt="A",
            qual=".", filter="PASS", info={"ASSEMBLY": "GRCh38"},
        )
        self.assertEqual(
            togovar_variant_url(variant),
            "https://grch38.togovar.org/variant/12-111803962-G-A",
        )

    def test_exact_allele_and_all_frequency_datasets_are_normalized(self):
        body = {"data": [{
            "id": "tgv123",
            "chromosome": "1",
            "position": 2,
            "reference": "C",
            "alternate": "T",
            "frequencies": [
                {"source": "jga_wgs", "ac": 2, "an": 1000, "af": 0.002,
                 "filter": ["PASS"]},
                {"source": "gnomad_exomes", "ac": 1, "an": 10000, "af": 0.0001,
                 "filter": ["PASS"]},
            ],
        }]}
        client = Client(body)
        result = TogoVarProvider(client).get_frequency(
            Variant("GRCh38", "1", 2, "C", "T")
        )

        self.assertEqual(client.url, API_URL)
        self.assertEqual(client.kwargs["dataset_version"], f"TogoVar API {API_VERSION}")
        self.assertEqual(client.kwargs["data"]["query"]["location"], {
            "chromosome": "1", "position": 2,
        })
        self.assertEqual(len(result), 2)
        self.assertEqual({item["population"] for item in result}, {
            "jga_wgs:global", "gnomad_exomes:global",
        })
        self.assertTrue(all(item["source"] == "TogoVar" for item in result))
        self.assertTrue(all(item["quality_status"] == "PASS" for item in result))
        self.assertEqual(result[0]["AF"], "0.002")
        self.assertEqual(result[0]["reported_AF"], "0.002")
        self.assertEqual(result[0]["upstream_dataset_version"],
                         "NOT_PROVIDED_BY_TOGOVAR_API")

    def test_configured_source_groups_keep_only_gnomad_and_tommo(self):
        body = {"data": [{
            "id": "tgv123", "chromosome": "1", "position": 2,
            "reference": "C", "alternate": "T",
            "frequencies": [
                {"source": "jga_wgs", "ac": 2, "an": 1000, "af": 0.002,
                 "filter": ["PASS"]},
                {"source": "tommo", "ac": 3, "an": 1000, "af": 0.003,
                 "filter": ["PASS"]},
                {"source": "gnomad_exomes", "ac": 1, "an": 10000, "af": 0.0001,
                 "filter": ["PASS"]},
                {"source": "gnomad_genomes", "ac": 1, "an": 20000, "af": 0.00005,
                 "filter": ["PASS"]},
            ],
        }]}
        result = TogoVarProvider(
            Client(body), frequency_sources=["gnomad", "tommo"]
        ).get_frequency(Variant("GRCh38", "1", 2, "C", "T"))
        self.assertEqual({item["upstream_dataset"] for item in result}, {
            "gnomad_exomes", "gnomad_genomes", "tommo",
        })

    def test_unknown_or_duplicate_source_groups_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown"):
            TogoVarProvider(Client({}), frequency_sources=["unknown"])
        with self.assertRaisesRegex(ValueError, "nonempty unique"):
            TogoVarProvider(Client({}), frequency_sources=["tommo", "tommo"])

    def test_same_position_with_a_different_allele_is_not_a_match(self):
        body = {"data": [{
            "id": "tgv-other", "chromosome": "1", "position": 2,
            "reference": "C", "alternate": "A", "frequencies": [],
        }]}
        result = TogoVarProvider(Client(body)).get_frequency(
            Variant("GRCh38", "1", 2, "C", "T")
        )
        self.assertIsNone(result)

    def test_missing_variant_or_frequency_is_not_zero_frequency(self):
        variant = Variant("GRCh38", "1", 2, "C", "T")
        self.assertIsNone(TogoVarProvider(Client({"data": []})).get_frequency(variant))
        body = {"data": [{
            "id": "tgv123", "chromosome": "1", "position": 2,
            "reference": "C", "alternate": "T",
        }]}
        self.assertIsNone(TogoVarProvider(Client(body)).get_frequency(variant))

    def test_non_pass_frequency_is_retained_but_marked_filtered(self):
        body = {"data": [{
            "id": "tgv123", "chromosome": "1", "position": 2,
            "reference": "C", "alternate": "T",
            "frequencies": [{
                "source": "tommo", "ac": 1, "an": 100, "af": 0.01,
                "filter": ["LowQual"],
            }],
        }]}
        result = TogoVarProvider(Client(body)).get_frequency(
            Variant("GRCh38", "1", 2, "C", "T")
        )
        self.assertEqual(result[0]["quality_status"], "FILTERED")

    def test_inconsistent_frequency_is_rejected(self):
        body = {"data": [{
            "id": "tgv123", "chromosome": "1", "position": 2,
            "reference": "C", "alternate": "T",
            "frequencies": [{
                "source": "tommo", "ac": 1, "an": 100, "af": 0.5,
                "filter": ["PASS"],
            }],
        }]}
        with self.assertRaisesRegex(ValueError, "Inconsistent"):
            TogoVarProvider(Client(body)).get_frequency(
                Variant("GRCh38", "1", 2, "C", "T")
            )

    def test_rounded_reported_af_is_recomputed_from_exact_counts(self):
        body = {"data": [{
            "id": "tgv123", "chromosome": "1", "position": 2,
            "reference": "C", "alternate": "T",
            "frequencies": [{
                "source": "gem_j_wga", "ac": 3536, "an": 15008,
                "af": 0.23600000143051147, "filter": ["PASS"],
            }],
        }]}
        result = TogoVarProvider(Client(body)).get_frequency(
            Variant("GRCh38", "1", 2, "C", "T")
        )
        self.assertEqual(result[0]["AF"], str(Decimal(3536) / Decimal(15008)))
        self.assertEqual(result[0]["reported_AF"], "0.23600000143051147")


if __name__ == "__main__":
    unittest.main()
