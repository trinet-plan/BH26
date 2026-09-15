import unittest

from acmg.providers.ensembl import EnsemblIdentityProvider


class FakeClient:
    def __init__(self, body):
        self.body = body

    def fetch(self, url, **kwargs):
        if "/sequence/" in url:
            if "100..101" in url:
                return {"body": "AG", "retrieved_at": "2026-01-01T00:00:00Z"}
            if "100..100" in url:
                return {"body": "A", "retrieved_at": "2026-01-01T00:00:00Z"}
            if "101..101" in url:
                return {"body": "G", "retrieved_at": "2026-01-01T00:00:00Z"}
        return {"body": self.body, "retrieved_at": "2026-01-01T00:00:00Z"}


class EnsemblIdentityTests(unittest.TestCase):
    def test_deletion_is_anchored_and_reference_validated(self):
        body = [{"input": "NM_1.2:c.1delG", "assembly_name": "GRCh38",
                 "seq_region_name": "1", "start": 101, "end": 101,
                 "allele_string": "G/-", "strand": 1}]
        provider = EnsemblIdentityProvider(FakeClient(body), "115")
        record = {"identity": {"TRANSCRIPT": "NM_1.2", "HGVSC": "c.1delG"}}
        candidate = provider.map_record(record)
        self.assertEqual(candidate["variant"],
                         {"assembly": "GRCh38", "chrom": "1", "pos": 100,
                          "ref": "AG", "alt": "A"})
        self.assertEqual(candidate["matched_identifiers"],
                         {"TRANSCRIPT": "NM_1.2", "HGVSC": "c.1delG"})

    def test_wrong_build_is_rejected(self):
        body = [{"input": "NM_1.2:c.1G>A", "assembly_name": "GRCh37",
                 "seq_region_name": "1", "start": 101, "end": 101,
                 "allele_string": "G/A", "strand": 1}]
        provider = EnsemblIdentityProvider(FakeClient(body), "115")
        with self.assertRaisesRegex(ValueError, "requested HGVS/build"):
            provider.map_record({"identity": {"TRANSCRIPT": "NM_1.2", "HGVSC": "c.1G>A"}})

    def test_json_wrapped_sequence_is_supported(self):
        client = FakeClient([])
        original = client.fetch

        def fetch(url, **kwargs):
            if "/sequence/" in url:
                return {"body": '{"seq":"G"}', "retrieved_at": "2026-01-01T00:00:00Z"}
            return original(url, **kwargs)

        client.fetch = fetch
        self.assertEqual(EnsemblIdentityProvider(client, "115").reference.sequence("1", 1, 1), "G")

    def test_reverse_strand_alleles_are_reverse_complemented(self):
        body = [{"input": "NM_1.2:c.1C>T", "assembly_name": "GRCh38",
                 "seq_region_name": "1", "start": 101, "end": 101,
                 "allele_string": "C/T", "strand": -1}]
        provider = EnsemblIdentityProvider(FakeClient(body), "115")
        candidate = provider.map_record(
            {"identity": {"TRANSCRIPT": "NM_1.2", "HGVSC": "c.1C>T"}})
        self.assertEqual(candidate["variant"]["ref"], "G")
        self.assertEqual(candidate["variant"]["alt"], "A")


if __name__ == "__main__":
    unittest.main()
