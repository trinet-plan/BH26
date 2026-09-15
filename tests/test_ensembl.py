import unittest
from pathlib import Path

from acmg.core.identity import reconcile
from acmg.core.input import audit_demo
from acmg.providers.http import CachedHttpClient
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
        candidate, annotation = provider.map_record_with_annotation(record)
        self.assertEqual(candidate["variant"],
                         {"assembly": "GRCh38", "chrom": "1", "pos": 100,
                          "ref": "AG", "alt": "A"})
        self.assertEqual(candidate["matched_identifiers"],
                         {"TRANSCRIPT": "NM_1.2", "HGVSC": "c.1delG"})
        self.assertEqual(annotation["transcript"], "NM_1.2")
        self.assertEqual(annotation["variant_key"], "GRCh38:1:100:AG:A")

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

    def test_demo_cache_resolves_all_records_offline(self):
        root = Path(__file__).resolve().parents[1]
        records = audit_demo(root / "demo-data")
        client = CachedHttpClient(root / "tests" / "fixtures" / "ensembl-cache", offline=True)
        provider = EnsemblIdentityProvider(client, "116")
        resolved = [reconcile(record, [provider.map_record(record)], provider.reference)
                    for record in records]
        counts = {status: sum(row["identity_status"] == status for row in resolved)
                  for status in ("VERIFIED", "CORRECTED", "PENDING")}
        self.assertEqual(counts, {"VERIFIED": 24, "CORRECTED": 4, "PENDING": 0})
        corrections = {row["source"]["variant_id"]: row["resolution"]["variant"]
                       for row in resolved if row["identity_status"] == "CORRECTED"}
        self.assertEqual(corrections["case1-var1"],
                         {"assembly": "GRCh38", "chrom": "11", "pos": 47351252,
                          "ref": "CT", "alt": "C"})
        self.assertEqual(corrections["case2-var2"],
                         {"assembly": "GRCh38", "chrom": "11", "pos": 47347665,
                          "ref": "AC", "alt": "A"})


if __name__ == "__main__":
    unittest.main()
