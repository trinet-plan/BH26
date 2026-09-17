import unittest
from pathlib import Path

from acmg_pipeline.automated_core.identity import reconcile
from acmg_pipeline.automated_core.input import audit_demo
from acmg_pipeline.providers.http import CachedHttpClient
from acmg_pipeline.providers.ensembl import EnsemblIdentityProvider


class FakeClient:
    def __init__(self, body):
        self.body = body
        self.urls = []

    def fetch(self, url, **kwargs):
        self.urls.append(url)
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

    def test_exact_refseq_protein_and_predictions_are_preserved(self):
        body = [{"input": "NM_1.2:c.1G>A", "assembly_name": "GRCh38",
                 "seq_region_name": "1", "start": 101, "end": 101,
                 "allele_string": "G/A", "strand": 1,
                 "transcript_consequences": [{
                     "transcript_id": "NM_1.2", "gene_symbol": "TEST",
                     "consequence_terms": ["missense_variant"],
                     "protein_id": "NP_1.1", "protein_start": 7, "protein_end": 7,
                     "amino_acids": "R/H", "alphamissense": {
                         "am_pathogenicity": 0.9, "am_class": "likely_pathogenic"},
                     "spliceai": {"DS_AG": 0.01, "DS_AL": 0.2, "DS_DG": 0, "DS_DL": 0.03},
                     "conservation": -1.2,
                 }]}]
        client = FakeClient(body)
        provider = EnsemblIdentityProvider(client, "116")
        record = {"identity": {"TRANSCRIPT": "NM_1.2", "HGVSC": "c.1G>A", "GENE": "TEST"}}
        _, annotation, predictions = provider.map_record_with_evidence(record)
        self.assertEqual(
            {key: annotation[key] for key in ("protein_id", "protein_start", "protein_end",
                                               "ref_aa", "alt_aa", "protein_length_change")},
            {"protein_id": "NP_1.1", "protein_start": 7, "protein_end": 7,
             "ref_aa": "R", "alt_aa": "H", "protein_length_change": 0})
        by_predictor = {item["predictor"]: item for item in predictions}
        self.assertEqual(by_predictor["AlphaMissense"]["score"], 0.9)
        self.assertEqual(by_predictor["SpliceAI"]["score"], 0.2)
        self.assertEqual(by_predictor["Ensembl Compara conservation"]["score"], -1.2)
        self.assertIn("refseq=1", client.urls[0])
        self.assertIn("protein=1", client.urls[0])

    def test_inframe_length_change_is_derived_from_changed_peptide(self):
        consequence = {"consequence_terms": ["inframe_deletion"]}
        self.assertEqual(__import__("acmg_pipeline.providers.ensembl", fromlist=["_protein_length_change"])
                         ._protein_length_change(consequence, "ABC", "-"), -3)

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
