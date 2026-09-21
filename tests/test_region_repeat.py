"""PM4/BP3's `region` evidence: a UniProt Repeat-feature first pass over the altered
protein interval - see region_repeat.py's own module docstring for why an unanswered
Repeat feature is read as repetitive=False rather than left unset (unlike
protein_region.py's critical_region_disrupted, which this shares its UniProt fetch
and overlap helpers with)."""

import unittest
from unittest.mock import patch

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.region_repeat import METHOD, RegionRepeatProvider


def uniprot_feature(feature_type, start, end, description=""):
    return {"type": feature_type, "description": description,
            "location": {"start": {"value": start}, "end": {"value": end}}}


class FakeClient:
    def __init__(self, uniprot_features=()):
        self.uniprot_features = list(uniprot_features)
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        return {"body": {"features": self.uniprot_features}, "body_sha256": "abc123",
                "retrieved_at": "2026-09-21T00:00:00Z"}


class RegionRepeatTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "3", 10146618, "G", "T")

    def region(self, start=100, end=105, *, gene="TEST", transcript="NM_000256.3",
               protein_id="NP_TEST.1", accession="P00000", uniprot_features=()):
        client = FakeClient(uniprot_features=uniprot_features)
        provider = RegionRepeatProvider(client, "116")
        with patch("acmg_pipeline.providers.region_repeat.gene_to_uniprot_accession",
                  return_value=accession):
            records = provider.get_region(self.variant, gene, transcript, protein_id, start, end)
        return records, client

    def test_no_repeat_feature_over_the_interval_is_repetitive_false_not_unset(self):
        """Unlike protein_region.py's critical_region_disrupted, absence IS the signal
        here - see the module docstring for why."""
        record = self.region(uniprot_features=[
            uniprot_feature("Domain", 1, 50, "Unrelated domain"),
        ])[0][0]
        self.assertFalse(record["repetitive"])
        self.assertFalse(record["nonfunctional_repeat"])
        self.assertTrue(record["functional_review_complete"])

    def test_a_compositional_bias_feature_also_sets_repetitive(self):
        """UniProt types a low-complexity/homopolymeric run (poly-Pro, poly-Gln...) as
        "Compositional bias", not "Repeat" - found empirically from FOXG1's real UniProt
        entry (P55316: "Compositional bias" 59-85 "Pro residues") after the 679-variant
        ground-truth run showed BP3 missing exactly this kind of real case."""
        record = self.region(uniprot_features=[
            uniprot_feature("Compositional bias", 95, 110, "Pro residues"),
        ])[0][0]
        self.assertTrue(record["repetitive"])
        self.assertFalse(record["functional_importance"])
        self.assertTrue(record["nonfunctional_repeat"])

    def test_a_repeat_feature_over_the_interval_sets_repetitive(self):
        record = self.region(uniprot_features=[
            uniprot_feature("Repeat", 95, 110, "ANK 3"),
        ])[0][0]
        self.assertTrue(record["repetitive"])
        self.assertFalse(record["functional_importance"])
        self.assertTrue(record["nonfunctional_repeat"])

    def test_a_repeat_with_an_overlapping_domain_is_not_nonfunctional(self):
        record = self.region(uniprot_features=[
            uniprot_feature("Repeat", 95, 110, "ANK 3"),
            uniprot_feature("Domain", 90, 150, "Ankyrin repeat domain"),
        ])[0][0]
        self.assertTrue(record["repetitive"])
        self.assertTrue(record["functional_importance"])
        self.assertFalse(record["nonfunctional_repeat"])

    def test_a_disordered_region_feature_does_not_count_as_functional(self):
        record = self.region(uniprot_features=[
            uniprot_feature("Repeat", 95, 110, "ANK 3"),
            uniprot_feature("Region", 90, 150, "Disordered"),
        ])[0][0]
        self.assertTrue(record["repetitive"])
        self.assertFalse(record["functional_importance"])

    def test_a_feature_with_an_unresolved_fuzzy_boundary_is_not_a_crash(self):
        """Real UniProt data: an unresolved/fuzzy terminus is recorded with value=null,
        which parses without KeyError/TypeError and must still be rejected."""
        record = self.region(uniprot_features=[
            {"type": "Repeat", "description": "fuzzy",
             "location": {"start": {"value": None}, "end": {"value": 110}}},
        ])[0][0]
        self.assertFalse(record["repetitive"])

    def test_a_repeat_outside_the_interval_does_not_count(self):
        record = self.region(start=100, end=105, uniprot_features=[
            uniprot_feature("Repeat", 500, 520, "Unrelated repeat"),
        ])[0][0]
        self.assertFalse(record["repetitive"])

    def test_no_uniprot_accession_yields_nothing(self):
        self.assertEqual(self.region(accession=None)[0], [])

    def test_no_feature_list_yields_nothing(self):
        client = FakeClient()
        client.fetch = lambda url, **kw: {"body": {}, "body_sha256": "x", "retrieved_at": "t"}
        provider = RegionRepeatProvider(client, "116")
        with patch("acmg_pipeline.providers.region_repeat.gene_to_uniprot_accession",
                  return_value="P00000"):
            records = provider.get_region(self.variant, "TEST", "NM_1", "NP_1", 100, 105)
        self.assertEqual(records, [])

    def test_missing_identifiers_are_not_an_error(self):
        self.assertEqual(self.region(gene="")[0], [])
        self.assertEqual(self.region(transcript="")[0], [])

    def test_a_non_positive_or_absent_start_yields_nothing(self):
        for start in (None, "100", 0, -1):
            with self.subTest(start=start):
                self.assertEqual(self.region(start=start)[0], [])

    def test_an_end_before_start_falls_back_to_start(self):
        record = self.region(start=100, end=50, uniprot_features=[
            uniprot_feature("Repeat", 100, 100, "single residue"),
        ])[0][0]
        self.assertEqual(record["end"], 100)
        self.assertTrue(record["repetitive"])

    def test_the_record_declares_itself_automated(self):
        record = self.region()[0][0]
        self.assertEqual(record["assessment_method"], "automated")
        self.assertEqual(record["method"], METHOD)
        self.assertNotIn("curator", record)
        self.assertTrue(reviewed_or_automated(record))

    def test_the_record_carries_the_protein_id_it_was_asked_about(self):
        """evaluate_region() rejects a region record whose protein_id doesn't match the
        annotation's - this must be the one passed in, not one read from UniProt."""
        record = self.region(protein_id="NP_ANNOTATED.2")[0][0]
        self.assertEqual(record["protein_id"], "NP_ANNOTATED.2")

    def test_two_variants_on_one_transcript_get_distinct_ids(self):
        first = self.region(start=100, end=100)[0][0]
        second = self.region(start=200, end=200)[0][0]
        self.assertNotEqual(first["evidence_id"], second["evidence_id"])

    def test_a_policy_version_is_required(self):
        with self.assertRaises(ValueError):
            RegionRepeatProvider(FakeClient(), "")


if __name__ == "__main__":
    unittest.main()
