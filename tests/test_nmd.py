"""NMD prediction from exon numbering - and the boundary case it declines to decide.

ClinGen's rule turns on the last exon and the last 50 nucleotides of the penultimate one.
Exon numbering settles everything upstream of that and nothing inside it, so the tests are
mostly about where the provider stops.
"""

import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.nmd import (
    METHOD, RULE_SOURCE, NmdPredictionProvider, parse_exon,
)


def consequence(**overrides):
    record = {
        "transcript_id": "ENST00000399249",
        "gene_symbol": "MYBPC3",
        "consequence_terms": ["frameshift_variant"],
        "exon": "2/34",
        "mane_select": "NM_000256.3",
    }
    record.update(overrides)
    return record


class FakeClient:
    def __init__(self, consequences):
        self.body = [{"transcript_consequences": consequences}]
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        return {"body": self.body, "body_sha256": "abc123",
                "retrieved_at": "2026-09-17T00:00:00Z"}


class NmdPredictionTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "11", 47351252, "CT", "C")

    def predict(self, consequences, gene="MYBPC3", transcript="NM_000256.3",
                hgvsc="NM_000256.3:c.278delA"):
        provider = NmdPredictionProvider(FakeClient(consequences), "116")
        return provider.get_nmd_prediction(self.variant, gene, transcript, hgvsc)

    def test_a_ptc_upstream_of_the_final_exons_is_predicted_to_undergo_nmd(self):
        records = self.predict([consequence(exon="2/34")])
        self.assertEqual(len(records), 1)
        self.assertIs(records[0]["predicted"], True)
        self.assertEqual(records[0]["rule_source"], RULE_SOURCE)
        self.assertEqual(records[0]["exon"], "2/34")

    def test_the_last_two_exons_are_left_to_a_curator(self):
        """"Within the last 50 nucleotides of the penultimate exon" needs that exon's
        length. Numbering does not carry it, and guessing would guess at exactly the
        distinction the rule draws."""
        for exon in ("34/34", "33/34", "8/8", "7/8"):
            with self.subTest(exon=exon):
                self.assertEqual(self.predict([consequence(exon=exon)]), [])

    def test_the_exon_just_before_the_boundary_is_still_decidable(self):
        records = self.predict([consequence(exon="32/34")])
        self.assertEqual(len(records), 1)

    def test_a_variant_spanning_a_junction_is_placed_by_its_first_exon(self):
        records = self.predict([consequence(exon="8-9/34")])
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["exon"], "8-9/34")

    def test_consequences_for_other_genes_are_ignored(self):
        records = self.predict([
            consequence(gene_symbol="SPI1", exon="1/5"),
            consequence(exon="2/34"),
        ])
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["exon"], "2/34")

    def test_a_non_truncating_consequence_is_not_an_nmd_question(self):
        self.assertEqual(
            self.predict([consequence(consequence_terms=["missense_variant"])]), [])

    def test_the_mane_consequence_wins_when_transcripts_disagree(self):
        """VEP answers about Ensembl transcripts; mane_select is what ties one to the
        RefSeq accession being evaluated."""
        records = self.predict([
            consequence(transcript_id="ENST00000544791", exon="30/31", mane_select=None),
            consequence(exon="2/34", mane_select="NM_000256.3"),
        ])
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["exon"], "2/34")

    def test_disagreeing_candidates_without_a_mane_match_decide_nothing(self):
        self.assertEqual(
            self.predict([
                consequence(exon="2/34", mane_select=None),
                consequence(transcript_id="ENST2", exon="9/27", mane_select=None),
            ]),
            [],
        )

    def test_candidates_agreeing_on_the_exon_are_usable_without_a_mane_match(self):
        records = self.predict([
            consequence(exon="2/34", mane_select=None),
            consequence(transcript_id="ENST2", exon="2/34", mane_select=None),
        ])
        self.assertEqual(len(records), 1)

    def test_no_exon_number_means_no_record(self):
        self.assertEqual(self.predict([consequence(exon=None)]), [])
        self.assertEqual(self.predict([]), [])

    def test_missing_inputs_are_not_an_error(self):
        self.assertEqual(self.predict([consequence()], gene=""), [])
        self.assertEqual(self.predict([consequence()], transcript=""), [])
        self.assertEqual(self.predict([consequence()], hgvsc=""), [])

    def test_it_asks_for_exon_numbers_rather_than_reusing_the_annotation_request(self):
        """The committed offline annotation cache is keyed on a URL without `numbers`."""
        client = FakeClient([consequence()])
        NmdPredictionProvider(client, "116").get_nmd_prediction(
            self.variant, "MYBPC3", "NM_000256.3", "NM_000256.3:c.278delA")
        self.assertIn("numbers=1", client.urls[0])

    def test_the_record_declares_itself_automated(self):
        record = self.predict([consequence()])[0]
        self.assertEqual(record["assessment_method"], "automated")
        self.assertEqual(record["method"], METHOD)
        self.assertNotIn("curator", record)
        self.assertTrue(reviewed_or_automated(record))

    def test_a_release_is_required(self):
        with self.assertRaises(ValueError):
            NmdPredictionProvider(FakeClient([]), "")

    def test_parse_exon_reads_the_shapes_vep_reports(self):
        self.assertEqual(parse_exon("8/34"), (8, 34))
        self.assertEqual(parse_exon("8-9/34"), (8, 34))
        for value in (None, "", "8", "x/y", 8):
            with self.subTest(value=value):
                self.assertIsNone(parse_exon(value))


if __name__ == "__main__":
    unittest.main()
