"""The protein-loss measurement, and the two gates it deliberately leaves unanswered.

NF07 splits strong from moderate on a fraction, which is arithmetic. NF04 and NF06 read the
same record and ask what the lost residues do, which is not. The provider supplies one and
not the others, so a curator answers one question with the number already in hand.
"""

import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.protein_region import METHOD, ProteinRegionProvider

from tests.test_nmd import consequence


class FakeClient:
    """Answers the VEP consequence request and the Ensembl lookup from one fake."""

    def __init__(self, consequences, length=213, protein_id="ENSP00000256474"):
        self.consequences = consequences
        self.length = length
        self.protein_id = protein_id
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        if "/lookup/id/" in url:
            translation = ({"length": self.length, "id": self.protein_id}
                           if self.length is not None else None)
            body = {"id": "ENST00000399249", "Translation": translation}
        else:
            body = [{"transcript_consequences": self.consequences}]
        return {"body": body, "body_sha256": "abc123",
                "retrieved_at": "2026-09-17T00:00:00Z"}


class ProteinRegionTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "3", 10146618, "G", "T")

    def region(self, protein_start=204, *, consequences=None, length=213,
               gene="MYBPC3", transcript="NM_000256.3", hgvsc="NM_000256.3:c.610G>T"):
        client = FakeClient(consequences or [consequence(exon="3/3")], length=length)
        provider = ProteinRegionProvider(client, "116")
        return provider.get_protein_region(self.variant, gene, transcript, hgvsc,
                                           protein_start), client

    def test_it_measures_what_is_lost_from_the_termination_codon_onwards(self):
        """VHL c.610G>T: p.Glu204Ter in a 213-residue protein, which the panel called
        PVS1_Moderate - 10 of 213 is under the rule set's 10% threshold."""
        records, _ = self.region(204, length=213)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["lost_residues"], 10)
        self.assertEqual(records[0]["total_protein_length"], 213)

    def test_the_codon_itself_counts_as_lost(self):
        records, _ = self.region(213, length=213)
        self.assertEqual(records[0]["lost_residues"], 1)

    def test_it_does_not_answer_what_the_lost_residues_do(self):
        record = self.region()[0][0]
        self.assertNotIn("critical_region_disrupted", record)
        self.assertNotIn("region_biologically_relevant", record)

    def test_a_position_past_the_protein_is_refused(self):
        self.assertEqual(self.region(500, length=213)[0], [])

    def test_a_missing_translation_length_yields_nothing(self):
        self.assertEqual(self.region(204, length=None)[0], [])

    def test_a_non_integer_or_absent_position_yields_nothing(self):
        for position in (None, "204", 0, -1):
            with self.subTest(protein_start=position):
                self.assertEqual(self.region(position)[0], [])

    def test_missing_identifiers_are_not_an_error(self):
        self.assertEqual(self.region(gene="")[0], [])
        self.assertEqual(self.region(transcript="")[0], [])
        self.assertEqual(self.region(hgvsc="")[0], [])

    def test_disagreeing_transcripts_decide_nothing(self):
        records, _ = self.region(consequences=[
            consequence(transcript_id="ENST1", mane_select=None),
            consequence(transcript_id="ENST2", mane_select=None),
        ])
        self.assertEqual(records, [])

    def test_it_reuses_the_consequence_request_the_nmd_rule_makes(self):
        """Same URL, so the client serves it from cache rather than asking VEP twice."""
        _records, client = self.region()
        self.assertIn("numbers=1", client.urls[0])
        self.assertIn("/lookup/id/ENST00000399249", client.urls[1])

    def test_the_record_declares_itself_automated(self):
        record = self.region()[0][0]
        self.assertEqual(record["assessment_method"], "automated")
        self.assertEqual(record["method"], METHOD)
        self.assertNotIn("curator", record)
        self.assertTrue(reviewed_or_automated(record))

    def test_two_variants_on_one_transcript_get_distinct_ids(self):
        first = self.region(204)[0][0]
        second = self.region(210)[0][0]
        self.assertNotEqual(first["evidence_id"], second["evidence_id"])

    def test_a_release_is_required(self):
        with self.assertRaises(ValueError):
            ProteinRegionProvider(FakeClient([]), "")


if __name__ == "__main__":
    unittest.main()
