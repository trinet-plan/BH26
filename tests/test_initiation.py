"""The downstream start codon, and the two gates on the same record it leaves alone.

IC02 is a property of the coding sequence and nothing else. IC01 and IC03 are not, and PVS1
reads IC01 first, so the path still stops for a curator - with this question answered.
"""

import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.initiation import (
    METHOD, InitiationProvider, downstream_in_frame_start,
)


def consequence(**overrides):
    record = {
        "transcript_id": "ENST00000382848",
        "gene_symbol": "GJB2",
        "consequence_terms": ["start_lost"],
        "mane_select": "NM_004004.6",
    }
    record.update(overrides)
    return record


class FakeClient:
    def __init__(self, consequences, cds):
        self.consequences = consequences
        self.cds = cds
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        body = ({"seq": self.cds} if "/sequence/id/" in url
                else [{"transcript_consequences": self.consequences}])
        return {"body": body, "body_sha256": "abc123",
                "retrieved_at": "2026-09-17T00:00:00Z"}


# ATG, two codons, then an in-frame ATG at codon 4.
CDS_WITH_RESTART = "ATG" + "AAA" + "CCC" + "ATG" + "TAA"
CDS_WITHOUT_RESTART = "ATG" + "AAA" + "CCC" + "TAA"


class DownstreamStartTests(unittest.TestCase):
    def test_it_finds_the_first_in_frame_start_after_the_initiation_codon(self):
        self.assertEqual(downstream_in_frame_start(CDS_WITH_RESTART), (4, "found"))

    def test_an_atg_out_of_frame_is_not_a_restart(self):
        """One base in, so the ATG never begins a codon."""
        self.assertEqual(downstream_in_frame_start("ATG" + "AAT" + "GAA" + "TAA"),
                         (None, "absent"))

    def test_the_initiation_codon_itself_is_not_the_answer(self):
        self.assertEqual(downstream_in_frame_start(CDS_WITHOUT_RESTART), (None, "absent"))

    def test_absence_is_an_answer_but_an_unusable_sequence_is_not(self):
        self.assertEqual(downstream_in_frame_start(CDS_WITHOUT_RESTART)[1], "absent")
        for value in (None, "", "ATGA", "ATG", 123):
            with self.subTest(value=value):
                self.assertEqual(downstream_in_frame_start(value), (None, "unusable"))

    def test_lower_case_sequence_is_read(self):
        self.assertEqual(downstream_in_frame_start(CDS_WITH_RESTART.lower()), (4, "found"))


class InitiationProviderTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "13", 20189546, "A", "G")

    def assess(self, cds=CDS_WITH_RESTART, consequences=None, gene="GJB2",
               transcript="NM_004004.6", hgvsc="NM_004004.6:c.2T>C"):
        client = FakeClient(consequences or [consequence()], cds)
        provider = InitiationProvider(client, "116")
        return provider.get_initiation_assessment(self.variant, gene, transcript, hgvsc), client

    def test_a_restart_codon_answers_ic02(self):
        records, _ = self.assess()
        self.assertEqual(len(records), 1)
        self.assertIs(records[0]["downstream_in_frame_start"], True)
        self.assertEqual(records[0]["downstream_start_codon"], 4)

    def test_no_restart_codon_is_also_an_answer(self):
        records, _ = self.assess(cds=CDS_WITHOUT_RESTART)
        self.assertEqual(len(records), 1)
        self.assertIs(records[0]["downstream_in_frame_start"], False)
        self.assertIsNone(records[0]["downstream_start_codon"])

    def test_an_unusable_sequence_yields_nothing(self):
        self.assertEqual(self.assess(cds="ATGA")[0], [])
        self.assertEqual(self.assess(cds=None)[0], [])

    def test_it_does_not_answer_the_other_two_gates(self):
        record = self.assess()[0][0]
        self.assertNotIn("intact_alternative_transcript", record)
        self.assertNotIn("upstream_pathogenic_evidence", record)

    def test_a_consequence_that_is_not_start_loss_is_not_this_question(self):
        self.assertEqual(
            self.assess(consequences=[consequence(consequence_terms=["stop_gained"])])[0], [])

    def test_consequences_for_other_genes_are_ignored(self):
        records, _ = self.assess(consequences=[
            consequence(gene_symbol="OTHER", transcript_id="ENST9"),
            consequence(),
        ])
        self.assertEqual(len(records), 1)

    def test_the_mane_consequence_wins_when_transcripts_disagree(self):
        records, _ = self.assess(consequences=[
            consequence(transcript_id="ENST9", mane_select=None),
            consequence(mane_select="NM_004004.6"),
        ])
        self.assertEqual(records[0]["ensembl_transcript"], "ENST00000382848")

    def test_disagreeing_candidates_without_a_mane_match_decide_nothing(self):
        self.assertEqual(self.assess(consequences=[
            consequence(transcript_id="ENST1", mane_select=None),
            consequence(transcript_id="ENST2", mane_select=None),
        ])[0], [])

    def test_missing_identifiers_are_not_an_error(self):
        self.assertEqual(self.assess(gene="")[0], [])
        self.assertEqual(self.assess(transcript="")[0], [])
        self.assertEqual(self.assess(hgvsc="")[0], [])

    def test_it_reads_the_coding_sequence_of_the_matched_transcript(self):
        _records, client = self.assess()
        self.assertIn("numbers=1", client.urls[0])
        self.assertIn("/sequence/id/ENST00000382848?type=cds", client.urls[1])

    def test_the_record_declares_itself_automated(self):
        record = self.assess()[0][0]
        self.assertEqual(record["assessment_method"], "automated")
        self.assertEqual(record["method"], METHOD)
        self.assertNotIn("curator", record)
        self.assertTrue(reviewed_or_automated(record))

    def test_a_release_is_required(self):
        with self.assertRaises(ValueError):
            InitiationProvider(FakeClient([], ""), "")


if __name__ == "__main__":
    unittest.main()
