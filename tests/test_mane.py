"""MANE Select as PVS1's transcript-relevance signal - and the case it refuses to decide.

The provider exists to unblock NF01, so the tests are about the boundary it holds while
doing that: only a match asserts relevance, a non-match asserts nothing, and the record
says it was derived rather than reviewed.
"""

import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.mane import METHOD, ManeTranscriptProvider, accession, parse


COLUMNS = ["#NCBI_GeneID", "Ensembl_Gene", "HGNC_ID", "symbol", "name", "RefSeq_nuc",
           "RefSeq_prot", "Ensembl_nuc", "Ensembl_prot", "MANE_status"]


def summary(rows):
    lines = ["\t".join(COLUMNS)]
    for symbol, nucleotide, status in rows:
        lines.append("\t".join(["1", "ENSG1", "HGNC:1", symbol, "name", nucleotide,
                                "NP_1.1", "ENST1.1", "ENSP1.1", status]))
    return "\n".join(lines) + "\n"


class FakeClient:
    def __init__(self, listing=None, text=None):
        self.listing = listing or "MANE.GRCh38.v1.4.summary.txt.gz MANE.GRCh38.v1.5.summary.txt.gz"
        self.text = text
        self.formats = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.formats.append(response_format)
        body = self.text if url.endswith(".gz") else self.listing
        return {"body": body, "retrieved_at": "2026-09-17T00:00:00Z"}


class ManeTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "11", 47351252, "CT", "C")

    def provider(self, rows, listing=None):
        client = FakeClient(listing=listing, text=summary(rows))
        return ManeTranscriptProvider.from_directory(client), client

    def assess(self, rows, gene, transcript):
        provider, _ = self.provider(rows)
        return provider.get_transcript_assessment(self.variant, gene, transcript)

    def test_the_evaluated_transcript_being_mane_select_asserts_relevance(self):
        records = self.assess([("MYBPC3", "NM_000256.3", "MANE Select")],
                              "MYBPC3", "NM_000256.3")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["relevance"], "RELEVANT")
        self.assertEqual(records[0]["mane_select_accession"], "NM_000256.3")

    def test_a_different_transcript_is_unresolved_rather_than_irrelevant(self):
        """resolve_erepo_transcripts.py found curation legitimately uses non-MANE
        transcripts - TNNT2's MANE pick is not the accession ERepo curates against - so
        NOT_RELEVANT here would call PVS1 inapplicable every time that happens."""
        self.assertEqual(
            self.assess([("TNNT2", "NM_001276345.2", "MANE Select")],
                        "TNNT2", "NM_001001430.3"),
            [],
        )

    def test_a_gene_without_a_mane_select_yields_nothing(self):
        self.assertEqual(self.assess([("OTHER", "NM_9.9", "MANE Select")],
                                     "MYBPC3", "NM_000256.3"), [])
        self.assertEqual(self.assess([("MYBPC3", "NM_000256.3", "MANE Plus Clinical")],
                                     "MYBPC3", "NM_000256.3"), [])

    def test_a_version_difference_does_not_withhold_the_signal(self):
        """The version moves when the sequence record is updated, not when the pick changes."""
        records = self.assess([("MYBPC3", "NM_000256.4", "MANE Select")],
                              "MYBPC3", "NM_000256.3")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["evaluated_accession"], "NM_000256.3")
        self.assertEqual(records[0]["mane_select_accession"], "NM_000256.4")

    def test_missing_gene_or_transcript_is_not_an_error(self):
        rows = [("MYBPC3", "NM_000256.3", "MANE Select")]
        self.assertEqual(self.assess(rows, "", "NM_000256.3"), [])
        self.assertEqual(self.assess(rows, "MYBPC3", ""), [])

    def test_an_exon_numbered_on_this_transcript_answers_nf03(self):
        """Once relevance is asserted of the transcript, an exon VEP can number *on that
        transcript* is in it by construction - no further source decides NF03."""
        provider, _ = self.provider([("MYBPC3", "NM_000256.3", "MANE Select")])
        record = provider.get_transcript_assessment(
            self.variant, "MYBPC3", "NM_000256.3", exon="2/34")[0]
        self.assertEqual(record["exon_relevance"], "RELEVANT")
        self.assertEqual(record["affected_exon"], "2/34")

    def test_without_an_exon_nf03_stays_unresolved_rather_than_denied(self):
        """The converse is never claimed: no exon number does not mean the exon is absent,
        which PVS1 would read as NOT_RELEVANT and call the criterion inapplicable."""
        record = self.assess([("MYBPC3", "NM_000256.3", "MANE Select")],
                             "MYBPC3", "NM_000256.3")[0]
        self.assertEqual(record["relevance"], "RELEVANT")
        self.assertNotIn("exon_relevance", record)
        self.assertNotIn("affected_exon", record)

    def test_the_record_declares_itself_automated(self):
        record = self.assess([("MYBPC3", "NM_000256.3", "MANE Select")],
                             "MYBPC3", "NM_000256.3")[0]
        self.assertEqual(record["assessment_method"], "automated")
        self.assertEqual(record["method"], METHOD)
        self.assertEqual(record["policy_version"], "MANE v1.5")
        self.assertNotIn("curator", record)
        self.assertTrue(reviewed_or_automated(record))

    def test_two_variants_on_one_transcript_get_distinct_ids(self):
        provider, _ = self.provider([("MYBPC3", "NM_000256.3", "MANE Select")])
        other = Variant("GRCh38", "11", 47335041, "C", "T")
        first = provider.get_transcript_assessment(self.variant, "MYBPC3", "NM_000256.3")[0]
        second = provider.get_transcript_assessment(other, "MYBPC3", "NM_000256.3")[0]
        self.assertNotEqual(first["evidence_id"], second["evidence_id"])

    def test_the_release_is_read_from_the_directory_not_assumed(self):
        provider, client = self.provider([("MYBPC3", "NM_000256.3", "MANE Select")])
        self.assertEqual(provider.release, "1.5")  # the newest of 1.4 and 1.5
        self.assertTrue(provider.summary_url.endswith("MANE.GRCh38.v1.5.summary.txt.gz"))

    def test_the_summary_is_fetched_as_compressed_text(self):
        provider, client = self.provider([("MYBPC3", "NM_000256.3", "MANE Select")])
        provider.get_transcript_assessment(self.variant, "MYBPC3", "NM_000256.3")
        self.assertIn("text-gz", client.formats)

    def test_a_directory_with_no_summary_is_an_error(self):
        with self.assertRaises(ValueError):
            ManeTranscriptProvider.from_directory(FakeClient(listing="index of /"))

    def test_a_summary_without_a_header_is_an_error(self):
        with self.assertRaises(ValueError):
            parse("no header here\n")

    def test_accession_drops_the_version(self):
        self.assertEqual(accession("NM_000256.3"), "NM_000256")
        self.assertIsNone(accession(None))


if __name__ == "__main__":
    unittest.main()
