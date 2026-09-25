"""The protein-loss measurement, and the NF04/NF06 UniProt-derived first pass.

NF07 splits strong from moderate on a fraction, which is arithmetic. NF04 and NF06 read the
same record and ask what the lost residues do - since 2026-09-17 this provider answers that
too, from UniProt feature overlap, always flagged in pvs1.py's own review_points as a
prediction rather than a curator's own review (see protein_region.py's module docstring for
why this default direction was picked, checked against the real VHL c.610G>T ground truth).
"""

import unittest
from unittest.mock import patch

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.protein_region import METHOD, UNIPROT_METHOD, ProteinRegionProvider

from tests.test_nmd import consequence


def uniprot_feature(feature_type, start, end, description=""):
    return {"type": feature_type, "description": description,
            "location": {"start": {"value": start}, "end": {"value": end}}}


class FakeClient:
    """Answers the VEP consequence request, the Ensembl lookup, and UniProt from one fake."""

    def __init__(self, consequences, length=213, protein_id="ENSP00000256474", uniprot_features=()):
        self.consequences = consequences
        self.length = length
        self.protein_id = protein_id
        self.uniprot_features = list(uniprot_features)
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        if "rest.uniprot.org" in url:
            body = {"features": self.uniprot_features}
        elif "/lookup/id/" in url:
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
               gene="MYBPC3", transcript="NM_000256.3", hgvsc="NM_000256.3:c.610G>T",
               accession=None, uniprot_features=()):
        client = FakeClient(consequences or [consequence(exon="3/3")], length=length,
                            uniprot_features=uniprot_features)
        provider = ProteinRegionProvider(client, "116")
        with patch("acmg_pipeline.providers.protein_region.gene_to_uniprot_accession",
                  return_value=accession):
            records = provider.get_protein_region(self.variant, gene, transcript, hgvsc,
                                                   protein_start)
        return records, client

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

    def test_no_uniprot_accession_leaves_nf04_nf06_unanswered(self):
        """Honest gap, not a fabricated default - accession=None simulates
        gene_to_uniprot_accession() failing to resolve."""
        record = self.region(accession=None)[0][0]
        self.assertIsNone(record.get("critical_region_disrupted"))
        self.assertNotIn("region_biologically_relevant", record)

    def test_a_named_domain_overlapping_the_lost_interval_sets_nf04(self):
        record = self.region(
            protein_start=100, length=213, accession="TEST1",
            uniprot_features=[uniprot_feature("Domain", 90, 150, "Test domain")],
        )[0][0]
        self.assertTrue(record["critical_region_disrupted"])
        self.assertTrue(record["region_biologically_relevant"])
        self.assertEqual(record["region_relevance_method"], UNIPROT_METHOD)

    def test_the_whole_protein_chain_feature_does_not_count_as_critical(self):
        """"Chain" always spans the full protein - it would make NF04 meaningless
        if it counted, so it is excluded (see module docstring)."""
        record = self.region(
            protein_start=204, length=213, accession="P40337",
            uniprot_features=[uniprot_feature("Chain", 1, 213, "Test protein")],
        )[0][0]
        self.assertIsNone(record["critical_region_disrupted"])

    def test_a_domain_spanning_most_of_the_protein_still_counts(self):
        """A real Domain feature counts regardless of how much of the protein it
        spans - MYOC's real "Olfactomedin-like" Domain (residues 244-503 of 504)
        is exactly this shape and is a real, legitimate functional domain; the
        actual MYOC false positive (2026-09-24) turned out to be driven by
        separate single-residue Binding site features, not this Domain - see
        test_a_single_residue_binding_site_does_not_count_as_critical."""
        record = self.region(
            protein_start=368, length=504, accession="Q99972",
            uniprot_features=[uniprot_feature("Domain", 244, 503, "Olfactomedin-like")],
        )[0][0]
        self.assertTrue(record["critical_region_disrupted"])

    def test_a_single_residue_binding_site_does_not_count_as_critical(self):
        """Real MYOC false positive (2026-09-24): 5 separate single-residue Ca2+
        binding sites (residues 380/428/429/477/478 of a 504-residue protein,
        PDB-backed) are scattered across nearly the whole back half of the
        protein, so almost any sufficiently-upstream truncation's lost interval
        overlaps at least one - not what PVS1's "critical functional domain"
        concept means (see uniprot_features.CRITICAL_FEATURE_TYPES's own
        comment)."""
        record = self.region(
            protein_start=300, length=504, accession="Q99972",
            uniprot_features=[uniprot_feature("Binding site", 380, 380, "")],
        )[0][0]
        self.assertIsNone(record["critical_region_disrupted"])

    def test_a_short_motif_does_not_count_as_critical(self):
        """Real PTEN false positive (2026-09-24): the "PDZ domain-binding" Motif
        (residues 401-403 of 403) sits at the very C-terminal tail; a "Motif" is
        UniProt's category for short linear sequence signals, not an identifiable
        functional domain in the ACMG PVS1 sense."""
        record = self.region(
            protein_start=378, length=403, accession="P60484",
            uniprot_features=[uniprot_feature("Motif", 401, 403, "PDZ domain-binding")],
        )[0][0]
        self.assertIsNone(record["critical_region_disrupted"])

    def test_an_active_site_still_counts_even_as_a_single_residue(self):
        """Unlike "Binding site"/"Motif", "Active site" marks the literal catalytic
        residue - exactly the kind of position PVS1's NF04 gate is asking about,
        so it is kept even though it is typically a single residue too (PTEN's
        own Active site 124: "Phosphocysteine intermediate")."""
        record = self.region(
            protein_start=124, length=403, accession="P60484",
            uniprot_features=[uniprot_feature("Active site", 124, 124, "Phosphocysteine intermediate")],
        )[0][0]
        self.assertTrue(record["critical_region_disrupted"])

    def test_vhl_c610gt_real_case_no_domain_but_still_relevant(self):
        """VHL c.610G>T ground truth (ClinGen: PVS1_Moderate): real UniProt data for
        P40337 has no Domain/Region over residues 204-213 (its domains stop at 192)
        and no "Disordered" annotation there either - region_biologically_relevant
        must still default True, or this real case would come out wrong."""
        record = self.region(
            protein_start=204, length=213, accession="P40337",
            uniprot_features=[
                uniprot_feature("Chain", 1, 213, "von Hippel-Lindau disease tumor suppressor"),
                uniprot_feature("Domain", 63, 154, "Beta domain"),
                uniprot_feature("Domain", 155, 192, "Alpha domain"),
            ],
        )[0][0]
        self.assertIsNone(record["critical_region_disrupted"])
        self.assertTrue(record["region_biologically_relevant"])

    def test_full_disordered_coverage_sets_nf06_false(self):
        record = self.region(
            protein_start=204, length=213, accession="TEST2",
            uniprot_features=[uniprot_feature("Region", 200, 213, "Disordered")],
        )[0][0]
        self.assertFalse(record["region_biologically_relevant"])

    def test_partial_disordered_coverage_still_defaults_relevant(self):
        """Only part of the lost interval is disordered - not proof the whole thing is
        unimportant, so the default (relevant) still holds."""
        record = self.region(
            protein_start=204, length=213, accession="TEST3",
            uniprot_features=[uniprot_feature("Region", 200, 208, "Disordered")],
        )[0][0]
        self.assertTrue(record["region_biologically_relevant"])

    def test_a_disordered_region_is_never_read_as_critical(self):
        record = self.region(
            protein_start=204, length=213, accession="TEST4",
            uniprot_features=[uniprot_feature("Region", 200, 213, "Disordered")],
        )[0][0]
        self.assertIsNone(record["critical_region_disrupted"])

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
