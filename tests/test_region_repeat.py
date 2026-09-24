"""PM4/BP3's `region` evidence: a UniProt Repeat/Compositional-bias feature first pass,
plus a Wootton-Federhen/SEG sequence-complexity fallback and isoform escalation, over
the altered protein interval - see region_repeat.py's own module docstring for why an
unanswered position is read as repetitive=False rather than left unset (unlike
protein_region.py's critical_region_disrupted, which this shares its UniProt fetch and
overlap helpers with)."""

import unittest
from unittest.mock import patch

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.region_repeat import METHOD, SEG_TRIGGER_K1, RegionRepeatProvider

RPGR_ORF15_WINDOW = "EGEGEEEEGEEEGEE"  # real UniProt Q92834-6 residues 1030-1045
DIVERSE_WINDOW = "MKTAYIAKQRQISFV"  # varied composition, not low-complexity
# A 20-residue rotating alphabet: any 12+ consecutive-residue slice of it is close to
# maximally diverse (SEG complexity well above the 2.2-bit trigger), so tests that
# only care about the feature-based signal can use it as filler without a stray
# low-complexity window (a real amino-acid sequence's own letter-frequency bias can
# dip under 2.2 in a short window by chance - found the hard way while writing these).
HIGH_COMPLEXITY_FILLER = ("ACDEFGHIKLMNPQRSTVWY" * 60)[:1000]


def uniprot_feature(feature_type, start, end, description=""):
    return {"type": feature_type, "description": description,
            "location": {"start": {"value": start}, "end": {"value": end}}}


class FakeClient:
    """Serves one UniProt entry per accession, keyed by the `.json` filename in the
    URL - `entries[accession]` is a dict with `features`/`sequence`/`comments` keys,
    any of which may be omitted. A `default_accession`'s entry gets a generous 1000
    "A" sequence when it supplies none, so tests that don't care about sequence
    complexity or isoform escalation don't have to construct one.
    """

    def __init__(self, default_accession="P00000", features=(), sequence=None, entries=None):
        self.entries = dict(entries or {})
        if default_accession not in self.entries:
            self.entries[default_accession] = {
                "features": list(features),
                "sequence": sequence if sequence is not None else "A" * 1000,
            }
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        accession = url.rsplit("/", 1)[-1].removesuffix(".json")
        entry = self.entries.get(accession)
        if entry is None:
            return {"body": {}, "body_sha256": f"missing-{accession}", "retrieved_at": "2026-09-21T00:00:00Z"}
        sequence = entry.get("sequence")
        body = {"features": entry.get("features", [])}
        if sequence is not None:
            body["sequence"] = {"length": len(sequence), "value": sequence}
        if entry.get("isoforms"):
            body["comments"] = [{
                "commentType": "ALTERNATIVE PRODUCTS",
                "isoforms": [{"isoformIds": [iso], "isoformSequenceStatus": "Described"}
                            for iso in entry["isoforms"]],
            }]
        return {"body": body, "body_sha256": f"sha-{accession}", "retrieved_at": "2026-09-21T00:00:00Z"}


class RegionRepeatTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "3", 10146618, "G", "T")

    def region(self, start=100, end=105, *, gene="TEST", transcript="NM_000256.3",
               protein_id="NP_TEST.1", accession="P00000", uniprot_features=(),
               sequence=None, entries=None):
        client = FakeClient(default_accession=accession or "P00000", features=uniprot_features,
                            sequence=sequence, entries=entries)
        provider = RegionRepeatProvider(client, "116")
        with patch("acmg_pipeline.providers.region_repeat.gene_to_uniprot_accession",
                  return_value=accession):
            records = provider.get_region(self.variant, gene, transcript, protein_id, start, end)
        return records, client

    def test_no_repeat_feature_over_the_interval_is_repetitive_false_not_unset(self):
        """Unlike protein_region.py's critical_region_disrupted, absence IS the signal
        here - see the module docstring for why."""
        record = self.region(sequence=HIGH_COMPLEXITY_FILLER,
                             uniprot_features=[uniprot_feature("Domain", 1, 50, "Unrelated domain")])[0][0]
        self.assertFalse(record["repetitive"])
        self.assertFalse(record["nonfunctional_repeat"])
        self.assertTrue(record["functional_review_complete"])

    def test_a_compositional_bias_feature_also_sets_repetitive(self):
        """UniProt types a low-complexity/homopolymeric run (poly-Pro, poly-Gln...) as
        "Compositional bias", not "Repeat" - found empirically from FOXG1's real UniProt
        entry (P55316: "Compositional bias" 59-85 "Pro residues") after the 679-variant
        ground-truth run showed BP3 missing exactly this kind of real case."""
        record = self.region(sequence=HIGH_COMPLEXITY_FILLER,
                             uniprot_features=[uniprot_feature("Compositional bias", 95, 110, "Pro residues")])[0][0]
        self.assertTrue(record["repetitive"])
        self.assertFalse(record["functional_importance"])
        self.assertTrue(record["nonfunctional_repeat"])

    def test_a_repeat_feature_over_the_interval_sets_repetitive(self):
        record = self.region(sequence=HIGH_COMPLEXITY_FILLER,
                             uniprot_features=[uniprot_feature("Repeat", 95, 110, "ANK 3")])[0][0]
        self.assertTrue(record["repetitive"])
        self.assertFalse(record["functional_importance"])
        self.assertTrue(record["nonfunctional_repeat"])

    def test_a_repeat_with_an_overlapping_domain_is_not_nonfunctional(self):
        record = self.region(sequence=HIGH_COMPLEXITY_FILLER, uniprot_features=[
            uniprot_feature("Repeat", 95, 110, "ANK 3"),
            uniprot_feature("Domain", 90, 150, "Ankyrin repeat domain"),
        ])[0][0]
        self.assertTrue(record["repetitive"])
        self.assertTrue(record["functional_importance"])
        self.assertFalse(record["nonfunctional_repeat"])

    def test_a_disordered_region_feature_does_not_count_as_functional(self):
        record = self.region(sequence=HIGH_COMPLEXITY_FILLER, uniprot_features=[
            uniprot_feature("Repeat", 95, 110, "ANK 3"),
            uniprot_feature("Region", 90, 150, "Disordered"),
        ])[0][0]
        self.assertTrue(record["repetitive"])
        self.assertFalse(record["functional_importance"])

    def test_a_feature_with_an_unresolved_fuzzy_boundary_is_not_a_crash(self):
        """Real UniProt data: an unresolved/fuzzy terminus is recorded with value=null,
        which parses without KeyError/TypeError and must still be rejected."""
        record = self.region(sequence=HIGH_COMPLEXITY_FILLER, uniprot_features=[
            {"type": "Repeat", "description": "fuzzy",
             "location": {"start": {"value": None}, "end": {"value": 110}}},
        ])[0][0]
        self.assertFalse(record["repetitive"])

    def test_a_repeat_outside_the_interval_does_not_count(self):
        record = self.region(start=100, end=105, sequence=HIGH_COMPLEXITY_FILLER,
                             uniprot_features=[uniprot_feature("Repeat", 500, 520, "Unrelated repeat")])[0][0]
        self.assertFalse(record["repetitive"])

    def test_no_uniprot_accession_yields_nothing(self):
        self.assertEqual(self.region(accession=None)[0], [])

    def test_no_feature_list_yields_nothing_but_not_an_error(self):
        """An entry that resolves but has no readable `features` list still produces a
        record, scored on sequence complexity alone (see module docstring)."""
        records, _ = self.region(sequence="A" * 1000, entries={"P00000": {}})
        self.assertEqual(len(records), 1)
        self.assertFalse(records[0]["repetitive"])

    def test_missing_identifiers_are_not_an_error(self):
        self.assertEqual(self.region(gene="")[0], [])
        self.assertEqual(self.region(transcript="")[0], [])

    def test_a_non_positive_or_absent_start_yields_nothing(self):
        for start in (None, "100", 0, -1):
            with self.subTest(start=start):
                self.assertEqual(self.region(start=start)[0], [])

    def test_an_end_before_start_falls_back_to_start(self):
        record = self.region(start=100, end=50, sequence=HIGH_COMPLEXITY_FILLER,
                             uniprot_features=[uniprot_feature("Repeat", 100, 100, "single residue")])[0][0]
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

    # -- Sequence-complexity fallback (RPGR ORF15 has no UniProt feature at all) --

    def test_a_low_complexity_window_sets_repetitive_with_no_feature_at_all(self):
        prefix = "M" * 1029
        record = self.region(start=1030, end=1044, sequence=prefix + RPGR_ORF15_WINDOW)[0][0]
        self.assertTrue(record["repetitive"])
        self.assertFalse(record["functional_importance"])
        self.assertIsNotNone(record["sequence_complexity"])
        self.assertLessEqual(record["sequence_complexity"], SEG_TRIGGER_K1)

    def test_a_diverse_window_does_not_set_repetitive(self):
        prefix = "M" * 99
        record = self.region(start=100, end=114, sequence=prefix + DIVERSE_WINDOW)[0][0]
        self.assertFalse(record["repetitive"])
        self.assertGreater(record["sequence_complexity"], SEG_TRIGGER_K1)

    def test_a_short_indel_with_ordinary_flanking_sequence_does_not_falsely_trigger(self):
        """A short (2-residue) interval whose immediate 12-mer neighborhood happens to
        dip under 2.2 bits by ordinary amino-acid frequency alone - not a real repeat -
        must not be flagged once the surrounding, more diverse sequence is also counted.
        Found on real GAA/HNF1A/PAH/LDLR/CDKL5 UniProt data (window=12 gave 1.53-2.15
        bits on real, non-repetitive PM4-met ground-truth variants after the first
        679-variant run with SEG_WINDOW=12 regressed PM4 66.7%->37.5% recall) - this is
        the same real GAA sequence (residues ~19-44) that produced that regression;
        see SEG_WINDOW's own comment."""
        real_gaa_fragment = (
            "MGVRHPPCSHRLLAVCALVSLATAALLGHILLHDFLLVPRELSGSSPVLE"
        )  # UniProt P10253 residues 1-51 (signal peptide + N-terminus)
        record = self.region(start=31, end=32, sequence=real_gaa_fragment.ljust(200, "A"))[0][0]
        self.assertFalse(record["repetitive"])

    def test_a_short_interval_is_padded_to_the_seg_window_before_scoring(self):
        """A 1-residue interval alone has zero complexity trivially - SEG's own minimum
        window (12) must be used, built from the surrounding real sequence."""
        sequence = "M" * 50 + RPGR_ORF15_WINDOW + "M" * 50
        record = self.region(start=53, end=53, sequence=sequence)[0][0]
        self.assertTrue(record["repetitive"])

    def test_no_sequence_leaves_complexity_unscored(self):
        records, _ = self.region(sequence=None, entries={"P00000": {"sequence": None}})
        self.assertIsNone(records[0]["sequence_complexity"])

    # -- Isoform escalation (RPGR-ORF15 numbering exceeds the canonical entry) --

    def test_a_position_past_the_canonical_length_escalates_to_a_covering_isoform(self):
        canonical_seq = "M" * 1020
        iso_seq = "M" * 1029 + RPGR_ORF15_WINDOW + "M" * 100
        entries = {
            "Q92834": {"sequence": canonical_seq, "isoforms": ["Q92834-6"]},
            "Q92834-6": {"sequence": iso_seq},
        }
        record = self.region(start=1030, end=1044, accession="Q92834", entries=entries)[0][0]
        self.assertEqual(record["source_version"], "Q92834-6")
        self.assertEqual(record["isoform_used"], "Q92834-6")
        self.assertTrue(record["repetitive"])

    def test_a_position_within_canonical_length_never_escalates(self):
        canonical_seq = HIGH_COMPLEXITY_FILLER
        entries = {"Q92834": {"sequence": canonical_seq, "isoforms": ["Q92834-6"]}}
        record = self.region(start=100, end=114, accession="Q92834", entries=entries)[0][0]
        self.assertEqual(record["source_version"], "Q92834")
        self.assertIsNone(record["isoform_used"])

    def test_no_isoform_covers_the_position_yields_nothing(self):
        entries = {
            "Q92834": {"sequence": "M" * 1020, "isoforms": ["Q92834-6"]},
            "Q92834-6": {"sequence": "M" * 1030},
        }
        records, _ = self.region(start=1100, end=1100, accession="Q92834", entries=entries)
        self.assertEqual(records, [])


if __name__ == "__main__":
    unittest.main()
