"""VCEP criteria specifications read as PVS1's mechanism gate, and the limits on that reading."""

import json
import tempfile
import unittest
from pathlib import Path

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.providers.cspec_applicability import (
    METHOD, POLICY_VERSION, CSpecApplicabilityProvider,
)

VARIANT = Variant("GRCh38", "2", 71590212, "TG", "AA")


def snapshot(entries, **overrides):
    document = {
        "schema_version": "1.0",
        "registry_status": "APPROVED",
        "source": "ClinGen Criteria Specification Registry (CSpec)",
        "source_url": "https://cspec.genome.network/cspec/api/svis",
        "source_version": "test-2026-09-18",
        "retrieved_at": "2026-09-18T00:00:00+00:00",
        "entries": entries,
    }
    document.update(overrides)
    return document


def entry(pvs1="Applicable", mondo=("MONDO:0015152",),
          inheritance=("Autosomal recessive inheritance",), svis=("GN180",), **criteria):
    return {
        "criteria": {"PVS1": pvs1, **criteria} if pvs1 is not None else dict(criteria),
        "criteria_missing": [],
        "mondo": list(mondo),
        "inheritance": list(inheritance),
        "svis": list(svis),
    }


class ProviderTests(unittest.TestCase):
    def provider(self, document):
        path = Path(self.tmp.name) / "snapshot.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        return CSpecApplicabilityProvider(path)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def records(self, document, gene="DYSF", variant=VARIANT):
        return self.provider(document).get_mechanism(variant, gene)

    # --- the mapping ----------------------------------------------------

    def test_applicable_becomes_an_established_mechanism_for_the_named_disease(self):
        record, = self.records(snapshot({"DYSF": entry("Applicable")}))
        self.assertIs(record["lof_mechanism_established"], True)
        self.assertEqual(record["gene"], "DYSF")
        self.assertEqual(record["condition"], "MONDO:0015152")
        self.assertEqual(record["inheritance"], "autosomal_recessive")
        self.assertEqual(record["category"], "gene_disease")
        self.assertEqual(record["variant_key"], VARIANT.key)
        self.assertEqual(record["quality_status"], "PASS")

    def test_not_applicable_becomes_false_rather_than_a_missing_value(self):
        """A VCEP that struck PVS1 out of its own specification has ruled on the premise.

        This is the MYH7 case. Reporting it as unknown would send a curator looking for a
        mechanism the panel already declined to accept.
        """
        record, = self.records(snapshot({"MYH7": entry("Not applicable")}), gene="MYH7")
        self.assertIs(record["lof_mechanism_established"], False)

    def test_every_record_says_the_mapping_was_a_decision_and_which_one(self):
        record, = self.records(snapshot({"DYSF": entry()}))
        self.assertEqual(record["method"], METHOD)
        self.assertEqual(record["policy_version"], POLICY_VERSION)
        self.assertEqual(record["cspec_pvs1_applicability"], "Applicable")
        self.assertIn("not a restatement", record["policy_note"])

    def test_the_record_is_automated_so_a_reviewed_assessment_still_outranks_it(self):
        record, = self.records(snapshot({"DYSF": entry()}))
        self.assertEqual(record["assessment_method"], "automated")
        self.assertNotIn("human_signoff", record)
        self.assertNotIn("reviewed_at", record)

    # --- what it refuses to answer --------------------------------------

    def test_only_pvs1_pp2_and_bp1_are_read_even_though_the_snapshot_carries_all_28(self):
        """Every other column (PM3 here) is reference for a curator only - see the
        provider's module docstring."""
        record, = self.records(snapshot({"DYSF": entry(
            "Applicable", PP2="Not applicable", BP1="Applicable", PM3="Applicable")}))
        self.assertNotIn("missense_mechanism_established", record)
        self.assertNotIn("spectrum_review_complete", record)
        self.assertNotIn("predominantly_truncating", record)

    def test_pp2_not_applicable_becomes_false_rather_than_a_missing_value(self):
        """Real false-positive source (2026-09-24, PP2 validation against the 679-variant
        ground truth): SCN2A/SCN1A (Epilepsy Sodium Channel VCEP) both mark PP2 "Not
        applicable" for every evidence strength, but mechanism.py's statistical fallback
        (gnomAD missense constraint) had no way to see that."""
        record, = self.records(snapshot({"DYSF": entry("Applicable", PP2="Not applicable")}))
        self.assertIs(record["pp2_applicable"], False)
        self.assertEqual(record["cspec_pp2_applicability"], "Not applicable")

    def test_bp1_not_applicable_becomes_false_rather_than_a_missing_value(self):
        record, = self.records(snapshot({"DYSF": entry("Applicable", BP1="Not applicable")}))
        self.assertIs(record["bp1_applicable"], False)
        self.assertEqual(record["cspec_bp1_applicability"], "Not applicable")

    def test_pp2_applicable_never_becomes_a_positive_claim(self):
        """A VCEP marking PP2 "Applicable" is not the same statement as a curator having
        reviewed missense_mechanism_established/spectrum_review_complete/
        low_benign_missense_variation for this gene - "Applicable" leaves the field unset
        rather than fabricating a True."""
        record, = self.records(snapshot({"DYSF": entry("Applicable", PP2="Applicable")}))
        self.assertNotIn("pp2_applicable", record)

    def test_bp1_applicable_never_becomes_a_positive_claim(self):
        record, = self.records(snapshot({"DYSF": entry("Applicable", BP1="Applicable")}))
        self.assertNotIn("bp1_applicable", record)

    def test_pp2_not_applicable_alone_is_enough_to_produce_a_record(self):
        """Even when PVS1 itself is never stated for this gene, a decided PP2/BP1 column
        must not be silently dropped."""
        record, = self.records(snapshot({"DYSF": entry(pvs1=None, PP2="Not applicable")}))
        self.assertNotIn("lof_mechanism_established", record)
        self.assertIs(record["pp2_applicable"], False)

    def test_a_gene_the_registry_does_not_cover_produces_nothing(self):
        self.assertEqual(self.records(snapshot({"DYSF": entry()}), gene="SCN5A"), [])

    def test_a_gene_whose_specification_never_states_pvs1_produces_nothing(self):
        document = snapshot({"DYSF": entry(pvs1=None, PM2="Applicable")})
        self.assertEqual(self.records(document), [])

    def test_an_unrecognised_applicability_wording_produces_nothing(self):
        """Neither true nor false: a value this provider does not understand is not evidence."""
        self.assertEqual(self.records(snapshot({"DYSF": entry("Applicable with caveats")})), [])

    def test_a_gene_with_no_mondo_term_produces_nothing(self):
        """A mechanism has to be about a disease; matching on a gene alone is what the
        disease gate exists to prevent."""
        self.assertEqual(self.records(snapshot({"DYSF": entry(mondo=())})), [])

    def test_no_gene_produces_nothing(self):
        self.assertEqual(self.records(snapshot({"DYSF": entry()}), gene=None), [])

    # --- scope and identity ---------------------------------------------

    def test_a_specification_covering_two_diseases_yields_one_record_each(self):
        records = self.records(
            snapshot({"MYH7": entry(mondo=("MONDO:0005045", "MONDO:0005021"))}), gene="MYH7")
        self.assertEqual(sorted(item["condition"] for item in records),
                         ["MONDO:0005021", "MONDO:0005045"])
        self.assertEqual(len({item["evidence_id"] for item in records}), 2,
                         "one id per gene-disease pair, or they collapse onto each other")

    def test_two_variants_in_one_gene_do_not_share_an_evidence_id(self):
        document = snapshot({"DYSF": entry()})
        other = Variant("GRCh38", "2", 71590300, "C", "T")
        first, = self.records(document)
        second, = self.records(document, variant=other)
        self.assertNotEqual(first["evidence_id"], second["evidence_id"])

    def test_more_than_one_inheritance_mode_is_left_unstated_rather_than_picked(self):
        record, = self.records(snapshot({"DYSF": entry(inheritance=(
            "Autosomal recessive inheritance", "Autosomal dominant inheritance"))}))
        self.assertIsNone(record["inheritance"])

    def test_an_inheritance_wording_the_provider_does_not_know_is_left_unstated(self):
        record, = self.records(snapshot({"DYSF": entry(inheritance=("Somatic mosaicism",))}))
        self.assertIsNone(record["inheritance"])

    # --- refusing to answer quietly --------------------------------------

    def test_a_snapshot_missing_its_provenance_raises_rather_than_answering_nothing(self):
        """Answering [] here is indistinguishable from a gene the registry does not cover,
        and PVS1 would report a missing mechanism for every variant in the run."""
        for field in ("registry_status", "source", "source_url", "retrieved_at", "entries"):
            with self.subTest(missing=field):
                document = snapshot({"DYSF": entry()})
                del document[field]
                with self.assertRaises(ValueError) as caught:
                    self.records(document)
                self.assertIn(field, str(caught.exception))

    def test_a_snapshot_whose_entries_are_not_an_object_raises(self):
        with self.assertRaises(ValueError):
            self.records(snapshot([]))

    def test_a_missing_file_raises(self):
        provider = CSpecApplicabilityProvider(Path(self.tmp.name) / "absent.json")
        with self.assertRaises(OSError):
            provider.get_mechanism(VARIANT, "DYSF")

    def test_the_file_is_read_once_per_provider_instance(self):
        provider = self.provider(snapshot({"DYSF": entry()}))
        provider.get_mechanism(VARIANT, "DYSF")
        provider.path.unlink()
        self.assertTrue(provider.get_mechanism(VARIANT, "DYSF"),
                        "a second lookup must not go back to disk")


class CommittedSnapshotTests(unittest.TestCase):
    """The provider against the snapshot the project actually ships."""

    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parent.parent
        cls.provider = CSpecApplicabilityProvider(root / "config" / "cspec_applicability.json")

    def test_it_answers_the_gene_the_other_two_producers_go_silent_for(self):
        record, = self.provider.get_mechanism(VARIANT, "DYSF")
        self.assertIs(record["lof_mechanism_established"], True)
        self.assertEqual(record["condition"], "MONDO:0015152")
        self.assertEqual(record["inheritance"], "autosomal_recessive")

    def test_it_reproduces_the_hand_written_review_on_myh7(self):
        records = self.provider.get_mechanism(VARIANT, "MYH7")
        self.assertTrue(records)
        for record in records:
            with self.subTest(condition=record["condition"]):
                self.assertIs(record["lof_mechanism_established"], False)


if __name__ == "__main__":
    unittest.main()
