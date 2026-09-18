"""CSpec applicability collection, and the boundaries the committed draft must keep.

Every test here is offline: collect() reads through a disk cache, so seeding that cache with
synthetic documents exercises the real fetch/parse path without touching the registry.
"""

import json
import tempfile
import unittest
from pathlib import Path

from test_data.collect_cspec_applicability import (
    ALL_CODES, applicability, collect, gene_entries,
)

DRAFT = Path(__file__).resolve().parent.parent / "config" / "cspec_applicability_draft.json"


def strength(value):
    return {"@type": "EvidenceStrength", "applicability": value}


def code(label, *values):
    return {"@type": "CriteriaCode", "label": label,
            "evidenceStrengths": [strength(value) for value in values]}


def gene(label="TEST", mondo=("MONDO:0000001",), inheritance=("Autosomal recessive inheritance",)):
    return {
        "@id": f"https://www.genenames.org/tools/search/#!/?query={label}",
        "label": label,
        "diseases": [{
            "label": item,
            "modeOfInheritance": [{"@label": mode} for mode in inheritance],
        } for item in mondo],
    }


def svi(gn_id, rule_sets):
    return {"@id": f"https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/{gn_id}",
            "ruleSets": rule_sets}


def seed(cache_dir, documents):
    """Write the exact filenames collect()'s fetch() looks for, so nothing hits the network."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / "svis.json").write_text(
        json.dumps({"data": [svi(gn, []) for gn in documents]}), encoding="utf-8")
    for gn_id, rule_sets in documents.items():
        (cache_dir / f"svi_{gn_id}.json").write_text(
            json.dumps(svi(gn_id, rule_sets)), encoding="utf-8")


class ApplicabilityFoldingTests(unittest.TestCase):
    """One `applicability` per criterion, folded from the strengths that carry it."""

    def test_any_applicable_strength_makes_the_criterion_applicable(self):
        self.assertEqual(
            applicability(code("PVS1", "Not applicable", "Applicable", "Not applicable")),
            "Applicable")

    def test_every_strength_saying_not_makes_it_not_applicable(self):
        self.assertEqual(applicability(code("PP2", "Not applicable", "Not applicable")),
                         "Not applicable")

    def test_the_registrys_mixed_casing_is_not_read_as_a_third_value(self):
        # The live registry emits both "Not Applicable" and "Not applicable" in one document.
        self.assertEqual(applicability(code("BP4", "Not Applicable", "Not applicable")),
                         "Not applicable")
        self.assertEqual(applicability(code("BP4", "Not Applicable", "APPLICABLE")), "Applicable")

    def test_a_qualified_affirmative_still_means_applicable(self):
        """The wordings a VCEP uses when it specified or adopted the criterion.

        Reading these as negatives is the bug that reported the Lysosomal Diseases VCEP as
        having struck out all 28 criteria for GAA, PVS1 included.
        """
        for wording in ("Applicable with VCEP specification", "Applicable as originally described"):
            with self.subTest(wording=wording):
                self.assertEqual(applicability(code("PVS1", "Not applicable", wording)),
                                 "Applicable")

    def test_not_applicable_for_this_vcep_is_a_negative_despite_containing_the_word(self):
        self.assertEqual(applicability(code("BP6", "Not Applicable for this VCEP")),
                         "Not applicable")

    def test_the_six_wordings_the_registry_actually_uses(self):
        """Counted across all 208 documents on 2026-09-18, not read off one of them."""
        affirmative = ("Applicable", "Applicable with VCEP specification",
                       "Applicable as originally described")
        negative = ("Not applicable", "Not Applicable", "Not Applicable for this VCEP")
        for wording in affirmative:
            self.assertEqual(applicability(code("PVS1", wording)), "Applicable", wording)
        for wording in negative:
            self.assertEqual(applicability(code("PVS1", wording)), "Not applicable", wording)

    def test_a_criterion_stating_nothing_is_unknown_rather_than_negative(self):
        self.assertIsNone(applicability(code("PS2")))
        self.assertIsNone(applicability({"label": "PS2"}))
        self.assertIsNone(applicability(code("PS2", "", "")))


class GeneScopeTests(unittest.TestCase):
    """The disease scope a mechanism record would be anchored to."""

    def test_symbol_disease_and_inheritance_are_read_together(self):
        entries = list(gene_entries({"genes": [gene("DYSF", ("MONDO:0015152",))]}))
        self.assertEqual(entries, [("DYSF", ["MONDO:0015152"], ["Autosomal recessive inheritance"])])

    def test_a_mondo_iri_is_accepted_where_the_label_is_not_an_identifier(self):
        entries = list(gene_entries({"genes": [{
            "@id": "https://www.genenames.org/tools/search/#!/?query=PTEN",
            "label": "PTEN",
            "diseases": [{"@id": "http://purl.obolibrary.org/obo/MONDO_0017623",
                          "label": "PTEN hamartoma tumor syndrome"}],
        }]}))
        self.assertEqual(entries, [("PTEN", ["MONDO:0017623"], [])])

    def test_the_symbol_falls_back_to_the_query_in_the_iri(self):
        entries = list(gene_entries({"genes": [{
            "@id": "https://www.genenames.org/tools/search/#!/?query=BRCA1", "diseases": []}]}))
        self.assertEqual(entries, [("BRCA1", [], [])])

    def test_a_gene_with_no_readable_symbol_is_skipped_not_guessed(self):
        self.assertEqual(list(gene_entries({"genes": [{"@id": "urn:opaque", "diseases": []}]})), [])

    def test_a_disease_with_no_identifier_contributes_no_scope(self):
        entries = list(gene_entries({"genes": [{
            "@id": "https://www.genenames.org/tools/search/#!/?query=TEST", "label": "TEST",
            "diseases": [{"label": "some disease with no identifier"}]}]}))
        self.assertEqual(entries, [("TEST", [], [])])


class CollectTests(unittest.TestCase):
    """End to end over a seeded cache."""

    def collect(self, documents):
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp) / "cspec"
            seed(cache, documents)
            return collect(cache)

    def test_a_gene_carries_its_criteria_scope_and_source_document(self):
        entries, stats = self.collect({
            "GN180": [{"genes": [gene("DYSF", ("MONDO:0015152",))],
                       "criteriaCodes": [code("PVS1", "Applicable"),
                                         code("PP2", "Not applicable")]}]})
        self.assertEqual(entries["DYSF"]["criteria"], {"PVS1": "Applicable", "PP2": "Not applicable"})
        self.assertEqual(entries["DYSF"]["mondo"], ["MONDO:0015152"])
        self.assertEqual(entries["DYSF"]["inheritance"], ["Autosomal recessive inheritance"])
        self.assertEqual(entries["DYSF"]["svis"], ["GN180"])
        self.assertEqual(stats["genes_with_criteria"], 1)
        self.assertEqual(stats["failures"], [])

    def test_criteria_the_document_never_states_are_listed_not_filled_in(self):
        entries, _ = self.collect({"GN001": [{"genes": [gene()],
                                              "criteriaCodes": [code("PVS1", "Applicable")]}]})
        missing = entries["TEST"]["criteria_missing"]
        self.assertNotIn("PVS1", missing)
        self.assertIn("PP2", missing)
        self.assertEqual(len(missing), len(ALL_CODES) - 1)

    def test_criteria_are_ordered_by_the_acmg_sequence_not_the_registrys(self):
        entries, _ = self.collect({"GN001": [{"genes": [gene()], "criteriaCodes": [
            code("BP7", "Applicable"), code("PVS1", "Applicable"), code("PM2", "Applicable")]}]})
        self.assertEqual(list(entries["TEST"]["criteria"]), ["PVS1", "PM2", "BP7"])

    def test_a_gene_specified_by_two_vceps_keeps_both_documents_and_diseases(self):
        entries, _ = self.collect({
            "GN002": [{"genes": [gene("MYH7", ("MONDO:0005045",))],
                       "criteriaCodes": [code("PVS1", "Not applicable")]}],
            "GN003": [{"genes": [gene("MYH7", ("MONDO:0005021",))],
                       "criteriaCodes": [code("PM2", "Applicable")]}],
        })
        self.assertEqual(entries["MYH7"]["svis"], ["GN002", "GN003"])
        self.assertEqual(entries["MYH7"]["mondo"], ["MONDO:0005021", "MONDO:0005045"])
        self.assertEqual(entries["MYH7"]["criteria"],
                         {"PVS1": "Not applicable", "PM2": "Applicable"})

    def test_every_rule_set_in_one_document_is_read(self):
        entries, _ = self.collect({"GN004": [
            {"genes": [gene("SGCA")], "criteriaCodes": [code("PVS1", "Applicable")]},
            {"genes": [gene("SGCB")], "criteriaCodes": [code("PVS1", "Not applicable")]},
        ]})
        self.assertEqual(entries["SGCA"]["criteria"]["PVS1"], "Applicable")
        self.assertEqual(entries["SGCB"]["criteria"]["PVS1"], "Not applicable")

    def test_a_document_with_no_criteria_is_counted_and_yields_no_gene_entry(self):
        entries, stats = self.collect({
            "GN005": [{"genes": [gene("EMPTY")], "criteriaCodes": []}],
            "GN006": [{"genes": [gene("REAL")], "criteriaCodes": [code("PVS1", "Applicable")]}],
        })
        self.assertNotIn("EMPTY", entries)
        self.assertEqual(stats["documents_without_criteria"], 1)
        # The gene is still named by the registry even though nothing was stated about it.
        self.assertEqual(stats["genes_named"], 2)
        self.assertEqual(stats["genes_with_criteria"], 1)

    def test_a_criterion_stating_no_strength_is_left_out_of_the_entry(self):
        entries, _ = self.collect({"GN007": [{"genes": [gene()], "criteriaCodes": [
            code("PVS1", "Applicable"), code("PS2")]}]})
        self.assertNotIn("PS2", entries["TEST"]["criteria"])
        self.assertIn("PS2", entries["TEST"]["criteria_missing"])


class CommittedDraftTests(unittest.TestCase):
    """What config/cspec_applicability_draft.json must keep being."""

    @classmethod
    def setUpClass(cls):
        cls.document = json.loads(DRAFT.read_text(encoding="utf-8"))
        cls.entries = cls.document["entries"]

    def test_it_stays_a_draft_that_names_its_source_and_its_rule(self):
        self.assertEqual(self.document["registry_status"], "DRAFT")
        self.assertEqual(self.document["entry_method"], "machine_collected")
        self.assertTrue(self.document["retrieved_at"])
        self.assertIn("cspec.genome.network", self.document["source_url"])
        self.assertIn("MUST NOT", self.document["purpose"])
        self.assertTrue(self.document["not_verified_by_a_second_reader"])

    def test_it_never_states_a_mechanism_in_the_review_files_vocabulary(self):
        """`PVS1: Applicable` is not `lof_mechanism_established`; the draft must not blur that.

        Deciding whether the two are the same sentence is the curator's call (see the
        collector's docstring), so a machine-collected file asserting the mechanism directly
        would have made that call silently.
        """
        raw = json.dumps(self.entries)
        for field in ("lof_mechanism_established", "pp2_applicable", "bp1_applicable",
                      "missense_mechanism_established", "human_signoff", "reviewed_at"):
            self.assertNotIn(field, raw, f"entries must not speak in {field}")
        # The prose is where the distinction belongs, and it has to keep naming it.
        self.assertIn("lof_mechanism_established", self.document["purpose"])

    def test_every_entry_carries_criteria_and_a_disease_scope(self):
        self.assertTrue(self.entries)
        for symbol, entry in self.entries.items():
            with self.subTest(gene=symbol):
                self.assertTrue(entry["criteria"], "an entry with no criteria is not usable")
                self.assertTrue(entry["mondo"], "a mechanism record needs a disease to be about")
                self.assertTrue(entry["svis"])
                for value in entry["criteria"].values():
                    self.assertIn(value, ("Applicable", "Not applicable"))
                stated = set(entry["criteria"]) | set(entry["criteria_missing"])
                self.assertEqual(stated, set(ALL_CODES))

    def test_the_statistics_describe_the_entries_beside_them(self):
        stats = self.document["statistics"]
        self.assertEqual(stats["genes_with_criteria"], len(self.entries))
        self.assertEqual(stats["genes_with_all_28"],
                         sum(1 for e in self.entries.values() if not e["criteria_missing"]))
        self.assertEqual(stats["failures"], [], "a partial crawl must not be published as a draft")

    def test_it_answers_the_gene_the_three_wired_up_routes_go_silent_for(self):
        # DYSF c.3498_3499delinsAA: ClinGen Dosage emits nothing for a recessive gene, G2P has
        # no record, and the hand-written review covers five genes - see doc section 1-7.
        dysf = self.entries["DYSF"]
        self.assertEqual(dysf["criteria"]["PVS1"], "Applicable")
        self.assertIn("MONDO:0015152", dysf["mondo"])
        self.assertEqual(dysf["inheritance"], ["Autosomal recessive inheritance"])

    def test_it_agrees_with_the_hand_written_review_on_pvs1(self):
        """The three genes where both sources have a disease-scoped answer.

        This is the evidence that reading VCEP applicability reproduces the judgment a curator
        already made by hand - not a licence to map the two fields onto each other.
        """
        review = json.loads(
            (DRAFT.parent / "gene-disease-review-decisions.json").read_text(encoding="utf-8"))
        compared = 0
        for group in review["groups"]:
            entry = self.entries.get(group["gene"])
            if not entry or group["condition"] not in entry["mondo"]:
                continue
            compared += 1
            with self.subTest(gene=group["gene"], condition=group["condition"]):
                self.assertEqual(entry["criteria"]["PVS1"] == "Applicable",
                                 group["lof_mechanism_established"])
        self.assertEqual(compared, 3, "the overlap with the hand-written review changed")


if __name__ == "__main__":
    unittest.main()
