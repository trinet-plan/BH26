"""MONDO's mapping set read for equivalence only, and never for a guess."""

import unittest

from acmg_pipeline.providers.mondo import MondoMappingProvider, METHOD, parse, version


HEADER = "\t".join(["subject_id", "subject_label", "predicate_id", "object_id",
                    "object_label", "mapping_justification"])
PREAMBLE = ["# curie_map:", "#   MONDO: http://purl.obolibrary.org/obo/MONDO_",
            "# mapping_set_id: https://w3id.org/mondo.sssom.tsv"]


def row(subject, predicate, obj):
    return "\t".join([subject, "label", predicate, obj, "label",
                      "semapv:UnspecifiedMatching"])


def mapping_set(rows, preamble=PREAMBLE):
    return "\n".join([*preamble, HEADER, *rows]) + "\n"


class FakeClient:
    def __init__(self, text):
        self.text = text
        self.calls = 0

    def fetch(self, url, *, response_format="json", **kwargs):
        self.calls += 1
        return {"body": self.text, "retrieved_at": "2026-09-17T00:00:00Z"}


class MondoMappingTests(unittest.TestCase):
    def provider(self, rows, preamble=PREAMBLE):
        return MondoMappingProvider(FakeClient(mapping_set(rows, preamble)))

    def test_an_exact_match_resolves_to_the_mondo_term(self):
        mapping = self.provider(
            [row("MONDO:0007750", "skos:exactMatch", "OMIM:143890")]).normalize("OMIM:143890")
        self.assertEqual(mapping["normalized_condition"], "MONDO:0007750")
        self.assertEqual(mapping["mapping_type"], "equivalent")
        self.assertEqual(mapping["input_condition"], "OMIM:143890")
        self.assertEqual(mapping["method"], METHOD)
        self.assertTrue(mapping["source_version"])
        self.assertEqual(mapping["assessment_method"], "automated")

    def test_a_mondo_identifier_resolves_to_itself_without_reading_the_file(self):
        """It is already in the target vocabulary, and the two are different match levels:
        an identifier match is EXACT, a mapped one is only EQUIVALENT."""
        provider = self.provider([])
        mapping = provider.normalize("MONDO:0007750")
        self.assertEqual(mapping["normalized_condition"], "MONDO:0007750")
        self.assertEqual(mapping["mapping_type"], "identity")
        self.assertEqual(provider.client.calls, 0)

    def test_a_related_concept_is_not_an_equivalence(self):
        """broadMatch says the two are related, which is the distinction the match levels
        exist to keep, so it is not read as the same disease."""
        mapping = self.provider(
            [row("MONDO:0000001", "skos:broadMatch", "OMIM:143890")]).normalize("OMIM:143890")
        self.assertIsNone(mapping)

    def test_an_identifier_resolving_to_two_terms_resolves_to_neither(self):
        """Picking one would choose a disease on the curator's behalf."""
        mapping = self.provider([
            row("MONDO:0007750", "skos:exactMatch", "OMIM:143890"),
            row("MONDO:0007751", "skos:exactMatch", "OMIM:143890"),
        ]).normalize("OMIM:143890")
        self.assertIsNone(mapping)

    def test_an_unmapped_or_unusable_identifier_resolves_to_nothing(self):
        provider = self.provider([row("MONDO:0007750", "skos:exactMatch", "OMIM:143890")])
        for condition in ("OMIM:999999", "Phenylketonuria", "", None, 143890):
            with self.subTest(condition=condition):
                self.assertIsNone(provider.normalize(condition))

    def test_orphanet_resolves_the_same_way(self):
        mapping = self.provider(
            [row("MONDO:0000044", "skos:exactMatch", "Orphanet:437")]).normalize("Orphanet:437")
        self.assertEqual(mapping["normalized_condition"], "MONDO:0000044")

    def test_a_changed_preamble_is_a_different_policy_version(self):
        """The set carries no dated release line, so the preamble is what pins a replay."""
        first = self.provider([]).normalize("OMIM:143890")
        changed = self.provider([row("MONDO:0007750", "skos:exactMatch", "OMIM:143890")],
                                preamble=[*PREAMBLE, "# mapping_set_version: 2026-09-16"])
        self.assertIsNone(first)
        self.assertTrue(changed.normalize("OMIM:143890")["source_version"])
        self.assertNotEqual(
            version(PREAMBLE),
            version([*PREAMBLE, "# mapping_set_version: 2026-09-16"]))

    def test_a_file_without_a_header_is_an_error_not_an_empty_index(self):
        with self.assertRaises(ValueError):
            parse("# only a preamble\n")


if __name__ == "__main__":
    unittest.main()
