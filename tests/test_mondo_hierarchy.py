"""MONDO ancestry read for relatedness only, and refused when it cannot be settled."""

import unittest

from acmg_pipeline.criteria.common import ontology_related
from acmg_pipeline.providers.mondo_hierarchy import (
    MAX_ANCESTORS, METHOD, MondoHierarchyProvider,
)


def page(obo_ids, total=None):
    return {"_embedded": {"terms": [{"obo_id": item} for item in obo_ids]},
            "page": {"totalElements": len(obo_ids) if total is None else total}}


class FakeClient:
    def __init__(self, body):
        self.body = body
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        return {"body": self.body, "retrieved_at": "2026-09-17T00:00:00Z"}


class MondoHierarchyTests(unittest.TestCase):
    def provider(self, body):
        return MondoHierarchyProvider(FakeClient(body))

    def test_ancestors_come_back_as_mondo_terms(self):
        provider = self.provider(page(["MONDO:0005045", "MONDO:0004994"]))
        result = provider.ancestors("MONDO:0007268")
        self.assertEqual(result["ancestors"], ["MONDO:0004994", "MONDO:0005045"])
        self.assertEqual(result["term"], "MONDO:0007268")
        self.assertEqual(result["method"], METHOD)
        self.assertIn("MONDO_0007268", provider.client.urls[0])

    def test_upper_level_classes_are_not_ancestors_for_this_purpose(self):
        """BFO's classes sit above every disease, so relating two diseases through them
        would make every pair of diseases related."""
        result = self.provider(
            page(["BFO:0000001", "BFO:0000002", "MONDO:0005045"])).ancestors("MONDO:0007268")
        self.assertEqual(result["ancestors"], ["MONDO:0005045"])

    def test_a_truncated_ancestry_resolves_to_nothing(self):
        """A short ancestry silently standing in for a long one would make a real relation
        look like an unrelated disease."""
        self.assertIsNone(
            self.provider(page(["MONDO:0005045"], total=12)).ancestors("MONDO:0007268"))
        self.assertIsNone(
            self.provider(page(["MONDO:0005045"], total=MAX_ANCESTORS + 1))
            .ancestors("MONDO:0007268"))

    def test_only_a_mondo_term_is_looked_up(self):
        provider = self.provider(page([]))
        for term in ("OMIM:115197", "", None, 7268):
            with self.subTest(term=term):
                self.assertIsNone(provider.ancestors(term))
        self.assertEqual(provider.client.urls, [])


class OntologyRelatedTests(unittest.TestCase):
    def test_related_in_either_direction(self):
        case = {"condition_ancestors": ["MONDO:0005045"]}
        record = {}
        self.assertTrue(ontology_related(case, "MONDO:0007268", record, "MONDO:0005045"))
        self.assertTrue(ontology_related(
            {}, "MONDO:0005045", {"condition_ancestors": ["MONDO:0005045"]}, "MONDO:0007268"))

    def test_the_same_term_is_a_match_not_a_relation(self):
        self.assertFalse(ontology_related(
            {"condition_ancestors": ["MONDO:0005045"]}, "MONDO:0005045", {}, "MONDO:0005045"))

    def test_unrelated_terms_and_missing_ancestry_are_not_related(self):
        self.assertFalse(ontology_related(
            {"condition_ancestors": ["MONDO:0005045"]}, "MONDO:0007268", {}, "MONDO:0009861"))
        self.assertFalse(ontology_related({}, "MONDO:0007268", {}, "MONDO:0005045"))
        self.assertFalse(ontology_related({}, None, {}, "MONDO:0005045"))


if __name__ == "__main__":
    unittest.main()
