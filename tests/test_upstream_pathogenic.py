"""IC03's asymmetry: a yes needs one placed variant, a no needs every one of them placed.

ClinVar positions a protein change only when it reads as a plain substitution, and the
nonsense and frameshift changes it cannot place are the ones most likely to sit near the
N-terminus. So "nothing upstream" is a much stronger claim than "something upstream", and
the provider only makes it when the data supports it.
"""

import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.upstream_pathogenic import (
    METHOD, UpstreamPathogenicProvider, gene_pathogenic_search,
)


def document(uid, protein_change, *, gene="GJB2", classification="Pathogenic"):
    return {
        "uid": uid, "accession": f"VCV{uid}", "protein_change": protein_change,
        "gene_sort": gene, "genes": [{"symbol": gene}],
        "germline_classification": {"description": classification},
    }


class FakeClient:
    def __init__(self, documents, *, count=None, retmax=5000):
        self.documents = documents
        self.count = len(documents) if count is None else count
        self.retmax = retmax
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        if "esearch" in url:
            body = {"esearchresult": {"count": str(self.count), "retmax": str(self.retmax),
                                      "idlist": [str(d["uid"]) for d in self.documents]}}
        else:
            result = {str(d["uid"]): d for d in self.documents}
            result["uids"] = list(result)
            body = {"result": result}
        return {"body": body, "body_sha256": "abc123",
                "retrieved_at": "2026-09-17T00:00:00Z"}


class UpstreamPathogenicTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "13", 20189546, "A", "G")

    def evidence(self, documents, codon=34, *, count=None, gene="GJB2",
                 transcript="NM_004004.6"):
        """Returns (fields, client). The fields are merged into IC02's record, not a
        record of their own - two initiation_assessments are a conflict to PVS1."""
        client = FakeClient(documents, count=count)
        provider = UpstreamPathogenicProvider(client, "2026-09-15")
        return provider.get_upstream_evidence(self.variant, gene, transcript, codon), client

    def test_one_placed_variant_upstream_settles_it_as_yes(self):
        fields, _ = self.evidence([document(1, "R32W"), document(2, "L90P")])
        self.assertIs(fields["upstream_pathogenic_evidence"], True)
        self.assertEqual(fields["upstream_variants"][0]["positions"], [32])

    def test_a_variant_at_the_restart_codon_is_not_upstream_of_it(self):
        fields, _ = self.evidence([document(1, "M34V")], codon=34)
        self.assertIs(fields["upstream_pathogenic_evidence"], False)

    def test_no_when_everything_returned_could_be_placed(self):
        fields, _ = self.evidence([document(1, "L90P"), document(2, "R143W")])
        self.assertIs(fields["upstream_pathogenic_evidence"], False)
        self.assertEqual(fields["upstream_unplaced_count"], 0)

    def test_an_unplaceable_variant_withholds_the_no(self):
        """A frameshift near the N-terminus carries no position and could be the evidence."""
        self.assertEqual(self.evidence([document(1, "L90P"), document(2, "")])[0], {})
        self.assertEqual(self.evidence([document(1, "L90P"), document(2, None)])[0], {})

    def test_an_unplaceable_variant_does_not_withhold_the_yes(self):
        fields, _ = self.evidence([document(1, "R32W"), document(2, "")])
        self.assertIs(fields["upstream_pathogenic_evidence"], True)
        self.assertEqual(fields["upstream_unplaced_count"], 1)

    def test_a_truncated_search_supports_neither_answer(self):
        self.assertEqual(self.evidence([document(1, "R32W")], count=9999)[0], {})

    def test_records_for_other_genes_are_ignored(self):
        fields, _ = self.evidence([document(1, "R32W", gene="OTHER"),
                                   document(2, "L90P")])
        self.assertIs(fields["upstream_pathogenic_evidence"], False)

    def test_the_search_is_not_restricted_to_missense(self):
        """A nonsense variant answers the question as well as a substitution does, so the
        term must not carry the molecular-consequence filter PM1 and PM5 use."""
        _records, client = self.evidence([document(1, "R32W")])
        self.assertIn("esearch", client.urls[0])
        self.assertNotIn("missense", client.urls[0])

    def test_a_codon_of_one_or_less_has_no_upstream_to_search(self):
        self.assertEqual(self.evidence([document(1, "R32W")], codon=1)[0], {})
        self.assertEqual(self.evidence([document(1, "R32W")], codon=0)[0], {})

    def test_a_missing_codon_or_identifier_is_not_an_error(self):
        self.assertEqual(self.evidence([document(1, "R32W")], codon=None)[0], {})
        self.assertEqual(self.evidence([document(1, "R32W")], gene="")[0], {})
        self.assertEqual(self.evidence([document(1, "R32W")], transcript="")[0], {})

    def test_it_answers_only_ic03(self):
        fields = self.evidence([document(1, "R32W")])[0]
        self.assertNotIn("intact_alternative_transcript", fields)
        self.assertNotIn("downstream_in_frame_start", fields)

    def test_the_merged_fields_do_not_overwrite_the_initiation_record(self):
        """They are folded into IC02's record, so every key is namespaced apart from the
        answer itself."""
        fields = self.evidence([document(1, "R32W")])[0]
        collisions = set(fields) - {"upstream_pathogenic_evidence"}
        self.assertTrue(all(key.startswith("upstream_") for key in collisions), collisions)

    def test_the_source_of_the_answer_travels_with_it(self):
        fields = self.evidence([document(1, "R32W")])[0]
        self.assertEqual(fields["upstream_method"], METHOD)
        self.assertTrue(fields["upstream_source"] and fields["upstream_source_version"])

    def test_the_search_reports_whether_it_was_complete(self):
        client = FakeClient([document(1, "R32W")], count=9999)
        found = gene_pathogenic_search(client, "2026-09-15", "GJB2")
        self.assertFalse(found["complete"])
        self.assertEqual(found["count"], 9999)

    def test_a_release_is_required(self):
        with self.assertRaises(ValueError):
            UpstreamPathogenicProvider(FakeClient([]), "")


if __name__ == "__main__":
    unittest.main()
