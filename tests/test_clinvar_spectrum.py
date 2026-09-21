"""Gene-wide ClinVar missense/truncating classification counts for BP1's
statistical suggestion (see acmg_pipeline.providers.gene_disease_draft)."""

import unittest
from urllib.parse import parse_qs, urlparse

from acmg_pipeline.providers.clinvar_spectrum import ClinvarSpectrumProvider


def summary(uid, classification, gene="TEST"):
    return {uid: {"accession": f"VCV{int(uid):09d}", "gene_sort": gene,
                 "germline_classification": {"description": classification}}}


class FakeClient:
    """Simulates real ESearch paging: `count` is the true total, and each
    call returns the `retstart:retstart+retmax` slice of the full id list -
    unlike a fixed-response fake, this actually exercises
    gene_consequence_search_all()'s multi-page loop, not just its single-
    page fallback path. `missense_count`/`truncating_count` (independent of
    how many ids are actually supplied) simulate ClinVar reporting a total
    the fake's own id lists can never fully deliver, the same "server says
    more exist than pages ever produce" shape a truly incomplete real
    search would have.
    """

    def __init__(self, *, missense_ids=(), truncating_ids=(), summaries=None,
                missense_count=None, truncating_count=None):
        self.missense_ids = list(missense_ids)
        self.truncating_ids = list(truncating_ids)
        self.summaries = summaries or {}
        self.missense_count = missense_count if missense_count is not None else len(self.missense_ids)
        self.truncating_count = (truncating_count if truncating_count is not None
                                 else len(self.truncating_ids))
        self.urls = []

    def fetch(self, url, **kwargs):
        self.urls.append(url)
        if "esearch.fcgi" in url:
            query = parse_qs(urlparse(url).query)
            truncating = "nonsense" in url
            ids = self.truncating_ids if truncating else self.missense_ids
            count = self.truncating_count if truncating else self.missense_count
            retmax = int(query["retmax"][0])
            retstart = int(query.get("retstart", ["0"])[0])
            page = ids[retstart:retstart + retmax]
            return {"body": {"esearchresult": {"count": str(count), "retmax": str(retmax),
                                               "idlist": page}},
                    "retrieved_at": "2026-09-17T00:00:00Z"}
        return {"body": {"result": {"uids": list(self.summaries), **self.summaries}},
                "retrieved_at": "2026-09-17T00:00:01Z"}


class ClinvarSpectrumTests(unittest.TestCase):
    def test_counts_pathogenic_and_benign_missense_and_pathogenic_truncating(self):
        summaries = {}
        summaries.update(summary("1", "Pathogenic"))
        summaries.update(summary("2", "Pathogenic"))
        summaries.update(summary("3", "Benign"))
        summaries.update(summary("4", "Pathogenic"))
        client = FakeClient(missense_ids=["1", "2", "3"], truncating_ids=["4"],
                            summaries=summaries)
        spectrum = ClinvarSpectrumProvider(client, "2026-09-17").get_spectrum("TEST")
        self.assertEqual(spectrum["pathogenic_missense_count"], 2)
        self.assertEqual(spectrum["benign_missense_count"], 1)
        self.assertEqual(spectrum["pathogenic_truncating_count"], 1)
        self.assertTrue(spectrum["complete"])
        for field in ("source", "source_version", "retrieved_at", "evidence_id"):
            self.assertTrue(spectrum[field], field)

    def test_an_incomplete_missense_search_yields_no_spectrum(self):
        # ClinVar reports 10 total but the fake only ever has 1 id to hand back across
        # pages - the same shape a real search that can't deliver its own count has.
        client = FakeClient(missense_ids=["1"], truncating_ids=["2"],
                            summaries={**summary("1", "Pathogenic"), **summary("2", "Pathogenic")},
                            missense_count=10)
        self.assertIsNone(ClinvarSpectrumProvider(client, "2026-09-17").get_spectrum("TEST"))

    def test_an_incomplete_truncating_search_yields_no_spectrum(self):
        client = FakeClient(missense_ids=["1"], truncating_ids=["2"],
                            summaries={**summary("1", "Pathogenic"), **summary("2", "Pathogenic")},
                            truncating_count=10)
        self.assertIsNone(ClinvarSpectrumProvider(client, "2026-09-17").get_spectrum("TEST"))

    def test_pages_through_a_gene_with_more_hits_than_one_page(self):
        """A gene whose ClinVar hits exceed one page (previously: get_spectrum()
        returned None outright for any such gene, real examples being APC/
        BRCA1/BRCA2 - see clinvar_spectrum.py's own module comment) is now
        paged through in full, not just accepted up to a fixed ceiling."""
        missense_ids = [str(i) for i in range(1, 8)]  # 7 ids, page_size=3 -> 3 pages
        summaries = {}
        for uid in missense_ids[:5]:
            summaries.update(summary(uid, "Pathogenic"))
        for uid in missense_ids[5:]:
            summaries.update(summary(uid, "Benign"))
        client = FakeClient(missense_ids=missense_ids, truncating_ids=["100"],
                            summaries={**summaries, **summary("100", "Pathogenic")})
        spectrum = ClinvarSpectrumProvider(client, "2026-09-17", page_size=3).get_spectrum("TEST")
        self.assertIsNotNone(spectrum)
        self.assertEqual(spectrum["pathogenic_missense_count"], 5)
        self.assertEqual(spectrum["benign_missense_count"], 2)
        self.assertEqual(spectrum["pathogenic_truncating_count"], 1)
        # 3 pages for missense (3+3+1) + 1 page for truncating = 4 esearch calls.
        self.assertEqual(sum("esearch.fcgi" in url for url in client.urls), 4)

    def test_records_for_a_different_gene_are_excluded(self):
        summaries = {**summary("1", "Pathogenic", gene="TEST"),
                     **summary("2", "Pathogenic", gene="OTHER")}
        client = FakeClient(missense_ids=["1", "2"], summaries=summaries)
        spectrum = ClinvarSpectrumProvider(client, "2026-09-17").get_spectrum("TEST")
        self.assertEqual(spectrum["pathogenic_missense_count"], 1)

    def test_a_conflicting_classification_counts_as_neither(self):
        client = FakeClient(missense_ids=["1"],
                            summaries=summary("1", "Conflicting classifications of pathogenicity"))
        spectrum = ClinvarSpectrumProvider(client, "2026-09-17").get_spectrum("TEST")
        self.assertEqual(spectrum["pathogenic_missense_count"], 0)
        self.assertEqual(spectrum["benign_missense_count"], 0)

    def test_no_gene_yields_none(self):
        self.assertIsNone(ClinvarSpectrumProvider(FakeClient(), "2026-09-17").get_spectrum(""))

    def test_a_release_is_required(self):
        with self.assertRaises(ValueError):
            ClinvarSpectrumProvider(FakeClient(), "")


if __name__ == "__main__":
    unittest.main()
