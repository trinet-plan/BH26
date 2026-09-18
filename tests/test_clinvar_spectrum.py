"""Gene-wide ClinVar missense/truncating classification counts for BP1's
statistical suggestion (see acmg_pipeline.providers.gene_disease_draft)."""

import unittest

from acmg_pipeline.providers.clinvar_spectrum import ClinvarSpectrumProvider


def summary(uid, classification, gene="TEST"):
    return {uid: {"accession": f"VCV{int(uid):09d}", "gene_sort": gene,
                 "germline_classification": {"description": classification}}}


class FakeClient:
    def __init__(self, *, missense_ids=(), truncating_ids=(), summaries=None,
                missense_complete=True, truncating_complete=True):
        self.missense_ids = list(missense_ids)
        self.truncating_ids = list(truncating_ids)
        self.summaries = summaries or {}
        self.missense_complete = missense_complete
        self.truncating_complete = truncating_complete
        self.urls = []

    def fetch(self, url, **kwargs):
        self.urls.append(url)
        if "esearch.fcgi" in url:
            truncating = "nonsense" in url
            ids = self.truncating_ids if truncating else self.missense_ids
            complete = self.truncating_complete if truncating else self.missense_complete
            count = len(ids) if complete else len(ids) + 10_000
            return {"body": {"esearchresult": {"count": str(count), "retmax": "5000",
                                               "idlist": list(ids)}},
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
        client = FakeClient(missense_ids=["1"], truncating_ids=["2"],
                            summaries={**summary("1", "Pathogenic"), **summary("2", "Pathogenic")},
                            missense_complete=False)
        self.assertIsNone(ClinvarSpectrumProvider(client, "2026-09-17").get_spectrum("TEST"))

    def test_an_incomplete_truncating_search_yields_no_spectrum(self):
        client = FakeClient(missense_ids=["1"], truncating_ids=["2"],
                            summaries={**summary("1", "Pathogenic"), **summary("2", "Pathogenic")},
                            truncating_complete=False)
        self.assertIsNone(ClinvarSpectrumProvider(client, "2026-09-17").get_spectrum("TEST"))

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
