import unittest
from unittest.mock import AsyncMock, patch

from acmg_pipeline import pp1_segregation_search


class PreferredPmidsTests(unittest.IsolatedAsyncioTestCase):
    """search_family_segregation() must try this variant's own ERepo evidence_pmids
    (preferred_pmids) before/alongside a fresh keyword search - confirmed empirically
    (2026-09-24) that a live search alone returns papers with zero overlap against real
    ERepo citations for the same variant (MYH7 c.1594T>C, MSH2 c.1012G>A, RUNX1
    c.601C>T all had 0/N overlap)."""

    async def test_a_preferred_pmid_is_judged_even_when_the_live_search_finds_nothing(self):
        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids", new=AsyncMock(return_value=[])),
            patch("acmg_pipeline.pipeline.fetch_full_text",
                  new=AsyncMock(return_value=("full text", "ok"))),
            patch(
                "acmg_pipeline.pp1_segregation_search._judge_paper_for_segregation",
                new=AsyncMock(return_value={
                    "family_segregation_described": True,
                    "inheritance_pattern": "autosomal dominant",
                    "relatives": [{"relationship": "father", "affected_status": True,
                                   "variant_status": True}],
                }),
            ),
        ):
            result = await pp1_segregation_search.search_family_segregation(
                "MYH7", "c.1594T>C", preferred_pmids=("11106718",))

        self.assertTrue(result.found)
        self.assertEqual(result.pmid, "11106718")

    async def test_preferred_pmids_are_tried_before_the_live_search_result(self):
        judged_pmids = []

        async def fake_judge(pipeline, gene, hgvsc, text, pmid):
            judged_pmids.append(pmid)
            return {"no_segregation_data_found": True}

        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids",
                  new=AsyncMock(return_value=["99999999"])),
            patch("acmg_pipeline.pipeline.fetch_full_text",
                  new=AsyncMock(return_value=("full text", "ok"))),
            patch("acmg_pipeline.pp1_segregation_search._judge_paper_for_segregation",
                  new=fake_judge),
        ):
            await pp1_segregation_search.search_family_segregation(
                "MYH7", "c.1594T>C", preferred_pmids=("11106718",), max_candidates=5)

        self.assertEqual(judged_pmids[0], "11106718")

    async def test_no_preferred_pmids_falls_back_to_the_live_search_alone(self):
        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids",
                  new=AsyncMock(return_value=["12345678"])) as mock_search,
            patch("acmg_pipeline.pipeline.fetch_full_text",
                  new=AsyncMock(return_value=(None, "not in PMC"))),
        ):
            result = await pp1_segregation_search.search_family_segregation("MYH7", "c.1594T>C")

        self.assertFalse(result.found)
        mock_search.assert_called()

    async def test_query_never_embeds_the_raw_hgvsc_notation(self):
        """Regression guard: search_candidate_pmids's own docstring says exact c./p.
        notation in a query returns ZERO PubMed results - the disease= text passed
        here must never contain the raw hgvsc string."""
        seen_disease_args = []

        async def fake_search(mcp, gene, hgvsp=None, disease=None, max_results=10):
            seen_disease_args.append(disease)
            return []

        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids", new=fake_search),
        ):
            await pp1_segregation_search.search_family_segregation("MYH7", "c.1594T>C")

        for disease in seen_disease_args:
            self.assertNotIn("c.1594T>C", disease or "")


if __name__ == "__main__":
    unittest.main()
