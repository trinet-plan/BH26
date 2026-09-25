import unittest
from unittest.mock import AsyncMock, patch

from acmg_pipeline import pp4_literature_search


class PreferredPmidsTests(unittest.IsolatedAsyncioTestCase):
    """search_diagnostic_yield() must try this variant's own ERepo evidence_pmids
    (preferred_pmids) before/alongside a fresh keyword search - confirmed empirically
    (2026-09-25) that a live search alone surfaces broad multi-gene-panel cohort papers
    (e.g. RPGR's yield diluted across a 266-gene retinal-dystrophy panel, 4.5%) instead
    of the narrower phenotype-matched paper a real curator would cite, pushing real
    PP4-MET cases (BRCA1, MSH2, RPGR) below PP4_MIN_POSTERIOR."""

    async def test_a_preferred_pmid_is_judged_even_when_the_live_search_finds_nothing(self):
        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids", new=AsyncMock(return_value=[])),
            patch("acmg_pipeline.pipeline.fetch_full_text",
                  new=AsyncMock(return_value=("full text", "ok"))),
            patch(
                "acmg_pipeline.pp4_literature_search._judge_paper_for_yield",
                new=AsyncMock(return_value={
                    "yield_percent_stated": 42.0,
                    "is_overall_yield": True,
                    "denominator_description": "all patients tested",
                    "sample_size": 100,
                    "quote": "42 of 100 patients had a pathogenic variant",
                }),
            ),
        ):
            result = await pp4_literature_search.search_diagnostic_yield(
                "RPGR", "RPGR-related retinopathy", preferred_pmids=("30209399",))

        self.assertTrue(result.found)
        self.assertEqual(result.pmid, "30209399")
        self.assertEqual(result.yield_fraction, 0.42)

    async def test_preferred_pmids_are_tried_before_the_live_search_result(self):
        judged_pmids = []

        async def fake_judge(pipeline, gene, phenotype_description, text, pmid):
            judged_pmids.append(pmid)
            return {"no_yield_statistic_found": True}

        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids",
                  new=AsyncMock(return_value=["99999999"])),
            patch("acmg_pipeline.pipeline.fetch_full_text",
                  new=AsyncMock(return_value=("full text", "ok"))),
            patch("acmg_pipeline.pp4_literature_search._judge_paper_for_yield",
                  new=fake_judge),
        ):
            await pp4_literature_search.search_diagnostic_yield(
                "RPGR", "RPGR-related retinopathy", preferred_pmids=("30209399",), max_candidates=5)

        self.assertEqual(judged_pmids[0], "30209399")

    async def test_no_preferred_pmids_falls_back_to_the_live_search_alone(self):
        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids",
                  new=AsyncMock(return_value=["12345678"])) as mock_search,
            patch("acmg_pipeline.pipeline.fetch_full_text",
                  new=AsyncMock(return_value=(None, "not in PMC"))),
        ):
            result = await pp4_literature_search.search_diagnostic_yield(
                "RPGR", "RPGR-related retinopathy")

        self.assertFalse(result.found)
        mock_search.assert_called()

    async def test_no_phenotype_description_short_circuits_without_a_search(self):
        result = await pp4_literature_search.search_diagnostic_yield(
            "RPGR", "", preferred_pmids=("30209399",))

        self.assertFalse(result.found)
        self.assertEqual(result.reason, "no_phenotype_description_available")


class CohortAscertainmentPreferenceTests(unittest.IsolatedAsyncioTestCase):
    """Real ground-truth misses (2026-09-25: BRCA1 c.191G>A, MSH2 c.1012G>A, RPGR
    c.492G>T) all found a genuine, correctly-extracted gene-specific yield statistic
    that was diluted by a broad multigene/NGS panel denominator (e.g. RPGR's yield
    diluted across a 266-gene retinal-dystrophy panel: 4.5%). A gene_specific candidate
    must win over an earlier broad_panel one, not just whichever is found first."""

    async def test_a_later_gene_specific_hit_is_preferred_over_an_earlier_broad_panel_one(self):
        judgments = {
            "11111111": {
                "yield_percent_stated": 4.5, "is_overall_yield": True,
                "cohort_ascertainment": "broad_panel", "sample_size": 5201,
            },
            "22222222": {
                "yield_percent_stated": 72.0, "is_overall_yield": True,
                "cohort_ascertainment": "gene_specific", "sample_size": 40,
            },
        }

        async def fake_judge(pipeline, gene, phenotype_description, text, pmid):
            return judgments[pmid]

        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids",
                  new=AsyncMock(return_value=["11111111", "22222222"])),
            patch("acmg_pipeline.pipeline.fetch_full_text",
                  new=AsyncMock(return_value=("full text", "ok"))),
            patch("acmg_pipeline.pp4_literature_search._judge_paper_for_yield", new=fake_judge),
        ):
            result = await pp4_literature_search.search_diagnostic_yield(
                "RPGR", "RPGR-related retinopathy")

        self.assertTrue(result.found)
        self.assertEqual(result.pmid, "22222222")
        self.assertEqual(result.yield_fraction, 0.72)
        self.assertEqual(result.cohort_ascertainment, "gene_specific")

    async def test_a_broad_panel_hit_is_kept_as_a_fallback_when_nothing_narrower_is_found(self):
        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids",
                  new=AsyncMock(return_value=["11111111"])),
            patch("acmg_pipeline.pipeline.fetch_full_text",
                  new=AsyncMock(return_value=("full text", "ok"))),
            patch("acmg_pipeline.pp4_literature_search._judge_paper_for_yield", new=AsyncMock(return_value={
                "yield_percent_stated": 4.5, "is_overall_yield": True,
                "cohort_ascertainment": "broad_panel", "sample_size": 5201,
            })),
        ):
            result = await pp4_literature_search.search_diagnostic_yield(
                "RPGR", "RPGR-related retinopathy")

        self.assertTrue(result.found)
        self.assertEqual(result.cohort_ascertainment, "broad_panel")
        self.assertEqual(result.yield_fraction, 0.045)

    async def test_query_construction_tries_a_panel_excluding_phrasing_first(self):
        seen_disease_args = []

        async def fake_search(mcp, gene, hgvsp=None, disease=None, max_results=10):
            seen_disease_args.append(disease)
            return []

        with (
            patch("acmg_pipeline.pipeline.connect_pubmed", new=AsyncMock(return_value=object())),
            patch("acmg_pipeline.pipeline.search_candidate_pmids", new=fake_search),
        ):
            await pp4_literature_search.search_diagnostic_yield("RPGR", "RPGR-related retinopathy")

        self.assertIn("multigene panel", seen_disease_args[0])
        self.assertIn("NOT", seen_disease_args[0])


if __name__ == "__main__":
    unittest.main()
