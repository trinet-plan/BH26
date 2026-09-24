import json
import unittest
from unittest.mock import AsyncMock, MagicMock, patch

import acmg_pipeline.pipeline as pl


def _mock_result(pmids):
    content = MagicMock()
    content.text = json.dumps({"pmids": pmids})
    result = MagicMock()
    result.content = [content]
    return result


class DualNotationQueryTests(unittest.IsolatedAsyncioTestCase):
    """search_candidate_pmids() must try both the 3-letter (standard HGVS) and
    1-letter amino-acid forms of a protein change - confirmed empirically
    (2026-09-24) that neither convention reliably wins across real variants
    (MYH7 Ser532Pro/S532P found 0/1, PTEN Pro38Ser/P38S found 1/0, RUNX1
    Arg201Ter/R201* found 0/1, BRCA1 Cys61Gly/C61G found 17/68)."""

    async def test_both_amino_acid_notations_are_queried(self):
        queries_seen = []

        async def fake_call_tool_safe(mcp, name, arguments):
            queries_seen.append(arguments["query"])
            return _mock_result([])

        with patch("acmg_pipeline.pipeline.call_tool_safe", new=fake_call_tool_safe):
            await pl.search_candidate_pmids(object(), "MYH7", hgvsp="p.Ser532Pro")

        self.assertIn("MYH7 AND Ser532Pro", queries_seen)
        self.assertIn("MYH7 AND S532P", queries_seen)

    async def test_a_stop_codon_change_gets_the_asterisk_one_letter_form(self):
        queries_seen = []

        async def fake_call_tool_safe(mcp, name, arguments):
            queries_seen.append(arguments["query"])
            return _mock_result([])

        with patch("acmg_pipeline.pipeline.call_tool_safe", new=fake_call_tool_safe):
            await pl.search_candidate_pmids(object(), "RUNX1", hgvsp="p.Arg201Ter")

        self.assertIn("RUNX1 AND Arg201Ter", queries_seen)
        self.assertIn("RUNX1 AND R201*", queries_seen)

    async def test_no_hgvsp_skips_the_amino_acid_queries_entirely(self):
        queries_seen = []

        async def fake_call_tool_safe(mcp, name, arguments):
            queries_seen.append(arguments["query"])
            return _mock_result([])

        with patch("acmg_pipeline.pipeline.call_tool_safe", new=fake_call_tool_safe):
            await pl.search_candidate_pmids(object(), "MYH7", hgvsp="N/A", disease="cardiomyopathy")

        self.assertTrue(all("AND" in q and "cardiomyopathy" in q or "novel variant" in q
                             for q in queries_seen))


if __name__ == "__main__":
    unittest.main()
