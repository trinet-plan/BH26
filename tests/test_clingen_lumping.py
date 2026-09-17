"""ClinGen's lumping decisions read strictly, and refused rather than read as empty."""

import unittest

from acmg_pipeline.providers.clingen_lumping import (
    METHOD, ClinGenLumpingProvider, parse_index, parse_report,
)


INDEX_HEADER = ('"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI",'
                '"SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"')
REPORT_URL = "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_1"


def index(rows=None):
    default = ['"AARS1","HGNC:20","CMT axonal type 2N","MONDO:0013212","AD","SOP10",'
               f'"Definitive","{REPORT_URL}","2024-03-14T16:00:00.000Z","CMT GCEP"']
    return "\n".join([
        '"CLINGEN GENE DISEASE VALIDITY CURATIONS","","","","","","","","",""',
        '"FILE CREATED: 2026-09-17","","","","","","","","",""',
        INDEX_HEADER,
        '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++",'
        '"+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","++++++++++"',
        *(default if rows is None else rows),
    ]) + "\n"


def report(included=("613287",), excluded=("616339", "619691"), date="07/11/2023",
           headings=("Included MIM Phenotypes", "Excluded MIM Phenotypes", "Evaluation Date"),
           pane='id="gdvt5"'):
    included_at, excluded_at, date_at = headings
    body = [f"<div {pane}>", "<h3>Lumping &amp; Splitting</h3>", f"<dt>{included_at}</dt>"]
    body += [f"<li>MIM:{item} - a phenotype</li>" for item in included]
    body += [f"<dt>{excluded_at}</dt>"]
    body += [f"<li>MIM:{item} - another phenotype</li>" for item in excluded]
    body += [f"<dt>{date_at}</dt><dd>{date}</dd>", "</div>"]
    return "<html><body><nav>Lumping and Splitting</nav>" + "".join(body) + "</body></html>"


class FakeResolver:
    """Stands in for the MONDO mapping set; only the terms it knows resolve."""

    MAP = {"OMIM:613287": "MONDO:0013212", "OMIM:616339": "MONDO:0014593",
           "OMIM:619691": "MONDO:0030517"}

    def normalize(self, condition):
        term = self.MAP.get(condition)
        return {"normalized_condition": term} if term else None


class FakeClient:
    def __init__(self, index_text, pages):
        self.index_text = index_text
        self.pages = pages
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        body = self.index_text if url.endswith("/download") else self.pages[url]
        return {"body": body, "retrieved_at": "2026-09-17T00:00:00Z"}


class ParseIndexTests(unittest.TestCase):
    def test_the_separator_row_is_not_a_curation(self):
        entries = parse_index(index())
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]["gene"], "AARS1")
        self.assertEqual(entries[0]["condition"], "MONDO:0013212")
        self.assertEqual(entries[0]["classification"], "Definitive")

    def test_a_curation_without_a_mondo_disease_is_skipped(self):
        entries = parse_index(index(rows=[
            '"AARS1","HGNC:20","label","","AD","SOP10","Definitive","url","2024","GCEP"']))
        self.assertEqual(entries, [])

    def test_a_file_without_the_header_is_an_error(self):
        with self.assertRaises(ValueError):
            parse_index('"CLINGEN GENE DISEASE VALIDITY CURATIONS","",""\n')


class ParseReportTests(unittest.TestCase):
    def test_the_two_phenotype_lists_are_kept_apart(self):
        parsed = parse_report(report())
        self.assertEqual(parsed["included"], ["OMIM:613287"])
        self.assertEqual(parsed["excluded"], ["OMIM:616339", "OMIM:619691"])
        self.assertEqual(parsed["evaluated_at"], "07/11/2023")

    def test_navigation_text_outside_the_pane_is_not_the_report(self):
        """The site's own menu says "Lumping and Splitting" on every page."""
        self.assertIsNone(parse_report("<nav>Lumping and Splitting</nav>"))

    def test_every_landmark_is_required(self):
        cases = {
            "no pane": report(pane='id="other"'),
            "no included heading": report(headings=("Lumped In", "Excluded MIM Phenotypes",
                                                    "Evaluation Date")),
            "no excluded heading": report(headings=("Included MIM Phenotypes", "Kept Out",
                                                    "Evaluation Date")),
            "no date": report(headings=("Included MIM Phenotypes", "Excluded MIM Phenotypes",
                                        "Reviewed On")),
            "no included phenotype": report(included=()),
        }
        for label, markup in cases.items():
            with self.subTest(label):
                self.assertIsNone(parse_report(markup))

    def test_an_empty_excluded_list_is_a_real_answer(self):
        """Nothing kept out is a decision a panel can make; nothing lumped in is not."""
        parsed = parse_report(report(excluded=()))
        self.assertEqual(parsed["excluded"], [])
        self.assertEqual(parsed["included"], ["OMIM:613287"])


class GetScopeTests(unittest.TestCase):
    def provider(self, index_text=None, page=None):
        client = FakeClient(index_text or index(), {REPORT_URL: page or report()})
        return ClinGenLumpingProvider(client, FakeResolver()), client

    def test_the_scope_comes_back_in_mondo_terms(self):
        provider, _ = self.provider()
        scope = provider.get_scope("AARS1", "MONDO:0013212")
        self.assertEqual(scope["included"], ["MONDO:0013212"])
        self.assertEqual(scope["excluded"], ["MONDO:0014593", "MONDO:0030517"])
        self.assertEqual(scope["source_version"], "07/11/2023")
        self.assertEqual(scope["method"], METHOD)
        self.assertEqual(scope["gene_disease_validity"], "Definitive")
        self.assertEqual(scope["unresolved_phenotypes"], [])

    def test_the_curated_disease_is_always_in_its_own_scope(self):
        """Even where its MIM entry does not resolve, the curation is about that disease."""
        provider, _ = self.provider(page=report(included=("999999",)))
        scope = provider.get_scope("AARS1", "MONDO:0013212")
        self.assertIn("MONDO:0013212", scope["included"])
        self.assertEqual(scope["unresolved_phenotypes"], ["OMIM:999999"])

    def test_an_unreadable_page_yields_no_scope_rather_than_an_empty_one(self):
        provider, _ = self.provider(page="<html><body>redesigned</body></html>")
        self.assertIsNone(provider.get_scope("AARS1", "MONDO:0013212"))

    def test_a_gene_disease_pair_with_no_curation_yields_nothing(self):
        provider, client = self.provider()
        self.assertIsNone(provider.get_scope("AARS1", "MONDO:9999999"))
        self.assertIsNone(provider.get_scope("BRCA1", "MONDO:0013212"))
        self.assertEqual(client.urls, [client.urls[0]])  # only the index was fetched

    def test_the_index_is_fetched_once_for_many_lookups(self):
        provider, client = self.provider()
        provider.get_scope("AARS1", "MONDO:0013212")
        provider.get_scope("AARS1", "MONDO:0013212")
        self.assertEqual(sum(1 for url in client.urls if url.endswith("/download")), 1)


if __name__ == "__main__":
    unittest.main()
