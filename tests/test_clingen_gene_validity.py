"""ClinGen Gene-Disease Validity bulk CSV parsing and per-gene filtering."""

import unittest

from acmg_pipeline.providers.clingen_gene_validity import ClinGenGeneValidityProvider


def csv_text(rows, created="2026-09-17"):
    header = ('"CLINGEN GENE DISEASE VALIDITY CURATIONS","","","","","","","","",""\n'
              f'"FILE CREATED: {created}","","","","","","","","",""\n'
              '"WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity","","","","","","","","",""\n'
              '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++",'
              '"+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","+++++++++++++++++++"\n'
              '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP",'
              '"CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n')
    body = "".join(
        (f'"{gene}","{hgnc}","{label}","{mondo}","{moi}","SOP8","{classification}",'
         f'"https://search.clinicalgenome.org/kb/gene-validity/{report}","{date}","{gcep}"\n')
        for gene, hgnc, label, mondo, moi, classification, report, date, gcep in rows
    )
    return header + body


class FakeClient:
    def __init__(self, text):
        self.text = text
        self.calls = 0

    def fetch(self, url, *, response_format="json", **kwargs):
        self.calls += 1
        return {"body": self.text, "retrieved_at": "2026-09-17T00:00:00Z"}


class ClinGenGeneValidityTests(unittest.TestCase):
    def provider(self, rows, created="2026-09-17"):
        return ClinGenGeneValidityProvider(FakeClient(csv_text(rows, created)))

    def test_matches_are_filtered_by_gene_symbol(self):
        provider = self.provider([
            ("MYH7", "HGNC:7577", "hypertrophic cardiomyopathy", "MONDO:0005045",
             "AD", "Definitive", "CGGV:1", "2021-10-07T16:00:00.000Z", "Test GCEP"),
            ("MYBPC3", "HGNC:7551", "hypertrophic cardiomyopathy", "MONDO:0005045",
             "AD", "Definitive", "CGGV:2", "2021-10-07T16:00:00.000Z", "Test GCEP"),
        ])
        rows = provider.get_validity("MYH7")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["gene"], "MYH7")
        self.assertEqual(rows[0]["condition"], "MONDO:0005045")
        self.assertEqual(rows[0]["classification"], "Definitive")
        self.assertEqual(rows[0]["moi"], "AD")

    def test_an_unknown_gene_returns_no_rows(self):
        provider = self.provider([
            ("MYH7", "HGNC:7577", "hypertrophic cardiomyopathy", "MONDO:0005045",
             "AD", "Definitive", "CGGV:1", "2021-10-07T16:00:00.000Z", "Test GCEP"),
        ])
        self.assertEqual(provider.get_validity("NOT_A_GENE"), [])

    def test_a_gene_with_multiple_diseases_returns_every_row(self):
        provider = self.provider([
            ("MYH7", "HGNC:7577", "hypertrophic cardiomyopathy", "MONDO:0005045",
             "AD", "Definitive", "CGGV:1", "2021-10-07T16:00:00.000Z", "Test GCEP"),
            ("MYH7", "HGNC:7577", "dilated cardiomyopathy", "MONDO:0013262",
             "AD", "Limited", "CGGV:2", "2020-01-01T16:00:00.000Z", "Test GCEP"),
        ])
        self.assertEqual(len(provider.get_validity("MYH7")), 2)

    def test_every_row_carries_full_provenance(self):
        provider = self.provider([
            ("MYH7", "HGNC:7577", "hypertrophic cardiomyopathy", "MONDO:0005045",
             "AD", "Definitive", "CGGV:1", "2021-10-07T16:00:00.000Z", "Test GCEP"),
        ])
        row = provider.get_validity("MYH7")[0]
        for field in ("source", "source_version", "retrieved_at", "evidence_id"):
            self.assertTrue(row[field], field)
        self.assertEqual(row["source_version"], "2026-09-17")

    def test_the_snapshot_is_fetched_only_once_across_repeated_lookups(self):
        client = FakeClient(csv_text([
            ("MYH7", "HGNC:7577", "hypertrophic cardiomyopathy", "MONDO:0005045",
             "AD", "Definitive", "CGGV:1", "2021-10-07T16:00:00.000Z", "Test GCEP"),
        ]))
        provider = ClinGenGeneValidityProvider(client)
        provider.get_validity("MYH7")
        provider.get_validity("MYBPC3")
        self.assertEqual(client.calls, 1)

    def test_a_response_missing_the_file_created_header_is_rejected(self):
        client = FakeClient('"GENE SYMBOL","GENE ID (HGNC)"\n"MYH7","HGNC:7577"\n')
        with self.assertRaises(ValueError):
            ClinGenGeneValidityProvider(client).get_validity("MYH7")


if __name__ == "__main__":
    unittest.main()
