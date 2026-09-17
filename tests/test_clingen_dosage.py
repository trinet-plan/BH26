"""ClinGen dosage as an automated LoF-mechanism signal - what it may and may not conclude.

The provider exists to unblock PVS1's G01 gate, so the tests are about the boundary it has
to hold while doing that: it is an automated signal that says so, it never decides the
recessive case, and it refuses to emit anything it cannot version.
"""

import unittest
from types import SimpleNamespace

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import reviewed_or_automated
from acmg_pipeline.providers.clingen_dosage import METHOD, ClinGenDosageProvider, parse


HEADER = (
    "#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHaploinsufficiency Score\t"
    "Haploinsufficiency Description\tTriplosensitivity Score\tDate Last Evaluated\t"
    "Haploinsufficiency Disease ID"
)


def curation_list(rows, release="16 Sep,2026"):
    lines = ["#ClinGen Gene Curation Results"]
    if release:
        lines.append(f"#{release}")
    lines.append("#Genomic Locations are reported on GRCh38")
    lines.append(HEADER)
    lines.extend(rows)
    return "\n".join(lines) + "\n"


def row(symbol, score, *, gene_id="1", description="", evaluated="2020-01-01", disease=""):
    return "\t".join([symbol, gene_id, "1p1", "chr1:1-2", score, description, "0",
                      evaluated, disease])


class FakeClient:
    def __init__(self, text):
        self.text = text
        self.calls = 0

    def fetch(self, url, *, response_format="json", **kwargs):
        self.calls += 1
        return {"body": self.text, "retrieved_at": "2026-09-17T00:00:00Z"}


class ClinGenDosageTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")

    def provider(self, rows, release="16 Sep,2026"):
        return ClinGenDosageProvider(FakeClient(curation_list(rows, release)))

    def mechanism(self, rows, gene, release="16 Sep,2026"):
        return self.provider(rows, release).get_mechanism(self.variant, gene)

    def test_sufficient_dosage_evidence_establishes_the_mechanism(self):
        records = self.mechanism([row("RUNX1", "3", disease="MONDO:0100083")], "RUNX1")
        self.assertEqual(len(records), 1)
        self.assertIs(records[0]["lof_mechanism_established"], True)
        self.assertEqual(records[0]["haploinsufficiency_score"], "3")

    def test_a_recessive_score_decides_nothing(self):
        """PAH scores 30, and the Phenylketonuria VCEP applied PVS1 to PAH c.806delT.

        Reading 30 as "not a LoF mechanism" would strike out recessive PVS1 wholesale, so
        no record is emitted and PVS1 keeps asking a curator for the mechanism.
        """
        self.assertEqual(self.mechanism([row("PAH", "30")], "PAH"), [])

    def test_scores_about_the_dominant_mechanism_are_negative(self):
        for score in ("0", "1", "2", "40"):
            with self.subTest(score=score):
                records = self.mechanism([row("GENE1", score)], "GENE1")
                self.assertEqual(len(records), 1)
                self.assertIs(records[0]["lof_mechanism_established"], False)

    def test_a_gene_absent_from_the_list_yields_nothing(self):
        self.assertEqual(self.mechanism([row("RUNX1", "3")], "RPE65"), [])
        self.assertEqual(self.mechanism([row("RUNX1", "3")], ""), [])

    def test_an_unversioned_file_is_not_usable_evidence(self):
        """Without the file's own release there is nothing to pin the policy to."""
        self.assertEqual(self.mechanism([row("RUNX1", "3")], "RUNX1", release=None), [])

    def test_the_record_declares_itself_automated_rather_than_curated(self):
        record = self.mechanism([row("RUNX1", "3")], "RUNX1")[0]
        self.assertEqual(record["assessment_method"], "automated")
        self.assertEqual(record["method"], METHOD)
        self.assertEqual(record["policy_version"], "16 Sep,2026")
        self.assertNotIn("curator", record)
        self.assertNotIn("reviewed_at", record)
        # The shape PVS1's _candidates() accepts for an automated assessment.
        self.assertTrue(reviewed_or_automated(record))

    def test_the_record_carries_the_provenance_every_evidence_item_needs(self):
        record = self.mechanism([row("RUNX1", "3")], "RUNX1")[0]
        for field in ("evidence_id", "source", "source_version", "retrieved_at",
                      "quality_status"):
            self.assertTrue(record[field], field)
        self.assertTrue(record["evidence_id"].startswith("urn:sha256:"))

    def test_the_curated_disease_scopes_the_record(self):
        """ClinGen curates haploinsufficiency against a named disease, so the record says
        which one. Scope and weight are separate: this stays an automated assessment."""
        record = self.mechanism([row("MYBPC3", "3", disease="MONDO:0005045")], "MYBPC3")[0]
        self.assertEqual(record["condition"], "MONDO:0005045")
        self.assertEqual(record["haploinsufficiency_disease_id"], "MONDO:0005045")
        self.assertEqual(record["assessment_method"], "automated")
        self.assertNotIn("curator", record)

    def test_a_row_naming_no_disease_stays_unscoped(self):
        """The score still says something about the gene, so the record is still emitted -
        it simply cannot settle a mechanism for a disease the list does not name."""
        record = self.mechanism([row("MYBPC3", "3")], "MYBPC3")[0]
        self.assertNotIn("condition", record)
        self.assertIsNone(record["haploinsufficiency_disease_id"])
        self.assertIs(record["lof_mechanism_established"], True)

    def test_two_variants_in_one_gene_get_distinct_evidence_ids(self):
        """automated_cli deduplicates evidence by evidence_id, so a gene-only id collapsed
        every variant in a gene down to one and left the rest without a mechanism."""
        provider = self.provider([row("MYBPC3", "3")])
        other = Variant("GRCh38", "11", 47335041, "C", "T")
        first = provider.get_mechanism(self.variant, "MYBPC3")[0]
        second = provider.get_mechanism(other, "MYBPC3")[0]
        self.assertNotEqual(first["evidence_id"], second["evidence_id"])
        self.assertIn(self.variant.key, first["evidence_id"])
        self.assertIn(other.key, second["evidence_id"])

    def test_parse_reports_a_file_missing_its_columns(self):
        for text in ("no header at all\n", "#Gene Symbol\tGene ID\nX\t1\n"):
            with self.subTest(text=text):
                with self.assertRaises(ValueError):
                    parse(text)

    def test_a_short_row_does_not_break_the_parse(self):
        genes, release = parse(curation_list([row("RUNX1", "3"), "TRUNCATED\t2"]))
        self.assertIn("RUNX1", genes)
        self.assertNotIn("TRUNCATED", genes)
        self.assertEqual(release, "16 Sep,2026")


if __name__ == "__main__":
    unittest.main()
