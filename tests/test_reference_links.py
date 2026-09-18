"""autopvs1_variant_url() - AutoPVS1's real positional deep link.

Added 2026-09-18: a live browser check (curator-suggested) found AutoPVS1
answers /variant/{build}/{chrom}-{pos}-{ref}-{alt} with the actual per-variant
PVS1 flowchart (confirmed against MYH9 hg19 22-36678800-G-A and this
project's own MYBPC3 c.278delA, hg38 11-47351252-CT-C, matching that
variant's own PVS1_very_strong result) - contradicting this module's earlier
"fixed URL, not deep-linkable" conclusion, which had never actually been
verified against this exact URL shape.
"""

import unittest
from urllib.parse import unquote

from acmg_pipeline.criteria.reference_links import (
    AUTOPVS1_URL, autopvs1_variant_url, clinvar_position_url,
    franklin_variant_url, reference_urls_for_criterion,
)
from acmg_pipeline.vcf_record import VariantRecord


def variant(**overrides):
    fields = dict(chrom="11", pos=47351252, id="v1", ref="CT", alt="C",
                  qual="", filter="", info={})
    fields.update(overrides)
    return VariantRecord(**fields)


class AutoPvs1VariantUrlTests(unittest.TestCase):
    def test_a_concrete_allele_builds_the_real_deep_link(self):
        self.assertEqual(
            autopvs1_variant_url(variant()),
            "https://autopvs1.bgi.com/variant/hg38/11-47351252-CT-C")

    def test_build_defaults_to_hg38_but_hg19_can_be_requested(self):
        self.assertEqual(
            autopvs1_variant_url(variant(chrom="22", pos=36678800, ref="G", alt="A"), build="hg19"),
            "https://autopvs1.bgi.com/variant/hg19/22-36678800-G-A")

    def test_a_placeholder_allele_falls_back_to_the_homepage_not_none(self):
        for alt in (".", ""):
            with self.subTest(alt=alt):
                self.assertEqual(autopvs1_variant_url(variant(alt=alt)), AUTOPVS1_URL)


class ClinvarPositionUrlTests(unittest.TestCase):
    """clinvar_position_url() - PS1's genomic-position window, confirmed live
    2026-09-18 (curator-suggested) against this project's own MYBPC3
    c.278delA (11:47351252): a +-10bp window returns 11 real ClinVar
    variants including a nearby pathogenic frameshift."""

    def test_default_window_is_ten_base_pairs_either_side(self):
        url = clinvar_position_url(VariantRecord(
            chrom="11", pos=47351252, id="v1", ref="CT", alt="C",
            qual="", filter="", info={}))
        term = unquote(url.split("term=")[1])
        self.assertEqual(term, "11[chr] AND 47351242:47351262[chrpos38]")

    def test_window_is_configurable(self):
        url = clinvar_position_url(VariantRecord(
            chrom="11", pos=47351252, id="v1", ref="CT", alt="C",
            qual="", filter="", info={}), window=3)
        term = unquote(url.split("term=")[1])
        self.assertEqual(term, "11[chr] AND 47351249:47351255[chrpos38]")

    def test_missing_position_returns_none(self):
        self.assertIsNone(clinvar_position_url(VariantRecord(
            chrom="11", pos=0, id="v1", ref="CT", alt="C",
            qual="", filter="", info={})))


class MultipleReferenceUrlTests(unittest.TestCase):
    def test_franklin_builds_a_grch38_small_variant_link(self):
        self.assertEqual(
            franklin_variant_url(variant(chrom="chr14", pos=23425971, ref="G", alt="A")),
            "https://franklin.genoox.com/clinical-db/variant/snp/"
            "chr14-23425971-G-A-hg38",
        )

    def test_population_criteria_include_togovar_gnomad_and_franklin(self):
        urls = reference_urls_for_criterion(
            "PM2", variant(chrom="14", pos=23425971, ref="G", alt="A")
        )
        self.assertEqual(urls, [
            "https://grch38.togovar.org/variant/14-23425971-G-A",
            "https://gnomad.broadinstitute.org/variant/"
            "14-23425971-G-A?dataset=gnomad_r4",
            "https://franklin.genoox.com/clinical-db/variant/snp/"
            "chr14-23425971-G-A-hg38",
        ])

    def test_a_criterion_without_a_specific_page_still_gets_franklin(self):
        self.assertEqual(
            reference_urls_for_criterion("PP3", variant()),
            [
                "https://franklin.genoox.com/clinical-db/variant/snp/"
                "chr11-47351252-CT-C-hg38"
            ],
        )

    def test_placeholder_allele_does_not_create_variant_specific_links(self):
        self.assertEqual(reference_urls_for_criterion("PM2", variant(alt=".")), [])


if __name__ == "__main__":
    unittest.main()
