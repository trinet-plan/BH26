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
    AUTOPVS1_URL, CLINVAR_URL, FRANKLIN_URL, GNOMAD_URL,
    PP3_BP4_CALIBRATION_URL, TOGOVAR_URL, UNIPROT_URL, autopvs1_variant_url,
    clinvar_codon_url, clinvar_position_url, franklin_variant_url,
    positional_allele, protein_residue, reference_urls_for_criterion,
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

    def test_missing_position_falls_back_to_the_unqueried_page(self):
        """A window around position 0 would be a search for the wrong place."""
        self.assertEqual(clinvar_position_url(VariantRecord(
            chrom="11", pos=0, id="v1", ref="CT", alt="C",
            qual="", filter="", info={})), CLINVAR_URL)


class MultipleReferenceUrlTests(unittest.TestCase):
    def test_franklin_builds_a_grch38_small_variant_link(self):
        self.assertEqual(
            franklin_variant_url(variant(chrom="chr14", pos=23425971, ref="G", alt="A")),
            "https://franklin.genoox.com/clinical-db/variant/snp/"
            "chr14-23425971-G-A-hg38",
        )

    def test_population_criteria_include_togovar_and_gnomad(self):
        """The document asks for gnomAD; TogoVar is the provider actually
        queried, so both are offered - see _URL_BUILDERS' own note."""
        urls = reference_urls_for_criterion(
            "PM2", variant(chrom="14", pos=23425971, ref="G", alt="A")
        )
        self.assertEqual(urls, [
            "https://grch38.togovar.org/variant/14-23425971-G-A",
            "https://gnomad.broadinstitute.org/variant/"
            "14-23425971-G-A?dataset=gnomad_r4",
        ])

    def test_pm1_and_pm5_lead_with_franklin_then_uniprot(self):
        """Franklin first: the document's own PM1/PM5 screenshots are its
        Region Viewer and its PM1 verdict."""
        urls = reference_urls_for_criterion(
            "PM1", variant(chrom="14", pos=23425971, ref="G", alt="A"))
        self.assertEqual(urls[0], "https://franklin.genoox.com/clinical-db/"
                                  "variant/snp/chr14-23425971-G-A-hg38")

    def test_pp3_and_bp4_point_at_the_calibration_the_pipeline_applies(self):
        for code in ("PP3", "BP4"):
            with self.subTest(code=code):
                self.assertEqual(reference_urls_for_criterion(code, variant()),
                                 [PP3_BP4_CALIBRATION_URL])

    def test_a_criterion_the_document_does_not_name_gets_no_link(self):
        """PS2/PS3/PM4/BP2/... are absent from doc/recs for expert board.docx.
        Franklin used to be appended to every code; it no longer is."""
        for code in ("PS2", "PS3", "PS4", "PM4", "PM6", "PP4", "PP5",
                     "BS3", "BS4", "BP1", "BP2", "BP3", "BP5", "BP6", "BP7"):
            with self.subTest(code=code):
                self.assertEqual(reference_urls_for_criterion(code, variant()), [])

    def test_placeholder_allele_falls_back_to_the_unqueried_pages(self):
        """Not an empty list: the curator still gets the right sites, just
        without this variant encoded into the URL."""
        self.assertEqual(reference_urls_for_criterion("PM2", variant(alt=".")),
                         [TOGOVAR_URL, GNOMAD_URL])


class ClinvarCodonUrlTests(unittest.TestCase):
    """PS1 asks for "the clinvar page/summary for that codon"."""

    def test_residue_is_read_from_hgvsp(self):
        self.assertEqual(
            protein_residue(variant(info={"HGVSP": "p.(Arg719Trp)"})), "Arg719")

    def test_residue_without_parentheses_is_read_too(self):
        self.assertEqual(
            protein_residue(variant(info={"HGVSP": "p.Arg719Trp"})), "Arg719")

    def test_unparseable_hgvsp_yields_no_residue(self):
        self.assertIsNone(protein_residue(variant(info={"HGVSP": "p.?"})))

    def test_the_search_is_by_gene_and_residue(self):
        url = clinvar_codon_url(variant(
            info={"GENE": "MYH7", "HGVSP": "p.(Arg719Trp)"}))
        self.assertEqual(unquote(url.split("term=")[1]), "MYH7[gene] AND Arg719")

    def test_without_a_residue_it_falls_back_to_the_genomic_window(self):
        url = clinvar_codon_url(variant(info={"GENE": "MYH7"}))
        self.assertEqual(unquote(url.split("term=")[1]),
                         "11[chr] AND 47351242:47351262[chrpos38]")


if __name__ == "__main__":
    unittest.main()


class PositionalAlleleGuardTests(unittest.TestCase):
    """A URL that resolves to nothing must not be emitted at all.

    None of these sites report a bad allele: checked live 2026-09-18, gnomAD,
    Franklin and AutoPVS1 all answer HTTP 200 (gnomAD byte-identically to a
    real variant) and TogoVar answers 403 for real and malformed alike - so
    the record's shape is the only thing that can be checked.
    """

    def concrete(self, **overrides):
        fields = dict(chrom="14", pos=23425971, ref="G", alt="A")
        fields.update(overrides)
        return positional_allele(variant(**fields))

    def test_a_concrete_allele_passes(self):
        self.assertEqual(self.concrete(), ("14", 23425971, "G", "A"))

    def test_multi_alt_is_rejected(self):
        self.assertIsNone(self.concrete(alt="A,T"))

    def test_symbolic_alt_is_rejected(self):
        self.assertIsNone(self.concrete(alt="<DEL>"))

    def test_upstream_deletion_alt_is_rejected(self):
        self.assertIsNone(self.concrete(alt="*"))

    def test_missing_alt_is_rejected(self):
        self.assertIsNone(self.concrete(alt="."))

    def test_position_zero_is_rejected(self):
        self.assertIsNone(self.concrete(pos=0))

    def test_a_chr_prefix_is_normalized_away(self):
        """TogoVar, gnomAD and AutoPVS1 all want a bare CHROM."""
        self.assertEqual(self.concrete(chrom="chr14"), ("14", 23425971, "G", "A"))

    def test_lowercase_bases_are_upper_cased(self):
        self.assertEqual(self.concrete(ref="g", alt="a"), ("14", 23425971, "G", "A"))

    def test_a_non_grch38_assembly_is_rejected(self):
        """Every positional builder here is GRCh38-only; a GRCh37 record used to
        be linked to whatever sits at those coordinates in GRCh38."""
        self.assertIsNone(self.concrete(info={"ASSEMBLY": "GRCh37"}))

    def test_a_declared_grch38_assembly_passes(self):
        self.assertEqual(self.concrete(info={"ASSEMBLY": "GRCh38"}),
                         ("14", 23425971, "G", "A"))

    def test_an_absent_assembly_still_means_grch38(self):
        self.assertEqual(self.concrete(info={}), ("14", 23425971, "G", "A"))

    def test_a_rejected_record_gets_landing_pages_not_variant_urls(self):
        """No URL may encode the variant wrongly; none may be missing either."""
        rejected = variant(chrom="14", pos=23425971, ref="G", alt="<DEL>",
                           info={"ASSEMBLY": "GRCh37"})
        expected = {
            "PM2": [TOGOVAR_URL, GNOMAD_URL],
            "BA1": [TOGOVAR_URL, GNOMAD_URL],
            "BS1": [TOGOVAR_URL, GNOMAD_URL],
            "BS2": [TOGOVAR_URL, GNOMAD_URL],
            "PM1": [FRANKLIN_URL, UNIPROT_URL],
            "PM5": [FRANKLIN_URL, UNIPROT_URL],
            "PS1": [CLINVAR_URL],
            "PVS1": [AUTOPVS1_URL],
        }
        for code, urls in expected.items():
            with self.subTest(code=code):
                self.assertEqual(reference_urls_for_criterion(code, rejected), urls)

    def test_every_documented_criterion_keeps_a_link_for_an_empty_record(self):
        empty = VariantRecord(chrom="14", pos=23425971, id="v1", ref=".",
                              alt=".", qual="", filter="", info={})
        for code in ("PVS1", "PS1", "PM1", "PM2", "PM3", "PM5", "PP1", "PP2",
                     "PP3", "BA1", "BS1", "BS2", "BP4"):
            with self.subTest(code=code):
                self.assertTrue(reference_urls_for_criterion(code, empty))

    def test_pvs1_still_offers_the_autopvs1_homepage(self):
        """PVS1 has always carried some link; the homepage is not broken, it
        just needs the curator to paste the variant in."""
        self.assertEqual(
            reference_urls_for_criterion("PVS1", variant(alt="<DEL>")),
            [AUTOPVS1_URL],
        )
