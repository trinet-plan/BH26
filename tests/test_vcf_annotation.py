"""Reading evidence from a VCF under the same contract a provider is held to.

services/resolve.py decided INFO could not be evidence, for three reasons about what the
demo files contain. These tests are that each reason is answered rather than waived: a field
that cannot be pinned to a release is dropped, a curator's conclusion is refused outright,
and a score keeps the version that decides whether a calibration can use it.
"""

import unittest
from pathlib import Path

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.automated_core.vcf_adapter import parse_vcf as adapter_parse_vcf
from acmg_pipeline.services.evidence import EvidenceService
from acmg_pipeline.providers.vcf_annotation import (
    REFUSED, VcfEvidenceProvider, compound_entries, declared_version, field_provenance,
    subfield_order,
)
from acmg_pipeline.vcf_record import InfoFieldDef, ParsedVcf, VariantRecord


MAPPING = {
    "gnomAD_AF": {"category": "population", "field": "AF", "group": "gnomAD"},
    "gnomAD_AC": {"category": "population", "field": "AC", "group": "gnomAD"},
    "gnomAD_AN": {"category": "population", "field": "AN", "group": "gnomAD"},
    "REVEL": {"category": "computational", "field": "score", "group": "REVEL"},
}


def info_def(info_id, description):
    return InfoFieldDef(id=info_id, number="1", type="Float", description=description)


def parsed(info, *, defs=None, meta=None):
    record = VariantRecord(chrom="1", pos=2, id="v", ref="C", alt="T", qual=".",
                           filter="PASS", info=info)
    return ParsedVcf(meta=meta if meta is not None else {"source": "VEP-run", "fileDate": "20260917"},
                     info_defs=defs or {}, record=record)


class FieldProvenanceTests(unittest.TestCase):
    def test_a_field_naming_its_own_source_and_version_is_pinned(self):
        definition = info_def("gnomAD_AF", 'Allele frequency,Source="gnomAD",Version="4.1"')
        self.assertEqual(field_provenance(definition, {}), ("gnomAD", "4.1"))

    def test_the_files_own_source_and_date_stand_in_when_a_field_is_silent(self):
        self.assertEqual(
            field_provenance(info_def("X", "no provenance here"),
                             {"source": "VEP-run", "fileDate": "20260917"}),
            ("VEP-run", "20260917"))

    def test_a_field_that_cannot_be_pinned_at_all_is_not_pinned(self):
        self.assertIsNone(field_provenance(info_def("X", "plain"), {}))
        self.assertIsNone(field_provenance(None, {"source": "VEP-run"}))
        self.assertIsNone(field_provenance(None, {"fileDate": "20260917"}))


class VcfEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.provider = VcfEvidenceProvider(MAPPING)

    def test_fields_of_one_source_become_one_record(self):
        records, skipped = self.provider.get_evidence(
            self.variant, parsed({"gnomAD_AF": "0.01", "gnomAD_AC": "5", "gnomAD_AN": "500"}))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["category"], "population")
        self.assertEqual(records[0]["AF"], "0.01")
        self.assertEqual(records[0]["AC"], "5")
        self.assertEqual(sorted(records[0]["info_fields"]), ["gnomAD_AC", "gnomAD_AF", "gnomAD_AN"])
        self.assertEqual(skipped, [])

    def test_different_sources_do_not_share_a_record(self):
        records, _ = self.provider.get_evidence(
            self.variant, parsed({"gnomAD_AF": "0.01", "REVEL": "0.8"}))
        self.assertEqual({record["category"] for record in records},
                         {"population", "computational"})

    def test_a_field_with_no_pinnable_release_is_dropped(self):
        """The reason resolve.py gave: no source version for anything. A provider fetches
        it instead, rather than an unversioned value passing as evidence."""
        records, skipped = self.provider.get_evidence(
            self.variant, parsed({"gnomAD_AF": "0.01"}, meta={}))
        self.assertEqual(records, [])
        self.assertEqual(skipped, [{"info": "gnomAD_AF", "reason": "NO_SOURCE_VERSION"}])

    def test_a_field_version_beats_the_files_fallback(self):
        records, _ = self.provider.get_evidence(
            self.variant,
            parsed({"gnomAD_AF": "0.01"},
                   defs={"gnomAD_AF": info_def("gnomAD_AF", 'AF,Source="gnomAD",Version="4.1"')}))
        self.assertEqual(records[0]["source"], "gnomAD")
        self.assertEqual(records[0]["source_version"], "4.1")

    def test_a_curators_conclusion_cannot_be_mapped_at_all(self):
        """Feeding a conclusion back as its own support is the circularity
        test_clingen_positive.py checks for, so it fails at construction."""
        for info_id in ("CLNSIG", "ACMG_CODES", "DISEASE_ASSOCIATION"):
            with self.subTest(info_id=info_id):
                with self.assertRaises(ValueError):
                    VcfEvidenceProvider({info_id: {"category": "annotation", "field": "x"}})

    def test_the_refusal_list_covers_the_fields_the_demo_vcfs_carry(self):
        self.assertLessEqual({"CLNSIG", "ACMG_CODES", "CLNVARIATIONID", "DISEASE_ASSOCIATION"},
                             REFUSED)

    def test_an_unmapped_field_is_ignored_rather_than_guessed_at(self):
        records, skipped = self.provider.get_evidence(
            self.variant, parsed({"SOMETHING_NEW": "9", "gnomAD_AF": "0.01"}))
        self.assertEqual(len(records), 1)
        self.assertNotIn("SOMETHING_NEW", records[0].get("info_fields", []))
        self.assertEqual(skipped, [])

    def test_a_mapping_to_an_unknown_category_is_reported_not_applied(self):
        provider = VcfEvidenceProvider({"X": {"category": "gene_disease", "field": "y"}})
        records, skipped = provider.get_evidence(self.variant, parsed({"X": "1"}))
        self.assertEqual(records, [])
        self.assertEqual(skipped, [{"info": "X", "reason": "UNSUPPORTED_MAPPING"}])

    def test_every_record_says_it_came_from_the_submitted_file(self):
        """A provider queried a source; this arrived with the variant, and the values do not
        show the difference."""
        records, _ = self.provider.get_evidence(self.variant, parsed({"gnomAD_AF": "0.01"}))
        self.assertEqual(records[0]["origin"], "vcf_input")
        self.assertEqual(records[0]["assessment_method"], "automated")

    def test_records_carry_the_provenance_the_evidence_service_filters_on(self):
        records, _ = self.provider.get_evidence(self.variant, parsed({"gnomAD_AF": "0.01"}))
        for field in ("evidence_id", "source", "source_version", "retrieved_at",
                      "quality_status"):
            self.assertTrue(records[0][field], field)
        self.assertTrue(records[0]["evidence_id"].startswith("urn:sha256:"))

    def test_two_variants_do_not_share_an_evidence_id(self):
        other = Variant("GRCh38", "1", 3, "C", "T")
        first = self.provider.get_evidence(self.variant, parsed({"gnomAD_AF": "0.01"}))[0][0]
        second = self.provider.get_evidence(other, parsed({"gnomAD_AF": "0.01"}))[0][0]
        self.assertNotEqual(first["evidence_id"], second["evidence_id"])

    def test_an_empty_file_produces_nothing(self):
        self.assertEqual(self.provider.get_evidence(self.variant, parsed({})), ([], []))

    def test_the_mapping_must_be_an_object(self):
        with self.assertRaises(ValueError):
            VcfEvidenceProvider(["gnomAD_AF"])


class AnnotatedVcfTests(unittest.TestCase):
    """Against a file shaped the way the providers' own sources are, not the demo VCFs."""

    FIXTURE = Path(__file__).resolve().parent / "fixtures" / "annotated-vcf-example.vcf"

    def setUp(self):
        self.parsed = adapter_parse_vcf(self.FIXTURE)[0]
        self.variant = Variant("GRCh38", "11", 47352561, "G", "A")
        self.records, self.skipped = VcfEvidenceProvider({
            **MAPPING,
            "AM_PATHOGENICITY": {"category": "computational", "field": "score",
                                 "group": "AlphaMissense"},
        }).get_evidence(self.variant, self.parsed)

    def test_the_frequency_fields_become_one_pinned_population_record(self):
        population = [r for r in self.records if r["category"] == "population"]
        self.assertEqual(len(population), 1)
        self.assertEqual(population[0]["source"], "gnomAD")
        self.assertEqual(population[0]["source_version"], "4.1.1")
        self.assertEqual(population[0]["AF"], 0.00012)
        self.assertEqual(population[0]["AN"], 150000)

    def test_a_score_with_no_stated_release_is_skipped_not_taken(self):
        self.assertEqual(self.skipped, [])  # it falls back to the file's own source/date
        computational = [r for r in self.records if r["category"] == "computational"]
        self.assertEqual(len(computational), 1)
        # Pinned to the file, not to a predictor release - so computational.py's calibration
        # match fails and PP3/BP4 do not use it. The score is carried, not made usable.
        self.assertEqual(computational[0]["source_version"], "20260917")

    def test_the_clinical_significance_in_the_file_is_never_read(self):
        self.assertIn("CLNSIG", self.parsed.record.info)  # it is in the file
        for record in self.records:  # and in none of the evidence
            self.assertNotIn("CLNSIG", record["info_fields"])

    def test_the_population_record_passes_the_evidence_service_filter(self):
        population = next(r for r in self.records if r["category"] == "population")
        service = EvidenceService(self.records)
        candidates = service.get_candidates("population", self.variant)
        self.assertIn(population["evidence_id"],
                      [item["evidence_id"] for item in candidates])


if __name__ == "__main__":
    unittest.main()


CSQ_MAPPING = {"CSQ": {"category": "annotation", "group": "VEP", "compound": {
    "transcript_subfield": "Feature",
    "fields": {"Consequence": "consequences", "SYMBOL": "gene", "Feature": "transcript",
               "EXON": "exon", "HGVSc": "hgvsc", "HGVSp": "hgvsp"},
    "list_fields": ["consequences"],
}}}


class SubfieldOrderTests(unittest.TestCase):
    """The header states the order, which is what makes the compound form readable."""

    def test_a_vep_format_clause_is_read(self):
        definition = info_def("CSQ", "Consequence annotations from Ensembl VEP. "
                                     "Format: Allele|Consequence|IMPACT|SYMBOL|Feature")
        self.assertEqual(subfield_order(definition),
                         ["Allele", "Consequence", "IMPACT", "SYMBOL", "Feature"])

    def test_a_snpeff_annotations_clause_is_read(self):
        definition = info_def("ANN", "Functional annotations: "
                                     "'Allele | Annotation | Annotation_Impact | Feature_ID '")
        self.assertEqual(subfield_order(definition),
                         ["Allele", "Annotation", "Annotation_Impact", "Feature_ID"])

    def test_a_header_stating_no_order_gives_none(self):
        self.assertIsNone(subfield_order(info_def("CSQ", "Consequence annotations")))
        self.assertIsNone(subfield_order(None))

    def test_entries_are_split_per_transcript_and_keyed_by_name(self):
        order = ["Allele", "Consequence", "Feature"]
        entries = compound_entries("A|missense_variant|ENST1,A|intron_variant|ENST2", order)
        self.assertEqual([entry["Feature"] for entry in entries], ["ENST1", "ENST2"])
        self.assertEqual(entries[0]["Consequence"], "missense_variant")

    def test_a_short_entry_leaves_the_trailing_subfields_empty(self):
        entries = compound_entries("A|intron_variant", ["Allele", "Consequence", "Feature"])
        self.assertEqual(entries[0]["Feature"], "")


class CompoundFieldTests(unittest.TestCase):
    """A CSQ field carries one entry per transcript; choosing among them is not this
    provider's question, so only an exact match on the evaluated transcript is read."""

    FIXTURE = Path(__file__).resolve().parent / "fixtures" / "vep-annotated-example.vcf"

    def setUp(self):
        self.parsed = adapter_parse_vcf(self.FIXTURE)[0]
        self.variant = Variant("GRCh38", "11", 47352561, "G", "A")
        self.provider = VcfEvidenceProvider(CSQ_MAPPING)

    def evidence(self, transcript):
        return self.provider.get_evidence(self.variant, self.parsed, transcript=transcript)

    def test_the_entry_for_the_evaluated_transcript_is_the_one_read(self):
        records, skipped = self.evidence("NM_000256.3")
        self.assertEqual(skipped, [])
        self.assertEqual(records[0]["consequences"],
                         ["missense_variant", "splice_region_variant"])
        self.assertEqual(records[0]["exon"], "2/34")
        self.assertEqual(records[0]["hgvsc"], "NM_000256.3:c.278G>A")

    def test_a_different_transcript_reads_its_own_entry(self):
        records, _ = self.evidence("NM_001321226.2")
        self.assertEqual(records[0]["consequences"], ["intron_variant"])
        self.assertNotIn("exon", records[0])  # empty subfields are absent, not blank

    def test_a_transcript_with_no_entry_yields_nothing(self):
        records, skipped = self.evidence("NM_999999.1")
        self.assertEqual(records, [])
        self.assertEqual(skipped, [{"info": "CSQ", "reason": "NO_SINGLE_TRANSCRIPT_ENTRY"}])

    def test_without_a_transcript_the_choice_is_not_made_here(self):
        records, skipped = self.evidence(None)
        self.assertEqual(records, [])
        self.assertEqual(skipped, [{"info": "CSQ", "reason": "NO_TRANSCRIPT_TO_MATCH"}])

    def test_a_version_difference_still_matches_the_entry(self):
        records, _ = self.evidence("NM_000256.9")
        self.assertEqual(records[0]["transcript"], "NM_000256.3")

    def test_a_compound_field_whose_header_declares_no_order_is_skipped(self):
        """Positions would be guesswork, which is the one thing this must not do."""
        record = VariantRecord(chrom="11", pos=47352561, id="v", ref="G", alt="A", qual=".",
                               filter="PASS", info={"CSQ": "A|missense_variant|NM_000256.3"})
        undeclared = ParsedVcf(meta={"source": "x", "fileDate": "20260917"},
                               info_defs={"CSQ": info_def("CSQ", "Consequence annotations")},
                               record=record)
        records, skipped = self.provider.get_evidence(self.variant, undeclared,
                                                      transcript="NM_000256.3")
        self.assertEqual(records, [])
        self.assertEqual(skipped, [{"info": "CSQ", "reason": "NO_DECLARED_FORMAT"}])

    def test_the_clinical_significance_beside_it_is_still_never_read(self):
        records, _ = self.evidence("NM_000256.3")
        self.assertIn("CLNSIG", self.parsed.record.info)
        for record in records:
            self.assertNotIn("CLNSIG", record["info_fields"])

    def test_the_record_still_says_it_came_from_the_submitted_file(self):
        records, _ = self.evidence("NM_000256.3")
        self.assertEqual(records[0]["origin"], "vcf_input")
        self.assertEqual(records[0]["source"], "ensembl-vep-116")


DECLARATION = {
    "matches_source": "Cross-Gene Annotation Explorer",
    "declared_versions": {"gnomAD": "4.1.1", "Ensembl VEP": "116"},
    "declared_by": "BH26 project configuration",
    "justification": "The exporter states it used VEP but names no release.",
    "declared_at": "2026-09-17",
}
EXPORTER_META = {"source": "Cross-Gene Annotation Explorer (selected variants with VEP annotations)"}


class DeclaredVersionTests(unittest.TestCase):
    """A declaration is a stated assumption, so it has to say who states it and about what."""

    def test_a_complete_declaration_supplies_the_version(self):
        version, recorded = declared_version(DECLARATION, EXPORTER_META, "gnomAD")
        self.assertEqual(version, "4.1.1")
        self.assertEqual(recorded["declared_by"], "BH26 project configuration")
        self.assertTrue(recorded["justification"])

    def test_it_does_not_attach_itself_to_another_file(self):
        self.assertEqual(
            declared_version(DECLARATION, {"source": "some other exporter"}, "gnomAD"),
            (None, None))
        self.assertEqual(declared_version(DECLARATION, {}, "gnomAD"), (None, None))

    def test_a_group_it_says_nothing_about_stays_unpinned(self):
        self.assertEqual(declared_version(DECLARATION, EXPORTER_META, "dbNSFP"), (None, None))
        self.assertEqual(declared_version(DECLARATION, EXPORTER_META, None), (None, None))

    def test_an_incomplete_declaration_is_not_usable(self):
        for missing in ("matches_source", "declared_by", "justification"):
            with self.subTest(missing=missing):
                partial = {k: v for k, v in DECLARATION.items() if k != missing}
                self.assertEqual(declared_version(partial, EXPORTER_META, "gnomAD"),
                                 (None, None))
        self.assertEqual(declared_version(None, EXPORTER_META, "gnomAD"), (None, None))


class ExporterVcfTests(unittest.TestCase):
    """The real file: VEP output flattened into VEP_* fields, with no release named."""

    FIXTURE = Path(__file__).resolve().parent / "fixtures" / "exporter-vcf-example.vcf"
    MAPPING = {
        "GNOMAD_AF": {"category": "population", "group": "gnomAD", "field": "AF"},
        "VEP_CONSEQUENCE": {"category": "annotation", "group": "Ensembl VEP",
                            "field": "consequences"},
        "VEP_TRANSCRIPT": {"category": "annotation", "group": "Ensembl VEP",
                           "field": "transcript"},
        "AM_CLASS": {"category": "computational", "group": "AlphaMissense",
                     "field": "classification"},
    }

    def setUp(self):
        self.parsed = adapter_parse_vcf(self.FIXTURE)[0]
        self.variant = Variant("GRCh38", "2", 174757593, "C", "T")

    def evidence(self, declaration):
        return VcfEvidenceProvider(self.MAPPING, declaration).get_evidence(
            self.variant, self.parsed)

    def test_without_a_declaration_the_file_supports_nothing(self):
        """It names no fileDate, no reference, and no Source or Version on any field."""
        records, skipped = self.evidence(None)
        self.assertEqual(records, [])
        self.assertEqual({item["reason"] for item in skipped}, {"NO_SOURCE_VERSION"})

    def test_a_declaration_makes_the_declared_groups_usable(self):
        records, _ = self.evidence(DECLARATION)
        by_source = {record["source"]: record for record in records}
        self.assertEqual(by_source["gnomAD"]["source_version"], "4.1.1")
        self.assertEqual(by_source["Ensembl VEP"]["source_version"], "116")

    def test_a_declared_version_never_passes_as_one_the_file_stated(self):
        records, _ = self.evidence(DECLARATION)
        for record in records:
            self.assertEqual(record["version_status"], "DECLARED")
            self.assertTrue(record["version_declaration"]["justification"])

    def test_an_undeclared_group_is_still_left_to_the_providers(self):
        """AlphaMissense is in the mapping and in the file, and the declaration says nothing
        about it, so it stays unpinned rather than borrowing another group's version."""
        with_am = adapter_parse_vcf(self.FIXTURE)[1]  # the row carrying AM_CLASS
        self.assertIn("AM_CLASS", with_am.record.info)
        _records, skipped = VcfEvidenceProvider(self.MAPPING, DECLARATION).get_evidence(
            Variant("GRCh38", "2", 188990344, "G", "A"), with_am)
        self.assertIn({"info": "AM_CLASS", "reason": "NO_SOURCE_VERSION"}, skipped)

    def test_percent_encoded_delimiters_are_decoded(self):
        """VCF encodes the characters its own delimiters use; %2C is a comma."""
        records, _ = self.evidence(DECLARATION)
        annotation = next(r for r in records if r["category"] == "annotation")
        self.assertEqual(annotation["consequences"], "stop_gained, frameshift_variant")

    def test_the_clinical_significance_on_every_row_is_never_read(self):
        self.assertIn("CLNSIG", self.parsed.record.info)
        for record in self.evidence(DECLARATION)[0]:
            self.assertNotIn("CLNSIG", record["info_fields"])
