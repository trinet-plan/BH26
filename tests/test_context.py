import unittest

from acmg_pipeline.automated_core.context import apply_context, context_summary, load_context


VARIANT = {"assembly": "GRCh38", "chrom": "1", "pos": 2, "ref": "C", "alt": "T"}
KEY = "GRCh38:1:2:C:T"
ENTRY = {"variant_key": KEY, "gene": "TEST", "hgvs_c": "c.1A>T", "caid": "CA000000",
         "resolved_by": "ClinGen Allele Registry"}
EXCEPTIONS = {"source": "ClinGen SVI BA1 exception list", "source_version": "2018",
              "reviewed_at": "2026-09-15", "complete": True,
              "entry_method": "manual_transcription", "variants": [ENTRY]}


def document(**overrides):
    return {"schema_version": "1.0", "context_version": "2026-09-15", **overrides}


class CuratedContextTests(unittest.TestCase):
    def record(self, variant=None):
        return {"record_id": "test", "variant": variant or VARIANT,
                "identity_provenance": [{"source": "test"}]}

    def test_disease_context_reaches_the_record(self):
        context = load_context(document(records={KEY: {
            "condition": "MONDO:0000001", "condition_label": "Test disease",
            "inheritance": "autosomal_dominant",
            "disease_frequency_threshold": {"max_credible_af": 0.001}}}))
        updated = apply_context(self.record(), context)
        self.assertEqual(updated["condition"], "MONDO:0000001")
        self.assertEqual(updated["condition_label"], "Test disease")
        self.assertEqual(updated["disease_frequency_threshold"]["max_credible_af"], 0.001)
        # Records without curated context are returned unchanged.
        other = apply_context(self.record({**VARIANT, "pos": 9}), context)
        self.assertNotIn("condition", other)

    def test_an_incomplete_exception_list_resolves_nothing(self):
        """An untranscribed list cannot claim a variant is absent from it."""
        context = load_context(document(ba1_exceptions={**EXCEPTIONS, "complete": False,
                                                        "variants": []}))
        self.assertNotIn("ba1_exception_assessment", apply_context(self.record(), context))

    def test_a_complete_list_resolves_both_answers(self):
        context = load_context(document(ba1_exceptions=EXCEPTIONS))
        listed = apply_context(self.record(), context)["ba1_exception_assessment"]
        self.assertTrue(listed["is_exception"])
        self.assertEqual(listed["source_version"], "2018")
        absent = apply_context(self.record({**VARIANT, "pos": 9}), context)
        self.assertFalse(absent["ba1_exception_assessment"]["is_exception"])

    def test_an_empty_but_complete_list_is_a_negative_answer(self):
        context = load_context(document(ba1_exceptions={**EXCEPTIONS, "variants": []}))
        assessment = apply_context(self.record(), context)["ba1_exception_assessment"]
        self.assertFalse(assessment["is_exception"])
        self.assertEqual(assessment["list_size"], 0)

    def test_each_entry_records_how_it_was_resolved(self):
        """A wrong key would silently exempt the wrong variant, so the mapping is recorded."""
        entry = {key: value for key, value in ENTRY.items() if key != "resolved_by"}
        with self.assertRaisesRegex(ValueError, "resolved_by"):
            load_context(document(ba1_exceptions={**EXCEPTIONS, "variants": [entry]}))
        with self.assertRaisesRegex(ValueError, "variant_key"):
            load_context(document(ba1_exceptions={**EXCEPTIONS, "variants": [{"gene": "TEST"}]}))

    def test_how_the_list_was_entered_must_be_stated(self):
        exceptions = {key: value for key, value in EXCEPTIONS.items() if key != "entry_method"}
        with self.assertRaisesRegex(ValueError, "entry_method"):
            load_context(document(ba1_exceptions=exceptions))
        with self.assertRaisesRegex(ValueError, "entry_method"):
            load_context(document(ba1_exceptions={**EXCEPTIONS, "entry_method": "guessed"}))

    def test_completeness_must_be_stated(self):
        exceptions = {key: value for key, value in EXCEPTIONS.items() if key != "complete"}
        with self.assertRaisesRegex(ValueError, "complete"):
            load_context(document(ba1_exceptions=exceptions))

    def test_provenance_is_required(self):
        exceptions = {key: value for key, value in EXCEPTIONS.items() if key != "source_version"}
        with self.assertRaisesRegex(ValueError, "source_version"):
            load_context(document(ba1_exceptions=exceptions))

    def test_context_version_is_required(self):
        with self.assertRaisesRegex(ValueError, "context_version"):
            load_context({"schema_version": "1.0"})

    def test_evidence_cannot_be_smuggled_in_as_context(self):
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            load_context(document(records={KEY: {"category": "population", "AF": "0.9"}}))

    def test_keys_must_be_variant_keys(self):
        with self.assertRaisesRegex(ValueError, "variant key"):
            load_context(document(records={"MYBPC3": {"condition": "MONDO:0000001"}}))

    def test_record_context_overrides_reused_variant_context(self):
        context = load_context(document(
            records={KEY: {"condition": "MONDO:variant"}},
            record_contexts={"case2:1:1": {"condition": "MONDO:record"}},
        ))
        first = apply_context({**self.record(), "record_id": "case1:1:1"}, context)
        second = apply_context({**self.record(), "record_id": "case2:1:1"}, context)
        self.assertEqual(first["condition"], "MONDO:variant")
        self.assertEqual(second["condition"], "MONDO:record")

    def test_unknown_record_context_field_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unsupported curated record context"):
            load_context(document(record_contexts={"case1:1:1": {"unknown": True}}))

    # --- gene_frequency_thresholds (BA1/BS1's own VCEP-specific numbers) --------------

    GENE_ENTRY = {
        "ba1": {"max_af": 0.001, "frequency_statistic": "faf95", "comparison": ">=",
               "source": "ClinGen CSpec test VCEP", "source_version": "1.0",
               "reviewed_at": "2026-09-22"},
        "bs1": {"max_credible_af": 0.0001, "frequency_statistic": "faf95", "comparison": ">=",
               "condition": "MONDO:0000001", "condition_scope": "gene_wide",
               "inheritance": "autosomal_dominant", "source": "ClinGen CSpec test VCEP",
               "source_version": "1.0", "reviewed_at": "2026-09-22"},
    }

    def test_gene_threshold_feeds_ba1_and_bs1_by_gene_not_variant_key(self):
        context = load_context(document(gene_frequency_thresholds={"TEST": self.GENE_ENTRY}))
        record = {**self.record(), "gene": "TEST"}
        updated = apply_context(record, context)
        self.assertEqual(updated["ba1_threshold_override"], self.GENE_ENTRY["ba1"])
        self.assertEqual(updated["disease_frequency_threshold"], self.GENE_ENTRY["bs1"])

    def test_an_unmatched_gene_resolves_nothing(self):
        context = load_context(document(gene_frequency_thresholds={"TEST": self.GENE_ENTRY}))
        updated = apply_context({**self.record(), "gene": "OTHER"}, context)
        self.assertNotIn("ba1_threshold_override", updated)
        self.assertNotIn("disease_frequency_threshold", updated)

    def test_a_more_specific_record_context_is_not_overwritten_by_the_gene_default(self):
        """A per-record curated disease_frequency_threshold is a more specific curation than
        a gene-wide one, and must win - the gene-wide entry only fills a gap."""
        context = load_context(document(
            gene_frequency_thresholds={"TEST": self.GENE_ENTRY},
            record_contexts={"test": {"disease_frequency_threshold": {"max_credible_af": 0.5}}},
        ))
        updated = apply_context({**self.record(), "gene": "TEST"}, context)
        self.assertEqual(updated["disease_frequency_threshold"]["max_credible_af"], 0.5)
        # ba1_threshold_override has no per-record entry here, so the gene-wide one still fills it.
        self.assertEqual(updated["ba1_threshold_override"], self.GENE_ENTRY["ba1"])

    def test_gene_frequency_thresholds_requires_its_own_provenance(self):
        incomplete = {"ba1": {"max_af": 0.001, "frequency_statistic": "faf95"}}  # no source/etc
        with self.assertRaisesRegex(ValueError, "requires"):
            load_context(document(gene_frequency_thresholds={"TEST": incomplete}))

    def test_gene_frequency_thresholds_rejects_an_unknown_field(self):
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            load_context(document(gene_frequency_thresholds={"TEST": {"pm2": {}}}))

    # --- gene_critical_domains (PM1's own VCEP-published codon ranges) --------------

    DOMAIN_ENTRY = {
        "source": "ClinGen CSpec test VCEP", "source_version": "1.0",
        "reviewed_at": "2026-09-22",
        "ranges": [{"start": 90, "end": 162, "label": "test domain"}],
    }

    def test_gene_critical_domains_feeds_the_record_by_gene(self):
        context = load_context(document(gene_critical_domains={"TEST": self.DOMAIN_ENTRY}))
        updated = apply_context({**self.record(), "gene": "TEST"}, context)
        self.assertEqual(updated["gene_critical_domains"], self.DOMAIN_ENTRY)

    def test_an_unmatched_gene_resolves_no_critical_domains(self):
        context = load_context(document(gene_critical_domains={"TEST": self.DOMAIN_ENTRY}))
        updated = apply_context({**self.record(), "gene": "OTHER"}, context)
        self.assertNotIn("gene_critical_domains", updated)

    def test_gene_critical_domains_requires_its_own_provenance(self):
        incomplete = {"source": "x", "source_version": "1", "reviewed_at": "2026-09-22"}  # no ranges
        with self.assertRaisesRegex(ValueError, "requires"):
            load_context(document(gene_critical_domains={"TEST": incomplete}))

    def test_gene_critical_domains_rejects_a_reversed_range(self):
        bad = {**self.DOMAIN_ENTRY, "ranges": [{"start": 200, "end": 100}]}
        with self.assertRaisesRegex(ValueError, "start<=end"):
            load_context(document(gene_critical_domains={"TEST": bad}))

    def test_summary_reports_what_was_loaded(self):
        summary = context_summary(load_context(document(ba1_exceptions=EXCEPTIONS,
                                                        records={KEY: {"condition": "x"}})))
        self.assertEqual(summary["records"], 1)
        self.assertEqual(summary["record_contexts"], 0)
        self.assertEqual(summary["ba1_exceptions"]["variants"], 1)
        self.assertTrue(summary["ba1_exceptions"]["complete"])
        # Every run records that the committed list was typed in by hand.
        self.assertEqual(summary["ba1_exceptions"]["entry_method"], "manual_transcription")
        self.assertIsNone(context_summary(None))


if __name__ == "__main__":
    unittest.main()
