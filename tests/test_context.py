import unittest

from acmg.core.context import apply_context, context_summary, load_context


VARIANT = {"assembly": "GRCh38", "chrom": "1", "pos": 2, "ref": "C", "alt": "T"}
KEY = "GRCh38:1:2:C:T"
ENTRY = {"variant_key": KEY, "gene": "TEST", "hgvs_c": "c.1A>T", "caid": "CA000000",
         "resolved_by": "ClinGen Allele Registry"}
EXCEPTIONS = {"source": "ClinGen SVI BA1 exception list", "source_version": "2018",
              "reviewed_at": "2026-09-15", "complete": True, "variants": [ENTRY]}


def document(**overrides):
    return {"schema_version": "1.0", "context_version": "2026-09-15", **overrides}


class CuratedContextTests(unittest.TestCase):
    def record(self, variant=None):
        return {"record_id": "test", "variant": variant or VARIANT,
                "identity_provenance": [{"source": "test"}]}

    def test_disease_context_reaches_the_record(self):
        context = load_context(document(records={KEY: {
            "condition": "MONDO:0000001", "inheritance": "autosomal_dominant",
            "disease_frequency_threshold": {"max_credible_af": 0.001}}}))
        updated = apply_context(self.record(), context)
        self.assertEqual(updated["condition"], "MONDO:0000001")
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

    def test_summary_reports_what_was_loaded(self):
        summary = context_summary(load_context(document(ba1_exceptions=EXCEPTIONS,
                                                        records={KEY: {"condition": "x"}})))
        self.assertEqual(summary["records"], 1)
        self.assertEqual(summary["ba1_exceptions"]["variants"], 1)
        self.assertTrue(summary["ba1_exceptions"]["complete"])
        self.assertIsNone(context_summary(None))


if __name__ == "__main__":
    unittest.main()
