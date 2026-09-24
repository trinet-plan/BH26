import unittest

from acmg_pipeline.criteria import bp5
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


VARIANT = VariantRecord(
    chrom="1", pos=1, id="test", ref="A", alt="T", qual=".", filter=".",
    info={"GENE": "MYH7", "HGVSC": "c.1594T>C", "HGVSP": "p.Ser532Pro"},
)

MATCHED_JSON = {
    "variant_matching": {"match_status": "matched", "match_type": "protein_notation",
                         "confidence": "high", "notes": "exact match"},
    "alternate_diagnosis_data": {
        "alternate_gene": "MYBPC3", "alternate_variant": "c.1224-52G>A",
        "phenotype_explained_by_alternate": True,
    },
    "single_study_only": True,
    "overall_evidence": {"direction": "BP5", "strength_hint": "not_clear",
                         "rationale": "Phenotype re-attributed to MYBPC3."},
}


class BP5FromJsonTests(unittest.TestCase):
    def test_round_trips_a_full_bp5_judgment(self):
        judgment = bp5.BP5Judgment.from_json(MATCHED_JSON)
        self.assertEqual(judgment.overall_evidence.direction, bp5.AlternateBasisDirection.BP5)
        self.assertEqual(judgment.alternate_diagnosis_data.alternate_gene, "MYBPC3")
        self.assertEqual(judgment.alternate_diagnosis_data.alternate_variant, "c.1224-52G>A")
        self.assertTrue(judgment.alternate_diagnosis_data.phenotype_explained_by_alternate)

    def test_tolerates_a_null_overall_evidence_block(self):
        """The LLM can emit "overall_evidence": null outright - same quirk ps3_bs3/ps4 guard against."""
        data = {**MATCHED_JSON, "overall_evidence": None}
        judgment = bp5.BP5Judgment.from_json(data)
        self.assertEqual(judgment.overall_evidence.direction, bp5.AlternateBasisDirection.NOT_CLEAR)

    def test_missing_alternate_diagnosis_data_defaults_to_all_none(self):
        data = {**MATCHED_JSON, "alternate_diagnosis_data": {}}
        judgment = bp5.BP5Judgment.from_json(data)
        self.assertIsNone(judgment.alternate_diagnosis_data.alternate_gene)
        self.assertIsNone(judgment.alternate_diagnosis_data.phenotype_explained_by_alternate)

    def test_string_shorthand_for_variant_matching_is_accepted(self):
        data = {**MATCHED_JSON, "variant_matching": "unsuccessful"}
        judgment = bp5.BP5Judgment.from_json(data)
        self.assertEqual(judgment.variant_matching.match_status.value, "unsuccessful")


class BP5FinalizeSafetyNetTests(unittest.TestCase):
    def test_bp5_with_alternate_gene_and_phenotype_explained_survives(self):
        judgment = bp5.BP5Judgment.from_json(MATCHED_JSON)
        result = bp5.finalize(judgment, pmid="11112222")
        self.assertEqual(result.effective_direction, bp5.AlternateBasisDirection.BP5)

    def test_bp5_with_no_alternate_gene_is_forced_to_not_clear(self):
        data = {**MATCHED_JSON, "alternate_diagnosis_data": {"phenotype_explained_by_alternate": True}}
        judgment = bp5.BP5Judgment.from_json(data)
        result = bp5.finalize(judgment, pmid="33334444")
        self.assertEqual(result.effective_direction, bp5.AlternateBasisDirection.NOT_CLEAR)
        self.assertTrue(any(h.severity == "warning" for h in result.curator_hints))

    def test_bp5_whose_alternate_finding_is_not_stated_to_explain_phenotype_is_forced_to_not_clear(self):
        data = {**MATCHED_JSON,
                "alternate_diagnosis_data": {"alternate_gene": "TTN",
                                             "phenotype_explained_by_alternate": False}}
        judgment = bp5.BP5Judgment.from_json(data)
        result = bp5.finalize(judgment, pmid="55556666")
        self.assertEqual(result.effective_direction, bp5.AlternateBasisDirection.NOT_CLEAR)

    def test_bp5_whose_alternate_finding_has_a_null_phenotype_explained_is_also_forced_to_not_clear(self):
        """null (not extracted) must not be silently treated as truthy."""
        data = {**MATCHED_JSON, "alternate_diagnosis_data": {"alternate_gene": "TTN"}}
        judgment = bp5.BP5Judgment.from_json(data)
        result = bp5.finalize(judgment, pmid="77778888")
        self.assertEqual(result.effective_direction, bp5.AlternateBasisDirection.NOT_CLEAR)

    def test_not_clear_direction_is_left_alone(self):
        data = {**MATCHED_JSON, "overall_evidence": {"direction": "not_clear", "rationale": "no evidence"}}
        judgment = bp5.BP5Judgment.from_json(data)
        result = bp5.finalize(judgment, pmid="99990000")
        self.assertEqual(result.effective_direction, bp5.AlternateBasisDirection.NOT_CLEAR)
        self.assertEqual(result.curator_hints, [])

    def test_unmatched_variant_generates_no_hints_even_with_bp5_direction(self):
        """An UNSUCCESSFUL match means the rest of the fields are noise - see ps4's own analogous guard."""
        data = {**MATCHED_JSON, "variant_matching": {"match_status": "unsuccessful"}}
        judgment = bp5.BP5Judgment.from_json(data)
        result = bp5.finalize(judgment, pmid="12121212")
        self.assertEqual(result.curator_hints, [])


class BP5PromptTests(unittest.TestCase):
    def test_prompt_names_the_target_variant_and_gene(self):
        prompt = bp5.build_prompt(VARIANT, ClinicalNoteExtraction(), "full text of the paper")
        self.assertIn("MYH7", prompt)
        self.assertIn("c.1594T>C", prompt)
        self.assertIn("full text of the paper", prompt)
        self.assertIn("BP5", prompt)

    def test_prompt_instructs_against_inferring_bp5_from_the_target_variant_alone(self):
        """The key false-positive mode this criterion must avoid: doubt about the
        target variant on its own is not the same as a positively identified
        alternate cause - see the module's own PROMPT_TEMPLATE step 4."""
        prompt = bp5.build_prompt(VARIANT, ClinicalNoteExtraction(), "text")
        self.assertIn("Do NOT infer", prompt)


class BP5StructuredEvidenceItemsTests(unittest.TestCase):
    def test_lists_the_extracted_alternate_finding(self):
        judgment = bp5.BP5Judgment.from_json(MATCHED_JSON)
        items = judgment.structured_evidence_items()
        labels = {item["label"] for item in items}
        self.assertIn("Alternate gene named", labels)
        self.assertIn("Alternate variant named", labels)
        phenotype_item = next(i for i in items if "explains the phenotype" in i["label"])
        self.assertTrue(phenotype_item["checked"])

    def test_omits_gene_and_variant_items_when_nothing_was_extracted(self):
        data = {**MATCHED_JSON, "alternate_diagnosis_data": {}}
        judgment = bp5.BP5Judgment.from_json(data)
        items = judgment.structured_evidence_items()
        labels = {item["label"] for item in items}
        self.assertNotIn("Alternate gene named", labels)
        self.assertNotIn("Alternate variant named", labels)


class BP5AggregationTests(unittest.TestCase):
    def test_empty_contributions_aggregate_to_not_clear(self):
        aggregated = bp5.aggregate_multi_paper_results([])
        self.assertEqual(aggregated.aggregated_direction, bp5.AlternateBasisDirection.NOT_CLEAR)


class BP5PipelineWiringTests(unittest.TestCase):
    """The two "forgot to wire the fast path" bugs this session already hit twice for
    BA1/BS1 and PM1 were each a default/registration that silently kept excluding the
    new code from an actual evaluation path. Guard the three analogous spots for BP5."""

    def test_bp5_is_a_literature_code(self):
        from acmg_pipeline.constants import LITERATURE_CODES
        self.assertIn("BP5", LITERATURE_CODES)

    def test_bp5_has_a_registered_engine(self):
        from acmg_pipeline.pipeline import ENGINE_BY_CRITERION, BP5_ENGINE
        self.assertIs(ENGINE_BY_CRITERION["BP5"], BP5_ENGINE)

    def test_bp5_is_in_both_default_literature_criteria_tuples(self):
        import inspect
        from acmg_pipeline.pipeline import (
            judge_variant_from_shared_input, judge_variant_from_structured_input,
        )
        for fn in (judge_variant_from_shared_input, judge_variant_from_structured_input):
            default = inspect.signature(fn).parameters["criteria"].default
            self.assertIn("BP5", default, f"{fn.__name__}'s default criteria tuple is missing BP5")


if __name__ == "__main__":
    unittest.main()
