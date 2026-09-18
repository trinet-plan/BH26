from acmg_pipeline.constants import CriterionStatus
import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.automated_engine import evaluate_prepared_record, make_services


RULES = {"PVS1": {
    "ruleset": {
        "name": "clingen_general_pvs1",
        "version": "2018+2023-splicing",
        "sources": [
            {"id": "acmg_amp_2015", "pmid": "25741868"},
            {"id": "clingen_pvs1_2018", "pmid": "30192042"},
            {"id": "clingen_splicing_2023", "pmid": "37352859"},
        ],
    },
    "rules": {"protein_loss_threshold": 0.10,
              "mechanism_source_precedence": ["better_source", "worse_source"]},
}}


class PVS1DecisionTreeTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.condition = "MONDO:TEST"
        self.input = {"variant": self.variant.to_dict(), "transcript": "NM_TEST.1",
                      "condition": self.condition}
        self.serial = 0

    def item(self, category, **values):
        self.serial += 1
        return {
            "category": category,
            "variant_key": self.variant.key,
            "evidence_id": f"test:{category}:{self.serial}",
            "source": "synthetic",
            "source_version": "1",
            "retrieved_at": "2026-09-16",
            "quality_status": "PASS",
            "curator": "test",
            "reviewed_at": "2026-09-16",
            **values,
        }

    def annotation(self, consequence="frameshift_variant"):
        return self.item("annotation", transcript="NM_TEST.1", gene="TEST",
                         consequences=[consequence], protein_id="NP_TEST.1")

    def mechanism(self, established=True, **values):
        """Curated for this case's disease unless a test says otherwise.

        A mechanism curated for no disease in particular no longer carries PVS1 to a verdict,
        so leaving this unscoped would stop every decision-tree test at the disease gate
        instead of testing the tree.
        """
        values.setdefault("condition", self.condition)
        return self.item("gene_disease", gene="TEST",
                         lof_mechanism_established=established, **values)

    def transcript(self, relevance="RELEVANT", exon_relevance="RELEVANT", **values):
        return self.item("transcript_assessment", transcript="NM_TEST.1",
                         relevance=relevance, exon_relevance=exon_relevance, **values)

    def nmd(self, predicted=True, **values):
        return self.item("nmd_prediction", transcript="NM_TEST.1", predicted=predicted,
                         exon="2/5", distance_to_final_junction=100,
                         rule_source="ClinGen PVS1 2018", **values)

    def region(self, **values):
        return self.item("protein_region", transcript="NM_TEST.1", **values)

    def population_lof(self, frequent=False, **values):
        return self.item("population_lof", transcript="NM_TEST.1",
                         lof_variants_frequent=frequent, **values)

    def evaluate_result(self, *items, input_data=None, config=RULES):
        return evaluate_prepared_record(input_data or self.input, make_services(list(items)),
                                        config, ["PVS1"])[0]

    def truncating_evidence(self, consequence="frameshift_variant", *, condition=None):
        scoped = {"condition": condition} if condition else {}
        return [
            self.annotation(consequence),
            self.mechanism(**scoped),
            self.transcript(**scoped),
            self.nmd(**scoped),
        ]

    def gene_level_evidence(self, consequence="frameshift_variant"):
        """Truncating evidence whose mechanism names no disease."""
        return [self.annotation(consequence),
                self.item("gene_disease", gene="TEST", lof_mechanism_established=True),
                self.transcript(), self.nmd()]

    def test_condition_missing_stops_before_a_verdict_but_keeps_the_variant_work(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        value = self.evaluate_result(*self.truncating_evidence(), input_data=input_data)
        self.assertEqual((value.status, value.strength), (CriterionStatus.UNKNOWN, None))
        self.assertEqual(value.evaluation_context["condition_status"], "NOT_PROVIDED")
        self.assertEqual(value.evaluation_context["applicability"], "MANUAL_REVIEW")
        self.assertIn("condition", value.missing_inputs)
        self.assertTrue(value.warnings)
        self.assertEqual([node["node_id"] for node in value.decision_trace], ["C01", "D01"])
        self.assertEqual(len(value.rules_used), 3)

        # The variant-level tree still ran, and says what PVS1 would have concluded.
        preliminary = value.provenance["preliminary_assessment"]
        self.assertTrue(preliminary["eligible_lof_variant"])
        self.assertEqual(preliminary["candidate_strength"], "very_strong")
        self.assertEqual(preliminary["decision_path"], "NF03")
        self.assertEqual([node["node_id"] for node in preliminary["decision_trace"]],
                         ["NF01", "NF02", "NF03"])

    def test_the_genes_curated_diseases_are_offered_when_none_was_supplied(self):
        """Named so a curator can supply the disease context, not so PVS1 can pick one."""
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:1", inheritance="AD"),
                 self.mechanism(False, condition="MONDO:2", evidence_id="second")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIsNone(value.strength)
        offered = value.provenance["candidate_conditions"]
        self.assertEqual([item["condition"] for item in offered], ["MONDO:1", "MONDO:2"])
        # The mechanism each one carries goes with it, because that is what differs between
        # a gene's diseases and what the curator is choosing between.
        self.assertEqual([item["lof_mechanism_established"] for item in offered], [True, False])
        self.assertEqual(offered[0]["inheritance"], "autosomal_dominant")
        self.assertEqual(offered[0]["inheritance_as_recorded"], "AD")
        self.assertIn("MONDO:1", " ".join(value.review_points))
        self.assertIn("MONDO:2", " ".join(value.review_points))

    def test_a_single_candidate_is_offered_and_still_not_taken(self):
        """85% of curated genes have one disease, which is exactly the shape that invites
        picking it - and a curated disease is not evidence that this case is about it."""
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        value = self.evaluate_result(*self.truncating_evidence(), input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual([item["condition"] for item in
                          value.provenance["candidate_conditions"]], [self.condition])
        self.assertEqual(value.evaluation_context["disease_match"], "UNKNOWN")
        self.assertIsNone(value.evaluation_context["mechanism_source"])

    def validity(self, condition, classification="Definitive", label="A disease", moi="AD",
                 **values):
        """A ClinGen gene-disease validity row: an association, never a mechanism."""
        return self.item("gene_disease_validity", gene="TEST", condition=condition,
                         condition_label=label, classification=classification, moi=moi,
                         **values)

    def test_a_validity_curation_adds_a_disease_the_mechanism_sources_do_not_name(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:1"),
                 self.validity("MONDO:2", classification="Limited", label="Another disease")]
        value = self.evaluate_result(*items, input_data=input_data)
        offered = value.provenance["candidate_conditions"]
        self.assertEqual([item["condition"] for item in offered], ["MONDO:1", "MONDO:2"])
        # The established mechanism leads; the association follows as context.
        self.assertIs(offered[0]["lof_mechanism_established"], True)
        self.assertIsNone(offered[1]["lof_mechanism_established"])
        self.assertEqual(offered[1]["gene_disease_validity"], "Limited")
        self.assertEqual(offered[1]["condition_label"], "Another disease")

    def test_a_validity_curation_never_reaches_the_mechanism_gate(self):
        """It says the gene and the disease are related, which is not a mechanism. Its own
        category keeps it out of the lookup rather than a rule that could be forgotten."""
        input_data = {**self.input, "condition": "MONDO:1"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.validity("MONDO:1", classification="Definitive")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIsNone(value.strength)
        self.assertIn("loss-of-function disease mechanism", value.missing_inputs)

    def test_one_disease_named_by_both_sources_is_offered_once(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:1"),
                 self.validity("MONDO:1", label="The disease")]
        value = self.evaluate_result(*items, input_data=input_data)
        offered = value.provenance["candidate_conditions"]
        self.assertEqual(len(offered), 1)
        self.assertEqual(offered[0]["condition_label"], "The disease")
        self.assertIs(offered[0]["lof_mechanism_established"], True)
        self.assertEqual(offered[0]["gene_disease_validity"], "Definitive")
        self.assertEqual(len(offered[0]["sources"]), 2)

    def test_mechanism_sources_that_disagree_leave_it_unstated(self):
        """Resolved here it would look settled; both readings stay visible instead."""
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:1"),
                 self.mechanism(False, condition="MONDO:1", evidence_id="second")]
        value = self.evaluate_result(*items, input_data=input_data)
        offered = value.provenance["candidate_conditions"]
        self.assertIsNone(offered[0]["lof_mechanism_established"])
        self.assertEqual(
            sorted(item["lof_mechanism_established"] for item in offered[0]["sources"]),
            [False, True])

    def test_nothing_is_offered_when_the_gene_has_no_curated_disease(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        value = self.evaluate_result(*self.gene_level_evidence(), input_data=input_data)
        self.assertNotIn("candidate_conditions", value.provenance)
        self.assertEqual(value.evaluation_context["applicability"], "NOT_EVALUATED")
        self.assertEqual(value.review_points, [])

    def test_another_genes_diseases_are_not_offered(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.item("gene_disease", gene="OTHER", condition="MONDO:9",
                           lof_mechanism_established=True)]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertNotIn("candidate_conditions", value.provenance)

    def test_a_preliminary_strength_never_becomes_the_criterion_strength(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        value = self.evaluate_result(*self.truncating_evidence(), input_data=input_data)
        self.assertIsNone(value.strength)
        self.assertIsNone(value.direction)
        self.assertIsNone(value.evidence_outcome)

    def test_condition_specific_mechanism_precedes_gene_level_fallback(self):
        input_data = {**self.input, "condition": "MONDO:1", "condition_label": "Disease"}
        items = self.truncating_evidence()
        items.append(self.mechanism(False, condition="MONDO:1"))
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertTrue(value.evaluation_context["condition_specific"])
        self.assertEqual(value.evaluation_context["mechanism_scope"], "CONDITION_SPECIFIC")

        fallback = self.evaluate_result(*self.gene_level_evidence(), input_data=input_data)
        self.assertEqual(fallback.status, CriterionStatus.UNKNOWN)
        self.assertEqual(fallback.evaluation_context["condition_status"], "PROVIDED")
        self.assertFalse(fallback.evaluation_context["condition_specific"])
        self.assertEqual(fallback.evaluation_context["applicability"], "NOT_EVALUATED")
        self.assertIn("disease-specific loss-of-function mechanism", fallback.missing_inputs)
        self.assertEqual(fallback.provenance["preliminary_assessment"]["candidate_strength"],
                         "very_strong")

    def test_other_condition_is_not_borrowed_and_unknown_is_not_negative(self):
        input_data = {**self.input, "condition": "MONDO:1"}
        other = self.mechanism(True, condition="MONDO:2")
        value = self.evaluate_result(self.annotation(), other, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIsNone(value.strength)
        # Not borrowed - but not discarded either: the curation is named and handed over.
        self.assertIn("disease-specific loss-of-function mechanism",
                      value.unresolved_requirements)
        self.assertEqual(value.evaluation_context["applicability"], "MANUAL_REVIEW")
        self.assertIn("MONDO:2", " ".join(value.review_points))

        unknown = self.mechanism(None, condition="MONDO:1")
        value = self.evaluate_result(self.annotation(), unknown, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("lof_mechanism_established", value.missing_inputs)

    def test_disease_match_is_exact_only_on_an_identifier_match(self):
        input_data = {**self.input, "condition": "MONDO:1"}
        items = self.truncating_evidence()
        items.append(self.mechanism(True, condition="MONDO:1"))
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["disease_match"], "EXACT")

        # A gene-level record carries no disease identifier, so nothing was matched and the
        # criterion must not be reported as if the disease had been confirmed.
        fallback = self.evaluate_result(*self.gene_level_evidence(), input_data=input_data)
        self.assertEqual(fallback.status, CriterionStatus.UNKNOWN)
        self.assertEqual(fallback.evaluation_context["disease_match"], "UNKNOWN")

    def test_a_mapped_condition_matches_but_is_reported_as_equivalent(self):
        """The case is recorded in OMIM and the curation in MONDO. They are the same disease
        once mapped, and the result says that a mapping was what lined them up."""
        mapping = {"input_condition": "OMIM:143890",
                   "normalized_condition": self.condition,
                   "mapping_type": "equivalent"}
        input_data = {**self.input, "condition": "OMIM:143890", "condition_mapping": mapping}
        value = self.evaluate_result(*self.truncating_evidence(), input_data=input_data)
        self.assertEqual((value.status, value.strength), (CriterionStatus.MET, "very_strong"))
        self.assertEqual(value.evaluation_context["disease_match"], "EQUIVALENT")
        self.assertTrue(value.evaluation_context["condition_specific"])

    def test_a_mapping_on_the_record_side_matches_the_same_way(self):
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="OMIM:143890", condition_mapping={
                     "input_condition": "OMIM:143890",
                     "normalized_condition": self.condition,
                     "mapping_type": "equivalent"})]
        value = self.evaluate_result(*items)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["disease_match"], "EQUIVALENT")

    def test_an_unmapped_identifier_from_another_vocabulary_does_not_match(self):
        """Without a mapping the two identifiers are different strings, and PVS1 must not
        read a resemblance into them - it stops and says the mechanism is not for this
        disease."""
        input_data = {**self.input, "condition": "OMIM:143890"}
        value = self.evaluate_result(*self.truncating_evidence(), input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.evaluation_context["disease_match"], "UNKNOWN")

    def test_inheritance_is_normalized_across_vocabularies(self):
        for declared, supplied in (("AR", "autosomal recessive"),
                                   ("autosomal_recessive", "AR"),
                                   ("Autosomal-Recessive", "autosomal_recessive")):
            with self.subTest(declared=declared, supplied=supplied):
                input_data = {**self.input, "inheritance": supplied}
                items = self.truncating_evidence()
                items.append(self.mechanism(True, inheritance=declared))
                value = self.evaluate_result(*items, input_data=input_data)
                self.assertEqual(value.status, CriterionStatus.MET)
                self.assertEqual(value.evaluation_context["inheritance"],
                                 "autosomal_recessive")
                self.assertEqual(value.evaluation_context["moi_match"], "MATCHED")

    def test_a_phenotype_lumped_into_the_curated_disease_is_used(self):
        """The panel decided this phenotype is the disease it curated, so the mechanism is
        about this case - an explicit decision, unlike an ontology relation."""
        input_data = {**self.input, "condition": "MONDO:0014593"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:0013212", condition_scope={
                     "included": ["MONDO:0013212", "MONDO:0014593"], "excluded": []})]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual((value.status, value.strength), (CriterionStatus.MET, "very_strong"))
        self.assertEqual(value.evaluation_context["disease_match"], "INCLUDED")
        self.assertTrue(value.evaluation_context["condition_specific"])

    def test_an_excluded_phenotype_is_answered_rather_than_asked_about(self):
        """ClinGen looked at this phenotype and kept it out, so a curator is told the answer
        instead of being asked the question again."""
        input_data = {**self.input, "condition": "MONDO:0030517"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:0013212", condition_scope={
                     "included": ["MONDO:0013212"], "excluded": ["MONDO:0030517"]})]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.evaluation_context["disease_match"], "EXCLUDED")
        # The curation is inapplicable, which is not the same as PVS1 being inapplicable: a
        # mechanism for this disease could still come from somewhere else, so the criterion
        # is unevaluated rather than ruled out.
        self.assertEqual(value.evaluation_context["applicability"], "NOT_EVALUATED")
        self.assertEqual(value.review_points, [])
        self.assertEqual(
            {node["node_id"]: node["result"] for node in value.decision_trace},
            {"C01": "PASS", "D01": "NOT_APPLICABLE"})

    def test_an_exclusion_outranks_an_ontology_relation(self):
        """A resemblance cannot reinstate a decision the panel already made."""
        input_data = {**self.input, "condition": "MONDO:0030517",
                      "condition_ancestors": ["MONDO:0013212"]}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:0013212", condition_scope={
                     "included": ["MONDO:0013212"], "excluded": ["MONDO:0030517"]})]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.evaluation_context["disease_match"], "EXCLUDED")

    def test_an_exclusion_on_one_curation_does_not_bury_a_match_on_another(self):
        input_data = {**self.input, "condition": "MONDO:0030517"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:0013212", condition_scope={
                     "included": ["MONDO:0013212"], "excluded": ["MONDO:0030517"]}),
                 self.mechanism(True, condition="MONDO:0030517", evidence_id="other")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["disease_match"], "EXACT")

    def test_an_exact_match_is_preferred_over_an_inclusion(self):
        input_data = {**self.input, "condition": "MONDO:0014593"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(False, condition="MONDO:0013212", condition_scope={
                     "included": ["MONDO:0013212", "MONDO:0014593"], "excluded": []}),
                 self.mechanism(True, condition="MONDO:0014593", evidence_id="exact")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["disease_match"], "EXACT")

    def test_a_parent_disease_curation_is_handed_to_a_curator_not_used(self):
        """The case names a subtype and the curation names the disease above it. That is a
        reason to ask a person, never a reason to decide: the same gene can lose function in
        one subtype and gain it in another."""
        input_data = {**self.input, "condition": "MONDO:0007268",
                      "condition_ancestors": ["MONDO:0005045", "MONDO:0004994"]}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:0005045")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIsNone(value.strength)
        self.assertEqual(value.evaluation_context["disease_match"], "PARENT_CHILD")
        self.assertEqual(value.evaluation_context["applicability"], "MANUAL_REVIEW")
        self.assertEqual(
            {node["node_id"]: node["result"] for node in value.decision_trace},
            {"C01": "PASS", "D01": "MANUAL_REVIEW"})
        self.assertTrue(value.review_points)
        # The related curation is attached, and the variant-level work is kept.
        self.assertIn("MONDO:0005045", {item.get("condition") for item in value.evidence})
        self.assertEqual(value.provenance["preliminary_assessment"]["candidate_strength"],
                         "very_strong")

    def test_a_subtype_curation_is_related_in_the_other_direction_too(self):
        """The case names the disease and the curation names a subtype of it. The ancestry
        sits on the record here rather than on the case."""
        input_data = {**self.input, "condition": "MONDO:0005045"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:0007268",
                                condition_ancestors=["MONDO:0005045"])]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.evaluation_context["disease_match"], "PARENT_CHILD")
        self.assertEqual(value.evaluation_context["applicability"], "MANUAL_REVIEW")

    def test_an_exact_match_is_never_downgraded_to_a_relation(self):
        input_data = {**self.input, "condition": "MONDO:0005045",
                      "condition_ancestors": ["MONDO:0004994"]}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:0005045")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["disease_match"], "EXACT")

    def test_unrelated_diseases_stay_unrelated_and_reach_a_curator(self):
        """No relation is invented, and the curation on file is still named: an expert panel
        has said what loss of function does in this gene, and whether that carries to another
        disease entity is a question for a person."""
        input_data = {**self.input, "condition": "MONDO:0007268",
                      "condition_ancestors": ["MONDO:0005045"]}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, condition="MONDO:0009861")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.evaluation_context["disease_match"], "UNKNOWN")
        self.assertEqual(value.evaluation_context["applicability"], "MANUAL_REVIEW")
        self.assertEqual(
            {node["node_id"]: node["value"] for node in value.decision_trace},
            {"C01": "PROVIDED", "D01": "OTHER_DISEASE_CURATED"})
        self.assertIn("MONDO:0009861", {item.get("condition") for item in value.evidence})
        self.assertEqual(value.provenance["preliminary_assessment"]["candidate_strength"],
                         "very_strong")

    def test_a_gene_with_no_curation_at_all_is_not_a_review(self):
        """Nothing on file is a different situation from something on file that may not
        apply, and a curator can only act on the second."""
        input_data = {**self.input, "condition": "MONDO:0007268"}
        items = [self.annotation(), self.transcript(), self.nmd()]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.evaluation_context["applicability"], "NOT_EVALUATED")
        self.assertEqual(value.review_points, [])

    def test_an_unscoped_record_alone_is_still_not_a_review(self):
        value = self.evaluate_result(*self.gene_level_evidence(),
                                     input_data={**self.input, "condition": "MONDO:0007268"})
        self.assertEqual(value.evaluation_context["applicability"], "NOT_EVALUATED")
        self.assertEqual(value.review_points, [])

    def test_an_unqualified_x_linked_mode_is_compatible_with_either_qualification(self):
        """A source that records "X-linked" without saying which zygosity is affected has
        not contradicted a case assessed as X-linked recessive. ACMG asks for compatible
        inheritance, not identical inheritance."""
        for supplied in ("XLR", "x-linked dominant", "x linked"):
            with self.subTest(supplied=supplied):
                input_data = {**self.input, "inheritance": supplied}
                items = [self.annotation(), self.transcript(), self.nmd(),
                         self.mechanism(True, inheritance="x_linked")]
                value = self.evaluate_result(*items, input_data=input_data)
                self.assertEqual(value.status, CriterionStatus.MET)
                self.assertEqual(value.evaluation_context["moi_match"], "MATCHED")

    def test_autosomal_modes_never_stand_in_for_one_another(self):
        input_data = {**self.input, "inheritance": "AR"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, inheritance="autosomal_dominant")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.evaluation_context["moi_match"], "MISMATCH")

    def test_mechanism_curated_for_another_inheritance_mode_is_not_borrowed(self):
        input_data = {**self.input, "inheritance": "autosomal recessive"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, inheritance="AD")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.evaluation_context["moi_match"], "MISMATCH")
        self.assertEqual(value.evaluation_context["applicability"], "MANUAL_REVIEW")
        self.assertIn("loss-of-function disease mechanism for this inheritance mode",
                      value.missing_inputs)
        # The rejected record is still reported, so a curator can tell a mechanism that is
        # absent from one that exists under a different mode.
        self.assertIn("AD", value.summary)
        self.assertIn("'autosomal_recessive'", value.summary)
        self.assertTrue(any(item.get("inheritance") == "AD" for item in value.evidence))

    def test_mode_scoped_mechanism_is_not_used_when_the_case_declares_no_mode(self):
        """Not borrowed, and not buried either: an unrecorded mode is a field somebody can
        fill in, so the curation is named and the question goes to a person."""
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, inheritance="AR")]
        value = self.evaluate_result(*items)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIsNone(value.strength)
        self.assertEqual(value.evaluation_context["moi_match"], "UNSTATED")
        self.assertEqual(value.evaluation_context["applicability"], "MANUAL_REVIEW")
        self.assertIn("Record the inheritance mode", " ".join(value.review_points))
        self.assertEqual(value.provenance["preliminary_assessment"]["candidate_strength"],
                         "very_strong")

    def test_a_contradicted_mode_asks_a_different_question_from_an_unstated_one(self):
        """Both fail to match; a curator fixes them differently, so they are not merged."""
        input_data = {**self.input, "inheritance": "AR"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, inheritance="AD")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.evaluation_context["moi_match"], "MISMATCH")
        self.assertEqual(value.evaluation_context["applicability"], "MANUAL_REVIEW")
        self.assertNotIn("Record the inheritance mode", " ".join(value.review_points))
        self.assertIn("autosomal_recessive", " ".join(value.review_points))

    def test_unscoped_mechanism_stays_usable_under_any_declared_mode(self):
        input_data = {**self.input, "inheritance": "AD"}
        value = self.evaluate_result(*self.truncating_evidence(), input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["moi_match"], "NOT_SCOPED")

    def test_unrecognized_inheritance_mode_is_not_treated_as_a_match(self):
        input_data = {**self.input, "inheritance": "AR"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, inheritance="recessive-ish")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.evaluation_context["moi_match"], "MISMATCH")

    def derived(self, established, method, **values):
        """A mechanism record from an automated source, as a provider emits one."""
        values.setdefault("condition", self.condition)
        return self.item("gene_disease", gene="TEST", lof_mechanism_established=established,
                         assessment_method="automated", method=method,
                         policy_version="1", source=method, **values)

    def test_the_higher_ranked_source_decides_and_the_other_stays_visible(self):
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.derived(True, "better_source"), self.derived(False, "worse_source")]
        value = self.evaluate_result(*items)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["mechanism_source"], "better_source")
        self.assertNotIn("Conflicting", value.summary)
        self.assertEqual({item.get("method") for item in value.evidence
                          if item["category"] == "gene_disease"},
                         {"better_source", "worse_source"})

    def test_a_reviewed_record_outranks_every_derived_source(self):
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.derived(False, "better_source"), self.mechanism(True)]
        value = self.evaluate_result(*items)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["mechanism_source"], "synthetic")

    def test_an_unranked_source_sits_below_every_ranked_one(self):
        """Not knowing where a source belongs must not promote it above the ones we placed."""
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.derived(True, "worse_source"), self.derived(False, "unlisted_source")]
        value = self.evaluate_result(*items)
        self.assertEqual(value.status, CriterionStatus.MET)
        self.assertEqual(value.evaluation_context["mechanism_source"], "worse_source")

    def test_two_records_of_equal_rank_are_still_a_conflict(self):
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.derived(True, "worse_source", evidence_id="a"),
                 self.derived(False, "worse_source", evidence_id="b")]
        value = self.evaluate_result(*items)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("Conflicting", value.summary)

    def test_a_non_list_precedence_is_a_rejected_rule_set(self):
        config = {"PVS1": {**RULES["PVS1"], "rules": {
            "protein_loss_threshold": 0.10, "mechanism_source_precedence": "better_source"}}}
        value = self.evaluate_result(*self.truncating_evidence(), config=config)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("PVS1.ruleset", value.missing_inputs)

    def test_conflicting_mechanisms_and_gene_mismatch_require_review(self):
        value = self.evaluate_result(self.annotation(), self.mechanism(True), self.mechanism(False))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        mismatch = self.item("gene_disease", gene="OTHER", lof_mechanism_established=True)
        value = self.evaluate_result(self.annotation(), mismatch)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)

    def test_rule_set_is_required_for_an_eligible_variant(self):
        value = self.evaluate_result(self.annotation(), config={})
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.missing_inputs, ["PVS1.ruleset"])

    def test_transcript_and_exon_relevance_stop_without_negative_inference(self):
        items = [self.annotation(), self.mechanism(), self.transcript("NOT_RELEVANT"), self.nmd()]
        self.assertEqual(self.evaluate_result(*items).status, CriterionStatus.UNKNOWN)
        items = [self.annotation(), self.mechanism(), self.transcript(exon_relevance="UNKNOWN"), self.nmd()]
        value = self.evaluate_result(*items)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("exon_relevance", value.unresolved_requirements)

    def test_nmd_escape_critical_region_is_strong(self):
        items = [self.annotation(), self.mechanism(), self.transcript(), self.nmd(False),
                 self.region(critical_region_disrupted=True)]
        value = self.evaluate_result(*items)
        self.assertEqual((value.status, value.strength), (CriterionStatus.MET, "strong"))

    def test_nmd_escape_population_and_region_can_make_pvs1_not_applicable(self):
        base = [self.annotation(), self.mechanism(), self.transcript(), self.nmd(False)]
        value = self.evaluate_result(*base, self.region(critical_region_disrupted=False,
                                                       region_biologically_relevant=True),
                                     self.population_lof(True))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        value = self.evaluate_result(*base, self.region(critical_region_disrupted=False,
                                                       region_biologically_relevant=False),
                                     self.population_lof(False))
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)

    def test_protein_loss_threshold_is_strictly_greater_than_ten_percent(self):
        base = [self.annotation(), self.mechanism(), self.transcript(), self.nmd(False),
                self.population_lof(False)]
        moderate = self.evaluate_result(
            *base, self.region(critical_region_disrupted=False,
                               region_biologically_relevant=True,
                               lost_residues=10, total_protein_length=100))
        strong = self.evaluate_result(
            *base, self.region(critical_region_disrupted=False,
                               region_biologically_relevant=True,
                               lost_residues=11, total_protein_length=100))
        self.assertEqual((moderate.status, moderate.strength), (CriterionStatus.MET, "moderate"))
        self.assertEqual((strong.status, strong.strength), (CriterionStatus.MET, "strong"))

        invalid = self.evaluate_result(
            *base, self.region(critical_region_disrupted=False,
                               region_biologically_relevant=True,
                               lost_residues=101, total_protein_length=100))
        self.assertEqual(invalid.status, CriterionStatus.UNKNOWN)

    def splice(self, outcome="OUT_OF_FRAME_PTC", disrupted=True, rescue=False, **values):
        return self.item("splice_assessment", transcript="NM_TEST.1", splice_outcome=outcome,
                         reading_frame_disrupted=disrupted, alternative_rescue=rescue, **values)

    def test_canonical_splice_routes_to_nmd_and_in_frame_paths(self):
        out_of_frame = [self.annotation("splice_acceptor_variant"), self.mechanism(),
                        self.splice(), self.transcript(), self.nmd()]
        self.assertEqual(self.evaluate_result(*out_of_frame).strength, "very_strong")

        in_frame = [self.annotation("splice_donor_variant"), self.mechanism(),
                    self.splice("IN_FRAME", False), self.transcript(),
                    self.region(critical_region_disrupted=False,
                                region_biologically_relevant=True,
                                lost_residues=5, total_protein_length=100),
                    self.population_lof(False)]
        self.assertEqual(self.evaluate_result(*in_frame).strength, "moderate")

    def test_splice_rescue_uncertainty_and_frame_conflict_are_distinct(self):
        base = [self.annotation("splice_donor_variant"), self.mechanism()]
        self.assertEqual(self.evaluate_result(*base, self.splice(rescue=True)).status,
                         CriterionStatus.UNKNOWN)
        uncertain = self.splice(outcome="UNCERTAIN", disrupted=False)
        self.assertEqual(self.evaluate_result(*base, uncertain).status, CriterionStatus.UNKNOWN)
        conflict = self.splice(outcome="IN_FRAME", disrupted=True)
        self.assertEqual(self.evaluate_result(*base, conflict).status, CriterionStatus.UNKNOWN)

    def test_noncanonical_rna_confirmed_splice_lof_is_eligible_and_tracks_use(self):
        rna = self.item("rna_assay", transcript="NM_TEST.1", lof_effect_confirmed=True,
                        splice_outcome="OUT_OF_FRAME_PTC", reading_frame_disrupted=True,
                        alternative_rescue=False)
        items = [self.annotation("synonymous_variant"), self.mechanism(), rna,
                 self.transcript(), self.nmd()]
        value = self.evaluate_result(*items)
        self.assertEqual((value.status, value.strength), (CriterionStatus.MET, "very_strong"))
        used = next(item for item in value.evidence if item["category"] == "rna_assay")
        self.assertEqual(used["used_by"], ["PVS1"])
        self.assertEqual(next(node for node in value.decision_trace
                              if node["node_id"] == "V01")["value"], "SPLICE_LOF_CONFIRMED")

    def initiation(self, alternative=False, downstream=True, pathogenic=False, **values):
        return self.item("initiation_assessment", transcript="NM_TEST.1",
                         intact_alternative_transcript=alternative,
                         downstream_in_frame_start=downstream,
                         upstream_pathogenic_evidence=pathogenic, **values)

    def test_start_loss_supporting_moderate_and_stopping_branches(self):
        base = [self.annotation("start_lost"), self.mechanism(), self.transcript()]
        supporting = self.evaluate_result(*base, self.initiation(pathogenic=False))
        moderate = self.evaluate_result(*base, self.initiation(pathogenic=True))
        self.assertEqual((supporting.status, supporting.strength), (CriterionStatus.MET, "supporting"))
        self.assertEqual((moderate.status, moderate.strength), (CriterionStatus.MET, "moderate"))
        self.assertEqual(self.evaluate_result(*base, self.initiation(alternative=True)).status,
                         CriterionStatus.UNKNOWN)
        unknown = self.initiation()
        unknown["downstream_in_frame_start"] = None
        self.assertEqual(self.evaluate_result(*base, unknown).status, CriterionStatus.UNKNOWN)

    def test_other_lof_is_unresolved_but_unrelated_consequence_is_not_applicable(self):
        other_lof = self.evaluate_result(self.annotation("transcript_ablation"))
        unrelated = self.evaluate_result(self.annotation("missense_variant"))
        self.assertEqual(other_lof.status, CriterionStatus.UNKNOWN)
        self.assertEqual(unrelated.status, CriterionStatus.UNKNOWN)


if __name__ == "__main__":
    unittest.main()
