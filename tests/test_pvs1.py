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
    "rules": {"protein_loss_threshold": 0.10},
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
        self.assertEqual(value.evaluation_context["applicability"], "NOT_EVALUATED")
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
        self.assertIn("loss-of-function disease mechanism", value.unresolved_requirements)

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

    def test_mechanism_curated_for_another_inheritance_mode_is_not_borrowed(self):
        input_data = {**self.input, "inheritance": "autosomal recessive"}
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, inheritance="AD")]
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.evaluation_context["moi_match"], "MISMATCH")
        self.assertIn("loss-of-function disease mechanism for this inheritance mode",
                      value.missing_inputs)
        # The rejected record is still reported, so a curator can tell a mechanism that is
        # absent from one that exists under a different mode.
        self.assertIn("AD", value.summary)
        self.assertIn("'autosomal_recessive'", value.summary)
        self.assertTrue(any(item.get("inheritance") == "AD" for item in value.evidence))

    def test_mode_scoped_mechanism_is_not_used_when_the_case_declares_no_mode(self):
        items = [self.annotation(), self.transcript(), self.nmd(),
                 self.mechanism(True, inheritance="AR")]
        value = self.evaluate_result(*items)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertEqual(value.evaluation_context["moi_match"], "MISMATCH")

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
