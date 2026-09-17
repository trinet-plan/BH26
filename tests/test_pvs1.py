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
        self.input = {"variant": self.variant.to_dict(), "transcript": "NM_TEST.1"}
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

    def test_condition_missing_uses_gene_level_mechanism_and_can_be_met(self):
        value = self.evaluate_result(*self.truncating_evidence())
        self.assertEqual((value.status, value.strength), (CriterionStatus.MET, "very_strong"))
        self.assertEqual(value.evaluation_context["condition_status"], "NOT_PROVIDED")
        self.assertEqual(value.evaluation_context["mechanism_scope"], "GENE_LEVEL")
        self.assertFalse(value.evaluation_context["condition_specific"])
        self.assertTrue(value.warnings)
        self.assertEqual(value.missing_inputs, value.unresolved_requirements)
        self.assertEqual([node["node_id"] for node in value.decision_trace],
                         ["C01", "G01", "G02", "V01", "NF01", "NF02", "NF03"])
        self.assertEqual(len(value.rules_used), 3)

    def test_condition_specific_mechanism_precedes_gene_level_fallback(self):
        input_data = {**self.input, "condition": "MONDO:1", "condition_label": "Disease"}
        items = self.truncating_evidence()
        items.append(self.mechanism(False, condition="MONDO:1"))
        value = self.evaluate_result(*items, input_data=input_data)
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertTrue(value.evaluation_context["condition_specific"])
        self.assertEqual(value.evaluation_context["mechanism_scope"], "CONDITION_SPECIFIC")

        fallback = self.evaluate_result(*self.truncating_evidence(), input_data=input_data)
        self.assertEqual(fallback.status, CriterionStatus.MET)
        self.assertEqual(fallback.evaluation_context["condition_status"], "PROVIDED")
        self.assertFalse(fallback.evaluation_context["condition_specific"])
        self.assertEqual(fallback.evaluation_context["mechanism_scope"], "GENE_LEVEL")

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
