import copy
import unittest

from acmg.core.models import CRITERIA, Status, Variant
from acmg.engine import evaluate_record, make_services


class CuratedCriteriaTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")
        self.input = {"variant": self.variant.to_dict(), "transcript": "NM_TEST.1", "condition": "test:disease"}
        self.base = {"variant_key": self.variant.key, "transcript": "NM_TEST.1", "condition": "test:disease",
                     "source": "synthetic", "source_version": "1", "retrieved_at": "2026-09-15",
                     "quality_status": "PASS", "curator": "test", "reviewed_at": "2026-09-15"}
        self.annotation = self.item("annotation", consequences=["missense_variant"], gene="TEST",
                                    protein_id="NP_TEST.1", protein_start=10, protein_end=10,
                                    ref_aa="R", alt_aa="W")

    def item(self, category, **values):
        return {**self.base, "category": category, "evidence_id": f"test:{category}", **values}

    def run_rule(self, code, *items):
        return evaluate_record(self.input, make_services([self.annotation, *items]), {}, [code])[0]

    def test_all_sixteen_without_evidence(self):
        values = evaluate_record(self.input, make_services([]), {})
        self.assertEqual([v.criterion for v in values], list(CRITERIA))
        self.assertEqual(sum(v.status == Status.DEPRECATED for v in values), 2)
        self.assertTrue(all(v.status in {Status.NOT_EVALUATED, Status.DEPRECATED} for v in values))

    def test_mechanism_is_disease_specific(self):
        item = self.item("gene_disease", gene="TEST", missense_mechanism_established=True,
                         spectrum_review_complete=True, low_benign_missense_variation=True,
                         predominantly_truncating=False)
        self.assertEqual(self.run_rule("PP2", item).status, Status.MET)
        self.assertEqual(self.run_rule("BP1", item).status, Status.NOT_MET)
        item["condition"] = "test:other-disease"
        self.assertEqual(self.run_rule("PP2", item).status, Status.NOT_EVALUATED)

    def test_bp1_requires_reviewed_spectrum(self):
        item = self.item("gene_disease", gene="TEST", missense_mechanism_established=False,
                         spectrum_review_complete=True, predominantly_truncating=True)
        self.assertEqual(self.run_rule("BP1", item).status, Status.MET)
        del item["spectrum_review_complete"]
        self.assertEqual(self.run_rule("BP1", item).status, Status.NOT_EVALUATED)

    def region(self, **values):
        return self.item("region", protein_id="NP_TEST.1", start=5, end=20, **values)

    HOTSPOT_POLICY = {"PM1": {"hotspot": {
        "window_aa": 5, "min_pathogenic": 3, "max_benign": 0,
        "method": "clinvar_local_density",
        "policy_source": "test policy", "policy_version": "PM1-hotspot-test"}}}

    def hotspot(self, **values):
        return self.region(region_type="mutational_hotspot", method="clinvar_local_density",
                           policy_version="PM1-hotspot-test", **values)

    def run_pm1(self, item, config=None, input_data=None, annotation=None):
        services = make_services([annotation or self.annotation, item])
        return evaluate_record(input_data or self.input, services, config or {}, ["PM1"])[0]

    def test_pm1_requires_a_declared_route(self):
        value = self.run_pm1(self.region(critical_functional_region=True, benign_depletion=True))
        self.assertEqual(value.status, Status.NOT_EVALUATED)
        self.assertIn("region_type", value.missing_inputs)

    def test_domain_overlap_alone_is_insufficient(self):
        item = self.region(region_type="critical_functional_domain", critical_functional_region=True)
        self.assertEqual(self.run_pm1(item).status, Status.NOT_EVALUATED)
        item["benign_depletion"] = True
        self.assertEqual(self.run_pm1(item).status, Status.MET)
        item["benign_depletion"] = False
        self.assertEqual(self.run_pm1(item).status, Status.NOT_MET)

    def test_critical_domain_does_not_require_pathogenic_enrichment(self):
        """A well-established active site stays PM1 even with few reported cases."""
        item = self.region(region_type="critical_functional_domain", critical_functional_region=True,
                           benign_depletion=True, pathogenic_enrichment=False)
        value = self.run_pm1(item)
        self.assertEqual(value.status, Status.MET)
        self.assertEqual(value.provenance["pm1_route"], "critical_functional_domain")

    def test_automated_evidence_cannot_claim_critical_domain(self):
        item = self.region(region_type="critical_functional_domain", assessment_method="automated",
                           method="clinvar_local_density", policy_version="PM1-hotspot-test",
                           critical_functional_region=True, benign_depletion=True)
        value = self.run_pm1(item)
        self.assertEqual(value.status, Status.NOT_EVALUATED)
        self.assertIn("curated_criticality", value.missing_inputs)

    def test_hotspot_route_accepts_automated_policy_evidence(self):
        item = self.hotspot(assessment_method="automated", pathogenic_count=3, benign_count=0)
        del item["curator"]
        del item["reviewed_at"]
        value = self.run_pm1(item, self.HOTSPOT_POLICY)
        self.assertEqual(value.status, Status.MET)
        self.assertEqual(value.provenance["assessment_method"], "automated")

    def test_hotspot_route_uses_configured_thresholds(self):
        item = self.hotspot(pathogenic_count=3, benign_count=0)
        value = self.run_pm1(item, self.HOTSPOT_POLICY)
        self.assertEqual(value.status, Status.MET)
        self.assertEqual(value.provenance["pm1_route"], "mutational_hotspot")
        self.assertEqual(value.provenance["pathogenic_count"], 3)

    def test_hotspot_benign_variation_is_a_real_negative(self):
        item = self.hotspot(pathogenic_count=4, benign_count=1)
        self.assertEqual(self.run_pm1(item, self.HOTSPOT_POLICY).status, Status.NOT_MET)

    def test_sparse_pathogenic_reports_are_missing_evidence(self):
        """Few ClinVar reports must not be read as proof that no hotspot exists."""
        item = self.hotspot(pathogenic_count=2, benign_count=0)
        value = self.run_pm1(item, self.HOTSPOT_POLICY)
        self.assertEqual(value.status, Status.NOT_EVALUATED)
        self.assertIn("pathogenic_density", value.missing_inputs)

    def test_hotspot_requires_versioned_policy(self):
        item = self.hotspot(pathogenic_count=4, benign_count=0)
        value = self.run_pm1(item)
        self.assertEqual(value.status, Status.NOT_EVALUATED)
        self.assertIn("PM1.hotspot.policy_version", value.missing_inputs)
        item["policy_version"] = "PM1-hotspot-other"
        value = self.run_pm1(item, self.HOTSPOT_POLICY)
        self.assertEqual(value.status, Status.NOT_EVALUATED)
        self.assertIn("policy_version", value.missing_inputs)

    def test_hotspot_counts_are_required(self):
        value = self.run_pm1(self.hotspot(mutational_hotspot=True), self.HOTSPOT_POLICY)
        self.assertEqual(value.status, Status.NOT_EVALUATED)
        self.assertIn("pathogenic_count", value.missing_inputs)

    def test_hotspot_assertion_must_match_its_counts(self):
        item = self.hotspot(pathogenic_count=4, benign_count=1, benign_depletion=True)
        value = self.run_pm1(item, self.HOTSPOT_POLICY)
        self.assertEqual(value.status, Status.MANUAL_REVIEW)
        self.assertTrue(value.review_points)

    def test_pm1_protein_level_without_condition(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        annotation = {key: value for key, value in self.annotation.items() if key != "condition"}
        item = {key: value for key, value in
                self.region(region_type="critical_functional_domain", critical_functional_region=True,
                            benign_depletion=True).items() if key != "condition"}
        value = self.run_pm1(item, input_data=input_data, annotation=annotation)
        self.assertEqual(value.status, Status.MET)
        self.assertEqual(value.provenance["assessment_scope"], "protein_level")
        self.assertEqual(value.provenance["condition_assessment"], "NOT_EVALUATED")
        self.assertTrue(value.review_points)
        item["benign_depletion"] = False
        value = self.run_pm1(item, input_data=input_data, annotation=annotation)
        self.assertEqual(value.status, Status.NOT_MET)
        self.assertEqual(value.review_points, [])

    def test_pm1_records_matched_condition(self):
        item = self.region(region_type="critical_functional_domain", critical_functional_region=True,
                           benign_depletion=True)
        value = self.run_pm1(item)
        self.assertEqual(value.status, Status.MET)
        self.assertEqual(value.provenance["condition_assessment"], "MATCHED")
        self.assertEqual(value.review_points, [])

    def test_pm1_does_not_borrow_other_disease_assessment(self):
        item = self.region(region_type="critical_functional_domain", critical_functional_region=True,
                           benign_depletion=True, condition="test:other-disease")
        value = self.run_pm1(item)
        self.assertEqual(value.status, Status.NOT_EVALUATED)
        self.assertIn("region", value.missing_inputs)

    def test_length_repeat_and_missing_function(self):
        self.annotation.update(consequences=["inframe_deletion"], protein_length_change=-1)
        item = self.region(nonfunctional_repeat=False, repetitive=False, functional_importance=True,
                           functional_review_complete=True)
        self.assertEqual(self.run_rule("PM4", item).status, Status.MET)
        self.assertEqual(self.run_rule("BP3", item).status, Status.NOT_MET)
        item.update(nonfunctional_repeat=True, repetitive=True, functional_importance=False)
        self.assertEqual(self.run_rule("PM4", item).status, Status.NOT_MET)
        self.assertEqual(self.run_rule("BP3", item).status, Status.MET)
        del item["functional_importance"]
        self.assertEqual(self.run_rule("BP3", item).status, Status.NOT_EVALUATED)

    def comparator(self, **values):
        return self.item("comparator", protein_id="NP_TEST.1", protein_start=10, ref_aa="R", alt_aa="W",
                         comparator_variant={**self.variant.to_dict(), "pos": 3, "ref": "G", "alt": "A"},
                         classification="Pathogenic", pathogenic_evidence_reviewed=True, independent_evidence=True,
                         mechanism_matches=True, splice_effect_checked=True, different_splice_mechanism=False,
                         primary_evidence=["test:primary-study"], **values)

    def test_same_vs_different_amino_acid(self):
        item = self.comparator()
        self.assertEqual(self.run_rule("PS1", item).status, Status.MET)
        self.assertEqual(self.run_rule("PM5", item).status, Status.NOT_EVALUATED)
        item["alt_aa"] = "Q"
        self.assertEqual(self.run_rule("PM5", item).status, Status.MET)

    def test_comparator_label_is_not_enough(self):
        item = self.comparator()
        item["independent_evidence"] = False
        self.assertEqual(self.run_rule("PS1", item).status, Status.MANUAL_REVIEW)
        item["independent_evidence"] = True
        item["comparator_variant"] = self.variant.to_dict()
        self.assertNotEqual(self.run_rule("PS1", item).status, Status.MET)

    def test_ps1_does_not_require_condition(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        annotation = {key: value for key, value in self.annotation.items() if key != "condition"}
        value = evaluate_record(input_data, make_services([annotation]), {}, ["PS1"])[0]
        self.assertEqual(value.missing_inputs, ["complete_comparator_search"])

    def test_automated_ps1_protein_match_can_be_met_without_condition(self):
        input_data = {key: value for key, value in self.input.items() if key != "condition"}
        annotation = {key: value for key, value in self.annotation.items() if key != "condition"}
        comparator = {key: value for key, value in self.comparator().items()
                      if key != "condition"}
        comparator.update(exact_protein_match=True, different_nucleotide_variant=True,
                          review_status_eligible=True, splice_effect_checked=True,
                          splice_conflict=False, conditions=["MONDO:0000001"])
        value = evaluate_record(input_data, make_services([annotation, comparator]), {}, ["PS1"])[0]
        self.assertEqual(value.status, Status.MET)
        self.assertEqual(value.provenance["assessment_scope"], "protein_level")
        self.assertEqual(value.provenance["condition_assessment"], "NOT_EVALUATED")
        self.assertTrue(value.review_points)

    def test_search_absence_requires_completeness(self):
        search = self.item("comparator_search", protein_id="NP_TEST.1", protein_start=10, complete=True)
        self.assertEqual(self.run_rule("PM5", search).status, Status.NOT_MET)
        search["complete"] = False
        self.assertEqual(self.run_rule("PM5", search).status, Status.NOT_EVALUATED)

    def test_bp7_and_rna_contradiction(self):
        self.annotation["consequences"] = ["synonymous_variant"]
        item = self.item("synonymous_assessment", outside_splice_critical_region=True,
                         no_predicted_splice_impact=True, not_conserved=True, contradictory_rna_evidence=False,
                         splice_prediction_evidence="test:splice", calibration_source="test:calibration",
                         conservation_evidence="test:conservation", position_rule_version="test:1")
        self.assertEqual(self.run_rule("BP7", item).status, Status.MET)
        item["contradictory_rna_evidence"] = True
        self.assertEqual(self.run_rule("BP7", item).status, Status.MANUAL_REVIEW)

    def test_pvs1_remains_provisional(self):
        self.annotation["consequences"] = ["frameshift_variant"]
        mechanism = self.item("gene_disease", gene="TEST", lof_mechanism_established=True)
        lof = self.item("lof_assessment", biologically_relevant_transcript=True, exon_relevant=True,
                        nmd_predicted=True, nmd_rule_source="test:nmd")
        value = self.run_rule("PVS1", mechanism, lof)
        self.assertEqual(value.status, Status.MANUAL_REVIEW)
        self.assertIsNone(value.strength)
        self.assertEqual(value.provenance["recommended_strength"], "very_strong")
        lof["nmd_predicted"] = False
        self.assertIsNone(self.run_rule("PVS1", mechanism, lof).provenance["recommended_strength"])
        mechanism["lof_mechanism_established"] = False
        self.assertEqual(self.run_rule("PVS1", mechanism, lof).status, Status.NOT_MET)

    def test_source_labels_do_not_affect_all_results(self):
        services = make_services([self.annotation])
        first = [r.to_dict() for r in evaluate_record(self.input, services, {})]
        altered = copy.deepcopy(self.input)
        altered.update(CLNSIG="Pathogenic", ACMG_CODES="PVS1,PS1,PM2", NOTE="ground truth")
        self.assertEqual(first, [r.to_dict() for r in evaluate_record(altered, services, {})])
