"""The records the automated providers emit have to satisfy PVS1's gates, not merely exist.

test_pvs1.py already covers the decision tree against hand-written evidence. What was never
checked is the seam between the two: whether ClinGen dosage, MANE and the NMD rule produce
records the tree actually accepts. Three separate failures lived in exactly that seam - a
gene-level evidence_id that collapsed every variant in one gene into one record, an exon the
annotation request never asked for, and two providers competing to supply one category - and
each looked fine from either side on its own.

So these tests drive real provider output, through the real evaluator, to a verdict.
"""

import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.automated_engine import evaluate_prepared_record, make_services
from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.providers.clingen_dosage import ClinGenDosageProvider
from acmg_pipeline.providers.mane import ManeTranscriptProvider
from acmg_pipeline.providers.nmd import NmdPredictionProvider

from tests.test_clingen_dosage import FakeClient as DosageClient, curation_list, row
from tests.test_mane import FakeClient as ManeClient, summary
from tests.test_nmd import FakeClient as VepClient, consequence

RULES = {"PVS1": {
    "ruleset": {
        "name": "clingen_general_pvs1",
        "version": "2018+2023-splicing",
        "sources": [{"id": "clingen_pvs1_2018", "pmid": "30192042"}],
    },
    "rules": {"protein_loss_threshold": 0.10},
}}

GENE = "MYBPC3"
TRANSCRIPT = "NM_000256.3"
HGVSC = "NM_000256.3:c.278delA"
CONDITION = "MONDO:0005045"


class Pvs1AutomatedGateTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "11", 47351252, "CT", "C")
        self.input = {"variant": self.variant.to_dict(), "transcript": TRANSCRIPT,
                      "condition": CONDITION}

    def annotation(self, consequence_term="frameshift_variant"):
        """Shaped like the record the Ensembl provider emits - note it has no exon."""
        return {
            "category": "annotation", "variant_key": self.variant.key,
            "evidence_id": f"test:annotation:{self.variant.key}",
            "source": "Ensembl VEP HGVS", "source_version": "116",
            "retrieved_at": "2026-09-17", "quality_status": "PASS",
            "transcript": TRANSCRIPT, "gene": GENE, "hgvsc": HGVSC,
            "consequences": [consequence_term], "protein_id": "NP_000247.2",
            "exon": None,
        }

    def dosage(self, score="3", gene=GENE):
        provider = ClinGenDosageProvider(DosageClient(curation_list([row(gene, score)])))
        return provider.get_mechanism(self.variant, gene)

    def curated_mechanism(self, established=True):
        """A disease-scoped mechanism, as an expert-panel specification supplies one.

        ClinGen dosage sensitivity scores a gene, not a gene-and-disease, so its records can
        no longer carry PVS1 past the disease gate on their own. The seams these tests are
        about are the MANE, NMD, region and initiation providers, so the mechanism is
        supplied here in the scoped form the gate requires, and the dosage provider has its
        own tests for what it can and cannot settle.
        """
        return [{
            "category": "gene_disease", "variant_key": self.variant.key,
            "evidence_id": f"curated:mechanism:{self.variant.key}",
            "source": "curator", "source_version": "1", "retrieved_at": "2026-09-17",
            "quality_status": "PASS", "curator": "test", "reviewed_at": "2026-09-17",
            "gene": GENE, "condition": CONDITION,
            "lof_mechanism_established": established,
        }]

    def mane(self, exon="2/34", accession=TRANSCRIPT):
        provider = ManeTranscriptProvider.from_directory(
            ManeClient(text=summary([(GENE, accession, "MANE Select")])))
        return provider.get_transcript_assessment(self.variant, GENE, TRANSCRIPT, exon=exon)

    def nmd(self, exon="2/34"):
        provider = NmdPredictionProvider(VepClient([consequence(exon=exon)]), "116")
        return provider.get_nmd_prediction(self.variant, GENE, TRANSCRIPT, HGVSC)

    def exon_from_vep(self, exon="2/34"):
        provider = NmdPredictionProvider(VepClient([consequence(exon=exon)]), "116")
        return provider.exon_on_transcript(GENE, TRANSCRIPT, HGVSC)

    def evaluate(self, *record_groups, consequence_term="frameshift_variant"):
        records = [self.annotation(consequence_term)]
        for group in record_groups:
            records.extend(group)
        return evaluate_prepared_record(self.input, make_services(records), RULES, ["PVS1"])[0]

    def nodes(self, result):
        return {node["node_id"]: node["result"] for node in result.decision_trace}

    def test_the_three_providers_together_carry_pvs1_to_very_strong(self):
        exon = self.exon_from_vep()
        result = self.evaluate(self.curated_mechanism(), self.mane(exon=exon), self.nmd())
        self.assertEqual(result.status, CriterionStatus.MET)
        self.assertEqual(result.strength, "very_strong")
        self.assertEqual(
            self.nodes(result),
            {"C01": "PASS", "D01": "PASS", "G01": "PASS", "G02": "PASS", "V01": "PASS",
             "NF01": "PASS", "NF02": "PASS", "NF03": "PASS"},
        )

    def test_dosage_sensitivity_alone_no_longer_reaches_a_verdict(self):
        """Dosage scores a gene, and PVS1 needs a mechanism for the disease being assessed.

        The finding is not discarded: the record is attached and the variant-level tree still
        runs, so a curator sees both the haploinsufficiency score and what PVS1 would have
        concluded once a disease-scoped mechanism is supplied.
        """
        exon = self.exon_from_vep()
        result = self.evaluate(self.dosage(), self.mane(exon=exon), self.nmd())
        self.assertEqual(result.status, CriterionStatus.UNKNOWN)
        self.assertEqual(result.evaluation_context["applicability"], "NOT_EVALUATED")
        self.assertEqual(self.nodes(result)["D01"], "UNKNOWN")
        self.assertIn("disease-specific loss-of-function mechanism", result.missing_inputs)
        self.assertTrue(any(item["source"].startswith("ClinGen") for item in result.evidence))
        preliminary = result.provenance["preliminary_assessment"]
        self.assertEqual(preliminary["candidate_strength"], "very_strong")
        self.assertEqual(preliminary["decision_path"], "NF03")

    def test_one_provider_supplies_the_transcript_assessment(self):
        """Two records in that category make _select_context_record() report a conflict, so
        the exon reaches NF03 through the MANE record rather than a second one."""
        records = self.mane() + self.nmd()
        categories = [item["category"] for item in records]
        self.assertEqual(categories.count("transcript_assessment"), 1)

    def test_each_variant_in_a_gene_gets_its_own_mechanism_record(self):
        """A gene-level evidence_id collapsed these into one, and the CLI deduplicates by
        that id, so every variant but the first lost its mechanism and stopped at G01."""
        other = Variant("GRCh38", "11", 47335041, "C", "T")
        provider = ClinGenDosageProvider(DosageClient(curation_list([row(GENE, "3")])))
        ids = {provider.get_mechanism(v, GENE)[0]["evidence_id"]
               for v in (self.variant, other)}
        self.assertEqual(len(ids), 2)

    def test_without_the_dosage_record_the_tree_stops_at_g01(self):
        result = self.evaluate(self.mane(), self.nmd())
        self.assertEqual(result.status, CriterionStatus.UNKNOWN)
        self.assertEqual(self.nodes(result)["G01"], "UNKNOWN")
        self.assertIn("loss-of-function disease mechanism", result.missing_inputs)

    def test_a_recessive_gene_stops_at_g01_rather_than_being_denied(self):
        """Score 30 emits nothing, so PVS1 asks for a mechanism instead of ruling it out -
        the PAH case, where the VCEP did apply PVS1."""
        result = self.evaluate(self.dosage(score="30"), self.mane(), self.nmd())
        self.assertEqual(self.nodes(result)["G01"], "UNKNOWN")

    def test_a_gene_scored_against_dosage_sensitivity_stops_before_the_mechanism_gate(self):
        """Score 40 says loss of function is not this *gene's* mechanism, which is a weaker
        statement than the gate needs, so the path stops at the disease match rather than at
        G02. PVS1 is withheld either way; only the reason a curator is given changes."""
        result = self.evaluate(self.dosage(score="40"), self.mane(), self.nmd())
        self.assertEqual(self.nodes(result)["D01"], "UNKNOWN")
        self.assertNotIn("G02", self.nodes(result))
        self.assertEqual(result.status, CriterionStatus.UNKNOWN)

    def test_a_scoped_mechanism_that_denies_loss_of_function_is_read_at_g02(self):
        result = self.evaluate(self.curated_mechanism(established=False), self.mane(),
                               self.nmd())
        self.assertEqual(self.nodes(result)["G02"], "NOT_APPLICABLE")
        self.assertEqual(result.evaluation_context["applicability"], "NOT_APPLICABLE")
        self.assertEqual(result.status, CriterionStatus.UNKNOWN)

    def test_without_the_mane_record_the_tree_stops_at_nf01(self):
        result = self.evaluate(self.curated_mechanism(), self.nmd())
        self.assertEqual(self.nodes(result)["NF01"], "UNKNOWN")
        self.assertIn("transcript_assessment", result.missing_inputs)

    def test_without_the_nmd_record_the_tree_stops_at_nf02(self):
        result = self.evaluate(self.curated_mechanism(), self.mane())
        self.assertEqual(self.nodes(result)["NF02"], "UNKNOWN")
        self.assertIn("nmd_prediction", result.missing_inputs)

    def test_without_an_exon_the_tree_stops_at_nf03_undenied(self):
        """NF03 must read this as unresolved, not NOT_RELEVANT: the exon being unnumbered
        says nothing about whether it is in the transcript."""
        result = self.evaluate(self.curated_mechanism(), self.mane(exon=None), self.nmd())
        self.assertEqual(self.nodes(result)["NF03"], "UNKNOWN")
        self.assertEqual(result.status, CriterionStatus.UNKNOWN)
        self.assertIn("exon_relevance", result.missing_inputs)

    def test_the_penultimate_exon_stops_at_nf02_with_the_exon_known(self):
        """The NMD rule declines the boundary case while NF03's question stays answerable."""
        exon = self.exon_from_vep("33/34")
        self.assertEqual(exon, "33/34")
        result = self.evaluate(self.curated_mechanism(), self.mane(exon=exon), self.nmd("33/34"))
        self.assertEqual(self.nodes(result)["NF02"], "UNKNOWN")

    def test_a_ptc_in_the_last_exon_takes_the_region_path_and_stops_at_nf06(self):
        """NMD is escaped, so PVS1 weighs the region instead - and the two gates that ask
        what the lost residues *do* are judgments no coordinate answers."""
        exon = self.exon_from_vep("34/34")
        result = self.evaluate(self.curated_mechanism(), self.mane(exon=exon), self.nmd("34/34"))
        nodes = self.nodes(result)
        self.assertEqual(nodes["NF02"], "PASS")
        self.assertEqual(nodes["NF04"], "UNKNOWN")   # critical region - curator
        self.assertEqual(nodes["NF05"], "PASS")      # no frequent LoF observed
        self.assertEqual(nodes["NF06"], "UNKNOWN")   # biological relevance - curator
        self.assertNotIn("NF07", nodes)
        self.assertIn("region_biological_relevance", result.missing_inputs)

    def test_a_non_mane_transcript_stops_at_nf01_undenied(self):
        result = self.evaluate(self.curated_mechanism(), self.mane(accession="NM_999999.1"), self.nmd())
        self.assertEqual(self.nodes(result)["NF01"], "UNKNOWN")

    def test_a_stop_gained_variant_takes_the_same_route(self):
        exon = self.exon_from_vep()
        result = self.evaluate(self.curated_mechanism(), self.mane(exon=exon), self.nmd(),
                               consequence_term="stop_gained")
        self.assertEqual(result.status, CriterionStatus.MET)
        self.assertEqual(self.nodes(result)["V01"], "PASS")

    def region(self, protein_start=204, length=213):
        """The measurement provider's record, with NF04/NF06 left unanswered as it does."""
        from tests.test_protein_region import FakeClient as RegionClient
        from acmg_pipeline.providers.protein_region import ProteinRegionProvider
        client = RegionClient([consequence(exon="3/3")], length=length)
        return ProteinRegionProvider(client, "116").get_protein_region(
            self.variant, GENE, TRANSCRIPT, HGVSC, protein_start)

    def relevance(self, **values):
        """The one judgment NF06 needs, as a curator would supply it."""
        return [{
            "category": "protein_region", "variant_key": self.variant.key,
            "evidence_id": f"curated:region:{self.variant.key}",
            "source": "curator", "source_version": "1", "retrieved_at": "2026-09-17",
            "quality_status": "PASS", "curator": "test", "reviewed_at": "2026-09-17",
            "transcript": TRANSCRIPT, "region_biologically_relevant": True, **values,
        }]

    def test_the_measurement_alone_still_stops_at_the_curator_judgment(self):
        result = self.evaluate(self.curated_mechanism(), self.mane(exon="3/3"), self.nmd("3/3"),
                               self.region())
        nodes = self.nodes(result)
        self.assertEqual(nodes["NF02"], "PASS")
        self.assertEqual(nodes["NF06"], "UNKNOWN")
        self.assertNotIn("NF07", nodes)
        self.assertIn("region_biological_relevance", result.missing_inputs)

    def test_measurement_plus_the_judgment_reaches_nf07_and_sets_the_strength(self):
        """VHL c.610G>T: 10 residues of 213 is 4.7%, under the 10% the rule set splits on,
        which is the PVS1_Moderate its expert panel recorded."""
        curated = self.relevance(lost_residues=10, total_protein_length=213)
        result = self.evaluate(self.curated_mechanism(), self.mane(exon="3/3"), self.nmd("3/3"), curated)
        self.assertEqual(result.status, CriterionStatus.MET)
        self.assertEqual(result.strength, "moderate")
        self.assertEqual(self.nodes(result)["NF07"], "PASS")

    def test_losing_more_than_a_tenth_of_the_protein_is_strong(self):
        curated = self.relevance(lost_residues=100, total_protein_length=213)
        result = self.evaluate(self.curated_mechanism(), self.mane(exon="3/3"), self.nmd("3/3"), curated)
        self.assertEqual(result.strength, "strong")

    def initiation(self, cds=None):
        """The IC02 answer, with IC01 and IC03 left unset as the provider leaves them."""
        from tests.test_initiation import CDS_WITH_RESTART, FakeClient as InitClient
        from tests.test_initiation import consequence as start_lost_consequence
        from acmg_pipeline.providers.initiation import InitiationProvider
        client = InitClient([start_lost_consequence(gene_symbol=GENE,
                                                    mane_select=TRANSCRIPT)],
                            cds or CDS_WITH_RESTART)
        return InitiationProvider(client, "116").get_initiation_assessment(
            self.variant, GENE, TRANSCRIPT, HGVSC)

    def curated_initiation(self, **values):
        """The judgments IC01 and IC03 need, as a curator would supply them."""
        return [{
            "category": "initiation_assessment", "variant_key": self.variant.key,
            "evidence_id": f"curated:initiation:{self.variant.key}",
            "source": "curator", "source_version": "1", "retrieved_at": "2026-09-17",
            "quality_status": "PASS", "curator": "test", "reviewed_at": "2026-09-17",
            "transcript": TRANSCRIPT, "intact_alternative_transcript": False,
            "downstream_in_frame_start": True, **values,
        }]

    def test_the_downstream_start_alone_still_stops_at_the_curator_judgment(self):
        """IC01 is read first, so answering IC02 does not advance the path by itself."""
        result = self.evaluate(self.curated_mechanism(), self.mane(exon="1/2"), self.initiation(),
                               consequence_term="start_lost")
        nodes = self.nodes(result)
        self.assertEqual(nodes["V01"], "PASS")
        self.assertEqual(nodes["IC00"], "PASS")
        self.assertEqual(nodes["IC01"], "UNKNOWN")
        self.assertIn("intact_alternative_transcript", result.missing_inputs)

    def test_an_intact_alternative_transcript_makes_pvs1_inapplicable(self):
        result = self.evaluate(
            self.curated_mechanism(), self.mane(exon="1/2"),
            self.curated_initiation(intact_alternative_transcript=True),
            consequence_term="start_lost")
        self.assertEqual(self.nodes(result)["IC01"], "NOT_APPLICABLE")
        self.assertEqual(result.status, CriterionStatus.UNKNOWN)

    def test_upstream_pathogenic_evidence_decides_moderate_against_supporting(self):
        for pathogenic, strength in ((True, "moderate"), (False, "supporting")):
            with self.subTest(upstream_pathogenic_evidence=pathogenic):
                result = self.evaluate(
                    self.curated_mechanism(), self.mane(exon="1/2"),
                    self.curated_initiation(upstream_pathogenic_evidence=pathogenic),
                    consequence_term="start_lost")
                self.assertEqual(result.status, CriterionStatus.MET)
                self.assertEqual(result.strength, strength)
                self.assertEqual(self.nodes(result)["IC03"], "PASS")

    def initiation_with_upstream(self, upstream=True, **overrides):
        """IC02's record with IC03's answer merged in, as the CLI assembles it."""
        from tests.test_upstream_pathogenic import FakeClient as ClinVarClient, document
        from acmg_pipeline.providers.upstream_pathogenic import UpstreamPathogenicProvider
        records = self.initiation()
        codon = records[0]["downstream_start_codon"]
        change = f"R{codon - 1}W" if upstream else f"L{codon + 5}P"
        reports = [document(1, change, gene=GENE)]
        fields = UpstreamPathogenicProvider(
            ClinVarClient(reports), "2026-09-15").get_upstream_evidence(
                self.variant, GENE, TRANSCRIPT, codon)
        records[0].update(fields)
        records[0].update(overrides)
        return records

    def test_both_derived_answers_live_on_one_initiation_record(self):
        """PVS1 selects one initiation_assessment and reads two as a conflict, so IC03 is
        merged into IC02's record rather than emitted beside it."""
        records = self.initiation_with_upstream()
        self.assertEqual(len(records), 1)
        self.assertIs(records[0]["downstream_in_frame_start"], True)
        self.assertIs(records[0]["upstream_pathogenic_evidence"], True)
        result = self.evaluate(self.curated_mechanism(), self.mane(exon="1/2"), records,
                               consequence_term="start_lost")
        self.assertNotIn("Conflicting", result.summary)

    def test_the_derived_answers_carry_the_path_to_ic03(self):
        """Only IC01 is left, and it is the judgment - the two derivable gates are done."""
        result = self.evaluate(self.curated_mechanism(), self.mane(exon="1/2"),
                               self.initiation_with_upstream(),
                               consequence_term="start_lost")
        self.assertEqual(self.nodes(result)["IC01"], "UNKNOWN")
        self.assertEqual(result.missing_inputs, ["intact_alternative_transcript"])

    def test_with_ic01_supplied_the_derived_answers_reach_met(self):
        for upstream, strength in ((True, "moderate"), (False, "supporting")):
            with self.subTest(upstream_pathogenic_evidence=upstream):
                records = self.initiation_with_upstream(
                    upstream=upstream, intact_alternative_transcript=False,
                    curator="test", reviewed_at="2026-09-17")
                result = self.evaluate(self.curated_mechanism(), self.mane(exon="1/2"), records,
                                       consequence_term="start_lost")
                self.assertEqual(result.status, CriterionStatus.MET)
                self.assertEqual(result.strength, strength)
                self.assertEqual(self.nodes(result)["IC02"], "PASS")
                self.assertEqual(self.nodes(result)["IC03"], "PASS")

    def test_automated_records_are_marked_as_such_in_the_evidence(self):
        """A curator reading the result has to see which gates were answered by derivation."""
        result = self.evaluate(self.curated_mechanism(), self.mane(exon="2/34"), self.nmd())
        derived = [item for item in result.evidence
                   if item.get("assessment_method") == "automated"]
        # The mechanism is curated here, because a derived one no longer passes the disease
        # gate; the two gates that are still answered by derivation have to say so.
        self.assertEqual(
            {item["category"] for item in derived},
            {"transcript_assessment", "nmd_prediction"},
        )
        for item in derived:
            self.assertTrue(item["method"] and item["policy_version"])


if __name__ == "__main__":
    unittest.main()
