from acmg_pipeline.constants import CriterionStatus
import json
import unittest
from pathlib import Path

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.automated_engine import evaluate_prepared_record, make_services
from acmg_pipeline.gene_disease import build_assessment_document, build_draft_document
from acmg_pipeline.providers.gene_disease_draft import GeneDiseaseDraftProvider


POLICY = {
    "policy_version": "test-v1",
    "policy_source": "Synthetic fastVEP-style triage policy",
    "pp2_min_mis_z": 3.09,
    "bp1_min_p_li": 0.9,
    "bp1_max_mis_z": 1.0,
    "bp1_max_pathogenic_missense": 3,
    "pvs1_min_p_li": 0.9,
    "pvs1_max_loeuf": 0.35,
}


def source(**values):
    return {
        "source": "test source",
        "source_version": "1",
        "retrieved_at": "2026-09-16T00:00:00Z",
        **values,
    }


class GeneDiseaseDraftProviderTests(unittest.TestCase):
    def setUp(self):
        self.provider = GeneDiseaseDraftProvider(POLICY)
        self.variant = Variant("GRCh38", "1", 100, "A", "G")
        self.validity = [source(
            gene="TEST", condition="MONDO:0000001", classification="Definitive"
        )]

    def build(self, **updates):
        values = {
            "variant_keys": [self.variant.key],
            "gene": "TEST",
            "condition": "MONDO:0000001",
            "transcript": "NM_TEST.1",
            "generated_at": "2026-09-16T01:00:00Z",
            "validity": self.validity,
            "constraint": source(gene="TEST", mis_z=3.5, p_li=0.98, loeuf=0.2),
            "clinvar_spectrum": source(
                gene="TEST",
                complete=True,
                pathogenic_missense_count=10,
                pathogenic_truncating_count=20,
                benign_missense_count=1,
            ),
        }
        values.update(updates)
        return self.provider.build(**values)[0]

    def test_generates_review_only_candidates(self):
        draft = self.build()
        self.assertEqual(draft["category"], "gene_disease_draft")
        self.assertEqual(draft["quality_status"], "DRAFT")
        self.assertFalse(draft["criterion_eligible"])
        self.assertEqual(draft["validity_state"], "SUPPORTED")
        self.assertEqual(draft["suggestions"]["PP2"]["status"], "CANDIDATE")
        self.assertEqual(draft["suggestions"]["BP1"]["status"], "NOT_SUGGESTED")
        self.assertEqual(draft["suggestions"]["PVS1"]["status"], "CANDIDATE")
        self.assertIsNone(draft["suggestions"]["PP2"]["criterion_met"])

    def test_draft_cannot_drive_a_criterion(self):
        annotation = source(
            category="annotation",
            quality_status="PASS",
            evidence_id="test:annotation",
            variant_key=self.variant.key,
            transcript="NM_TEST.1",
            condition="MONDO:0000001",
            gene="TEST",
            consequences=["missense_variant"],
        )
        prepared = {
            "record_id": "draft-safety",
            "variant": self.variant.to_dict(),
            "transcript": "NM_TEST.1",
            "condition": "MONDO:0000001",
            "identity_provenance": [{"source": "test"}],
        }
        value = evaluate_prepared_record(
            prepared,
            make_services([annotation, self.build()]),
            {},
            ["PP2"],
        )[0]
        self.assertEqual(value.status, CriterionStatus.UNKNOWN)
        self.assertIn("gene_disease", value.missing_inputs)

    def test_incomplete_spectrum_is_not_read_as_absence(self):
        draft = self.build(clinvar_spectrum=source(
            gene="TEST",
            complete=False,
            pathogenic_missense_count=0,
            pathogenic_truncating_count=20,
            benign_missense_count=0,
        ))
        self.assertEqual(draft["suggestions"]["BP1"]["status"], "INSUFFICIENT")

    def test_curated_negative_validity_blocks_all_suggestions(self):
        draft = self.build(validity=[source(
            gene="TEST", condition="MONDO:0000001", classification="Refuted"
        )])
        self.assertEqual(draft["validity_state"], "CURATED_NEGATIVE")
        self.assertEqual(
            {item["status"] for item in draft["suggestions"].values()},
            {"CURATED_NEGATIVE"},
        )

    def test_limited_validity_remains_insufficient_not_negative(self):
        draft = self.build(validity=[source(
            gene="TEST", condition="MONDO:0000001", classification="Limited"
        )])
        self.assertEqual(draft["validity_state"], "LIMITED")
        self.assertEqual(
            {item["status"] for item in draft["suggestions"].values()},
            {"INSUFFICIENT"},
        )

    def test_gene_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Constraint gene"):
            self.build(constraint=source(gene="OTHER", mis_z=3.5, p_li=1.0, loeuf=0.1))

    def test_each_variant_gets_a_unique_evidence_id(self):
        records = self.provider.build(
            variant_keys=[self.variant.key, "GRCh38:1:101:A:T"],
            gene="TEST",
            condition="MONDO:0000001",
            transcript="NM_TEST.1",
            generated_at="2026-09-16T01:00:00Z",
            validity=self.validity,
        )
        self.assertEqual(len({item["evidence_id"] for item in records}), 2)

    def test_builds_a_document_from_versioned_groups(self):
        document = build_draft_document({
            "source_version": "test-sources-v1",
            "retrieved_at": "2026-09-16T01:00:00Z",
            "policy": POLICY,
            "groups": [{
                "variant_keys": [self.variant.key],
                "gene": "TEST",
                "condition": "MONDO:0000001",
                "transcript": "NM_TEST.1",
                "validity": self.validity,
            }],
        })
        self.assertEqual(document["source_version"], "test-sources-v1")
        self.assertEqual(len(document["evidence"]), 1)

    def test_expands_machine_authored_assessments_without_claiming_human_signoff(self):
        document = build_assessment_document({
            "source_version": "test-review-v1",
            "assessed_at": "2026-09-16T01:00:00Z",
            "method": "test source synthesis",
            "policy_version": "test-policy-v1",
            "groups": [{
                "variant_keys": [self.variant.key],
                "gene": "TEST",
                "condition": "MONDO:0000001",
                "transcript": "NM_TEST.1",
                "lof_mechanism_established": False,
                "pp2_applicable": False,
                "bp1_applicable": False,
                "sources": [{"url": "https://example.test/source"}],
            }],
        }, {"evidence": [source(evidence_id="test:base")]})
        self.assertEqual(len(document["evidence"]), 2)
        assessment = document["evidence"][1]
        self.assertEqual(assessment["assessment_method"], "automated")
        self.assertFalse(assessment["human_signoff"])

    def test_committed_demo_review_honors_vcep_and_disease_context(self):
        path = Path(__file__).parents[1] / "config" / "gene-disease-review-decisions.json"
        assessments = build_assessment_document(json.loads(path.read_text(encoding="utf-8")))
        self.assertEqual(len(assessments["evidence"]), 15)

        def evaluate(gene, transcript, variant, condition, code):
            annotation = source(
                category="annotation",
                quality_status="PASS",
                evidence_id=f"test:annotation:{gene}",
                variant_key=variant.key,
                transcript=transcript,
                gene=gene,
                consequences=["missense_variant"],
            )
            prepared = {
                "record_id": f"demo-{gene}",
                "variant": variant.to_dict(),
                "transcript": transcript,
                "condition": condition,
                "identity_provenance": [{"source": "test"}],
            }
            return evaluate_prepared_record(
                prepared,
                make_services([annotation, *assessments["evidence"]]),
                {},
                [code],
            )[0]

        myh7 = evaluate(
            "MYH7", "NM_000257.4", Variant("GRCh38", "14", 23425971, "G", "A"),
            "MONDO:0005045", "PP2",
        )
        tnni3_hcm = evaluate(
            "TNNI3", "NM_000363.5", Variant("GRCh38", "19", 55156239, "G", "A"),
            "MONDO:0005045", "PP2",
        )
        tnni3_arvc = evaluate(
            "TNNI3", "NM_000363.5", Variant("GRCh38", "19", 55156239, "G", "A"),
            "MONDO:0016587", "PP2",
        )
        self.assertEqual(myh7.status, CriterionStatus.UNKNOWN)
        self.assertEqual(tnni3_hcm.status, CriterionStatus.MET)
        self.assertEqual(tnni3_arvc.status, CriterionStatus.UNKNOWN)


if __name__ == "__main__":
    unittest.main()
