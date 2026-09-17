import copy
import unittest

from acmg_pipeline.automated_core.identity import evaluation_inputs, reconcile
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.automated_core.reference import normalize


class MemoryReference:
    def __init__(self, sequence):
        self.bases = sequence

    def sequence(self, chrom, start, end):
        if not 1 <= start <= end <= len(self.bases):
            raise ValueError("out of bounds")
        return self.bases[start - 1:end]


class NormalizationTests(unittest.TestCase):
    def test_deletion_left_aligns_in_repeat(self):
        ref = MemoryReference("CAAAAAG")
        variant = Variant("GRCh38", "1", 4, "AA", "A")
        self.assertEqual(normalize(variant, ref), Variant("GRCh38", "1", 1, "CA", "C"))

    def test_insertion_left_aligns(self):
        ref = MemoryReference("CAAAAAG")
        variant = Variant("GRCh38", "1", 5, "A", "AA")
        self.assertEqual(normalize(variant, ref), Variant("GRCh38", "1", 1, "C", "CA"))

    def test_mnv_trims_without_changing_effect(self):
        ref = MemoryReference("ACGT")
        self.assertEqual(normalize(Variant("GRCh38", "1", 1, "ACG", "ATG"), ref),
                         Variant("GRCh38", "1", 2, "C", "T"))

    def test_reference_mismatch(self):
        with self.assertRaisesRegex(ValueError, "REF_MISMATCH"):
            normalize(Variant("GRCh38", "1", 1, "T", "G"), MemoryReference("A"))

    def test_first_base_anchor(self):
        ref = MemoryReference("AAAAC")
        variant = normalize(Variant("GRCh38", "1", 2, "AA", "A"), ref)
        self.assertEqual(variant.pos, 1)
        self.assertEqual(variant.ref, "AA")


class IdentityTests(unittest.TestCase):
    def setUp(self):
        self.ref = MemoryReference("ACGT")
        self.variant = Variant("GRCh38", "1", 2, "C", "T").to_dict()
        self.record = {
            "record_id": "case1:1:1", "source": {"variant_id": "example"},
            "raw_variant": self.variant, "parsed_variant": self.variant,
            "identity": {"TRANSCRIPT": "NM_000001.1", "HGVSC": "c.2C>T"},
            "issues": [], "identity_status": "PENDING", "resolution": None,
            "source_record": "ACMG_CODES=PS1;CLNSIG=Pathogenic",
        }
        self.candidate = {
            "variant": self.variant, "source": "synthetic-test", "source_version": "1",
            "retrieved_at": "2026-09-14T00:00:00Z",
            "matched_identifiers": dict(self.record["identity"]),
        }

    def test_verified_mapping(self):
        result = reconcile(self.record, [self.candidate], self.ref)
        self.assertEqual(result["identity_status"], "VERIFIED")
        payload = evaluation_inputs([result])[0]
        self.assertNotIn("source_record", payload)

    def test_missing_alt_can_be_corrected_only_with_mapping(self):
        self.record["parsed_variant"] = None
        self.record["issues"] = ["INVALID_VARIANT: missing ALT"]
        self.assertEqual(reconcile(self.record, [], self.ref)["identity_status"], "PENDING")
        self.assertEqual(reconcile(self.record, [self.candidate], self.ref)["identity_status"],
                         "CORRECTED")

    def test_disagreement_cannot_be_overridden_by_order(self):
        other = {**self.candidate, "variant": {**self.variant, "alt": "G"}}
        for candidates in ([other, self.candidate], [self.candidate, other]):
            result = reconcile(self.record, candidates, self.ref)
            self.assertEqual(result["identity_status"], "PENDING")
            self.assertIn("CONFLICTING_IDENTITY_CANDIDATES", result["issues"])

    def test_transcript_version_mismatch(self):
        bad = copy.deepcopy(self.candidate)
        bad["matched_identifiers"]["TRANSCRIPT"] = "NM_000001.2"
        self.assertEqual(reconcile(self.record, [bad], self.ref)["identity_status"], "PENDING")

    def test_labels_do_not_change_evaluation_input(self):
        first = evaluation_inputs([reconcile(self.record, [self.candidate], self.ref)])
        self.record["source_record"] = "ACMG_CODES=BA1;CLNSIG=Benign"
        second = evaluation_inputs([reconcile(self.record, [self.candidate], self.ref)])
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
