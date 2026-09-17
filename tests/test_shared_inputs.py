import inspect
import unittest
from importlib import import_module

from acmg_pipeline.automated_core.models import CRITERIA
from acmg_pipeline.automated_engine import evaluate_record, make_services
from acmg_pipeline.automated_va_spec import export_document
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


class SharedInputInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.variant = VariantRecord(
            chrom="1",
            pos=2,
            id="shared-input-test",
            ref="C",
            alt="T",
            qual="99",
            filter="PASS",
            info={"GENE": "TEST", "HGVSC": "c.2C>T"},
        )
        self.clinical_note = ClinicalNoteExtraction.from_json({
            "proband": {
                "phenotype": {"affected_status": True},
                "genotype": {"variant_status": True, "zygosity": "heterozygous"},
            },
            "family": {"inheritance_pattern": "autosomal_dominant"},
        })

    def test_every_criterion_has_the_shared_public_signature(self):
        for code in CRITERIA:
            with self.subTest(code=code):
                evaluate = import_module(f"acmg_pipeline.criteria.{code.lower()}").evaluate
                parameters = list(inspect.signature(evaluate).parameters.values())
                self.assertEqual([value.name for value in parameters[:2]],
                                 ["variant", "clinical_note"])
                self.assertIs(parameters[0].annotation, VariantRecord)
                self.assertIs(parameters[1].annotation, ClinicalNoteExtraction)

    def test_engine_accepts_shared_objects_and_va_spec_remains_valid(self):
        results = evaluate_record(
            self.variant,
            self.clinical_note,
            make_services([]),
            {},
        )
        self.assertEqual({result.criterion for result in results}, set(CRITERIA))
        self.assertTrue(all(result.variant["assembly"] == "GRCh38" for result in results))
        document = export_document([{
            "record_id": self.variant.id,
            "variant": results[0].variant,
            "results": [result.to_dict() for result in results],
        }])
        self.assertEqual(document["validated_by"]["schema_version"], "1.0.1")
        self.assertEqual(len(document["records"][0]["criterion_assessments"]), len(CRITERIA))


if __name__ == "__main__":
    unittest.main()
