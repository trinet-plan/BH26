import unittest

from acmg.core.models import Variant
from acmg.providers.clinvar import ClinVarProvider


XML = """<ClinVarResult-Set>
<VariationArchive Accession="VCV000042506" Version="7" VariationID="42506">
  <ClassifiedRecord>
    <SimpleAllele>
      <Location Assembly="GRCh38" Chr="11" start="47343571"
                referenceAlleleVCF="G" alternateAlleleVCF="A" />
    </SimpleAllele>
    <GermlineClassification ReviewStatus="criteria provided, multiple submitters">
      <Description>Benign</Description>
    </GermlineClassification>
  </ClassifiedRecord>
</VariationArchive>
</ClinVarResult-Set>"""


class Client:
    def fetch(self, url, **kwargs):
        self.url = url
        self.kwargs = kwargs
        return {"body": XML, "retrieved_at": "2026-09-15T00:00:00Z"}


class ClinVarProviderTests(unittest.TestCase):
    def test_vcv_record_corroborates_exact_variant(self):
        variant = Variant("GRCh38", "11", 47343571, "G", "A")
        evidence, identity = ClinVarProvider(Client(), "2026-09-15").get_record(
            "VCV000042506", variant)
        self.assertEqual(evidence["accession"], "VCV000042506.7")
        self.assertEqual(evidence["classifications"], ["Benign"])
        self.assertIn("NOT_PP5_BP6", evidence["use_restriction"])
        self.assertEqual(identity["matched_identifiers"],
                         {"CLNVARIATIONID": "VCV000042506"})

    def test_variant_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "expected GRCh38"):
            ClinVarProvider(Client(), "2026-09-15").get_record(
                "VCV000042506", Variant("GRCh38", "11", 1, "G", "A"))

    def test_versioned_request_records_current_version_advance(self):
        evidence, identity = ClinVarProvider(Client(), "2026-09-15").get_record(
            "VCV000042506.6", Variant("GRCh38", "11", 47343571, "G", "A"))
        self.assertTrue(evidence["version_advanced"])
        self.assertEqual(evidence["requested_accession"], "VCV000042506.6")
        self.assertEqual(identity["matched_identifiers"],
                         {"CLNVARIATIONID": "VCV000042506.6"})


if __name__ == "__main__":
    unittest.main()
