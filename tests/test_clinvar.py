import unittest

from acmg.core.models import Variant
from acmg.providers.clinvar import ClinVarComparatorProvider, ClinVarProvider


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

COMPARATOR_XML = """<ClinVarResult-Set>
<VariationArchive Accession="VCV000000123" Version="4" VariationID="123">
  <ClassifiedRecord><SimpleAllele><Location><SequenceLocation Assembly="GRCh38" Chr="1"
    positionVCF="200" referenceAlleleVCF="C" alternateAlleleVCF="T" /></Location>
    <HGVSlist><HGVS><NucleotideExpression><Expression>NM_1.2:c.2C&gt;T</Expression>
    </NucleotideExpression></HGVS><HGVS><ProteinExpression>
    <Expression>NP_1.1:p.Arg7His</Expression></ProteinExpression></HGVS></HGVSlist></SimpleAllele>
    <Classifications><GermlineClassification><ReviewStatus>reviewed by expert panel</ReviewStatus>
    <Description>Pathogenic</Description><ConditionList><Trait><Name>
    <ElementValue Type="Preferred">Test disease</ElementValue></Name>
    <XRef DB="MONDO" ID="0000001" /></Trait></ConditionList>
    </GermlineClassification></Classifications></ClassifiedRecord>
</VariationArchive></ClinVarResult-Set>"""


class Client:
    def fetch(self, url, **kwargs):
        self.url = url
        self.kwargs = kwargs
        return {"body": XML, "retrieved_at": "2026-09-15T00:00:00Z"}


class ComparatorClient:
    def fetch(self, url, **kwargs):
        if "esearch.fcgi" in url:
            return {"body": {"esearchresult": {"count": "1", "retmax": "1",
                                                 "idlist": ["123"]}},
                    "retrieved_at": "2026-09-15T00:00:00Z"}
        return {"body": COMPARATOR_XML, "retrieved_at": "2026-09-15T00:00:00Z"}


class ComparatorEnsembl:
    def map_record_with_evidence(self, record):
        variant = {"assembly": "GRCh38", "chrom": "1", "pos": 200,
                   "ref": "C", "alt": "T"}
        annotation = {"protein_id": "NP_1.1", "protein_start": 7,
                      "ref_aa": "R", "alt_aa": "H", "hgvsp": "NP_1.1:p.Arg7His"}
        predictions = [{"predictor": "SpliceAI", "score": 0.02}]
        return {"variant": variant}, annotation, predictions


class IntervalComparatorEnsembl:
    def map_record_with_evidence(self, record):
        variant = {"assembly": "GRCh38", "chrom": "1", "pos": 199,
                   "ref": "AC", "alt": "GT"}
        annotation = {"protein_id": "NP_1.1", "protein_start": 6, "protein_end": 7,
                      "ref_aa": "GR", "alt_aa": "GH", "hgvsp": "NP_1.1:p.Arg7His"}
        return {"variant": variant}, annotation, []


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

    def test_ps1_comparator_search_builds_protein_level_evidence(self):
        query = Variant("GRCh38", "1", 100, "G", "A")
        annotation = {"transcript": "NM_1.2", "gene": "TEST", "protein_id": "NP_1.1",
                      "protein_start": 7, "ref_aa": "R", "alt_aa": "H",
                      "hgvsp": "NP_1.1:p.Arg7His"}
        search, matches = ClinVarComparatorProvider(
            ComparatorClient(), "2026-09-15", ComparatorEnsembl()
        ).search_ps1(annotation, query, 0.01)
        self.assertTrue(search["complete"])
        self.assertEqual(len(matches), 1)
        match = matches[0]
        self.assertEqual(match["classification"], "Pathogenic")
        self.assertTrue(match["review_status_eligible"])
        self.assertFalse(match["splice_conflict"])
        self.assertIn("MONDO:0000001", match["conditions"])

    def test_ps1_accepts_mnv_interval_with_one_matching_residue_change(self):
        query = Variant("GRCh38", "1", 100, "G", "A")
        annotation = {"transcript": "NM_1.2", "gene": "TEST", "protein_id": "NP_1.1",
                      "protein_start": 7, "ref_aa": "R", "alt_aa": "H",
                      "hgvsp": "NP_1.1:p.Arg7His"}
        _, matches = ClinVarComparatorProvider(
            ComparatorClient(), "2026-09-15", IntervalComparatorEnsembl()
        ).search_ps1(annotation, query, 0.01)
        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["comparator_protein_interval"]["ref_aa"], "GR")
        self.assertFalse(matches[0]["splice_effect_checked"])


if __name__ == "__main__":
    unittest.main()
