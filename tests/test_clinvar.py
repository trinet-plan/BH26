import unittest

from acmg.core.models import Variant
from acmg.providers.clinvar import (
    ClinVarComparatorProvider, ClinVarHotspotProvider, ClinVarProvider,
)


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


HOTSPOT_POLICY = {"window_aa": 5, "min_pathogenic": 3, "max_benign": 0,
                  "method": "clinvar_local_density", "policy_source": "test policy",
                  "policy_version": "PM1-hotspot-test"}

SUMMARIES = {
    "1": {"accession": "VCV000000001", "gene_sort": "TEST", "protein_change": "R248W",
          "germline_classification": {"description": "Pathogenic",
                                      "review_status": "criteria provided, multiple submitters",
                                      "trait_set": [{"trait_name": "Test disease"}]}},
    "2": {"accession": "VCV000000002", "gene_sort": "TEST", "protein_change": "R248Q",
          "germline_classification": {"description": "Likely pathogenic", "trait_set": []}},
    "3": {"accession": "VCV000000003", "gene_sort": "TEST", "protein_change": "G245S",
          "germline_classification": {"description": "Pathogenic", "trait_set": []}},
    # Outside the +/-5 window.
    "4": {"accession": "VCV000000004", "gene_sort": "TEST", "protein_change": "P300L",
          "germline_classification": {"description": "Pathogenic", "trait_set": []}},
    # Counted neither way.
    "5": {"accession": "VCV000000005", "gene_sort": "TEST", "protein_change": "R249T",
          "germline_classification": {"description": "Uncertain significance", "trait_set": []}},
    "6": {"accession": "VCV000000006", "gene_sort": "TEST", "protein_change": "R249K",
          "germline_classification": {"description": "Conflicting classifications of pathogenicity",
                                      "trait_set": []}},
    # Not positionable.
    "7": {"accession": "VCV000000007", "gene_sort": "TEST", "protein_change": "R248fs",
          "germline_classification": {"description": "Pathogenic", "trait_set": []}},
    # Another gene sharing the locus.
    "8": {"accession": "VCV000000008", "gene_sort": "OTHER", "protein_change": "R248W",
          "germline_classification": {"description": "Pathogenic", "trait_set": []}},
}

BENIGN_SUMMARY = {"accession": "VCV000000009", "gene_sort": "TEST", "protein_change": "S246N",
                  "germline_classification": {"description": "Likely benign", "trait_set": []}}


# What ClinVar reports on the target protein NP_1.1, which is what the counts must use.
PROTEIN_EXPRESSIONS = {
    "1": "NP_1.1:p.Arg248Trp", "2": "NP_1.1:p.Arg248Gln", "3": "NP_1.1:p.Gly245Ser",
    "6": "NP_1.1:p.Arg249Lys", "9": "NP_1.1:p.Ser246Asn", "10": "NP_1.1:p.Arg248His",
    # An alternate isoform puts this one in the window; on NP_1.1 it is outside.
    "11": "NP_1.1:p.Ala275Thr",
}

ALTERNATE_ISOFORM_SUMMARY = {
    "accession": "VCV000000011", "gene_sort": "TEST", "protein_change": "A248T, A275T",
    "germline_classification": {"description": "Pathogenic", "trait_set": []}}


def candidate_xml(uid):
    expression = PROTEIN_EXPRESSIONS.get(uid)
    protein = (f"<HGVS><ProteinExpression><Expression>{expression}</Expression>"
               "</ProteinExpression></HGVS>") if expression else ""
    return ("<ClinVarResult-Set><VariationArchive "
            f'Accession="VCV{int(uid):09d}" Version="1" VariationID="{uid}">'
            f"<ClassifiedRecord><SimpleAllele><HGVSlist>{protein}</HGVSlist>"
            "</SimpleAllele></ClassifiedRecord></VariationArchive></ClinVarResult-Set>")


class HotspotClient:
    def __init__(self, uids=None, summaries=None, count=None):
        self.uids = uids or sorted(SUMMARIES)
        self.summaries = summaries or SUMMARIES
        self.count = len(self.uids) if count is None else count
        self.urls = []

    def fetch(self, url, **kwargs):
        self.urls.append(url)
        if "esearch.fcgi" in url:
            return {"body": {"esearchresult": {"count": str(self.count), "retmax": "500",
                                               "idlist": list(self.uids)}},
                    "retrieved_at": "2026-09-15T00:00:00Z"}
        if "efetch.fcgi" in url:
            uid = url.split("&id=", 1)[1].split("&", 1)[0]
            return {"body": candidate_xml(uid), "retrieved_at": "2026-09-15T00:00:02Z"}
        return {"body": {"result": {"uids": list(self.uids), **self.summaries}},
                "retrieved_at": "2026-09-15T00:00:01Z"}


class ClinVarHotspotTests(unittest.TestCase):
    annotation = {"transcript": "NM_1.2", "gene": "TEST", "protein_id": "NP_1.1",
                  "protein_start": 248, "protein_end": 248, "ref_aa": "R", "alt_aa": "H"}

    def search(self, client=None, policy=None):
        variant = Variant("GRCh38", "1", 100, "G", "A")
        provider = ClinVarHotspotProvider(client or HotspotClient(), "2026-09-15",
                                          policy or HOTSPOT_POLICY)
        return provider.search_hotspot(self.annotation, variant)

    def test_counts_only_positioned_in_window_missense_of_the_same_gene(self):
        search, region = self.search()
        self.assertTrue(search["complete"])
        self.assertEqual(region["region_type"], "mutational_hotspot")
        self.assertEqual(region["start"], 243)
        self.assertEqual(region["end"], 253)
        self.assertEqual(region["pathogenic_count"], 3)
        self.assertEqual(region["benign_count"], 0)
        self.assertEqual([item["variation_id"] for item in region["counted_variants"]
                          if item["bucket"] == "pathogenic"], ["1", "2", "3"])
        # Conflicting records are neither pathogenic nor benign, but stay visible for audit.
        self.assertEqual([item["variation_id"] for item in region["counted_variants"]
                          if item["bucket"] == "conflicting"], ["6"])
        self.assertIn("Test disease", region["counted_conditions"])

    def test_query_variant_does_not_support_its_own_density(self):
        summaries = {**SUMMARIES, "10": {
            "accession": "VCV000000010", "gene_sort": "TEST", "protein_change": "R248H",
            "variation_set": [{"canonical_spdi": "NC_000001.11:99:G:A"}],
            "germline_classification": {"description": "Pathogenic", "trait_set": []}}}
        client = HotspotClient(uids=sorted(summaries, key=int), summaries=summaries)
        _, region = self.search(client)
        self.assertEqual(region["self_excluded"], 1)
        self.assertEqual(region["pathogenic_count"], 3)
        self.assertEqual([item["bucket"] for item in region["counted_variants"]
                          if item["variation_id"] == "10"], ["query_variant"])

    def test_other_alleles_at_the_same_position_still_count(self):
        summaries = {**SUMMARIES, "10": {
            "accession": "VCV000000010", "gene_sort": "TEST", "protein_change": "R248P",
            "variation_set": [{"canonical_spdi": "NC_000001.11:99:G:C"}],
            "germline_classification": {"description": "Pathogenic", "trait_set": []}}}
        client = HotspotClient(uids=sorted(summaries, key=int), summaries=summaries)
        _, region = self.search(client)
        self.assertEqual(region["self_excluded"], 0)
        self.assertEqual(region["pathogenic_count"], 4)

    def test_benign_variation_in_window_is_counted(self):
        summaries = {**SUMMARIES, "9": BENIGN_SUMMARY}
        client = HotspotClient(uids=sorted(summaries), summaries=summaries)
        _, region = self.search(client)
        self.assertEqual(region["benign_count"], 1)
        self.assertEqual(region["pathogenic_count"], 3)

    def test_alternate_isoform_numbering_is_not_counted(self):
        """A window hit through another isoform's numbering must not become a count."""
        summaries = {**SUMMARIES, "11": ALTERNATE_ISOFORM_SUMMARY}
        client = HotspotClient(uids=sorted(summaries, key=int), summaries=summaries)
        _, region = self.search(client)
        self.assertEqual(region["pathogenic_count"], 3)
        self.assertEqual(region["outside_window_excluded"], 1)
        excluded = next(item for item in region["counted_variants"]
                        if item["variation_id"] == "11")
        self.assertEqual(excluded["bucket"], "outside_window_on_protein")
        self.assertEqual(excluded["reported_positions"], [248, 275])
        self.assertEqual(excluded["protein_positions"], [275])

    def test_candidate_without_this_protein_is_not_counted(self):
        summaries = {**SUMMARIES, "12": {
            "accession": "VCV000000012", "gene_sort": "TEST", "protein_change": "R248L",
            "germline_classification": {"description": "Pathogenic", "trait_set": []}}}
        client = HotspotClient(uids=sorted(summaries, key=int), summaries=summaries)
        _, region = self.search(client)
        self.assertEqual(region["pathogenic_count"], 3)
        self.assertEqual(region["unplaced_excluded"], 1)

    def test_region_evidence_is_automated_and_policy_versioned(self):
        _, region = self.search()
        self.assertEqual(region["assessment_method"], "automated")
        self.assertNotIn("curator", region)
        self.assertEqual(region["policy_version"], "PM1-hotspot-test")
        self.assertEqual(region["method"], "clinvar_local_density")
        self.assertIn("NOT_CRITICAL_DOMAIN", region["use_restriction"])

    def test_truncated_search_emits_no_density_evidence(self):
        client = HotspotClient(count=900)
        search, region = self.search(client)
        self.assertFalse(search["complete"])
        self.assertIsNone(region)
        self.assertEqual(search["quality_status"], "INCOMPLETE_SEARCH")

    def test_incomplete_policy_is_rejected(self):
        policy = {key: value for key, value in HOTSPOT_POLICY.items() if key != "policy_version"}
        with self.assertRaisesRegex(ValueError, "policy_version"):
            self.search(policy=policy)

    def test_missing_protein_annotation_is_rejected(self):
        variant = Variant("GRCh38", "1", 100, "G", "A")
        provider = ClinVarHotspotProvider(HotspotClient(), "2026-09-15", HOTSPOT_POLICY)
        annotation = {key: value for key, value in self.annotation.items() if key != "protein_id"}
        with self.assertRaises(ValueError):
            provider.search_hotspot(annotation, variant)


if __name__ == "__main__":
    unittest.main()
