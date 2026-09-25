"""G2P curation read for a disease-scoped mechanism, and refused where it does not settle one."""

import unittest

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.providers.gene2phenotype import (
    METHOD, Gene2PhenotypeProvider, mondo_accession,
)


def summary(records):
    return {"gene_symbol": "TEST", "records_summary": records}


def entry(stable_id="G2P0001", mechanism="loss of function", confidence="definitive",
          genotype="monoallelic_autosomal"):
    return {"stable_id": stable_id, "molecular_mechanism": mechanism,
            "confidence": confidence, "genotype": genotype, "panels": ["Cardiac"],
            "disease": "TEST-related disease"}


def detail(stable_id="G2P0001", accessions=("MONDO:0007268",), last_updated="2025-09-05",
           mechanism="loss of function", mechanism_support="evidence"):
    return {
        "stable_id": stable_id,
        "last_updated": last_updated,
        "molecular_mechanism": {"mechanism": mechanism, "mechanism_support": mechanism_support},
        "disease": {"name": "TEST-related disease", "ontology_terms": [
            {"accession": item, "source": "Mondo" if item.startswith("MONDO") else "OMIM"}
            for item in accessions]},
    }


class FakeClient:
    """Answers the summary URL and each record URL, and counts what was asked for."""

    def __init__(self, records, details=None):
        self.records = records
        self.details = details or {}
        self.urls = []

    def fetch(self, url, *, response_format="json", **kwargs):
        self.urls.append(url)
        if url.endswith("/summary/"):
            body = summary(self.records)
        else:
            body = self.details[url.rstrip("/").rsplit("/", 1)[-1]]
        return {"body": body, "retrieved_at": "2026-09-17T00:00:00Z", "body_sha256": "abc"}


class Gene2PhenotypeTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant("GRCh38", "1", 2, "C", "T")

    def mechanism(self, records, details=None):
        client = FakeClient(records, details)
        provider = Gene2PhenotypeProvider(client)
        return provider.get_mechanism(self.variant, "TEST"), client

    def test_a_definitive_loss_of_function_record_is_scoped_to_its_mondo_disease(self):
        records, _ = self.mechanism([entry()], {"G2P0001": detail()})
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertIs(record["lof_mechanism_established"], True)
        self.assertEqual(record["condition"], "MONDO:0007268")
        self.assertEqual(record["inheritance"], "autosomal_dominant")
        self.assertEqual(record["method"], METHOD)
        self.assertEqual(record["assessment_method"], "automated")
        self.assertEqual(record["source_version"], "2025-09-05")
        self.assertEqual(record["g2p_mechanism_support"], "evidence")
        self.assertNotIn("curator", record)

    def test_a_non_loss_of_function_mechanism_is_read_as_a_denial(self):
        for mechanism in ("gain of function", "dominant negative",
                          "undetermined non-loss-of-function"):
            with self.subTest(mechanism=mechanism):
                records, _ = self.mechanism(
                    [entry(mechanism=mechanism)], {"G2P0001": detail(mechanism=mechanism)})
                self.assertIs(records[0]["lof_mechanism_established"], False)

    def test_an_undetermined_mechanism_is_still_emitted_as_unresolved(self):
        """G2P has not settled it, which is not a statement that loss of function is ruled
        out - the same reading the dosage provider gives a score of 30. Unlike before,
        this is no longer dropped: it is emitted with lof_mechanism_established=None so
        pvs1.py's own mechanism-conflict detection can still see it (see the module
        docstring for the real MYOC case this was found from)."""
        records, client = self.mechanism(
            [entry(mechanism="undetermined")], {"G2P0001": detail(mechanism="undetermined")})
        self.assertEqual(len(records), 1)
        self.assertIsNone(records[0]["lof_mechanism_established"])
        self.assertEqual(len(client.urls), 2)  # the record detail IS now requested

    def test_a_loss_of_function_call_without_evidence_support_is_unresolved(self):
        """Real MYOC false positive (2026-09-24): an older G2P record claimed "loss of
        function" for one glaucoma subtype with mechanism_support="inferred" (inferred
        from the variant types curated, not a functional study) - not what ACMG PVS1
        means by an ESTABLISHED mechanism, so it is downgraded to unresolved rather than
        trusted at face value."""
        records, _ = self.mechanism(
            [entry()], {"G2P0001": detail(mechanism_support="inferred")})
        self.assertIsNone(records[0]["lof_mechanism_established"])
        self.assertEqual(records[0]["g2p_mechanism_support"], "inferred")

    def test_a_denial_without_evidence_support_is_also_unresolved(self):
        """The same downgrade applies symmetrically to a "gain of function"/"dominant
        negative" denial - an inferred denial should not out-rank a differently-resolved
        record for a related condition any more than an inferred affirmation should."""
        records, _ = self.mechanism(
            [entry(mechanism="gain of function")],
            {"G2P0001": detail(mechanism="gain of function", mechanism_support="inferred")})
        self.assertIsNone(records[0]["lof_mechanism_established"])

    def test_validity_below_moderate_emits_nothing_and_never_a_denial(self):
        """Confidence is a gene-disease validity statement; reading it as a mechanism would
        be the conflation gene_disease_draft.py exists to prevent."""
        for confidence in ("limited", "disputed", "refuted"):
            with self.subTest(confidence=confidence):
                records, _ = self.mechanism([entry(confidence=confidence)])
                self.assertEqual(records, [])

    def test_a_record_naming_no_mondo_term_is_not_emitted(self):
        """A disease name is not an identifier, and matching on names is the route the
        disease-match levels rule out."""
        records, _ = self.mechanism(
            [entry()], {"G2P0001": detail(accessions=("115197",))})
        self.assertEqual(records, [])

    def test_two_mondo_terms_are_not_a_tie_to_break_here(self):
        self.assertIsNone(mondo_accession(
            {"ontology_terms": [{"accession": "MONDO:1"}, {"accession": "MONDO:2"}]}))

    def test_each_gene_disease_pair_keeps_its_own_record(self):
        """One gene losing function in one disease and not in another is the case a
        gene-level source cannot separate, and the reason this provider exists."""
        records, _ = self.mechanism(
            [entry("G2P0001"), entry("G2P0002", mechanism="gain of function")],
            {"G2P0001": detail("G2P0001", ("MONDO:0007268",)),
             "G2P0002": detail("G2P0002", ("MONDO:0008647",), mechanism="gain of function")})
        self.assertEqual(
            {item["condition"]: item["lof_mechanism_established"] for item in records},
            {"MONDO:0007268": True, "MONDO:0008647": False})
        self.assertEqual(len({item["evidence_id"] for item in records}), 2)

    def test_the_allelic_requirement_becomes_the_inheritance_mode(self):
        for genotype, mode in (("biallelic_autosomal", "autosomal_recessive"),
                               ("monoallelic_X_hemizygous", "x_linked"),
                               ("monoallelic_X_heterozygous", "x_linked"),
                               ("mitochondrial", "mitochondrial")):
            with self.subTest(genotype=genotype):
                records, _ = self.mechanism(
                    [entry(genotype=genotype)], {"G2P0001": detail()})
                self.assertEqual(records[0]["inheritance"], mode)

    def test_an_unrecognized_allelic_requirement_leaves_the_mode_unset(self):
        records, _ = self.mechanism(
            [entry(genotype="something_new")], {"G2P0001": detail()})
        self.assertIsNone(records[0]["inheritance"])


if __name__ == "__main__":
    unittest.main()
