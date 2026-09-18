"""curatorInfo folded into hasEvidenceItems, not its own extension.

Added 2026-09-18 per the user's direction: PM1/PM5's UniProt domain facts
(and PM3's ClinVar mentions) are the structured result of a real lookup this
pipeline performed, the same kind of fact the automated engine already
reports via hasEvidenceItems for population/annotation evidence - unlike
referenceLink, a bare pointer to a page never actually read (see
export.build_reference_extensions()'s own docstring for why that one stays
an extension).

GA4GH's hasEvidenceItems has no generic/open "StudyResult" member (only
CohortAlleleFrequencyStudyResult/ExperimentalVariantFunctionalImpactStudy
Result, Statement, EvidenceLine, or a bare iriReference are valid - checked
directly against the installed ga4gh.va_spec models), so the facts are
nested as their own EvidenceLine, the same shape build_paper_evidence_line()
already uses for per-paper literature items.
"""

import unittest
from unittest.mock import patch

from acmg_pipeline.criteria.curator_info import DomainFeature, UniprotDomainContext
from acmg_pipeline.export import build_reference_extensions
from acmg_pipeline.vcf_record import VariantRecord


def variant(**overrides):
    fields = dict(chrom="14", pos=23884280, id="v1", ref="C", alt="T",
                  qual="", filter="", info={"GENE": "MYH7", "HGVSC": "c.2155C>T"})
    fields.update(overrides)
    return VariantRecord(**fields)


class CuratorInfoEvidenceItemTests(unittest.TestCase):
    def test_pm1_curator_info_becomes_a_nested_evidence_line(self):
        context = UniprotDomainContext(
            accession="P12883", position=719,
            covering_domains=[DomainFeature("Domain", 700, 750, "Myosin motor")],
        )
        with patch.dict("acmg_pipeline.export._CODE_CURATOR_INFO_FETCHERS",
                         {"PM1": lambda v: context}, clear=True):
            extensions, evidence_items = build_reference_extensions("PM1", variant())

        self.assertFalse(any(e.name == "curatorInfo" for e in extensions))
        self.assertEqual(len(evidence_items), 1)
        item = evidence_items[0]
        self.assertEqual(item["type"], "EvidenceLine")
        self.assertEqual(item["directionOfEvidenceProvided"], "neutral")
        facts_ext = next(e for e in item["extensions"] if e["name"] == "curatorInfo")
        self.assertEqual(facts_ext["value"]["uniprotAccession"], "P12883")
        self.assertEqual(facts_ext["value"]["proteinPosition"], 719)

    def test_no_fetcher_for_the_code_yields_no_evidence_items(self):
        extensions, evidence_items = build_reference_extensions("PS1", variant())
        self.assertEqual(evidence_items, [])

    def test_a_fetcher_that_finds_nothing_yields_no_evidence_items(self):
        with patch.dict("acmg_pipeline.export._CODE_CURATOR_INFO_FETCHERS",
                         {"PM1": lambda v: None}, clear=True):
            extensions, evidence_items = build_reference_extensions("PM1", variant())
        self.assertEqual(evidence_items, [])

    def test_the_item_validates_inside_a_real_evidence_line(self):
        """The whole point: GA4GH's own EvidenceLine model must accept this
        nested item inside hasEvidenceItems without complaint."""
        from ga4gh.va_spec.base.core import EvidenceLine, Method

        context = UniprotDomainContext(accession="P12883", position=719)
        with patch.dict("acmg_pipeline.export._CODE_CURATOR_INFO_FETCHERS",
                         {"PM1": lambda v: context}, clear=True):
            _, evidence_items = build_reference_extensions("PM1", variant())

        line = EvidenceLine(
            id="evline:test", directionOfEvidenceProvided="neutral",
            hasEvidenceItems=evidence_items,
            specifiedBy=Method(methodType="PM1", name="test"),
        )
        self.assertEqual(len(line.hasEvidenceItems), 1)


if __name__ == "__main__":
    unittest.main()
