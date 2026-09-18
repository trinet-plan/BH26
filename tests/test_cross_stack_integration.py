"""Both halves of the merged stack running in one pass.

`tests/test_main_integration.py` already covers main's 28-line contract, but it
stubs the literature layer out entirely (`judge_variant_from_shared_input` ->
`{}`) and feeds an empty evidence list, so every automated code lands on
NOT_EVALUATED. That shape cannot catch a regression where one side's real
output stops reaching the shared export.

The tests here drive the seam in `evaluate_variant_evidence_lines` with both
halves live at once:

  * evidence-cli (`acmg_pipeline.automated_engine.evaluate_record`) runs unpatched against real
    population and computational evidence chosen to make PM2 and PP3 MET;
  * main's literature path runs its real finalize/aggregate/export chain, with
    only the two genuinely external calls replaced - PMID resolution (ERepo)
    and the per-paper LLM judgment.

Nothing here touches the network or a vLLM endpoint.
"""

import asyncio
import json
import os
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest
from ga4gh.va_spec.base.core import EvidenceLine
from jsonschema import Draft202012Validator

os.environ.setdefault("VLLM_BASE_URL", "http://127.0.0.1:8000/v1")
os.environ.setdefault("VLLM_API_KEY", "test-only")

from acmg_pipeline.services.resolve import StaticEvidenceResolver
from acmg_pipeline.classification import (
    ALL_ACMG_CODES,
    AUTOMATED_CODES,
    LITERATURE_CODES,
)
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.common import MatchStatus, PaperContribution, VariantMatchingResult
from acmg_pipeline.criteria import ps3_bs3
from acmg_pipeline.pipeline import evaluate_variant_evidence_lines
from acmg_pipeline.vcf_record import VariantRecord

GENE = "MYH7"
HGVSC = "c.1594T>C"
TRANSCRIPT = "NM_000257.4"
VARIANT_KEY = "GRCh38:14:23883114:G:A"
PMIDS = ["16983074", "17351073"]


def _variant() -> VariantRecord:
    return VariantRecord(
        chrom="14",
        pos=23883114,
        id="cross-stack:MYH7-S532P",
        ref="G",
        alt="A",
        qual=".",
        filter="PASS",
        info={
            "ASSEMBLY": "GRCh38",
            "GENE": GENE,
            "HGVSC": HGVSC,
            "HGVSP": "p.Ser532Pro",
            "transcript": TRANSCRIPT,
        },
    )


def _normalized_evidence() -> list[dict]:
    """Evidence-cli input sized to drive PM2 (absent from gnomAD) and PP3 (REVEL)."""
    common = {
        "variant_key": VARIANT_KEY,
        "source_version": "4.1.1",
        "retrieved_at": "2026-09-15T00:00:00Z",
        "quality_status": "PASS",
    }
    return [
        {
            **common,
            "category": "population",
            "evidence_id": "https://example.test/evidence/gnomad-absent",
            "source": "gnomAD",
            "population": "global",
            "AC": 0,
            "AN": 120000,
            "AF": 0,
            "callable": True,
        },
        {
            **common,
            "category": "annotation",
            "evidence_id": "https://example.test/evidence/annotation",
            "source": "Ensembl",
            "transcript": TRANSCRIPT,
            "consequences": ["missense_variant"],
        },
        {
            **common,
            "category": "computational",
            "evidence_id": "https://example.test/evidence/revel",
            "source": "dbNSFP",
            "transcript": TRANSCRIPT,
            "predictor": "REVEL",
            "predictor_version": "dbNSFP-4.8a",
            "mechanism": "protein",
            "score": 0.95,
        },
    ]


def _ps3_contribution(pmid: str) -> PaperContribution:
    """A per-paper PS3 judgment that survives main's own consistency guards.

    The numeric `key_findings` matter: `detect_no_quantitative_evidence` forces
    a purely qualitative extraction to not_clear, so a stub without numbers
    would silently test the guard rather than the export path.
    """
    judgment = ps3_bs3.PS3BS3Judgment(
        variant_matching=VariantMatchingResult(match_status=MatchStatus.MATCHED),
        experiments=[
            ps3_bs3.ExperimentExtraction(
                assay_type="load-clamped laser trap assay",
                experimental_system="purified mouse cardiac myosin",
                readout="maximal force-generating capacity F(max)",
                comparator="wild type",
                result_direction=ps3_bs3.ResultDirection.FUNCTIONALLY_ABNORMAL,
                key_findings=["F(max) was depressed to 65% of wild type"],
            )
        ],
        overall_evidence=ps3_bs3.OverallEvidence(
            direction=ps3_bs3.OverallDirection.PS3,
            rationale="Motor function is measurably impaired relative to wild type.",
        ),
    )
    return PaperContribution(
        pmid,
        ps3_bs3.finalize(judgment, GENE, "Cardiomyopathy VCEP", "PS3", pmid=pmid),
    )


async def _fake_judge_single_paper(engine, mcp, pmid, *args, criterion=None, **kwargs):
    """Stand in for the vLLM round trip only; aggregation and export stay real."""
    code = criterion if criterion is not None else args[5]
    return _ps3_contribution(pmid) if code == "PS3" else None


def _run_both_stacks() -> list[dict]:
    config = json.loads(Path("config/demo-rules.json").read_text(encoding="utf-8"))
    with (
        patch(
            "acmg_pipeline.pipeline.resolve_pmids_for_variant",
            new=AsyncMock(return_value=(PMIDS, "erepo")),
        ),
        patch(
            "acmg_pipeline.pipeline.judge_single_paper",
            new=_fake_judge_single_paper,
        ),
        patch(
            "acmg_pipeline.export.reference_links.reference_urls_for_criterion",
            side_effect=lambda code, _variant: [f"https://example.test/{code}"],
        ),
        patch.dict("acmg_pipeline.export._CODE_CURATOR_INFO_FETCHERS", {}, clear=True),
    ):
        return asyncio.run(
            evaluate_variant_evidence_lines(
                _variant(),
                ClinicalNoteExtraction(),
                automated_config=config,
                evidence_resolver=StaticEvidenceResolver(_normalized_evidence()),
                mcp=object(),
                erepo_client=object(),
                vcep_name="Cardiomyopathy VCEP",
            )
        )


@pytest.fixture(scope="module")
def lines_by_code() -> dict[str, dict]:
    lines = _run_both_stacks()
    assert tuple(line["specifiedBy"]["methodType"] for line in lines) == ALL_ACMG_CODES
    return {line["specifiedBy"]["methodType"]: line for line in lines}


def _details(line: dict) -> dict:
    """Each detail field (status/direction/...) is its own top-level
    extension now, not grouped under one bh26AssessmentDetails object
    (2026-09-18, per the user's direction)."""
    from acmg_pipeline.automated_va_spec import details_from_extensions
    return details_from_extensions(line["extensions"])


def test_evidence_cli_criteria_are_really_evaluated(lines_by_code):
    assert _details(lines_by_code["PM2"])["status"] == "met"
    assert _details(lines_by_code["PP3"])["status"] == "met"
    assert lines_by_code["PM2"]["directionOfEvidenceProvided"] == "supports"
    assert lines_by_code["PP3"]["directionOfEvidenceProvided"] == "supports"


def test_literature_criteria_are_really_evaluated(lines_by_code):
    ps3 = lines_by_code["PS3"]
    assert _details(ps3)["status"] == "met"
    assert ps3["directionOfEvidenceProvided"] == "supports"
    assert {doc["pmid"] for doc in ps3["reportedIn"]} == set(PMIDS)
    assert len(ps3["hasEvidenceItems"]) == len(PMIDS)


def test_both_stacks_appear_in_one_document(lines_by_code):
    """The point of the merge: neither side may shadow or duplicate the other."""
    automated_met = {
        code for code in AUTOMATED_CODES
        if _details(lines_by_code[code])["status"] == "met"
    }
    literature_met = {
        code for code in LITERATURE_CODES
        if _details(lines_by_code[code])["status"] == "met"
    }
    assert automated_met and literature_met
    assert automated_met.isdisjoint(literature_met)
    assert len(lines_by_code) == len(ALL_ACMG_CODES) == 28
    assert len({line["id"] for line in lines_by_code.values()}) == 28


def test_combined_document_is_valid_va_spec(lines_by_code):
    schema_validator = Draft202012Validator(EvidenceLine.model_json_schema())
    for code, line in lines_by_code.items():
        EvidenceLine.model_validate(line)
        schema_validator.validate(line)
        assert line["specifiedBy"]["methodType"] == code


def test_every_line_passes_the_1_0_1_validator(lines_by_code):
    """Not just the 16 automated codes' scored lines, as before the merge."""
    from acmg_pipeline.export import validate_integrated_line

    for code, line in lines_by_code.items():
        validate_integrated_line(line, code)


def test_all_28_ids_use_the_single_evline_scheme(lines_by_code):
    """One id form for the whole document, independent of side and status.

    evidence-cli's own `urn:bh26:evidence-line:<sha256>` form hashed the
    evidence outcome, so a threshold change moved the id of an unchanged
    variant/criterion. Every line is now addressable as
    evline:{gene}_{safe_hgvsc}_{code}.
    """
    from acmg_pipeline.export import evidence_line_id

    safe_hgvsc = HGVSC.replace(">", "_").replace(".", "_")
    for code, line in lines_by_code.items():
        assert line["id"] == evidence_line_id(code, GENE, safe_hgvsc)
    assert len({line["id"] for line in lines_by_code.values()}) == 28


def test_scored_lines_keep_direction_and_outcome_consistent(lines_by_code):
    """The ACMG cross-field rule now covers the literature codes too."""
    scored = {
        code: line for code, line in lines_by_code.items()
        if line.get("evidenceOutcome")
    }
    assert "PS3" in scored and "PM2" in scored  # one from each half
    for code, line in scored.items():
        outcome = line["evidenceOutcome"]["primaryCoding"]["code"]
        expected = "neutral" if outcome.endswith("_not_met") else (
            "disputes" if code.startswith("B") else "supports"
        )
        assert line["directionOfEvidenceProvided"] == expected


def test_default_strength_agrees_with_evidence_cli_where_they_overlap():
    """The prefix rule main derives must not drift from evidence-cli's map."""
    from acmg_pipeline.criteria.common import DEFAULT_STRENGTH
    from acmg_pipeline.export import default_strength

    for code, expected in DEFAULT_STRENGTH.items():
        assert default_strength(code) == expected, code
    for code in ALL_ACMG_CODES:
        default_strength(code)  # every one of the 28 is answerable


def test_outcome_code_omits_the_default_strength_suffix():
    """'PP1_supporting' is not a legal ACMG outcome code; plain 'PP1' is."""
    from acmg_pipeline.common import PaperContribution, MatchStatus, VariantMatchingResult
    from acmg_pipeline.criteria import segregation
    from acmg_pipeline.export import build_evidence_line

    judgment = segregation.SegregationJudgment(
        variant_matching=VariantMatchingResult(match_status=MatchStatus.MATCHED),
        families=[segregation.FamilySegregationData(
            family_id="F1", affected_with_variant=3, affected_without_variant=0,
            unaffected_with_variant=0, unaffected_without_variant=2,
        )],
        overall_evidence=segregation.SegregationEvidence(
            direction=segregation.SegregationDirection.PP1,
            rationale="Variant segregates with disease.",
        ),
    )
    aggregated = segregation.aggregate_multi_paper_results([
        PaperContribution("12345678", segregation.finalize(judgment, pmid="12345678"))
    ])
    with patch(
        "acmg_pipeline.export.reference_links.reference_urls_for_criterion",
        return_value=["https://example.test/PP1"],
    ):
        line = build_evidence_line(aggregated, "MYH7", "c.2155C>T", "PP1")

    assert line["evidenceOutcome"]["primaryCoding"]["code"] == "PP1"
    assert line["strengthOfEvidenceProvided"]["primaryCoding"]["code"] == "supporting"


def test_unimplemented_codes_stay_placeholders(lines_by_code):
    """PS2/PM3/BP2 and the rest belong to neither side and must not be faked."""
    for code in set(ALL_ACMG_CODES) - AUTOMATED_CODES - LITERATURE_CODES:
        assert _details(lines_by_code[code])["status"] != "met"


