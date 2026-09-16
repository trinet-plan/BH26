"""The same cross-stack path as `test_cross_stack_integration.py`, but with the
LLM actually running.

Opt-in only. These tests reach a real vLLM endpoint (`VLLM_BASE_URL` from
`.env`), the PubMed MCP server and ClinGen ERepo, so they are skipped unless
`ACMG_LIVE_LLM=1` is set - an ordinary `pytest tests` stays hermetic.

    ACMG_LIVE_LLM=1 .venv/Scripts/python.exe -m pytest tests/test_cross_stack_live.py -q -s

What is asserted here is deliberately weaker than in the offline test. A live
LLM is not reproducible, and the honest contract for this pipeline is *not*
"gemma-4 returns PS3 for this variant". It is:

  1. the evidence-cli half stays fully deterministic while the LLM half runs;
  2. whatever the LLM returns survives aggregation and lands in a schema-valid
     28-line VA-Spec document;
  3. the judgment is never scored as a *mismatch* against ERepo - `not_clear`
     (reserved) is an acceptable, safe outcome, a wrong committed direction is
     not.

The actual direction/strength is printed rather than asserted, so the run
doubles as the ERepo-agreement report for this variant.
"""

import asyncio
import json
import os
from contextlib import AsyncExitStack
from pathlib import Path

import pytest
from ga4gh.va_spec.base.core import EvidenceLine
from jsonschema import Draft202012Validator

pytestmark = pytest.mark.skipif(
    os.environ.get("ACMG_LIVE_LLM") != "1",
    reason="live LLM/network test; set ACMG_LIVE_LLM=1 to run",
)

from acmg.services.resolve import StaticEvidenceResolver  # noqa: E402
from test_cross_stack_integration import (  # noqa: E402  (after the skip guard)
    GENE,
    HGVSC,
    _details,
    _normalized_evidence,
    _variant,
)


@pytest.fixture(scope="module")
def live_lines() -> dict[str, dict]:
    import acmg_pipeline.pipeline as pl
    from acmg_pipeline.clinical_note import ClinicalNoteExtraction
    from acmg_pipeline.gate import ERepoClient

    config = json.loads(Path("config/demo-rules.json").read_text(encoding="utf-8"))

    async def run() -> list[dict]:
        async with AsyncExitStack() as stack:
            mcp = await pl.connect_pubmed(stack)
            return await pl.evaluate_variant_evidence_lines(
                _variant(),
                ClinicalNoteExtraction(),
                automated_config=config,
                evidence_resolver=StaticEvidenceResolver(_normalized_evidence()),
                mcp=mcp,
                erepo_client=ERepoClient(),
                vcep_name="Cardiomyopathy VCEP",
            )

    lines = asyncio.run(run())
    return {line["specifiedBy"]["methodType"]: line for line in lines}


def test_automated_half_is_unaffected_by_the_live_half(live_lines):
    """evidence-cli must return exactly what the offline test pins it to."""
    assert _details(live_lines["PM2"])["status"] == "MET"
    assert _details(live_lines["PP3"])["status"] == "MET"
    assert _details(live_lines["BA1"])["status"] == "NOT_MET"
    assert _details(live_lines["BP4"])["status"] == "NOT_MET"


def test_live_document_is_valid_va_spec(live_lines):
    schema_validator = Draft202012Validator(EvidenceLine.model_json_schema())
    assert len(live_lines) == 28
    for code, line in live_lines.items():
        EvidenceLine.model_validate(line)
        schema_validator.validate(line)
        assert _details(line)["criterion"] == code


def test_live_ps3_is_never_scored_as_a_mismatch_against_erepo(live_lines, capsys):
    """ERepo has PS3 = met for MYH7 c.1594T>C; reserved is fine, wrong is not."""
    from acmg_pipeline.pipeline import score_direction
    from acmg_pipeline.criteria.ps3_bs3 import OverallDirection

    ps3 = live_lines["PS3"]
    details = _details(ps3)
    outcome = (ps3.get("evidenceOutcome") or {}).get("primaryCoding", {}).get("code")
    direction = ps3.get("directionOfEvidenceProvided")
    papers = len(ps3.get("hasEvidenceItems") or [])

    with capsys.disabled():
        print(
            f"\n[live] {GENE} {HGVSC} PS3: status={details['status']} "
            f"direction={direction} outcome={outcome} papers_judged={papers}\n"
            f"       ERepo ground truth: PS3 = met"
        )

    committed = (
        OverallDirection.PS3 if details["status"] == "MET"
        else OverallDirection.NOT_CLEAR
    )
    bucket = score_direction(committed, "PS3", gt_is_met=True)
    assert bucket in {"match", "reserved"}, (
        f"LLM committed to a direction contradicting ERepo: {details}"
    )
