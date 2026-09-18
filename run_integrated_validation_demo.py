"""
run_integrated_validation_demo.py

Runs the REAL integrated pipeline (acmg_pipeline.pipeline.
evaluate_variant_evidence_lines() - live PubMed+LLM literature judgment for
PS3/BS3/PS4, plus the real automated engine for the 16 Layer-1 codes, plus
stubs for the rest) for the 4 demo-data primary variants, USING THEIR REAL
CLINICAL NOTES (unlike run_integrated_validation_64.py, which has no
clinical note text for the ERepo-only variants and passes an empty one) -
a smaller, faster first pass than the full 64-variant run, per the user's
request (2026-09-17) to start with just the 4 demo cases.

The 4 target variants (each demo case's primary/case-defining variant,
not the "noise" background variants also present in demo-data):
  case1-var1: MYBPC3 c.278delA
  case2-var1: MYBPC3 c.2905+1G>A
  case3-var1: MYH7 c.2155C>T
  case4-var1: DSG2 c.1592T>G

Real GRCh38 coordinates are reused from the already-committed offline
fixture cache (tests/fixtures/ensembl-cache/external-cache) via the same
prepare-demo-online resolution run_automated_validation_64.py already
does - see this script's _resolve_demo_coordinates().
"""

from __future__ import annotations

import asyncio
import json
import sys
import tempfile
from contextlib import AsyncExitStack
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import acmg_pipeline.automated_cli as automated_cli_module
import acmg_pipeline.pipeline as pl
from acmg_pipeline.classification import ALL_ACMG_CODES, classify
from acmg_pipeline.clinical_note import extract_clinical_note
from acmg_pipeline.fulltext_cache import DiskBackedFullTextCache
from acmg_pipeline.gate import ERepoClient
from acmg_pipeline.llm_cache import DiskBackedLLMCache
from acmg_pipeline.pipeline_interface import _evidence_from_line
from acmg_pipeline.vcf_record import VariantRecord
from test_data.full_criteria_ground_truth import entries_for

ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "va_spec_output" / "integrated"

CASES = [
    ("case1", "case1-var1", "MYBPC3", "c.278delA", "case1_clinical_note_v1.txt"),
    ("case2", "case2-var1", "MYBPC3", "c.2905+1G>A", "case2_clinical_note_v1.txt"),
    ("case3", "case3-var1", "MYH7", "c.2155C>T", "case3_clinical_note_v1.txt"),
    ("case4", "case4-var1", "DSG2", "c.1592T>G", "case4_clinical_note_v2.txt"),
]


def _safe_hgvsc(hgvsc: str) -> str:
    return hgvsc.replace(">", "_").replace(".", "_").replace("+", "p").replace("-", "m").replace("*", "s")


def _resolve_demo_coordinates() -> dict:
    """Real GRCh38 coordinates + transcript for every demo-data ALT record, offline/cached."""
    work = Path(tempfile.mkdtemp())
    prepared = work / "prepared"
    exit_code = automated_cli_module.main([
        "prepare-demo-online", "--input-dir", str(ROOT / "demo-data"),
        "--cache-dir", str(ROOT / "tests" / "fixtures" / "ensembl-cache"),
        "--evidence-cache-dir", str(ROOT / "tests" / "fixtures" / "external-cache"),
        "--output-dir", str(prepared), "--ensembl-release", "116",
        "--with-gnomad", "--gnomad-release", "4.1.1",
        "--with-clinvar", "--clinvar-release", "2026-09-15",
        "--with-pm1-hotspot", "--rules", str(ROOT / "config" / "demo-rules.json"),
        "--with-dbnsfp", "--offline",
    ])
    if exit_code:
        raise RuntimeError(f"prepare-demo-online exited {exit_code}")
    resolved = json.loads((prepared / "variants.json").read_text(encoding="utf-8"))
    by_variant_id = {}
    for rec in resolved["records"]:
        matched = rec["identity_provenance"][0]["matched_identifiers"]
        by_variant_id[rec["source"]["variant_id"]] = {
            "variant": rec["variant"],
            "transcript": matched.get("TRANSCRIPT"),
        }
    return by_variant_id


async def main() -> None:
    # One timestamp for the whole run, prefixed onto every output filename so
    # files from different runs sort together and never silently clobber an
    # earlier run's output for the same variant.
    run_ts = datetime.now().strftime("%Y%m%d%H%M%S")

    coords = _resolve_demo_coordinates()

    automated_config = json.loads((ROOT / "config" / "demo-rules.json").read_text(encoding="utf-8"))
    automated_config["evidence_cache_dir"] = str(ROOT / "tests" / "fixtures" / "external-cache")
    # Live for this run (2026-09-17), not the usual offline/cached demo
    # convention: population_sources just switched from gnomAD to TogoVar
    # (h.muroda, "Add configurable TogoVar population frequency providers"),
    # and tests/fixtures/external-cache has zero TogoVar entries yet. A live
    # fetch here also warms that cache (CachedHttpClient writes through),
    # so a later run can go back to offline=True and still see this data.
    automated_config["offline"] = False
    automated_config["ensembl_release"] = "116"

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    full_text_cache = DiskBackedFullTextCache("cache/pubmed_fulltext")
    llm_cache = DiskBackedLLMCache("cache/llm_judgments")

    match_n = 0
    compared_n = 0

    async with AsyncExitStack() as stack:
        mcp = await pl.connect_pubmed(stack)
        pl.show("[MCP] Connected to PubMed")
        erepo_client = ERepoClient()

        for case_id, variant_id, gene, hgvsc, note_file in CASES:
            pl.show(f"\n{'#'*70}\n# {case_id} ({variant_id}): {gene} {hgvsc}\n{'#'*70}")
            coord_entry = coords.get(variant_id)
            if not coord_entry:
                pl.show(f"  -> ERROR: no resolved coordinates for {variant_id}, skipping")
                continue

            note_text = (ROOT / "demo-data" / note_file).read_text(encoding="utf-8")
            clinical_note = extract_clinical_note(note_text)

            variant = VariantRecord(
                chrom=coord_entry["variant"]["chrom"], pos=coord_entry["variant"]["pos"],
                id=variant_id, ref=coord_entry["variant"]["ref"], alt=coord_entry["variant"]["alt"],
                qual="", filter="",
                info={"GENE": gene, "TRANSCRIPT": coord_entry["transcript"], "HGVSC": hgvsc},
            )

            lines = await pl.evaluate_variant_evidence_lines(
                variant, clinical_note,
                automated_config=automated_config,
                mcp=mcp, erepo_client=erepo_client,
                full_text_cache=full_text_cache, llm_cache=llm_cache,
            )

            safe = _safe_hgvsc(hgvsc)
            (OUTPUT_DIR / f"{run_ts}_{gene}_{safe}.json").write_text(
                json.dumps(lines, indent=2, ensure_ascii=False), encoding="utf-8",
            )

            evidence_lines = dict(zip(ALL_ACMG_CODES, lines))
            evidence = [_evidence_from_line(code, evidence_lines[code]) for code in ALL_ACMG_CODES]
            result = classify(evidence)

            gt_entries = entries_for(gene, hgvsc)
            outcomes = {e.variant_outcome for e in gt_entries if e.variant_outcome}
            met_str = ", ".join(f"{e.code}({e.strength.value})" for e in result.met) or "(none)"
            pl.show(f"\n[Integrated classify()] {gene} {hgvsc}: category={result.category.value} "
                    f"(score={result.score})")
            pl.show(f"  met: {met_str}")
            if len(outcomes) == 1:
                real_outcome = outcomes.pop()
                compared_n += 1
                is_match = result.category.value == real_outcome
                match_n += is_match
                pl.show(f"  vs. real={real_outcome} -> {'MATCH' if is_match else 'DIFFERS'}")
            else:
                pl.show(f"  (no single known ground-truth classification to compare against: {outcomes})")

    pl.show(f"\n{'='*70}\n[run_integrated_validation_demo] Summary\n{'='*70}")
    if compared_n:
        pl.show(f"{compared_n} variant(s) compared; {match_n}/{compared_n} matched ({match_n/compared_n:.1%})")


if __name__ == "__main__":
    asyncio.run(main())
