"""
run_integrated_validation_64.py

Runs the REAL integrated pipeline (acmg_pipeline.pipeline.
evaluate_variant_evidence_lines() - live PubMed+LLM literature judgment for
PS3/BS3/PS4, plus the real automated engine for the 16 Layer-1 codes, plus
stubs for the rest) across the 64 ERepo/democase variants in test_data/
full_criteria_ground_truth.py, and compares classification.classify()'s
resulting overall ACMG category against each variant's real, already-known
classification.

[Why this exists, and how it differs from the other two 64-variant scripts]
  run_validation_64.py validates the literature engine alone (PS3/BS3/PS4)
  against per-criterion ground truth. run_automated_validation_64.py
  validates the automated engine alone (16 codes) the same way. Neither
  produces the actual thing this project's own output is supposed to be:
  one 28-EvidenceLine VA-Spec document per variant, both engines combined,
  fed through classify() for a real overall category. This script is that
  - per the user's explicit direction (2026-09-17): "統合したva_specを出力
  するのが目的なので、それに合わせて欲しい" (the goal is producing the
  INTEGRATED va_spec output, so align the validation with that).

[Why this needs real, resolved genomic coordinates, unlike the automated-
 only script's ALT="." placeholder trick]
  run_automated_validation_64.py's synthetic VCF rows use ALT="." (parsing
  fails on purpose, forcing identity resolution purely from HGVS - see that
  script's docstring). acmg_pipeline.pipeline._automated_variant() does NOT
  go through that same reconcile()-based correction path -
  evaluate_variant_evidence_lines() constructs acmg_pipeline.automated_core.
  models.Variant directly from variant.chrom/pos/ref/alt, which validates
  ref/alt as real ACGT alleles immediately (see that class's __post_init__)
  - a placeholder ALT="." would raise ValueError there. So this script
  first resolves real GRCh38 coordinates for all 64 variants the same way
  run_automated_validation_64.py does (prepare-demo-online's Ensembl VEP
  HGVS resolution, monkeypatching audit_demo the same way), caches them to
  test_data/erepo_variant_coordinates.json, and uses those real coordinates
  to build each variant's VariantRecord here.

[Cost]
  Every variant makes REAL PubMed MCP + LLM calls for its PS3/BS3/PS4
  literature judgment (not offline/cacheable the way the automated side
  is) - full text IS cached (cache/pubmed_fulltext/, shared with run_
  validation_64.py) once fetched, but the LLM call itself runs every time.
  Expect several minutes per variant with multiple cited papers (see run_
  validation_64.py's own historical timings) - this is NOT a quick offline
  check like the other two 64-variant scripts. Use --limit for a smoke
  test before running the full 64.

[What gets saved]
  Each variant's real 28-line VA-Spec EvidenceLine array is written to
  va_spec_output/integrated/{gene}_{safe_hgvsc}.json - this project's own
  actual deliverable, not just a validation artifact.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import sys
import tempfile
from contextlib import AsyncExitStack
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import acmg_pipeline.automated_cli as automated_cli_module
import acmg_pipeline.pipeline as pl
from acmg_pipeline.automated_core.input import audit_vcf
from acmg_pipeline.classification import ALL_ACMG_CODES, classify
from acmg_pipeline.fulltext_cache import DiskBackedFullTextCache
from acmg_pipeline.gate import ERepoClient
from acmg_pipeline.inputs import empty_clinical_note
from acmg_pipeline.llm_cache import DiskBackedLLMCache
from acmg_pipeline.pipeline_interface import _evidence_from_line
from acmg_pipeline.vcf_record import VariantRecord
from test_data.full_criteria_ground_truth import entries_for, unique_variants

ROOT = Path(__file__).resolve().parent
TRANSCRIPTS_PATH = ROOT / "test_data" / "erepo_variant_transcripts.json"
COORDS_PATH = ROOT / "test_data" / "erepo_variant_coordinates.json"
CACHE_DIR = ROOT / "cache" / "erepo_automated_ensembl"
EVIDENCE_CACHE_DIR = ROOT / "cache" / "erepo_automated_evidence"
OUTPUT_DIR = ROOT / "va_spec_output" / "integrated"


def _safe_hgvsc(hgvsc: str) -> str:
    return hgvsc.replace(">", "_").replace(".", "_").replace("+", "p").replace("-", "m").replace("*", "s")


def _resolve_coordinates(variants: list[tuple[str, str]], transcripts: dict) -> dict:
    """Real GRCh38 coordinates per erepo64:<line>:1 record_id (cached to disk)."""
    if COORDS_PATH.is_file():
        return json.loads(COORDS_PATH.read_text(encoding="utf-8"))

    lines = ["##fileformat=VCFv4.2", "##reference=GRCh38", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for i, (gene, hgvsc) in enumerate(variants, 1):
        entry = transcripts.get(f"{gene}|{hgvsc}")
        if not entry or not entry.get("transcript"):
            continue
        info = f"GENE={gene};TRANSCRIPT={entry['transcript']};HGVSC={hgvsc}"
        lines.append(f"1\t1\terepo64-{i}\tA\t.\t80\tPASS\t{info}")
    vcf_path = Path(tempfile.mkdtemp()) / "synthetic.vcf"
    vcf_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def _audit_demo_override(_directory):
        return audit_vcf(vcf_path, case_id="erepo64", assembly="GRCh38")

    automated_cli_module.audit_demo = _audit_demo_override
    work = Path(tempfile.mkdtemp())
    prepared = work / "prepared"
    exit_code = automated_cli_module.main([
        "prepare-demo-online", "--input-dir", str(work),
        "--cache-dir", str(CACHE_DIR), "--evidence-cache-dir", str(EVIDENCE_CACHE_DIR),
        "--output-dir", str(prepared), "--ensembl-release", "116",
        "--with-gnomad", "--gnomad-release", "4.1.1",
        "--with-clinvar", "--clinvar-release", "2026-09-15",
        "--with-pm1-hotspot", "--rules", str(ROOT / "config" / "demo-rules.json"),
        "--with-dbnsfp", "--offline",
    ])
    if exit_code:
        raise RuntimeError(f"prepare-demo-online exited {exit_code}")
    resolved = json.loads((prepared / "variants.json").read_text(encoding="utf-8"))
    coords = {rec["record_id"]: rec["variant"] for rec in resolved["records"]}
    COORDS_PATH.write_text(json.dumps(coords, indent=2, sort_keys=True), encoding="utf-8")
    return coords


async def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="only process the first N variants (smoke test)")
    args = parser.parse_args()

    transcripts = json.loads(TRANSCRIPTS_PATH.read_text(encoding="utf-8"))
    variants = list(unique_variants())
    coords = _resolve_coordinates(variants, transcripts)
    if args.limit:
        variants = variants[: args.limit]

    automated_config = json.loads((ROOT / "config" / "demo-rules.json").read_text(encoding="utf-8"))
    automated_config["evidence_cache_dir"] = str(EVIDENCE_CACHE_DIR)
    automated_config["offline"] = True
    automated_config["ensembl_release"] = "116"

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    full_text_cache = DiskBackedFullTextCache("cache/pubmed_fulltext")
    llm_cache = DiskBackedLLMCache("cache/llm_judgments")

    match_n = 0
    compared_n = 0
    skipped = []

    async with AsyncExitStack() as stack:
        mcp = await pl.connect_pubmed(stack)
        pl.show("[MCP] Connected to PubMed")
        erepo_client = ERepoClient()

        for i, (gene, hgvsc) in enumerate(variants, 1):
            # _resolve_coordinates()'s synthetic VCF has 3 header/meta lines
            # (##fileformat, ##reference, #CHROM...) before the first data
            # row, and audit_vcf() builds record_id from the real file LINE
            # NUMBER (not this loop's 1-based variant index) - so variant i
            # is always at line i+3.
            record_id = f"erepo64:{i + 3}:1"
            coord = coords.get(record_id)
            entry = transcripts.get(f"{gene}|{hgvsc}")
            if not coord or not entry or not entry.get("transcript"):
                skipped.append(f"{gene} {hgvsc} (no resolved coordinates)")
                continue

            pl.show(f"\n{'#'*70}\n# [{i}/{len(variants)}] {gene} {hgvsc}\n{'#'*70}")
            gt_entries = entries_for(gene, hgvsc)
            outcomes = {e.variant_outcome for e in gt_entries if e.variant_outcome}
            if len(outcomes) != 1:
                skipped.append(f"{gene} {hgvsc} (no single known ground-truth classification)")
                continue
            real_outcome = outcomes.pop()

            variant = VariantRecord(
                chrom=coord["chrom"], pos=coord["pos"], id=record_id,
                ref=coord["ref"], alt=coord["alt"], qual="", filter="",
                info={"GENE": gene, "TRANSCRIPT": entry["transcript"], "HGVSC": hgvsc},
            )

            try:
                lines = await pl.evaluate_variant_evidence_lines(
                    variant, empty_clinical_note(),
                    automated_config=automated_config,
                    mcp=mcp, erepo_client=erepo_client,
                    full_text_cache=full_text_cache, llm_cache=llm_cache,
                )
            except Exception as exc:
                pl.show(f"  -> ERROR evaluating {gene} {hgvsc}: {exc!r}")
                skipped.append(f"{gene} {hgvsc} (error: {exc!r})")
                continue

            safe = _safe_hgvsc(hgvsc)
            (OUTPUT_DIR / f"{gene}_{safe}.json").write_text(
                json.dumps(lines, indent=2, ensure_ascii=False), encoding="utf-8",
            )

            evidence_lines = dict(zip(ALL_ACMG_CODES, lines))
            evidence = [_evidence_from_line(code, evidence_lines[code]) for code in ALL_ACMG_CODES]
            result = classify(evidence)
            compared_n += 1
            is_match = result.category.value == real_outcome
            match_n += is_match
            met_str = ", ".join(f"{e.code}({e.strength.value})" for e in result.met) or "(none)"
            pl.show(f"\n[Integrated classify()] {gene} {hgvsc}: category={result.category.value} "
                    f"(score={result.score}) vs. real={real_outcome} "
                    f"-> {'MATCH' if is_match else 'DIFFERS'}")
            pl.show(f"  met: {met_str}")

    pl.show(f"\n{'='*70}\n[run_integrated_validation_64] Summary\n{'='*70}")
    pl.show(f"{compared_n} variant(s) compared; {match_n}/{compared_n} matched "
            f"({match_n/compared_n:.1%})" if compared_n else "0 variant(s) compared")
    if skipped:
        pl.show(f"\n{len(skipped)} variant(s) skipped:")
        for s in skipped:
            pl.show(f"  - {s}")


if __name__ == "__main__":
    asyncio.run(main())
