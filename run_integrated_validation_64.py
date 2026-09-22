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
from datetime import datetime
from contextlib import AsyncExitStack
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import acmg_pipeline.automated_cli as automated_cli_module
import acmg_pipeline.pipeline as pl
from acmg_pipeline.automated_core.input import audit_vcf
from acmg_pipeline.classification import ALL_ACMG_CODES, classification_to_dict, classify
from acmg_pipeline.fulltext_cache import DiskBackedFullTextCache
from acmg_pipeline.gate import ERepoClient
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
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
    parser.add_argument(
        "--only-file", type=Path, default=None,
        help="JSON file with a list of [gene, hgvsc] pairs - restrict the run to just "
             "these variants. The full 64-variant list is still iterated internally "
             "(skipping everything not in this set) so each variant's erepo64:<line>:1 "
             "record_id still matches the cached coordinates in "
             "test_data/erepo_variant_coordinates.json.",
    )
    args = parser.parse_args()

    # One timestamp for the whole run, prefixed onto every output filename so
    # files from different runs sort together and never silently clobber an
    # earlier run's output for the same variant.
    run_ts = datetime.now().strftime("%Y%m%d%H%M%S")

    transcripts = json.loads(TRANSCRIPTS_PATH.read_text(encoding="utf-8"))
    variants = list(unique_variants())
    coords = _resolve_coordinates(variants, transcripts)
    if args.limit:
        variants = variants[: args.limit]

    only_set = None
    if args.only_file:
        only_set = {tuple(pair) for pair in json.loads(args.only_file.read_text(encoding="utf-8"))}
    total_to_process = len(only_set) if only_set is not None else len(variants)

    automated_config = json.loads((ROOT / "config" / "demo-rules.json").read_text(encoding="utf-8"))
    automated_config["evidence_cache_dir"] = str(EVIDENCE_CACHE_DIR)
    # Live (2026-09-18), not offline: EVIDENCE_CACHE_DIR (cache/erepo_automated_
    # evidence/) had never been warmed for most of these 64 variants - not just
    # the providers wired up this session (ClinGen dosage/MANE/NMD/protein-
    # region/splice-default/initiation/gene-disease-draft/clinvar-spectrum) but
    # even the base Ensembl VEP annotation PVS1 needs first. With offline=True
    # every one of those was a silent OFFLINE_CACHE_MISS, so nearly every
    # automated criterion fell through to unknown/not_met regardless of the
    # variant - a run that looked complete (58/64 compared, 13 matched) but
    # never actually exercised the automated engine. A live run writes through
    # to this same cache (CachedHttpClient), so a later run can go back to
    # offline=True and still see this data - same convention run_integrated_
    # validation_demo.py already uses for its own cache.
    automated_config["offline"] = False
    automated_config["ensembl_release"] = "116"
    # PM4/BP3 (UniProt Repeat/Compositional-bias + SEG complexity) and BP7 (VEP
    # splice_region_variant + SpliceAI) - see providers/region_repeat.py and
    # providers/synonymous_assessment.py. Off by default in ProviderEvidenceResolver
    # and, until now, never turned on by this script - pipeline.py's own resolver
    # construction only just started reading these two keys (2026-09-22).
    automated_config["with_region_assessment"] = True
    automated_config["with_bp7_splice_assessment"] = True

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    full_text_cache = DiskBackedFullTextCache("cache/pubmed_fulltext")
    llm_cache = DiskBackedLLMCache("cache/llm_judgments")

    match_n = 0
    compared_n = 0
    skipped = []

    erepo_client = ERepoClient()

    async def process_one(gene, hgvsc, record_id, coord, entry):
        """One variant's full evaluation, on its OWN fresh MCP connection (see the
        `for attempt` loop below for why this is a whole fresh AsyncExitStack per
        attempt, not a connection shared across the run or reused across retries).
        Returns (clinical_note_used, lines) - lines is None on unrecoverable failure.
        """
        try:
            try:
                erepo_lookup = erepo_client.lookup(gene, hgvsc)
            except Exception:
                await asyncio.sleep(5)
                erepo_lookup = erepo_client.lookup(gene, hgvsc)
        except Exception as exc:
            pl.show(f"  -> ERROR looking up {gene} {hgvsc} in ERepo: {exc!r}")
            return None, None, f"{gene} {hgvsc} (ERepo lookup error: {exc!r})"

        if erepo_lookup.condition_id:
            clinical_note = ClinicalNoteExtraction(
                diagnosis=erepo_lookup.condition_label,
                condition_id=erepo_lookup.condition_id,
            )
        else:
            # ERepo has nothing at all for some of this dataset's variants (the
            # democase-sourced entries, e.g. MYBPC3 c.278delA/c.2905+1G>A/c.836del -
            # real HCM variants ERepo simply has never curated). Per the user's
            # explicit direction (2026-09-19): try a live literature search for what
            # disease the variant is reported to cause before giving up to an empty
            # clinical note.
            try:
                from acmg_pipeline import condition_from_literature_search
                lit_condition = await condition_from_literature_search.search_condition_from_literature(
                    gene, hgvsc,
                )
                if lit_condition.found:
                    from acmg_pipeline import hpo_mondo_extraction
                    clinical_note = await hpo_mondo_extraction.resolve_diagnosis_mondo(
                        ClinicalNoteExtraction(diagnosis=lit_condition.condition_name)
                    )
                else:
                    clinical_note = empty_clinical_note()
            except (KeyboardInterrupt, SystemExit):
                raise
            except BaseException as exc:
                # BaseException, not Exception: this opens its OWN separate MCP
                # connection (search_condition_from_literature() -> pipeline.
                # connect_pubmed()) and can hit the same BaseExceptionGroup-not-
                # Exception failure mode the evaluate_variant_evidence_lines() call
                # below does - see that call's own comment. A live search/LLM call
                # failing transiently must not end the whole run - fall back to no
                # clinical note (this variant's PP1/BS4/PP4/PVS1-mechanism results
                # come back UNKNOWN, same as the "not found" case already handles)
                # rather than failing the variant outright, since the automated and
                # PS3/BS3/PS4 evidence lines still do not depend on this.
                pl.show(f"  -> Literature condition search failed for {gene} {hgvsc}: {exc!r}")
                clinical_note = empty_clinical_note()

        # A fresh MCP connection for THIS variant alone, not the one shared across the
        # whole run the first two versions of this fix used. Confirmed necessary,
        # not just cautious, 2026-09-22: reconnecting on the SAME outer AsyncExitStack
        # after a failure stopped working part-way through a run (every subsequent
        # variant failed identically with CancelledError('Cancelled via cancel scope
        # <same id>') - once anyio cancels a scope, everything nested under it stays
        # cancelled, and connect_pubmed() called again on that same stack inherits the
        # poisoned scope instead of getting a clean one). A wholly separate
        # AsyncExitStack per variant can never inherit another variant's cancellation.
        for attempt in (1, 2):
            try:
                async with AsyncExitStack() as variant_stack:
                    mcp = await pl.connect_pubmed(variant_stack)
                    lines = await pl.evaluate_variant_evidence_lines(
                        variant_from(gene, hgvsc, record_id, coord, entry), clinical_note,
                        automated_config=automated_config,
                        mcp=mcp, erepo_client=erepo_client,
                        full_text_cache=full_text_cache, llm_cache=llm_cache,
                    )
                return clinical_note, lines, None
            except (KeyboardInterrupt, SystemExit):
                raise
            except BaseException as exc:
                if attempt == 1:
                    pl.show(f"  -> ERROR evaluating {gene} {hgvsc}: {exc!r} - "
                           f"retrying once on a fresh connection")
                    continue
                pl.show(f"  -> ERROR evaluating {gene} {hgvsc} (after retry): {exc!r}")
                return clinical_note, None, f"{gene} {hgvsc} (error: {exc!r})"

    def variant_from(gene, hgvsc, record_id, coord, entry):
        return VariantRecord(
            chrom=coord["chrom"], pos=coord["pos"], id=record_id,
            ref=coord["ref"], alt=coord["alt"], qual="", filter="",
            info={"GENE": gene, "TRANSCRIPT": entry["transcript"], "HGVSC": hgvsc},
        )

    k = 0
    for i, (gene, hgvsc) in enumerate(variants, 1):
        if only_set is not None and (gene, hgvsc) not in only_set:
            continue
        k += 1
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

        pl.show(f"\n{'#'*70}\n# [{k}/{total_to_process}] {gene} {hgvsc}\n{'#'*70}")
        gt_entries = entries_for(gene, hgvsc)
        outcomes = {e.variant_outcome for e in gt_entries if e.variant_outcome}
        if len(outcomes) != 1:
            skipped.append(f"{gene} {hgvsc} (no single known ground-truth classification)")
            continue
        real_outcome = outcomes.pop()

        clinical_note, lines, failure = await process_one(gene, hgvsc, record_id, coord, entry)
        if failure:
            skipped.append(failure)
            continue

        evidence_lines = dict(zip(ALL_ACMG_CODES, lines))
        evidence = [_evidence_from_line(code, evidence_lines[code]) for code in ALL_ACMG_CODES]
        result = classify(evidence)

        safe = _safe_hgvsc(hgvsc)
        (OUTPUT_DIR / f"{run_ts}_{gene}_{safe}.json").write_text(
            json.dumps(
                {"classification": classification_to_dict(result), "evidence_lines": lines},
                indent=2, ensure_ascii=False,
            ),
            encoding="utf-8",
        )

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
