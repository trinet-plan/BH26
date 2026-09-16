"""
run_validation_64.py

VALIDATION-mode run (see acmg_pipeline/gate.py's GateMode) of the LLM
literature-judgment pipeline across every (variant, criterion) pair in
test_data/full_criteria_ground_truth.py that this project has a real
JudgmentEngine for (PS3/BS3/PS4/PP1/BS4 - see acmg_pipeline.classification.
IMPLEMENTED_CODES). Ground truth is never shown to the LLM; it is only
used afterward to score the result (see score_direction() in pipeline.py) -
that is what "validation mode" means here, matching gate.py's own
GateMode.VALIDATION docstring ("hides the ground truth, forces the
downstream pipeline to run regardless").

227 (variant, criterion) pairs across 60 of the 64 dataset variants (the
other 4 have no ground truth for any of the 5 implemented codes - e.g.
they're PP4/Layer-1-only entries).

[Full-text caching, revised 2026-09-16]
  A variant with ground truth for multiple implemented criteria (e.g.
  RUNX1 c.601C>T: PS4/PP1/BS4) cites the SAME PMIDs for every one of
  those criteria - only the judgment prompt/schema differs per criterion,
  not the paper. An in-memory dict recreated per variant only avoided
  re-fetching WITHIN one variant; it still refetched a PMID shared by two
  DIFFERENT variants, and forgot everything as soon as the script exited.
  Replaced with acmg_pipeline.fulltext_cache.DiskBackedFullTextCache - ONE
  instance for the entire run (not one per variant), persisted to
  cache/pubmed_fulltext/ across runs too, since a paper's full text (or
  the fact that it has none) never changes. See that module's docstring
  for the size-cap safety net (max_files/evict_batch).

Usage:
  python3 run_validation_64.py                  # full 227-pair run
  python3 run_validation_64.py --limit 5         # first 5 variants only (smoke test)
"""

import argparse
import asyncio
import sys
from contextlib import AsyncExitStack
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.classification import IMPLEMENTED_CODES
from acmg_pipeline.fulltext_cache import DiskBackedFullTextCache
from acmg_pipeline.gate import ERepoClient
from test_data.full_criteria_ground_truth import GROUND_TRUTH, GroundTruthEntry

import acmg_pipeline.pipeline as pl


def _variant_criterion_pairs() -> dict[tuple[str, str], list[GroundTruthEntry]]:
    by_variant: dict[tuple[str, str], list[GroundTruthEntry]] = {}
    for e in GROUND_TRUTH:
        if e.criterion in IMPLEMENTED_CODES:
            by_variant.setdefault((e.gene, e.hgvsc), []).append(e)
    return by_variant


def _ground_truth_str(e: GroundTruthEntry) -> str:
    strength_suffix = f"_{e.strength.value.replace('_', ' ').title().replace(' ', '')}" if e.strength else ""
    return f"{e.criterion}{strength_suffix} = {e.status.value}"


async def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="only process the first N variants (smoke test)")
    parser.add_argument("--cache-dir", default="cache/pubmed_fulltext", help="disk cache directory")
    parser.add_argument("--cache-max-files", type=int, default=5000)
    parser.add_argument("--cache-evict-batch", type=int, default=100)
    args = parser.parse_args()

    full_text_cache = DiskBackedFullTextCache(
        args.cache_dir, max_files=args.cache_max_files, evict_batch=args.cache_evict_batch,
    )
    pl.show(f"[cache] {args.cache_dir} ({len(full_text_cache)} entrie(s) already cached, "
            f"max_files={args.cache_max_files}, evict_batch={args.cache_evict_batch})")

    by_variant = _variant_criterion_pairs()
    variant_keys = list(by_variant.keys())
    if args.limit:
        variant_keys = variant_keys[: args.limit]

    total_pairs = sum(len(by_variant[k]) for k in variant_keys)
    pl.show(f"[run_validation_64] {len(variant_keys)} variant(s), {total_pairs} (variant, criterion) pair(s)")

    async with AsyncExitStack() as stack:
        mcp = await pl.connect_pubmed(stack)
        pl.show("[MCP] Connected to PubMed")
        erepo_client = ERepoClient()

        tally = {"match": 0, "mismatch": 0, "reserved": 0}
        mismatches = []
        pair_i = 0

        for gene, hgvsc in variant_keys:
            entries = by_variant[(gene, hgvsc)]
            hgvsp = entries[0].hgvsp
            equivalents = [hgvsp] if hgvsp and hgvsp != "N/A" else [hgvsc]

            # No direct erepo_client.lookup() here - resolve_pmids_for_variant()
            # (ERepo first, live PubMed search fallback) is the single shared
            # PMID resolver every caller should go through, per the user's
            # 2026-09-16 direction. `disease` is unavailable from a
            # GroundTruthEntry (unlike VariantRecord.info["DISEASE_ASSOCIATION"]
            # in judge_variant_from_structured_input()), so the search
            # fallback here only gets the protein-change and generic
            # "novel variant" queries, not the disease-keyword one - a real
            # but accepted degradation given what this script's own data
            # source (test_data/full_criteria_ground_truth.py) carries.
            pmids, _pmid_source = await pl.resolve_pmids_for_variant(mcp, erepo_client, gene, hgvsc, hgvsp)
            if not pmids:
                pl.show("  -> No PMIDs from ERepo or PubMed search; skipping this variant entirely")
                continue

            for e in entries:
                pair_i += 1
                pl.show(f"\n### Pair {pair_i}/{total_pairs}: {gene} {hgvsc} / {e.criterion} ###")
                engine = pl.ENGINE_BY_CRITERION[e.criterion]
                ground_truth_str = _ground_truth_str(e)

                aggregated = await pl.judge_variant(
                    engine, mcp, pmids=pmids, gene=gene, hgvsc=hgvsc, hgvsp=hgvsp,
                    equivalents=equivalents, vcep_name=None, criterion=e.criterion,
                    ground_truth=ground_truth_str, full_text_cache=full_text_cache,
                )
                direction = aggregated.aggregated_direction
                _, gt_is_met = pl.parse_ground_truth(ground_truth_str)
                bucket = pl.score_direction(direction, e.criterion, gt_is_met)
                tally[bucket] += 1
                if bucket == "mismatch":
                    mismatches.append(
                        f"{gene} {hgvsc} ({e.criterion}): ground_truth={ground_truth_str}, "
                        f"LLM direction={direction.value}"
                    )

        pl.show(f"\n{'='*70}\n[Final tally over {sum(tally.values())} pair(s)]\n{'='*70}")
        pl.show(f"match={tally['match']}  mismatch={tally['mismatch']}  reserved(not_clear)={tally['reserved']}")
        committed = tally["match"] + tally["mismatch"]
        if committed:
            pl.show(f"Accuracy among committed judgments: {tally['match']}/{committed} = {tally['match']/committed:.1%}")
        if mismatches:
            pl.show("\nMismatches:")
            for m in mismatches:
                pl.show(f"  - {m}")


if __name__ == "__main__":
    asyncio.run(main())
