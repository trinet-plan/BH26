"""
test_data/resolve_erepo_transcripts.py

One-time (re-runnable) data-preparation step for validating the 16
automated Layer-1 criteria (acmg_pipeline.classification.AUTOMATED_CODES)
against test_data/full_criteria_ground_truth.py's 64 ERepo/democase-derived
variants - NOT just the 4 democase cases test_automated_criteria_ground_
truth.py already covers.

[The problem this solves]
  full_criteria_ground_truth.py only carries (gene, hgvsc) - no RefSeq
  transcript accession, no genomic coordinates. The automated engine's
  identity resolution (acmg_pipeline.automated_core.identity.reconcile())
  needs a TRANSCRIPT to query Ensembl VEP's HGVS endpoint
  (providers/ensembl.py: f"{transcript}:{hgvsc}") - c. notation is
  transcript-relative, so pairing the right hgvsc with the WRONG
  transcript accession can silently resolve to a different genomic
  position or fail outright.

[Why not just use MANE Select for every gene]
  Tried first, via NCBI's MANE.GRCh38.v1.5.summary.txt (fetched fresh,
  https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/) - it does
  NOT always agree with the transcript ERepo's own curation actually used
  (e.g. TNNT2: MANE Select is NM_001276345.2, but the democase VCF's own
  curated ACMG_CODES were asserted against NM_000364.4 - a real, known
  clinical-transcript-nomenclature discrepancy for this gene). Silently
  using MANE Select everywhere risks resolving several variants against
  the wrong transcript without any error (VEP would happily accept a
  syntactically valid HGVS on a real but wrong transcript).

[The actual method - real ERepo data corroborates the transcript]
  ERepo's classifications API (the SAME endpoint acmg_pipeline.gate.
  ERepoClient already calls, live-verified 2026-09-15) returns a full
  "hgvs" array with EVERY equivalent representation ERepo itself computed
  for that variant, including transcript-qualified ones, e.g. for MYH7
  c.2155C>T: ["NM_000257.4:c.2155C>T", "NC_000014.9:g.23425971G>A", ...].
  For each of the 64 variants, this script:
    1. Queries ERepo live (same GET .../evrepo/api/classifications?gene=
       &hgvs= endpoint) and collects every "NM_...:c...." entry.
    2. Takes the MANE Select accession (from the fetched summary file,
       cached locally as gene_transcripts_mane.json) as the preferred
       transcript, and looks for an EXACT "{mane_accession}:{hgvsc}"
       match in ERepo's own list first (highest confidence - ERepo's own
       curation and MANE Select agree).
    3. Falls back to the same base accession at a different version if
       an exact version match isn't found (RefSeq transcript versions are
       usually coding-sequence-stable across minor version bumps).
    4. Falls back to MANE Select alone (not corroborated by ERepo) only
       if ERepo has no matching entry at all for this variant (true for
       some democase-sourced entries ERepo never curated).
  Each variant's resolution_method is recorded in the output JSON so a
  downstream consumer (run_automated_validation_64.py) can weight/flag
  results accordingly - a "mane_fallback" entry is a real, lower-
  confidence resolution, not equivalent to an "exact_erepo_match" one.

[Output]
  test_data/erepo_variant_transcripts.json: {f"{gene}|{hgvsc}": {
    "transcript": "NM_...", "resolution_method": "exact_erepo_match" |
    "erepo_other_version" | "mane_fallback" | "unresolved", "mane_accession":
    "NM_..."}}. Committed to the repo (small, ~64 entries) so downstream
  scripts don't need to re-query ERepo/NCBI every run - re-run this script
  directly if the ground-truth dataset's variant list changes.
"""

from __future__ import annotations

import gzip
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import requests

from test_data.collectors.full_criteria_ground_truth import unique_variants

ROOT = Path(__file__).resolve().parent.parent.parent
MANE_CACHE = ROOT / "test_data" / "fetched_data" / "gene_transcripts_mane.json"
OUTPUT_PATH = ROOT / "test_data" / "fetched_data" / "erepo_variant_transcripts.json"
MANE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"


def _latest_mane_summary_url() -> str:
    with urllib.request.urlopen(MANE_URL, timeout=30) as resp:
        listing = resp.read().decode("utf-8", errors="ignore")
    names = sorted(set(re.findall(r"MANE\.GRCh38\.v[0-9.]+\.summary\.txt\.gz", listing)))
    if not names:
        raise RuntimeError("Could not find a MANE summary file in the NCBI FTP listing")
    return MANE_URL + names[-1]


def _fetch_mane_transcripts(genes: set[str]) -> dict[str, str]:
    if MANE_CACHE.is_file():
        cached = json.loads(MANE_CACHE.read_text(encoding="utf-8"))
        if genes.issubset(cached):
            return cached
    url = _latest_mane_summary_url()
    print(f"[resolve_erepo_transcripts] Fetching MANE summary: {url}")
    with urllib.request.urlopen(url, timeout=60) as resp:
        raw = gzip.decompress(resp.read()).decode("utf-8")
    mane: dict[str, str] = {}
    lines = raw.splitlines()
    header = lines[0].lstrip("#").split("\t")
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        if row.get("symbol") in genes and row.get("MANE_status") == "MANE Select":
            mane[row["symbol"]] = row["RefSeq_nuc"]
    missing = genes - set(mane)
    if missing:
        print(f"[resolve_erepo_transcripts] WARNING: no MANE Select transcript for: {sorted(missing)}")
    MANE_CACHE.write_text(json.dumps(mane, indent=2, sort_keys=True), encoding="utf-8")
    return mane


def _erepo_hgvs_list(gene: str, hgvsc: str) -> list[str]:
    resp = requests.get(
        "https://erepo.clinicalgenome.org/evrepo/api/classifications",
        params={"gene": gene, "hgvs": hgvsc}, timeout=15,
    )
    resp.raise_for_status()
    data = resp.json()
    out: list[str] = []
    for interp in data.get("variantInterpretations", []):
        out.extend(interp.get("hgvs", []))
    return out


def _resolve_one(gene: str, hgvsc: str, mane_accession: str | None) -> dict:
    hgvs_list = _erepo_hgvs_list(gene, hgvsc)
    nm_entries = [h for h in hgvs_list if re.match(r"^NM_\d+\.\d+:c\.", h)]
    by_accession: dict[str, str] = {}  # base accession (no version) -> full "NM_x.y:c...."
    for entry in nm_entries:
        accession_part = entry.split(":", 1)[0]
        base = accession_part.split(".")[0]
        by_accession.setdefault(base, entry)

    if mane_accession:
        mane_base = mane_accession.split(".")[0]
        exact = f"{mane_accession}:{hgvsc}"
        if exact in nm_entries:
            return {"transcript": mane_accession, "resolution_method": "exact_erepo_match",
                    "mane_accession": mane_accession}
        if mane_base in by_accession:
            other = by_accession[mane_base]
            transcript = other.split(":", 1)[0]
            return {"transcript": transcript, "resolution_method": "erepo_other_version",
                    "mane_accession": mane_accession}
        return {"transcript": mane_accession, "resolution_method": "mane_fallback",
                "mane_accession": mane_accession}
    return {"transcript": None, "resolution_method": "unresolved", "mane_accession": None}


def main() -> None:
    variants = list(unique_variants())
    genes = {gene for gene, _ in variants}
    mane = _fetch_mane_transcripts(genes)

    resolved: dict[str, dict] = {}
    if OUTPUT_PATH.is_file():
        resolved = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))

    for i, (gene, hgvsc) in enumerate(variants, 1):
        key = f"{gene}|{hgvsc}"
        if key in resolved:
            continue
        print(f"[{i}/{len(variants)}] {gene} {hgvsc}...", end=" ")
        try:
            entry = _resolve_one(gene, hgvsc, mane.get(gene))
        except requests.RequestException as exc:
            entry = {"transcript": mane.get(gene), "resolution_method": "erepo_query_failed",
                      "mane_accession": mane.get(gene), "error": str(exc)}
        print(entry["resolution_method"], "->", entry["transcript"])
        resolved[key] = entry
        OUTPUT_PATH.write_text(json.dumps(resolved, indent=2, sort_keys=True), encoding="utf-8")
        time.sleep(0.2)  # polite pacing against ERepo's live API

    counts: dict[str, int] = {}
    for entry in resolved.values():
        counts[entry["resolution_method"]] = counts.get(entry["resolution_method"], 0) + 1
    print(f"\n[resolve_erepo_transcripts] {len(resolved)} variant(s) resolved -> {OUTPUT_PATH}")
    print(f"  by method: {counts}")


if __name__ == "__main__":
    main()
