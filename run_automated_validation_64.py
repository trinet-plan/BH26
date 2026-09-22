"""
run_automated_validation_64.py

Validates the 16 automated Layer-1 criteria (acmg_pipeline.classification.
AUTOMATED_CODES) against ALL 64 unique variants in test_data.
full_criteria_ground_truth.py (ERepo + democase sourced) - not just the 4
demo-data cases test_automated_criteria_ground_truth.py already covers.
Companion to run_validation_64.py (the equivalent full-dataset validation
for this project's own 3 literature codes).

[Why this needs its own synthetic VCF, not demo-data's]
  full_criteria_ground_truth.py only carries (gene, hgvsc) - no RefSeq
  transcript accession, no genomic coordinates. See test_data/
  resolve_erepo_transcripts.py's own docstring for how each variant's
  transcript was resolved (54 exact ERepo hgvs-list matches, 4 same-
  accession-different-version, 6 demo-data-derived MANE Select fallbacks -
  run that script first, or use its committed output, test_data/
  erepo_variant_transcripts.json).

  Every synthetic VCF row uses ALT="." (CHROM/POS/REF are harmless
  placeholders - "1"/1/"A") - the exact same "identity is unknown from the
  source row alone, resolve it via HGVS" convention demo-data/case1_
  variants_v2.vcf already uses for case1-var1/case2-var2 (see
  acmg_pipeline/automated_core/input.py: an unparseable ALT makes
  parsed_variant None, and reconcile() then trusts the corroborated
  Ensembl VEP HGVS candidate instead - see identity.py's own docstring,
  "Candidates must come from an identity adapter, not from source
  classification labels").

[Why audit_demo() is monkeypatched, not modified]
  acmg_pipeline.automated_cli's prepare-demo-online/evaluate command
  handlers (200+ lines: Ensembl identity resolution, gnomAD/ClinVar/
  dbNSFP/PM1-hotspot wiring, evidence/manifest assembly) are exactly what
  this script needs to reuse unchanged - but they all start from
  `records = audit_demo(args.input_dir)`, which is hardcoded to
  demo-data/case{1..4}_variants_v2.vcf specifically (acmg_pipeline.
  automated_core.input.audit_demo loops `for case in range(1, 5)`, not a
  generic "read this directory" function). Rather than duplicate the
  200+ lines of provider-wiring logic that follows, this script
  monkeypatches acmg_pipeline.automated_cli's own `audit_demo` name (only
  in this process, only for the duration of this script) to call the
  generic acmg_pipeline.automated_core.input.audit_vcf() against this
  script's own synthetic VCF instead - every line of CLI logic downstream
  of that call runs completely unmodified.

[Live network, cached for reproducibility]
  These 64 variants have no committed fixture cache (unlike demo-data's
  tests/fixtures/ensembl-cache, external-cache). The FIRST run of this
  script makes real gnomAD/ClinVar/Ensembl/dbNSFP network calls (via the
  same CachedHttpClient every other provider path uses) and writes them
  to cache/erepo_automated_{ensembl,evidence}/ (gitignored, same
  cache/ convention as fulltext_cache.py's PubMed cache) - every
  subsequent run of this exact 64-variant set replays from that cache with
  zero network calls, exactly as verified for the demo-data comparison
  (see test_automated_criteria_ground_truth.py's commit history / the
  online-vs-offline check discussed with the user 2026-09-17).

[Comparison scope and method]
  Only criteria in AUTOMATED_CODES are compared (this dataset also has
  ground truth for PS3/BS3/PS4/PP1/BS4/PP4 and clinical-record-only codes,
  none of which this automated engine implements). Unlike the demo-data
  comparison (which only had curator MET assertions to check against),
  full_criteria_ground_truth.py has explicit NOT_MET entries too, so this
  script compares both directions: does the engine's met/not_met status
  agree with the ground truth's status, for every (gene, hgvsc, code) pair
  where both exist. Report only - not a pass/fail gate, matching every
  other ground-truth comparison in this project (see test_full_criteria_
  ground_truth.py and test_automated_criteria_ground_truth.py's own
  docstrings for why).
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import acmg_pipeline.automated_cli as automated_cli_module
from acmg_pipeline.automated_core.input import audit_vcf
from acmg_pipeline.classification import AUTOMATED_CODES
from acmg_pipeline.providers.clingen_gene_validity import ClinGenGeneValidityProvider
from acmg_pipeline.providers.http import CachedHttpClient
from test_data.full_criteria_ground_truth import entries_for, unique_variants

ROOT = Path(__file__).resolve().parent
TRANSCRIPTS_PATH = ROOT / "test_data" / "erepo_variant_transcripts.json"
DISEASE_CONTEXTS_PATH = ROOT / "test_data" / "erepo_disease_contexts.json"
CURATED_CONTEXT_PATH = ROOT / "config" / "curated-context.json"
CACHE_DIR = ROOT / "cache" / "erepo_automated_ensembl"
EVIDENCE_CACHE_DIR = ROOT / "cache" / "erepo_automated_evidence"


def _gene_level_fallback_conditions(genes: set[str]) -> dict[str, str]:
    """{gene: MONDO condition} for genes with exactly one Definitive/Strong
    ClinGen Gene-Disease Validity row - a fallback for when ERepo has no
    variant-specific interpretation to name a condition at all.

    Confirmed empirically (2026-09-19): 71/87 of the automated 679-variant
    run's remaining PVS1 misses were "Condition was not provided" - ERepo's
    own per-variant interpretation (test_data/erepo_disease_contexts.json)
    only resolves 266/679 variants, and the other 413 are disproportionately
    in genes (PTEN, MYBPC3, CDH1, ...) where LOF variants and a PVS1 call
    are common. Checked before implementing: 68/99 of this dataset's genes
    have exactly one Definitive/Strong row, which would newly resolve a
    condition for ~331 currently-unresolved variants.

    Deliberately conservative: a gene with zero or with MULTIPLE distinct
    Definitive/Strong conditions (different modes of inheritance, or
    genuinely different diseases) gets no fallback - picking one of several
    plausible diseases would be a guess, not a suggestion grounded in a
    single clear answer. This mirrors ClinGenGeneValidityProvider.
    get_validity()'s own use_restriction, "CANDIDATE_CONDITION_SUGGESTION_
    ONLY": a real per-variant interpretation (when ERepo has one) is always
    preferred over this gene-level guess - see the caller.
    """
    client = CachedHttpClient(str(EVIDENCE_CACHE_DIR), offline=True)
    provider = ClinGenGeneValidityProvider(client)
    fallback = {}
    for gene in sorted(genes):
        try:
            rows = provider.get_validity(gene)
        except Exception:  # noqa: BLE001 - no cached snapshot yet; just skip
            continue
        strong = {row["condition"] for row in rows if row.get("classification") in ("Definitive", "Strong")}
        if len(strong) == 1:
            fallback[gene] = next(iter(strong))
    return fallback


def _write_context(variants: list[tuple[str, str]], prepared: Path, out_path: Path) -> int:
    """Build the --context document `evaluate` reads: config/curated-
    context.json's ba1_exceptions (unchanged) plus a per-variant
    record_contexts entry, preferring test_data/erepo_disease_contexts.json
    (ERepo's own interpreted condition for each variant, generated by
    test_data/collect_all_disease_contexts.py) and falling back to
    _gene_level_fallback_conditions() when ERepo has none - the same shape
    run_pvs1_validation.py's write_context() already uses for its own,
    PVS1-only variant subset.

    Without at least the ERepo path, curated-context.json's record_contexts
    (keyed by the 4 demo-data cases' own record ids, e.g. "case1:25:1") has
    no entry for any of this dataset's erepo64-N records, so PVS1's
    mechanism gate and every gene_disease-association-gated code (BP1, PP2)
    get no condition to evaluate against at all, for every variant,
    regardless of which provider flags are enabled - confirmed empirically
    2026-09-19 (see this script's own git history / commit message for the
    full before/after).

    Returns how many of `variants` got a resolved condition.
    """
    contexts = json.loads(DISEASE_CONTEXTS_PATH.read_text(encoding="utf-8"))["contexts"]
    records = json.loads((prepared / "variants.json").read_text(encoding="utf-8"))["records"]
    by_variant_id = {record["source"]["variant_id"]: record["record_id"] for record in records}
    gene_fallback = _gene_level_fallback_conditions({gene for gene, _ in variants})
    # Keyed by gene symbol - see curated-context.json's own gene_frequency_thresholds entries
    # and criteria/ba1.py's/bs1.py's own docstrings for what "ba1"/"bs1" feed. automated_core.
    # context.apply_context() already reads this same shape gene-keyed for the live pipeline
    # path (pipeline.py's _apply_curated_context()); this script's own evaluate step goes
    # through automated_output.run_internal() -> the SAME load_context()/apply_context(), but
    # record_contexts is keyed by record_id, not gene, so the gene-keyed table is resolved to
    # each matching record_id here instead of relying on apply_context()'s own gene lookup
    # (this script's prepared records carry no top-level "gene" field for it to read).
    curated = json.loads(CURATED_CONTEXT_PATH.read_text(encoding="utf-8"))
    gene_thresholds = curated.get("gene_frequency_thresholds") or {}
    # PM1's own gene-keyed table - unlike gene_frequency_thresholds (BA1/BS1), this one is
    # not gated on a resolved disease condition at all below: PM1 evaluate_region()'s curated
    # critical-domain check (criteria/regions.py) is disease-agnostic by design (the same
    # `disease_required=False` the pre-existing hotspot route already uses), so a variant
    # this script could not resolve ANY condition for (the `continue` a few lines down)
    # would otherwise never see its gene's own critical-domain data either, despite PM1
    # never having needed a condition for it.
    gene_domains = curated.get("gene_critical_domains") or {}

    document = {
        "schema_version": "1.0",
        "context_version": "erepo64-automated-validation",
        "source": f"{CURATED_CONTEXT_PATH.name}'s ba1_exceptions + {DISEASE_CONTEXTS_PATH.name}'s "
                  "per-variant ERepo-interpreted conditions, falling back to a gene-level "
                  "ClinGen Gene-Disease Validity condition where ERepo has none",
        "ba1_exceptions": curated["ba1_exceptions"],
        "records": {}, "record_contexts": {},
    }
    resolved = 0
    for index, (gene, hgvsc) in enumerate(variants, 1):
        record_id = by_variant_id.get(f"erepo64-{index}")
        if not record_id:
            continue
        context = contexts.get(f"{gene}|{hgvsc}")
        condition = context["condition"] if context and context.get("condition") else gene_fallback.get(gene)
        record_context = {}
        if gene in gene_domains:
            record_context["gene_critical_domains"] = gene_domains[gene]
        if not condition:
            if record_context:
                document["record_contexts"][record_id] = record_context
            continue
        record_context["condition"] = condition
        gene_entry = gene_thresholds.get(gene)
        if gene_entry:
            if gene_entry.get("bs1"):
                record_context["disease_frequency_threshold"] = gene_entry["bs1"]
            if gene_entry.get("ba1"):
                record_context["ba1_threshold_override"] = gene_entry["ba1"]
        document["record_contexts"][record_id] = record_context
        resolved += 1
    out_path.write_text(json.dumps(document, ensure_ascii=False, indent=2), encoding="utf-8")
    return resolved


def _build_synthetic_vcf(variants: list[tuple[str, str]], transcripts: dict[str, dict],
                          out_path: Path) -> int:
    lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##caveat=Synthetic records for test_data/full_criteria_ground_truth.py's 64 "
        "ERepo/democase variants - CHROM/POS/REF are placeholders, ALT=\".\" deliberately "
        "so identity is resolved purely from TRANSCRIPT+HGVSC via Ensembl VEP, same "
        "convention as demo-data/case1_variants_v2.vcf's case1-var1/case2-var2.",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    included = 0
    for i, (gene, hgvsc) in enumerate(variants, 1):
        entry = transcripts.get(f"{gene}|{hgvsc}")
        if not entry or not entry.get("transcript"):
            continue
        info = f"GENE={gene};TRANSCRIPT={entry['transcript']};HGVSC={hgvsc}"
        lines.append(f"1\t1\terepo64-{i}\tA\t.\t80\tPASS\t{info}")
        included += 1
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return included


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="only process the first N variants")
    parser.add_argument("--offline", action="store_true",
                        help="require the cache to already have everything (no network fallback)")
    args = parser.parse_args()

    transcripts = json.loads(TRANSCRIPTS_PATH.read_text(encoding="utf-8"))
    variants = list(unique_variants())
    if args.limit:
        variants = variants[: args.limit]

    work = Path(tempfile.mkdtemp(prefix="bh26_erepo64_automated_"))
    print(f"[run_automated_validation_64] work dir: {work}")
    vcf_path = work / "erepo64_variants.vcf"
    n_included = _build_synthetic_vcf(variants, transcripts, vcf_path)
    print(f"[run_automated_validation_64] {n_included}/{len(variants)} variant(s) have a resolved "
          f"transcript and are included in this run")

    def _audit_demo_override(_directory):
        return audit_vcf(vcf_path, case_id="erepo64", assembly="GRCh38")

    automated_cli_module.audit_demo = _audit_demo_override

    prepared = work / "prepared"
    prepared_pass1 = work / "prepared_pass1"
    evaluated = work / "evaluated"
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    _COMMON_PREPARE_ARGS = [
        "--input-dir", str(work),  # ignored by the override above
        "--cache-dir", str(CACHE_DIR),
        "--evidence-cache-dir", str(EVIDENCE_CACHE_DIR),
        "--ensembl-release", "116",
        # --with-togovar, not --with-gnomad (mutually exclusive - argparse
        # rejects both): this project moved population_sources from gnomAD
        # to TogoVar (see acmg_pipeline.automated_cli._population_sources()'s
        # own docstring), and BS1/PM2 read population_sources, not a
        # standalone gnomAD record - --with-gnomad alone left BS1 unable to
        # find "a reliable population observation" for any variant.
        "--with-togovar",
        "--with-clinvar", "--clinvar-release", "2026-09-15",
        "--with-pm1-hotspot", "--rules", str(ROOT / "config" / "demo-rules.json"),
        "--with-dbnsfp",
        # The rest were missing entirely (2026-09-19 bug found via the
        # 679-variant ground-truth run): without them, PVS1/BA1/BS1/BP1/
        # BP3/BP7/PM4/PP2 had no disease-mechanism, gene-validity, domain,
        # or NMD context to evaluate against at all, so every one of those
        # 8 codes came back UNKNOWN for every variant with zero exceptions
        # (0/65-0/76 true positives, but also 0 false positives - the
        # engine never guessed, it just never got the inputs it needed).
        # See run_pvs1_validation.py's own PROVIDERS list, which already
        # used most of these correctly for its narrower PVS1-only check.
        "--with-clingen-dosage",
        "--with-gene2phenotype",
        "--with-mane-transcript",
        "--with-nmd-prediction",
        "--with-splice-default",
        "--with-mondo-mapping",
        "--with-clingen-lumping",
        "--with-mondo-hierarchy",
        "--with-clingen-gene-validity",
        "--with-region-assessment",
        "--with-bp7-splice-assessment",
    ]
    if args.offline:
        _COMMON_PREPARE_ARGS.append("--offline")

    # Pass 1: resolve identities only, so context.json can be built (it needs
    # prepared/variants.json's variant_id -> record_id mapping - see
    # _write_context()). --with-gene-disease-draft is deliberately left off
    # this pass: identity['CONDITION'] is never set yet, so the suggestion
    # would always come back INSUFFICIENT anyway (see --context's own help
    # text in automated_cli.py) - no point paying for it twice.
    exit_code = automated_cli_module.main([
        "prepare-demo-online", "--output-dir", str(prepared_pass1), *_COMMON_PREPARE_ARGS,
    ])
    if exit_code:
        print(f"[run_automated_validation_64] prepare-demo-online (pass 1) exited {exit_code}; stopping")
        sys.exit(exit_code)

    context_path = work / "context.json"
    n_with_condition = _write_context(variants, prepared_pass1, context_path)
    print(f"[run_automated_validation_64] {n_with_condition}/{len(variants)} variant(s) "
          f"have an ERepo-interpreted disease condition to evaluate against")

    # Pass 2: re-run with --context now available, so identity['CONDITION'] is
    # set before --with-gene-disease-draft's evidence is generated - the
    # per-provider caches from pass 1 make this replay almost everything,
    # paying real cost only for what actually depends on the condition.
    exit_code = automated_cli_module.main([
        "prepare-demo-online", "--output-dir", str(prepared), "--context", str(context_path),
        "--with-clinvar-spectrum", "--with-gene-disease-draft", *_COMMON_PREPARE_ARGS,
    ])
    if exit_code:
        print(f"[run_automated_validation_64] prepare-demo-online (pass 2) exited {exit_code}; stopping")
        sys.exit(exit_code)

    evaluate_args = [
        "evaluate", "--input", str(prepared / "variants.json"),
        "--evidence", str(prepared / "evidence.json"),
        "--config", str(ROOT / "config" / "demo-rules.json"),
        "--context", str(context_path),
        "--output-dir", str(evaluated), "--internal-only",
    ]
    if args.offline:
        evaluate_args.append("--offline")
    exit_code = automated_cli_module.main(evaluate_args)
    if exit_code:
        print(f"[run_automated_validation_64] evaluate exited {exit_code}; stopping")
        sys.exit(exit_code)

    results = json.loads((evaluated / "results.json").read_text(encoding="utf-8"))
    by_variant_id: dict[str, dict] = {r["source"]["variant_id"]: r for r in results["records"]}

    match_n = 0
    compared_n = 0
    mismatches = []
    unresolved_variants = []

    print(f"\n{'='*70}\n[Engine vs. full_criteria_ground_truth.py, AUTOMATED_CODES only]\n{'='*70}")
    for i, (gene, hgvsc) in enumerate(variants, 1):
        entry = transcripts.get(f"{gene}|{hgvsc}")
        variant_id = f"erepo64-{i}"
        if not entry or not entry.get("transcript"):
            unresolved_variants.append(f"{gene} {hgvsc}")
            continue
        record = by_variant_id.get(variant_id)
        if record is None:
            continue
        by_code = {r["criterion"]: r for r in record["results"]}
        gt_entries = [e for e in entries_for(gene, hgvsc) if e.criterion in AUTOMATED_CODES]
        for gt in gt_entries:
            compared_n += 1
            result = by_code.get(gt.criterion)
            engine_status = result["status"] if result else "unknown"
            engine_summary = result.get("summary") if result else None
            expected = "met" if gt.status.value == "met" else "not_met"
            is_match = engine_status == expected
            match_n += is_match
            gt_strength = f"/{gt.strength.value}" if gt.strength else ""
            label = f"{gene} {hgvsc} {gt.criterion} (ground truth: {expected}{gt_strength}, {entry['resolution_method']})"
            if is_match:
                print(f"  MATCH    {label} -> engine: {engine_status}")
            else:
                reason = f" [{engine_summary}]" if engine_summary else ""
                line = f"{label} -> engine: {engine_status}{reason}"
                print(f"  MISMATCH {line}")
                mismatches.append(line)

    print(f"\n{compared_n} (variant, AUTOMATED_CODES criterion) ground-truth pair(s) checked; "
          f"{match_n}/{compared_n} matched.")
    if mismatches:
        print(f"\n{len(mismatches)} mismatch(es):")
        for m in mismatches:
            print(f"  - {m}")
    if unresolved_variants:
        print(f"\n{len(unresolved_variants)} variant(s) skipped (no resolved transcript): "
              f"{unresolved_variants}")


if __name__ == "__main__":
    main()
