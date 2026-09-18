"""Run PVS1 end to end over the ground-truth variants and compare with the expert panel.

Two questions, deliberately kept apart, because they have very different answers and only one
of them is a measurement:

  --contexts erepo     Each variant gets the disease ERepo interpreted it under. This is the
                       measurement: realistic input, and whatever the disease gate then does
                       is what it would do in use.

  --contexts curated   Each variant gets the disease its own gene's curated mechanism names,
                       and the inheritance mode curated with it. This is NOT a measurement -
                       the input is chosen so the gate opens - it checks that the disease gate
                       hands off to the decision tree and that the tree runs to a verdict.

The gap between the two is the point. If `erepo` withholds where `curated` reaches MET, the
disease context is what is missing, not the decision tree.

Neither mode supplies the expert panel's PVS1 outcome to the evaluator. The disease context is
input (ACMG's minimal input is a variant plus a condition); the outcome is only ever compared
against, and lives in test_data/full_criteria_ground_truth.py.

Both modes need network access on a cold cache. Everything else replays from cache/.
"""
import argparse
import json
import sys
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import run_automated_validation_64 as base  # noqa: E402
import acmg_pipeline.automated_cli as cli  # noqa: E402
from acmg_pipeline.automated_core.input import audit_vcf  # noqa: E402
from test_data.collectors.full_criteria_ground_truth import entries_for, unique_variants  # noqa: E402

CONTEXTS_PATH = ROOT / "test_data" / "fetched_data" / "pvs1_disease_contexts.json"
RULES_PATH = ROOT / "config" / "demo-rules.json"
CURATED_CONTEXT_PATH = ROOT / "config" / "curated-context.json"

# --contexts curated: the disease each gene's own curated mechanism names, and the mode it was
# curated under. Read off ClinGen dosage and G2P on 2026-09-18; refresh with
# providers/clingen_dosage.py and providers/gene2phenotype.py if a source moves.
CURATED = {
    "RUNX1|c.601C>T": ("MONDO:0100083", None),            # ClinGen dosage, mode-unscoped
    "RPE65|c.495+1dup": ("MONDO:0008765", "autosomal recessive"),
    "MYBPC3|c.278delA": ("MONDO:0007268", "autosomal dominant"),
    "MYBPC3|c.2905+1G>A": ("MONDO:0007268", "autosomal dominant"),
    "MYBPC3|c.836del": ("MONDO:0007268", "autosomal dominant"),
    "USH2A|c.8559-2A>G": ("MONDO:0010169", "autosomal recessive"),
    "CDH1|c.387+5G>A": ("MONDO:0100488", "autosomal dominant"),
}

PROVIDERS = [
    "--with-clingen-dosage", "--with-gene2phenotype", "--with-mane-transcript",
    "--with-nmd-prediction", "--with-splice-default", "--with-mondo-mapping",
    "--with-clingen-lumping", "--with-mondo-hierarchy", "--with-clingen-gene-validity",
]


def selected(mode):
    """(gene, hgvsc, condition, inheritance) for every variant this mode can evaluate."""
    transcripts = json.loads(base.TRANSCRIPTS_PATH.read_text(encoding="utf-8"))
    rows = []
    if mode == "erepo":
        contexts = json.loads(CONTEXTS_PATH.read_text(encoding="utf-8"))["contexts"]
        lookup = {key: (value["condition"], None) for key, value in contexts.items()}
    else:
        lookup = CURATED
    for gene, hgvsc in unique_variants():
        key = f"{gene}|{hgvsc}"
        if key not in lookup or not any(e.criterion == "PVS1" for e in entries_for(gene, hgvsc)):
            continue
        if not transcripts.get(key, {}).get("transcript"):
            continue
        rows.append((gene, hgvsc, *lookup[key]))
    return rows, transcripts


def write_vcf(rows, transcripts, path):
    lines = ["##fileformat=VCFv4.2", "##reference=GRCh38",
             "##caveat=Synthetic records; identity is resolved from TRANSCRIPT+HGVSC via "
             "Ensembl VEP, the same convention as run_automated_validation_64.py",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for index, (gene, hgvsc, _, _) in enumerate(rows, 1):
        transcript = transcripts[f"{gene}|{hgvsc}"]["transcript"]
        lines.append(f"1\t1\tpvs1-{index}\tA\t.\t80\tPASS\t"
                     f"GENE={gene};TRANSCRIPT={transcript};HGVSC={hgvsc}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_context(rows, prepared, path):
    """The disease context as a curated-context document, which is how the evaluator takes it.

    It does not travel through the VCF: prepared records are built from an allowlist that a
    synthetic INFO key does not reach, so a CONDITION written there is silently dropped.
    """
    records = json.loads((prepared / "variants.json").read_text(encoding="utf-8"))["records"]
    by_variant = {record["source"]["variant_id"]: record["record_id"] for record in records}
    document = {
        "schema_version": "1.0", "context_version": "pvs1-validation",
        "source": "test_data/pvs1_disease_contexts.json, or each gene's own curated mechanism",
        # Reuse the repository's reviewed BA1 exception block; only the disease contexts differ.
        "ba1_exceptions": json.loads(
            CURATED_CONTEXT_PATH.read_text(encoding="utf-8"))["ba1_exceptions"],
        "records": {}, "record_contexts": {},
    }
    for index, (_, _, condition, inheritance) in enumerate(rows, 1):
        record_id = by_variant.get(f"pvs1-{index}")
        if not record_id:
            continue
        entry = {"condition": condition}
        if inheritance:
            entry["inheritance"] = inheritance
        document["record_contexts"][record_id] = entry
    path.write_text(json.dumps(document, ensure_ascii=False, indent=1), encoding="utf-8")


def report(rows, evaluated, mode):
    results = json.loads((evaluated / "results.json").read_text(encoding="utf-8"))
    by_variant = {record["source"]["variant_id"]: record for record in results["records"]}
    tally, applicability = Counter(), Counter()
    print(f"\n{'':3}{'gene':9}{'hgvsc':17}{'panel':22}{'ours':22}"
          f"{'applicability':16}{'match':14}{'stop'}")
    for index, (gene, hgvsc, _, _) in enumerate(rows, 1):
        truth = next(e for e in entries_for(gene, hgvsc) if e.criterion == "PVS1")
        applied = truth.status.value == "met"
        want = f"MET {truth.strength.value}" if applied else "withheld"
        record = by_variant.get(f"pvs1-{index}")
        criterion = next((item for item in (record or {}).get("results", [])
                          if item["criterion"] == "PVS1"), None)
        if criterion is None:
            tally[(applied, "no record")] += 1
            print(f"  x{gene:8}{hgvsc:17}{want:22}{'no record'}")
            continue
        # PVS1 never returns NOT_MET here - withholding is UNKNOWN - so the comparison is
        # applied against withheld, with the strength checked where both applied.
        got_applied = criterion["status"] == "met"
        got = f"MET {criterion.get('strength')}" if got_applied else "withheld"
        context = criterion.get("evaluation_context") or {}
        stop = ([node["node_id"] for node in criterion.get("decision_trace", [])] or ["-"])[-1]
        tally[(applied, got_applied)] += 1
        applicability[context.get("applicability")] += 1
        agree = applied == got_applied and (
            not applied or criterion.get("strength") == truth.strength.value)
        print(f"{'OK ' if agree else '  x'}{gene:8}{hgvsc:17}{want:22}{got:22}"
              f"{str(context.get('applicability')):16}{str(context.get('disease_match')):14}{stop}")
        for point in criterion.get("review_points", []):
            print(f"{'':9}  review: {point[:104]}")
    print(f"\n[{mode}] {len(rows)} variants")
    print(f"  panel applied  & applied : {tally[(True, True)]}")
    print(f"  panel applied  & withheld: {tally[(True, False)] + tally[(True, 'no record')]}")
    print(f"  panel withheld & withheld: {tally[(False, False)]}")
    print(f"  panel withheld & applied : {tally[(False, True)]}   <- false application")
    print(f"  applicability: {dict(applicability)}")
    if mode == "curated":
        print("  NOTE: contexts were chosen to open the gate; this is a connection check, "
              "not an accuracy measurement.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contexts", choices=("erepo", "curated"), default="erepo")
    parser.add_argument("--keep", action="store_true", help="print the working directory")
    args = parser.parse_args()

    rows, transcripts = selected(args.contexts)
    if not rows:
        raise SystemExit(f"no variants for --contexts {args.contexts}")
    work = Path(tempfile.mkdtemp(prefix=f"bh26_pvs1_{args.contexts}_"))
    write_vcf(rows, transcripts, work / "cases.vcf")
    print(f"[{args.contexts}] {len(rows)} variants with a transcript and a disease context")

    cli.audit_demo = lambda _directory: audit_vcf(work / "cases.vcf", case_id="pvs1",
                                                  assembly="GRCh38")
    prepared, evaluated = work / "prepared", work / "evaluated"
    base.CACHE_DIR.mkdir(parents=True, exist_ok=True)
    base.EVIDENCE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if cli.main(["prepare-demo-online", "--input-dir", str(work),
                 "--cache-dir", str(base.CACHE_DIR),
                 "--evidence-cache-dir", str(base.EVIDENCE_CACHE_DIR),
                 "--output-dir", str(prepared), "--ensembl-release", "116",
                 "--rules", str(RULES_PATH), *PROVIDERS]):
        raise SystemExit(1)
    write_context(rows, prepared, work / "context.json")
    if cli.main(["evaluate", "--input", str(prepared / "variants.json"),
                 "--evidence", str(prepared / "evidence.json"),
                 "--config", str(RULES_PATH), "--context", str(work / "context.json"),
                 "--output-dir", str(evaluated), "--internal-only"]):
        raise SystemExit(1)
    report(rows, evaluated, args.contexts)
    if args.keep:
        print(f"\nworking directory: {work}")


if __name__ == "__main__":
    main()
