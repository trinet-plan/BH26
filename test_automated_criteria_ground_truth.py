"""
test_automated_criteria_ground_truth.py

Cross-checks the 16 automated Layer-1 criteria (acmg_pipeline.classification.
AUTOMATED_CODES - merged into main from h.muroda's codex/integrate-evidence-
cli branch, see commits a575059/6d05552/c5b303f) against real human-curator
ground truth, the same way test_full_criteria_ground_truth.py already does
for this project's own 3 literature codes.

[Where the ground truth comes from]
  Each democase/case*_variants_v2.vcf row's own INFO column carries an
  ACMG_CODES=... field - the human curator's actual per-variant assertion
  (e.g. "PVS1_VeryStrong,PM2_Moderate"), written by hand into these files
  well before evidence-cli existed (see NOTE= in the same INFO column for
  the curation history/citations). This is a DIFFERENT ground-truth source
  from test_data/full_criteria_ground_truth.py (ERepo/democase-doc-derived,
  no genomic coordinates) - the automated engine needs real CHROM/POS/REF/
  ALT to query gnomAD/ClinVar/Ensembl/dbNSFP, which only the demo VCFs
  provide, so this script is scoped to the 4 demo cases (28 ALT records)
  rather than the full 373-entry ERepo dataset.

[Why this uses the CLI's --internal-only path]
  This script only needs the internal per-criterion results, and
  --internal-only returns exactly those without also running VA-Spec
  export, which this comparison does not read.

  This paragraph used to say the export raised AttributeError:
  'VariantPathogenicityEvidenceLine' object has no attribute 'schema_id'.
  That is no longer true - as of 2026-09-18 the installed ga4gh.va_spec
  model defines schema_id(), tests/test_va_spec.py covers it, and
  tests/test_demo_pipeline.py runs the full export and asserts
  va_spec_export == "VALIDATED". The note is kept rather than deleted
  because it was quoted as evidence of an open bug months after it had
  been fixed.

[Why this only checks MET assertions, not NOT_MET]
  ACMG_CODES= lists only the codes the curator found MET for that variant -
  it is not an exhaustive "every code was considered and these are the only
  MET ones" statement (e.g. a curator who didn't bother checking BP1 for a
  clearly-pathogenic frameshift wouldn't write "BP1 not met"). So this
  script only asserts "if the curator wrote CODE here, does the automated
  engine also say MET for CODE" - it does NOT treat a code's absence from
  ACMG_CODES as a ground-truth NOT_MET, and does not fail on the engine
  finding MET for something the curator didn't mention (logged separately,
  informationally, as "engine-only MET" - could be curator omission or a
  genuine engine false positive, this script cannot tell which).

[A known, expected gap: PVS1]
  PVS1 needs gene-level loss-of-function mechanism evidence
  (acmg_pipeline.criteria.pvs1._resolve_mechanism(), via a "gene_disease"
  services provider) that neither this script's CLI invocation nor
  evidence-cli's own tests/test_demo_pipeline.py wires in (that needs the
  separate build-gene-disease-drafts/build-gene-disease-evidence CLI
  steps, which operate on a totally different config, not exercised here).
  So PVS1 always comes back "unknown" here even where the curator asserted
  PVS1_VeryStrong - a real, known limitation of this offline invocation,
  not a defect in this script or a real disagreement with the engine's
  actual PVS1 judgment.
"""

from __future__ import annotations

import json
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.automated_cli import main as automated_main
from acmg_pipeline.classification import AUTOMATED_CODES
from test_harness import Harness

ROOT = Path(__file__).resolve().parent
DEMO_VCFS = sorted(ROOT.glob("democase/case*_variants_v2.vcf"))

h = Harness()
check = h.check


# ============================================================================
# [1] Parse each demo VCF's own ACMG_CODES= INFO field as ground truth
# ============================================================================

def _norm_strength(raw: str) -> str | None:
    """'VeryStrong' -> 'very_strong'; 'present' or '' -> None (no assertion)."""
    raw = re.sub(r"\([^)]*\)", "", raw)  # strip "(low_confidence_mouse_model)" etc.
    if raw in ("", "present"):
        return None
    return re.sub(r"(?<!^)(?=[A-Z])", "_", raw).lower()


def _parse_acmg_codes(field: str) -> dict[str, str | None]:
    """'PVS1_VeryStrong,PM2_Moderate' -> {'PVS1': 'very_strong', 'PM2': 'moderate'}."""
    out: dict[str, str | None] = {}
    if not field:
        return out
    for token in field.split(","):
        m = re.match(r"^([A-Z]{2}\d)(?:_(.+))?$", token)
        if not m:
            continue
        code, strength_raw = m.group(1), m.group(2) or ""
        out[code] = _norm_strength(strength_raw)
    return out


def _load_ground_truth() -> dict[str, dict[str, str | None]]:
    """variant_id (e.g. 'case1-var1') -> {code: strength_or_None}."""
    ground_truth: dict[str, dict[str, str | None]] = {}
    for vcf_path in DEMO_VCFS:
        for line in vcf_path.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            variant_id = cols[2]
            info = cols[7]
            m = re.search(r"ACMG_CODES=([^;]*)", info)
            if m:
                ground_truth[variant_id] = _parse_acmg_codes(m.group(1))
    return ground_truth


print("[1] Ground truth extracted from democase VCFs' own ACMG_CODES= field")
ground_truth = _load_ground_truth()
n_asserted = sum(len(codes) for codes in ground_truth.values())
check(f"at least 10 curator MET assertions collected across the 4 demo cases (got {n_asserted})",
      n_asserted >= 10)
check("every parsed code is a recognized AUTOMATED_CODES member or a literature code (sanity check)",
      all(code in AUTOMATED_CODES or code in {"PS2", "PS3", "PS4", "PM1", "PM5", "PP1", "PP3", "BS1"}
          for codes in ground_truth.values() for code in codes))


# ============================================================================
# [2] Run the automated engine offline (same fixtures as tests/test_demo_
#     pipeline.py), --internal-only to skip the currently-broken VA-Spec
#     export step (see module docstring)
# ============================================================================
print("\n[2] Running the automated engine over all 28 democase ALT records (offline, fixture-cached)")

_work = Path(tempfile.mkdtemp(prefix="bh26_automated_gt_check_"))
_prepared = _work / "prepared"
_evaluated = _work / "evaluated"

_prepare_exit = automated_main([
    "prepare-demo-online", "--input-dir", str(ROOT / "democase"),
    "--cache-dir", str(ROOT / "tests" / "fixtures" / "ensembl-cache"),
    "--evidence-cache-dir", str(ROOT / "tests" / "fixtures" / "external-cache"),
    "--output-dir", str(_prepared), "--ensembl-release", "116",
    "--with-gnomad", "--gnomad-release", "4.1.1",
    "--with-clinvar", "--clinvar-release", "2026-09-15",
    "--with-pm1-hotspot", "--rules", str(ROOT / "config" / "demo-rules.json"),
    "--with-dbnsfp",
    "--offline",
])
check("prepare-demo-online exits 0 (offline, fully fixture-cached)", _prepare_exit == 0)

_evaluate_exit = automated_main([
    "evaluate", "--input", str(_prepared / "variants.json"),
    "--evidence", str(_prepared / "evidence.json"), "--offline",
    "--config", str(ROOT / "config" / "demo-rules.json"),
    "--context", str(ROOT / "config" / "curated-context.json"),
    "--output-dir", str(_evaluated), "--internal-only",
])
check("evaluate --internal-only exits 0", _evaluate_exit == 0)

results = json.loads((_evaluated / "results.json").read_text(encoding="utf-8"))
records = results["records"]
check("all 28 democase ALT records evaluated", len(records) == 28)


# ============================================================================
# [3] Compare: for every curator-asserted MET code, does the engine agree?
#     (Report only - see module docstring for why this isn't a hard gate.)
# ============================================================================
print("\n[3] Engine vs. curator ACMG_CODES= ground truth (report only, not a pass/fail gate)\n")

match_n = 0
compared_n = 0
mismatches = []
engine_only_met = []

for record in records:
    variant_id = record["source"]["variant_id"]
    truth = ground_truth.get(variant_id, {})
    if not truth:
        continue
    by_code = {r["criterion"]: r for r in record["results"]}

    for code, expected_strength in truth.items():
        if code not in AUTOMATED_CODES:
            continue  # literature-side code (e.g. PS3/PS4) - not this script's concern
        compared_n += 1
        result = by_code.get(code)
        engine_status = result["status"] if result else "unknown"
        engine_strength = result.get("strength") if result else None
        engine_summary = result.get("summary") if result else None
        is_match = engine_status == "met"
        match_n += is_match
        label = f"{variant_id} {code} (curator: MET{'/' + expected_strength if expected_strength else ''})"
        if is_match:
            print(f"  MATCH    {label} -> engine: met"
                  f"{'/' + engine_strength if engine_strength else ''}")
        else:
            reason = f" [{engine_summary}]" if engine_summary else ""
            print(f"  MISMATCH {label} -> engine: {engine_status}{reason}")
            mismatches.append(f"{label} -> engine: {engine_status}{reason}")

    for code, result in by_code.items():
        if code in AUTOMATED_CODES and result["status"] == "met" and code not in truth:
            engine_only_met.append(f"{variant_id} {code} (curator ACMG_CODES did not mention it)")

print(f"\n{compared_n} curator-asserted MET call(s) checked across {len([v for v in ground_truth.values() if v])} "
      f"variant(s) with ground truth; {match_n}/{compared_n} confirmed MET by the automated engine.")
if mismatches:
    print(f"\n{len(mismatches)} mismatch(es) (engine did not confirm a curator-asserted MET call):")
    for m in mismatches:
        print(f"  - {m}")
if engine_only_met:
    print(f"\n{len(engine_only_met)} engine-only MET call(s) (engine said MET, curator's ACMG_CODES didn't "
          f"list this code for this variant - may be curator omission, may be a real engine false positive; "
          f"not distinguishable from this script alone):")
    for m in engine_only_met:
        print(f"  - {m}")


h.report_and_exit()
