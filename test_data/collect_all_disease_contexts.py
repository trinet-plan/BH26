"""Build test_data/erepo_disease_contexts.json from the live ClinGen ERepo API,
for EVERY variant in full_criteria_ground_truth.py (679 as of 2026-09-19), not
just the ~31 carrying a PVS1 call - see test_data/collect_pvs1_disease_
contexts.py's own docstring for why the PVS1-only variant needed one at all,
and reuse of the same reasoning here.

[Why this exists as a separate script, not a generalization of the PVS1 one]
  collect_pvs1_disease_contexts.py's own docstring/purpose text and output
  schema are scoped to "the ground-truth variants carrying a PVS1 call" -
  run_pvs1_validation.py still reads that exact file. Widening its variant
  filter in place would silently change what that file means without
  changing its name or docstring. This script is the same fetch/normalize
  logic, applied to every variant, writing a separate, similarly-named
  output file that run_automated_validation_64.py's own context-building
  step (not run_pvs1_validation.py) reads instead.

[Why every automated code needs this, not just PVS1]
  Confirmed empirically (2026-09-19, 679-variant automated validation run):
  BP1/PP2 (gene_disease-association-gated) and PVS1 (mechanism-gated) all
  came back UNKNOWN for every one of 679 variants with zero exceptions -
  not because their provider flags were missing (a first fix attempt that
  measurably helped elsewhere), but because run_automated_validation_64.py
  was passing config/curated-context.json unmodified, whose record_contexts
  are keyed by the 4 demo-data cases' own record ids (case1:25:1, ...) and
  contain nothing for this dataset's erepo64-N records - every disease
  lookup had no condition to look up at all.
"""
import json
import re
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from test_data.full_criteria_ground_truth import unique_variants  # noqa: E402

API = "https://erepo.clinicalgenome.org/evrepo/api/classifications"
OUT = ROOT / "test_data" / "erepo_disease_contexts.json"
USER_AGENT = "BH26-ACMG/0.1"
_MONDO = re.compile(r"(MONDO[:_]\d+)")


def normalize(hgvsc):
    """"c.278delA" and "c.278del" are the same variant written two ways."""
    return re.sub(r"(del|dup|ins)[ACGT]+$", r"\1", hgvsc.strip())


def fetch(gene):
    url = f"{API}?{urllib.parse.urlencode({'gene': gene})}"
    request = urllib.request.Request(
        url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def main():
    targets = list(unique_variants())
    wanted = {(gene, normalize(hgvsc)): (gene, hgvsc) for gene, hgvsc in targets}
    found, errors = {}, []
    for gene in sorted({gene for gene, _ in targets}):
        try:
            document = fetch(gene)
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            errors.append(f"{gene}: {type(exc).__name__}: {exc}")
            continue
        for interpretation in document.get("variantInterpretations", []):
            match = _MONDO.search(interpretation.get("@id") or "")
            if not match:
                continue
            condition = match.group(1).replace("_", ":")
            label = (interpretation.get("condition") or {}).get("label")
            for hgvs in interpretation.get("hgvs", []):
                if ":c." not in hgvs:
                    continue
                key = (gene, normalize("c." + hgvs.split(":c.")[-1]))
                original = wanted.get(key)
                if original and original not in found:
                    found[original] = {"condition": condition, "condition_label": label,
                                       "erepo_url": interpretation.get("@id")}
        print(f"{gene}: {sum(1 for key in found if key[0] == gene)} resolved", flush=True)

    document = {
        "schema_version": "1.0",
        "purpose": ("Disease context for EVERY full_criteria_ground_truth.py variant, as "
                    "recorded by the ERepo interpretation itself. This is evaluator INPUT - "
                    "no ground-truth outcome for any criterion is here and must not be supplied."),
        "source": "ClinGen Evidence Repository API",
        "source_url": API,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "total_variants": len(targets),
        "resolved": len(found),
        "unresolved": [f"{gene}|{hgvsc}" for gene, hgvsc in targets if (gene, hgvsc) not in found],
        "errors": errors,
        "contexts": {f"{gene}|{hgvsc}": value for (gene, hgvsc), value in sorted(found.items())},
    }
    OUT.write_text(json.dumps(document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"\n{len(found)}/{len(targets)} resolved -> {OUT}")


if __name__ == "__main__":
    main()
