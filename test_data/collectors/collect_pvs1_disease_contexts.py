"""Rebuild test_data/pvs1_disease_contexts.json from the live ClinGen ERepo API.

PVS1 cannot be evaluated without a disease context, so measuring it needs one per variant.
This records the disease each ground-truth variant was actually interpreted under, taken from
the same ERepo interpretation the expert-panel call comes from.

The disease is an INPUT, not a label. ACMG's own minimal input is a variant plus a condition,
and supplying it is supplying what a clinician would. What must never be supplied is the PVS1
outcome, which lives in test_data/pvs1_erepo_cases.py's registry and is only ever compared
against - see that file's own `purpose` for the same rule.

Not every variant resolves. ERepo publishes a subset of its interpretations through this
endpoint, and a gene's response simply may not contain the variant the ground truth names.
Those are recorded as unresolved rather than filled in from elsewhere: a disease chosen by
this script would be a disease nobody interpreted the variant under.
"""
import json
import re
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from test_data.collectors.full_criteria_ground_truth import entries_for, unique_variants  # noqa: E402

API = "https://erepo.clinicalgenome.org/evrepo/api/classifications"
OUT = ROOT / "test_data" / "fetched_data" / "pvs1_disease_contexts.json"
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
    targets = [(gene, hgvsc) for gene, hgvsc in unique_variants()
               if any(entry.criterion == "PVS1" for entry in entries_for(gene, hgvsc))]
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
        "purpose": ("Disease context for the ground-truth variants carrying a PVS1 call, as "
                    "recorded by the ERepo interpretation itself. This is evaluator INPUT - "
                    "the expert panel's PVS1 outcome is not here and must not be supplied."),
        "source": "ClinGen Evidence Repository API",
        "source_url": API,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "variants_with_a_pvs1_call": len(targets),
        "resolved": len(found),
        "unresolved": [f"{gene}|{hgvsc}" for gene, hgvsc in targets if (gene, hgvsc) not in found],
        "errors": errors,
        "contexts": {f"{gene}|{hgvsc}": value for (gene, hgvsc), value in sorted(found.items())},
    }
    OUT.write_text(json.dumps(document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"\n{len(found)}/{len(targets)} resolved -> {OUT}")


if __name__ == "__main__":
    main()
