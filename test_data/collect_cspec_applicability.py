"""Rebuild config/cspec_applicability_draft.json from the live ClinGen CSpec registry.

[Why this exists]
  PVS1 stops at G01 unless something states that loss of function is an established disease
  mechanism. The three routes wired up today all go silent for a recessive, non-cardiomyopathy
  gene: ClinGen Dosage emits nothing for a score of 30 (see providers/clingen_dosage.py's own
  docstring for why that is deliberate), Gene2Phenotype returns an empty records_summary for a
  gene it has not curated, and config/gene-disease-review-decisions.json covers five genes.
  DYSF c.3498_3499delinsAA - Pathogenic with PVS1 at the LGMD VCEP - hits all three and comes
  out Uncertain Significance.

  A VCEP that wrote a criteria specification has already answered this, and CSpec publishes
  those specifications as machine-readable JSON. This collects them.

[What is read]
  /cspec/api/svis lists every specification with its GN id. Each document's ruleSets carry
  the genes the specification is scoped to - with their MONDO disease and mode of inheritance -
  and the criteriaCodes whose evidenceStrengths state applicability.

  A criterion is recorded "Applicable" when at least one of its evidence strengths is, and
  "Not applicable" when every strength says it is not. The registry's other rendering of the
  same document (/cspec/SequenceVariantInterpretation/id/<numeric>) states one `applicability`
  per criterion directly, and the two agree on all 28 of GN180 (DYSF).

  That comparison alone is NOT enough to trust the fold, and a first version of this file got
  it wrong by relying on it. GN180 writes every affirmative as the bare word "Applicable", so
  matching that word exactly passed there while silently turning GN010's "Applicable with VCEP
  specification" into "Not applicable" - which reported the Lysosomal Diseases VCEP as having
  struck out all 28 criteria for GAA, including the PVS1 it actually applies. The six wordings
  the registry really uses are listed at APPLICABLE_PREFIX below, counted across all 208
  documents rather than read off one.

[What this is NOT]
  A VCEP marking PVS1 "Applicable" is a statement about its own specification, not the
  sentence `lof_mechanism_established` makes. They lined up on all three genes where the
  hand-written review and CSpec both have a disease-scoped answer (MYBPC3, MYH7, TNNI3), but
  that is evidence for the mapping, not the mapping itself - which is why the output is a
  DRAFT nothing reads at runtime, and why it is not written as `lof_mechanism_established`.

  PP2/BP1 do NOT carry over the same way. `pp2_applicable` in the hand-written review means
  the criterion is in scope and leaves met/not-met to the mechanism fields beside it, while
  CSpec's "Not applicable" means the VCEP struck the criterion out. Mapping one onto the other
  would change what the field says. Only the output's own vocabulary is recorded here.

  Nothing here is a variant-level judgment, and no ERepo classification is read.

[Cost]
  209 requests, ~15MB, all cached under --cache-dir so a re-run is free. The run took well
  under a minute against the live registry with zero failures (2026-09-18).

Usage:
  python -m test_data.collect_cspec_applicability
  python -m test_data.collect_cspec_applicability --cache-dir cache/cspec --report
"""
from __future__ import annotations

import argparse
import json
import re
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = ROOT / "config" / "cspec_applicability_draft.json"
DEFAULT_CACHE = ROOT / "cache" / "cspec"

SVIS_URL = "https://cspec.genome.network/cspec/api/svis"
DOC_URL = "https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/{}"
HEADERS = {"Accept": "application/json", "User-Agent": "BH26-cspec-applicability-collector/1.0"}

# The registry states applicability in six wordings, and a VCEP that specified or adopted a
# criterion says so in the value rather than in a separate field:
#   Applicable / Applicable with VCEP specification / Applicable as originally described
#   Not applicable / Not Applicable / Not Applicable for this VCEP
# So a strength is applicable when its value BEGINS with "applicable" once case-folded -
# exact equality against "applicable" silently drops the two qualified affirmatives, and
# "Not Applicable for this VCEP" must not be caught by a substring test.
APPLICABLE_PREFIX = "applicable"
ALL_CODES = (
    "PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
)


def fetch(url: str, cache_dir: Path, name: str) -> dict:
    path = cache_dir / f"{name}.json"
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    body = None
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers=HEADERS), timeout=60) as response:
                body = response.read().decode("utf-8")
            break
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
    cache_dir.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    return json.loads(body)


def applicability(code: dict) -> str | None:
    """Applicable if any evidence strength is; Not applicable if every one says not."""
    values = [str(entry.get("applicability") or "") for entry in code.get("evidenceStrengths") or []]
    values = [value for value in values if value]
    if not values:
        return None
    return ("Applicable"
            if any(value.casefold().startswith(APPLICABLE_PREFIX) for value in values)
            else "Not applicable")


def gene_entries(rule_set: dict):
    """(symbol, [MONDO ids], [modes of inheritance]) for every gene a rule set scopes.

    A gene whose symbol cannot be read is skipped rather than guessed: the symbol is the key
    every consumer will look this up by.
    """
    for gene in rule_set.get("genes") or []:
        match = re.search(r"query=([A-Za-z0-9_.-]+)", gene.get("@id", "") or "")
        symbol = gene.get("label") or (match.group(1) if match else None)
        if not symbol:
            continue
        mondo, inheritance = [], []
        for disease in gene.get("diseases") or []:
            identifier = disease.get("label") or ""
            if not identifier.startswith("MONDO:"):
                obo = re.search(r"MONDO_(\d+)", disease.get("@id", "") or "")
                identifier = f"MONDO:{obo.group(1)}" if obo else ""
            if identifier:
                mondo.append(identifier)
            for mode in disease.get("modeOfInheritance") or []:
                if mode.get("@label"):
                    inheritance.append(mode["@label"])
        yield symbol, sorted(set(mondo)), sorted(set(inheritance))


def collect(cache_dir: Path) -> tuple[dict, dict]:
    svis = fetch(SVIS_URL, cache_dir, "svis")["data"]
    gn_ids = []
    for svi in svis:
        match = re.search(r"/id/(GN\d+)", svi.get("@id", "") or "")
        if match:
            gn_ids.append(match.group(1))

    criteria_by_gene: dict[str, dict[str, str]] = defaultdict(dict)
    scope_by_gene: dict[str, dict[str, set]] = defaultdict(
        lambda: {"mondo": set(), "inheritance": set(), "svis": set()})
    empty_documents, failures = [], []

    for gn_id in gn_ids:
        try:
            document = fetch(DOC_URL.format(gn_id), cache_dir, f"svi_{gn_id}")
        except Exception as exc:
            failures.append({"svi": gn_id, "error": repr(exc)})
            continue
        stated_anything = False
        for rule_set in document.get("ruleSets") or []:
            codes = {}
            for code in rule_set.get("criteriaCodes") or []:
                label, value = code.get("label"), applicability(code)
                if label and value:
                    codes[label] = value
            for symbol, mondo, inheritance in gene_entries(rule_set):
                if codes:
                    criteria_by_gene[symbol].update(codes)
                    stated_anything = True
                scope_by_gene[symbol]["mondo"].update(mondo)
                scope_by_gene[symbol]["inheritance"].update(inheritance)
                scope_by_gene[symbol]["svis"].add(gn_id)
        if not stated_anything:
            empty_documents.append(gn_id)

    entries = {}
    for symbol, codes in sorted(criteria_by_gene.items()):
        scope = scope_by_gene[symbol]
        entries[symbol] = {
            "criteria": {code: codes[code] for code in ALL_CODES if code in codes},
            # A criterion the document never mentions is absent rather than guessed.
            "criteria_missing": [code for code in ALL_CODES if code not in codes],
            "mondo": sorted(scope["mondo"]),
            "inheritance": sorted(scope["inheritance"]),
            "svis": sorted(scope["svis"]),
        }
    stats = {
        "svis_listed": len(svis),
        "documents_fetched": len(gn_ids) - len(failures),
        "documents_without_criteria": len(empty_documents),
        "genes_named": len(scope_by_gene),
        "genes_with_criteria": len(entries),
        "genes_with_all_28": sum(1 for entry in entries.values() if not entry["criteria_missing"]),
        "failures": failures,
    }
    return entries, stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE,
                        help="Raw registry responses; present ones are reused, so a re-run is free")
    parser.add_argument("--report", action="store_true", help="Print a coverage summary")
    args = parser.parse_args()

    entries, stats = collect(args.cache_dir)
    document = {
        "schema_version": "1.0",
        "registry_status": "DRAFT",
        "purpose": (
            "Per-gene ACMG criterion applicability transcribed from public ClinGen VCEP criteria "
            "specifications (CSpec). Entries are a curator-review queue and MUST NOT be read as "
            "runtime evidence until a provider is written for them and a curator records "
            "reviewed_at. In particular, 'PVS1: Applicable' is a VCEP statement about its own "
            "specification and is NOT the same sentence as gene-disease-review-decisions.json's "
            "lof_mechanism_established; see test_data/collect_cspec_applicability.py."
        ),
        "source": "ClinGen Criteria Specification Registry (CSpec)",
        "source_url": SVIS_URL,
        "retrieved_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "entry_method": "machine_collected",
        "applicability_rule": (
            "A criterion is Applicable when at least one of its evidenceStrengths is Applicable, "
            "and Not applicable when every strength says it is not. Verified against the "
            "registry's own per-criterion applicability on GN180 (DYSF): all 28 agree."
        ),
        "not_verified_by_a_second_reader": True,
        "statistics": stats,
        "entries": entries,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.output} ({len(entries)} genes)")

    if args.report:
        for key, value in stats.items():
            if key != "failures":
                print(f"  {key:32} {value}")
        print(f"  {'failures':32} {len(stats['failures'])}")
        pvs1 = Counter(entry["criteria"].get("PVS1") for entry in entries.values()
                       if "PVS1" in entry["criteria"])
        print("  PVS1: " + ", ".join(f"{key}={value}" for key, value in pvs1.most_common()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
