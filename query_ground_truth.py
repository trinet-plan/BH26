"""
query_ground_truth.py

Command-line interface to test_data/full_criteria_ground_truth.py's
GROUND_TRUTH dataset (1028 entries as of 2026-09-16, see test_case_ground_
truth.md) - filters it down to one or more criteria (plus optional
gene/source/tier filters) instead of the full 28-code dataset. Intended
for whoever is implementing a specific Layer-1 code or PP4 and wants just
that code's real (gene, hgvsc, status, strength) fixtures to validate
their own logic against, without writing any Python themselves.

Examples:
  python3 query_ground_truth.py --criterion PVS1
  python3 query_ground_truth.py --criterion PP4 --format json
  python3 query_ground_truth.py --criterion PM2 PP3 --gene PTEN
  python3 query_ground_truth.py --criterion BA1 --source erepo --format csv > ba1_cases.csv
  python3 query_ground_truth.py --list-criteria
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from acmg_pipeline.classification import ALL_ACMG_CODES
from test_data.collectors.full_criteria_ground_truth import GROUND_TRUTH, GroundTruthEntry


def _matches(e: GroundTruthEntry, args) -> bool:
    if args.criterion and e.criterion not in args.criterion:
        return False
    if args.gene and e.gene not in args.gene:
        return False
    if args.source and e.source != args.source:
        return False
    if args.tier and e.tier not in args.tier:
        return False
    if args.status and e.status.value != args.status:
        return False
    return True


def _to_dict(e: GroundTruthEntry) -> dict:
    return dict(
        gene=e.gene, hgvsc=e.hgvsc, hgvsp=e.hgvsp, criterion=e.criterion,
        status=e.status.value, strength=e.strength.value if e.strength else None,
        source=e.source, tier=e.tier, variant_outcome=e.variant_outcome, notes=e.notes,
    )


def _print_table(entries: list[GroundTruthEntry]) -> None:
    if not entries:
        print("(no matching entries)")
        return
    for e in entries:
        strength = e.strength.value if e.strength else "-"
        print(f"{e.criterion:5s} {e.status.value:8s} {strength:12s} "
              f"{e.gene:10s} {e.hgvsc:25s} [{e.source}/{e.tier}]"
              + (f"  outcome={e.variant_outcome}" if e.variant_outcome else ""))


def _print_json(entries: list[GroundTruthEntry]) -> None:
    print(json.dumps([_to_dict(e) for e in entries], ensure_ascii=False, indent=2))


def _print_csv(entries: list[GroundTruthEntry]) -> None:
    writer = csv.DictWriter(sys.stdout, fieldnames=[
        "gene", "hgvsc", "hgvsp", "criterion", "status", "strength",
        "source", "tier", "variant_outcome", "notes",
    ])
    writer.writeheader()
    for e in entries:
        writer.writerow(_to_dict(e))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--criterion", nargs="+", metavar="CODE",
                         help="one or more ACMG codes, e.g. --criterion PVS1 PM1")
    parser.add_argument("--gene", nargs="+", metavar="GENE", help="restrict to these genes")
    parser.add_argument("--source", choices=["erepo", "democase"], help="restrict to this data source")
    parser.add_argument("--tier", nargs="+", choices=["A", "B", "C"], help="restrict to these confidence tiers")
    parser.add_argument("--status", choices=["met", "not_met"], help="restrict to met or not_met entries")
    parser.add_argument("--format", choices=["table", "json", "csv"], default="table")
    parser.add_argument("--list-criteria", action="store_true",
                         help="list all 28 ACMG codes with their entry counts, then exit")
    args = parser.parse_args()

    if args.list_criteria:
        from collections import Counter
        counts = Counter(e.criterion for e in GROUND_TRUTH)
        for code in ALL_ACMG_CODES:
            print(f"{code:5s} {counts.get(code, 0)}")
        return

    if args.criterion:
        for code in args.criterion:
            if code not in ALL_ACMG_CODES:
                parser.error(f"{code!r} is not a recognized ACMG/AMP 2015 code")

    matched = [e for e in GROUND_TRUTH if _matches(e, args)]

    if args.format == "table":
        _print_table(matched)
        print(f"\n{len(matched)} entrie(s) matched", file=sys.stderr)
    elif args.format == "json":
        _print_json(matched)
    elif args.format == "csv":
        _print_csv(matched)


if __name__ == "__main__":
    main()
