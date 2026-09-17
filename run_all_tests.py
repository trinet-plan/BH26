"""
run_all_tests.py

Runs every one of this project's hand-rolled test_*.py scripts (see
test_harness.py - no pytest/unittest here) as its own subprocess, in the
same order they've always been run manually, and prints a combined
per-file pass/fail summary at the end.

Each test file is still a fully standalone, individually-runnable script
(`python3 test_classification.py` keeps working exactly as before) - this
is only a convenience wrapper around what was previously a manual
`for f in test_*.py; do python3 "$f"; done` shell loop, so a full
regression check is one command instead of copy-pasting that loop.

Usage:
  python3 run_all_tests.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

TEST_FILES = [
    "test_ps3_bs3_judgment.py",
    "test_ps3_bs3_ps4_gate.py",
    "test_ps3_bs3_ps4_gate_full.py",
    "test_classification.py",
    "test_full_criteria_ground_truth.py",
    "test_automated_criteria_ground_truth.py",
    "test_pp4_pp1_bs4.py",
    "test_pp1_bs4_pp4_engine.py",
]


def main() -> None:
    root = Path(__file__).resolve().parent
    results: list[tuple[str, int]] = []

    for name in TEST_FILES:
        print(f"\n{'#'*70}\n# {name}\n{'#'*70}")
        proc = subprocess.run([sys.executable, str(root / name)], cwd=root)
        results.append((name, proc.returncode))

    print(f"\n{'='*70}\n[run_all_tests] Summary\n{'='*70}")
    all_ok = True
    for name, code in results:
        status = "PASS" if code == 0 else f"FAIL (exit={code})"
        print(f"  {status:16} {name}")
        all_ok = all_ok and code == 0

    print(f"{'='*70}")
    print("[run_all_tests] ALL SUITES PASSED" if all_ok else "[run_all_tests] ONE OR MORE SUITES FAILED")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
