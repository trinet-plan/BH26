"""
test_harness.py

Shared check()/pass-fail-counting logic for this project's hand-rolled
test_*.py scripts (test_ps3_bs3_judgment.py, test_ps3_bs3_ps4_gate.py,
test_ps3_bs3_ps4_gate_full.py, test_classification.py,
test_full_criteria_ground_truth.py). No pytest/unittest here by design -
these scripts predate this module and were written as plain top-to-bottom
`check(label, cond)` calls with a running tally, one process per file
(`python3 test_X.py`), each ending in `sys.exit(1 if failed else 0)` so CI
or a shell `&&` chain can treat them as a pass/fail gate.

[Why this module exists, 2026-09-17]
  Every one of those 5 files used to carry its own byte-for-byte identical
  ~10-line `passed = 0 / failed = 0 / def check(...)` block. Extracted here
  so there is exactly one definition of what "a check" means across the
  whole test suite, instead of 5 copies that could silently drift.

Each test file creates its OWN Harness() instance (not a module-level
global here) so multiple test files each importing this module never
share passed/failed state with each other.
"""

from __future__ import annotations

import sys


class Harness:
    def __init__(self, verbose: bool = True) -> None:
        self.passed = 0
        self.failed = 0
        self.failures: list[str] = []
        self.verbose = verbose

    def check(self, label: str, cond: bool) -> None:
        if cond:
            self.passed += 1
            if self.verbose:
                print(f"  OK   {label}")
        else:
            self.failed += 1
            self.failures.append(label)
            print(f"  FAIL {label}")

    def report_and_exit(self, width: int = 40, show_total: bool = False, list_failures: bool = False) -> None:
        """
        Prints the final summary box and calls sys.exit() with the same
        convention every one of these scripts already used: 0 if nothing
        failed, 1 otherwise.
        """
        total_suffix = f" (total {self.passed + self.failed})" if show_total else ""
        print(f"\n{'='*width}\n{self.passed} passed, {self.failed} failed{total_suffix}\n{'='*width}")
        if list_failures and self.failures:
            print("Failed cases:")
            for label in self.failures:
                print(f"  - {label}")
        sys.exit(1 if self.failed else 0)
