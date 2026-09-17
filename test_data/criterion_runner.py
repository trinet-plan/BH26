"""
test_data/criterion_runner.py

A generic harness for testing ONE ACMG code's judgment function against
this project's ground-truth dataset (test_data/full_criteria_ground_truth.
py). This is for whoever implements a Layer-1 code (see acmg_pipeline.
criteria.stubs.AUTOMATED_RULE_BASED_OTHER_TEAM) or PP4: plug your function
in, get a pass/fail report against real ERepo/demo-case data, without
writing your own comparison loop.

[The contract your judgment function must satisfy]
  def judge(
      variant: VariantRecord,
      clinical_note: ClinicalNoteExtraction,
  ) -> tuple[CriterionStatus, Optional[Strength]]:
      ...
      return CriterionStatus.MET, Strength.MODERATE   # or:
      return CriterionStatus.NOT_MET, None            # strength must be None when NOT_MET

  This is exactly the (status, strength) pair acmg_pipeline.classification.
  CriterionEvidence itself needs (see classification.py's CriterionEvidence.
  __post_init__), so a passing judge() can be wired straight into
  classification.classify() later - no adapter layer needed.

[Usage]
  from test_data.criterion_runner import run_criterion_check

  def judge_pvs1(variant, clinical_note) -> tuple[CriterionStatus, Optional[Strength]]:
      ...  # your real PVS1 decision-tree logic
      return CriterionStatus.MET, Strength.VERY_STRONG

  report = run_criterion_check("PVS1", judge_pvs1)
  report.print_summary()
  # or, to fail a CI-style script on any mismatch:
  import sys; sys.exit(1 if report.mismatched else 0)

Run this file directly for a worked example against a deliberately naive
placeholder judge function (see __main__ below) - it demonstrates the
report format without requiring a real implementation yet.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from acmg_pipeline.classification import Strength
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.gate import CriterionStatus
from acmg_pipeline.vcf_record import VariantRecord
from test_data.full_criteria_ground_truth import GroundTruthEntry, entries_for_criterion

JudgeFn = Callable[
    [VariantRecord, ClinicalNoteExtraction],
    tuple[CriterionStatus, Optional[Strength]],
]


@dataclass
class MismatchDetail:
    entry: GroundTruthEntry
    got_status: CriterionStatus
    got_strength: Optional[Strength]


@dataclass
class CriterionCheckResult:
    criterion: str
    total: int
    matched: int
    mismatched: list[MismatchDetail] = field(default_factory=list)
    errored: list[tuple[GroundTruthEntry, Exception]] = field(default_factory=list)

    def print_summary(self) -> None:
        print(f"[{self.criterion}] {self.matched}/{self.total} matched "
              f"({len(self.mismatched)} mismatched, {len(self.errored)} raised an error)")
        for m in self.mismatched:
            e = m.entry
            expected = f"{e.status.value}" + (f"/{e.strength.value}" if e.strength else "")
            got = f"{m.got_status.value}" + (f"/{m.got_strength.value}" if m.got_strength else "")
            print(f"  MISMATCH {e.gene} {e.hgvsc}: expected {expected}, got {got} "
                  f"[{e.source}/{e.tier}]" + (f"  ({e.notes})" if e.notes else ""))
        for e, exc in self.errored:
            print(f"  ERROR    {e.gene} {e.hgvsc}: {exc!r}")


def run_criterion_check(criterion: str, judge: JudgeFn) -> CriterionCheckResult:
    """
    Calls judge(variant, clinical_note) once per ground-truth entry for
    `criterion` and
    compares (status, strength) against that entry's known real value.
    Entries for the same (gene, hgvsc) with a different criterion are never
    passed to `judge` - you only ever see the criterion you're testing.
    """
    entries = entries_for_criterion(criterion)
    result = CriterionCheckResult(criterion=criterion, total=len(entries), matched=0)
    for e in entries:
        variant = VariantRecord(
            chrom="", pos=0, id="", ref="", alt="", qual="", filter="",
            info={"GENE": e.gene, "HGVSC": e.hgvsc, "HGVSP": e.hgvsp},
        )
        clinical_note = ClinicalNoteExtraction()
        try:
            got_status, got_strength = judge(variant, clinical_note)
        except Exception as exc:
            result.errored.append((e, exc))
            continue
        is_match = got_status == e.status and (got_strength == e.strength if got_status == CriterionStatus.MET else True)
        if is_match:
            result.matched += 1
        else:
            result.mismatched.append(MismatchDetail(e, got_status, got_strength))
    return result


if __name__ == "__main__":
    # Worked example: a deliberately naive placeholder for PVS1 that always
    # says NOT_MET, just to show what a report looks like before any real
    # logic exists. Replace this with an import of the real implementation
    # once it exists (e.g. `from acmg_pipeline.criteria.pvs1 import judge`).
    def placeholder_judge(
        variant: VariantRecord,
        clinical_note: ClinicalNoteExtraction,
    ) -> tuple[CriterionStatus, Optional[Strength]]:
        del variant, clinical_note
        return CriterionStatus.NOT_MET, None

    report = run_criterion_check("PVS1", placeholder_judge)
    report.print_summary()
