"""
stubs.py
Placeholder entries for the seven ACMG/AMP 2015 codes that remain outside
the integrated literature and automated-evidence evaluators.

[Why these exist]
  classification.classify() needs a CriterionEvidence for every code it is
  asked about, or it has no way to distinguish "evaluated, not met" from
  "never evaluated at all" - and conflating those two would silently treat
  an unimplemented criterion as if it had been checked and found lacking,
  which overstates this project's actual coverage. stub_evidence() returns
  an honest NOT_EVALUATED placeholder instead.

[Why seven codes share one file instead of seven near-empty modules]
  These seven split into two kinds of remaining workflow gap:
    - PS2, PM3, PM6, BS2, BP2, BP5: literature-adjacent Layer-3 codes, but
      the evidence for them lives in the PATIENT's own clinical/genetic-
      testing records (parentage testing, phasing, the patient's own
      pedigree, a healthy-carrier record) - not in papers. Out of scope for
      a literature-reading LLM pipeline by design, not by omission.
    - PP4: a phenotype-specific criterion still awaiting its dedicated
      implementation.
  Neither kind needs its own bespoke prompt/schema/finalize() the way
  ps4.py or segregation.py do - they need exactly one thing, a
  NOT_EVALUATED marker - so one shared function serves all seven, rather than
  seven copies of the same four-line stub. Whichever team/person implements a
  real one later has one obvious interface to match: build a module shaped
  like ps4.py or segregation.py, then register it in registry.py in place
  of the stub.
"""

from __future__ import annotations

from acmg_pipeline.classification import ALL_ACMG_CODES, IMPLEMENTED_CODES, CriterionEvidence
from acmg_pipeline.gate import CriterionStatus

STUB_CODES = [code for code in ALL_ACMG_CODES if code not in IMPLEMENTED_CODES]

# These six criteria depend primarily on patient/family records. PP4 remains
# a phenotype-specific workflow not implemented by either integrated engine.
LITERATURE_ADJACENT_BUT_CLINICAL_RECORD = {"PS2", "PM3", "PM6", "BS2", "BP2", "BP5"}
PHENOTYPE_WORKFLOW = {"PP4"}


# CriterionEvidence's `status` field is acmg_pipeline.gate.CriterionStatus,
# which has MET/NOT_MET/UNKNOWN (not a literal NOT_EVALUATED member) -
# UNKNOWN is that enum's own "not evaluated" value (see gate.py: "the
# criterion code itself does not exist in ERepo for this variant" - the
# same concept, just phrased for ERepo lookups specifically). Reusing it
# here avoids introducing a fourth status value that classify() would have
# to treat identically to UNKNOWN anyway.
_NOT_EVALUATED = CriterionStatus.UNKNOWN


def stub_evidence(code: str, note: str = "") -> CriterionEvidence:
    """
    Returns a NOT_EVALUATED placeholder CriterionEvidence for `code` - no
    judgment has been attempted, by this pipeline or (as far as this
    project knows) anyone else yet.
    """
    if code not in STUB_CODES:
        raise ValueError(
            f"{code!r} is not a stub code - it either has a real "
            f"implementation ({sorted(IMPLEMENTED_CODES)}) or is not a "
            f"recognized ACMG/AMP 2015 code."
        )
    reason = (
        "depends on the patient's own clinical/genetic-testing records, not literature"
        if code in LITERATURE_ADJACENT_BUT_CLINICAL_RECORD
        else "phenotype-specific criterion not implemented by the integrated evaluators"
    )
    return CriterionEvidence(
        code=code,
        status=_NOT_EVALUATED,
        source=f"stub ({reason}){': ' + note if note else ''}",
    )


def all_stub_evidence() -> list[CriterionEvidence]:
    """Convenience: a NOT_EVALUATED CriterionEvidence for every stub code at once."""
    return [stub_evidence(code) for code in STUB_CODES]
