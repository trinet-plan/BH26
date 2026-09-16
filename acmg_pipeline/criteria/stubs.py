"""
stubs.py
Placeholder entries for every ACMG/AMP 2015 code this project does not yet
implement a real judgment for - 23 of the 28 codes. Real implementations
exist only for PS3/BS3 (ps3_bs3.py), PS4 (ps4.py), and PP1/BS4
(segregation.py) - see acmg_pipeline.classification.IMPLEMENTED_CODES.

[Why these exist]
  classification.classify() needs a CriterionEvidence for every code it is
  asked about, or it has no way to distinguish "evaluated, not met" from
  "never evaluated at all" - and conflating those two would silently treat
  an unimplemented criterion as if it had been checked and found lacking,
  which overstates this project's actual coverage. stub_evidence() returns
  an honest NOT_EVALUATED placeholder instead.

[Why 23 codes share one file instead of 23 near-empty modules]
  Per doc/BH26_participant_briefing_v3_en.md's own Layer classification and
  design doc section 15-5, these 23 split into two genuinely different
  kinds of "not our problem (yet)":
    - PS2, PM3, PM6, BS2, BP2, BP5: literature-adjacent Layer-3 codes, but
      the evidence for them lives in the PATIENT's own clinical/genetic-
      testing records (parentage testing, phasing, the patient's own
      pedigree, a healthy-carrier record) - not in papers. Out of scope for
      a literature-reading LLM pipeline by design, not by omission.
    - PVS1, PS1, PM1, PM2, PM4, PM5, PP2, PP3, PP4, PP5, BA1, BS1, BP1,
      BP3, BP4, BP6, BP7: Layer-1 "Automated Evidence" codes that
      doc/acmg_tool_comparison_v1_en.md found existing tools (InterVar,
      AutoPVS1, gnomAD/ClinVar lookups, in-silico predictor thresholds)
      already handle well by rule/lookup - a different team's
      responsibility per this project's scope split (see the "他メンバーの
      作成しているcriteria" integration this module is meant to make easy),
      not something an LLM pipeline should reimplement from scratch.
  Neither kind needs its own bespoke prompt/schema/finalize() the way
  ps4.py or segregation.py do - they need exactly one thing, a
  NOT_EVALUATED marker - so one shared function serves all 23, rather than
  23 copies of the same four-line stub. Whichever team/person implements a
  real one later has one obvious interface to match: build a module shaped
  like ps4.py or segregation.py, then register it in registry.py in place
  of the stub.
"""

from __future__ import annotations

from acmg_pipeline.classification import ALL_ACMG_CODES, IMPLEMENTED_CODES, CriterionEvidence
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.gate import CriterionStatus
from acmg_pipeline.vcf_record import VariantRecord

STUB_CODES = [code for code in ALL_ACMG_CODES if code not in IMPLEMENTED_CODES]

# Documents *why* each stub code is out of scope for this LLM pipeline
# specifically (as opposed to "nobody has gotten to it yet") - see the
# module docstring. Purely informational; stub_evidence() works the same
# for every STUB_CODE regardless of which reason applies.
LITERATURE_ADJACENT_BUT_CLINICAL_RECORD = {"PS2", "PM3", "PM6", "BS2", "BP2", "BP5"}
AUTOMATED_RULE_BASED_OTHER_TEAM = set(STUB_CODES) - LITERATURE_ADJACENT_BUT_CLINICAL_RECORD


# CriterionEvidence's `status` field is acmg_pipeline.gate.CriterionStatus,
# which has MET/NOT_MET/UNKNOWN (not a literal NOT_EVALUATED member) -
# UNKNOWN is that enum's own "not evaluated" value (see gate.py: "the
# criterion code itself does not exist in ERepo for this variant" - the
# same concept, just phrased for ERepo lookups specifically). Reusing it
# here avoids introducing a fourth status value that classify() would have
# to treat identically to UNKNOWN anyway.
_NOT_EVALUATED = CriterionStatus.UNKNOWN


def stub_evidence(
    code: str,
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    note: str = "",
) -> CriterionEvidence:
    """
    Returns a NOT_EVALUATED placeholder CriterionEvidence for `code` - no
    judgment has been attempted, by this pipeline or (as far as this
    project knows) anyone else yet.
    """
    # Stub criteria do not inspect either input yet, but accepting both here
    # keeps their public boundary identical to implemented criteria.
    del variant, clinical_note
    if code not in STUB_CODES:
        raise ValueError(
            f"{code!r} is not a stub code - it either has a real "
            f"implementation ({sorted(IMPLEMENTED_CODES)}) or is not a "
            f"recognized ACMG/AMP 2015 code."
        )
    reason = (
        "depends on the patient's own clinical/genetic-testing records, not literature"
        if code in LITERATURE_ADJACENT_BUT_CLINICAL_RECORD
        else "rule/lookup-based criterion, expected from another team's automated-evidence module"
    )
    return CriterionEvidence(
        code=code,
        status=_NOT_EVALUATED,
        source=f"stub ({reason}){': ' + note if note else ''}",
    )


def all_stub_evidence(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
) -> list[CriterionEvidence]:
    """Convenience: a NOT_EVALUATED CriterionEvidence for every stub code at once."""
    return [stub_evidence(code, variant, clinical_note) for code in STUB_CODES]
