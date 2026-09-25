"""
stubs.py
Placeholder entries for the four ACMG/AMP 2015 codes that remain outside
every integrated evaluator (literature, automated, the PP1/BS4/PP4
phenotype-segregation engine added 2026-09-17, and the PS2/PM6 de-novo
rule added 2026-09-19).

[Why these exist]
  classification.classify() needs a CriterionEvidence for every code it is
  asked about, or it has no way to distinguish "evaluated, not met" from
  "never evaluated at all" - and conflating those two would silently treat
  an unimplemented criterion as if it had been checked and found lacking,
  which overstates this project's actual coverage. stub_evidence() returns
  an honest NOT_EVALUATED placeholder instead.

[Why the remaining stub codes share one file instead of one near-empty module each]
  PM3 and BP2 are literature-adjacent Layer-3 codes, but the evidence for
  them lives in the PATIENT's own clinical/genetic-testing records (trans/
  cis phasing against another variant) - not in papers, and not in any
  population database either. Out of scope for this project's provider
  pipeline by design, not by omission: no provider here has access to a
  second variant's phase relative to this one. Neither needs its own
  bespoke prompt/schema/finalize() the way ps4.py or segregation.py do -
  they need exactly one thing, a NOT_EVALUATED marker - so one shared
  function serves both, rather than two copies of the same four-line stub.
  Whichever team/person implements a real one later has one obvious
  interface to match: build a module shaped like ps4.py or segregation.py,
  then register it in registry.py in place of the stub.

[BS2 and BP5 are no longer here]
  BS2 moved to acmg_pipeline.constants.AUTOMATED_CODES / IMPLEMENTED_CODES
  on 2026-09-25, once acmg_pipeline.criteria.bs2 connected gnomAD's own
  homozygote_count/hemizygote_count fields as the project's explicit,
  documented reinterpretation of "observed in a healthy adult" (see
  bs2.py's own module docstring for why this is a deliberate scope choice,
  not a reversal of the reasoning below - and why it only covers
  autosomal-recessive/X-linked inheritance, not autosomal-dominant). BP5
  moved to acmg_pipeline.constants.LITERATURE_CODES on 2026-09-22 (see
  acmg_pipeline.criteria.bp5's own module docstring).

[PP1/BS4/PP4 and PS2/PM6 are no longer here]
  PP1/BS4/PP4 moved to acmg_pipeline.constants.PHENOTYPE_SEGREGATION_CODES /
  IMPLEMENTED_CODES on 2026-09-17, once acmg_pipeline.criteria.
  pp1_bs4_pp4_engine connected the real ClinGen 2024 evaluator (pulled in
  from r-kobayashi's pp4_pp1_bs4 branch). PS2/PM6 moved to
  acmg_pipeline.constants.DE_NOVO_CODES / IMPLEMENTED_CODES on 2026-09-19,
  once acmg_pipeline.criteria.ps2_pm6 connected a rule-based de-novo
  evaluator over ClinicalNoteExtraction.de_novo. Both still come back
  UNKNOWN when the input needed to decide is missing (no diagnosis to
  search the literature for; father/mother variant_status not reported) -
  the same honest-gap behavior stub_evidence() below provides - but that
  now happens inside each engine itself, not here.
"""

from __future__ import annotations

from acmg_pipeline.classification import CriterionEvidence
from acmg_pipeline.constants import IMPLEMENTED_CODES, STUB_CODES, CriterionStatus

# These two criteria depend on a second variant's trans/cis phase relative
# to this one - data no provider in this project's pipeline has access to.
LITERATURE_ADJACENT_BUT_CLINICAL_RECORD = {"PM3", "BP2"}


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
    reason = "depends on the patient's own clinical/genetic-testing records, not literature"
    return CriterionEvidence(
        code=code,
        status=_NOT_EVALUATED,
        source=f"stub ({reason}){': ' + note if note else ''}",
    )


def all_stub_evidence() -> list[CriterionEvidence]:
    """Convenience: a NOT_EVALUATED CriterionEvidence for every stub code at once."""
    return [stub_evidence(code) for code in sorted(STUB_CODES)]
