"""Canonical ACMG criterion groups and the pipeline's three assessment states."""

from __future__ import annotations

from enum import Enum


PATHOGENIC_CODES = (
    "PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4",
    "PM5", "PM6", "PP1", "PP2", "PP3", "PP4", "PP5",
)
BENIGN_CODES = ("BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7")
ALL_ACMG_CODES = PATHOGENIC_CODES + BENIGN_CODES

# The literature workflow deliberately owns only these three criteria.  The
# remaining historical literature criteria remain explicit UNKNOWN stubs.
LITERATURE_CODES = frozenset({"PS3", "BS3", "PS4"})
AUTOMATED_CODES = frozenset({
    "PVS1", "PS1", "PM1", "PM2", "PM4", "PM5", "PP2", "PP3", "PP5",
    "BA1", "BS1", "BP1", "BP3", "BP4", "BP6", "BP7",
})
# Connected 2026-09-17 via acmg_pipeline.criteria.pp1_bs4_pp4_engine (the
# ClinGen 2024 PP4+PP1/BS4 Bayesian-points evaluator pulled in from
# r-kobayashi's pp4_pp1_bs4 branch, see that module's own docstring). Real
# evidence for these 3 requires a diagnostic-yield statistic that
# acmg_pipeline.pp4_literature_search finds via a live PubMed+LLM search -
# with no diagnosis to search for, or nothing confirmed found, they
# correctly evaluate to UNKNOWN, the same "honest gap" behavior every
# other unimplemented code already uses, so keeping them in
# IMPLEMENTED_CODES (as opposed to STUB_CODES) does not overstate today's
# real coverage.
PHENOTYPE_SEGREGATION_CODES = frozenset({"PP1", "BS4", "PP4"})
IMPLEMENTED_CODES = LITERATURE_CODES | AUTOMATED_CODES | PHENOTYPE_SEGREGATION_CODES
STUB_CODES = frozenset(ALL_ACMG_CODES) - IMPLEMENTED_CODES
AUTOMATED_CRITERIA = tuple(code for code in ALL_ACMG_CODES if code in AUTOMATED_CODES)


class CriterionStatus(str, Enum):
    """Every criterion is either applied, ruled out, or explicitly unknown."""

    MET = "met"
    NOT_MET = "not_met"
    UNKNOWN = "unknown"
