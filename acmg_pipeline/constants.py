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
IMPLEMENTED_CODES = LITERATURE_CODES | AUTOMATED_CODES
STUB_CODES = frozenset(ALL_ACMG_CODES) - IMPLEMENTED_CODES
AUTOMATED_CRITERIA = tuple(code for code in ALL_ACMG_CODES if code in AUTOMATED_CODES)


class CriterionStatus(str, Enum):
    """Every criterion is either applied, ruled out, or explicitly unknown."""

    MET = "met"
    NOT_MET = "not_met"
    UNKNOWN = "unknown"
