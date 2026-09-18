from acmg_pipeline.constants import CriterionStatus
"""Small explicit domain contracts; source labels never become rule decisions."""

from dataclasses import asdict, dataclass, field
import re
from typing import Optional

from acmg_pipeline.constants import AUTOMATED_CRITERIA, CriterionStatus


CRITERIA = AUTOMATED_CRITERIA


@dataclass(frozen=True)
class Variant:
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str

    def __post_init__(self):
        if self.assembly != "GRCh38":
            raise ValueError("Only GRCh38 is supported")
        if type(self.pos) is not int or self.pos < 1:
            raise ValueError("POS must be a positive 1-based integer")
        chrom = self.chrom.removeprefix("chr")
        if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
            raise ValueError("Unsupported chromosome")
        object.__setattr__(self, "chrom", chrom)
        for allele in (self.ref, self.alt):
            if not re.fullmatch("[ACGT]+", allele):
                raise ValueError("REF/ALT must be nonempty ACGT alleles; missing/SV ALT unsupported")
        if self.ref == self.alt:
            raise ValueError("REF and ALT must differ")

    @property
    def key(self):
        return f"{self.assembly}:{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    def to_dict(self):
        return asdict(self)


@dataclass
class CriterionResult:
    criterion: str
    status: CriterionStatus
    variant: dict
    summary: str
    strength: Optional[str] = None
    direction: Optional[str] = None
    evidence_outcome: Optional[str] = None
    evidence: list = field(default_factory=list)
    # No separate `unresolved_requirements` field: PVS1's _finish() (the only
    # place that ever populated it) always set it to this exact same list -
    # tests/test_clingen_positive.py had an assertion proving the invariant
    # rather than a reason for the second field to exist. Removed 2026-09-18
    # per the user's direction to fold it into missing_inputs.
    missing_inputs: list[str] = field(default_factory=list)
    review_points: list[str] = field(default_factory=list)
    conflict_flags: list[str] = field(default_factory=list)
    provenance: dict = field(default_factory=dict)
    evaluation_context: Optional[dict] = None
    decision_trace: list[dict] = field(default_factory=list)
    rules_used: list[dict] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def __post_init__(self):
        if self.criterion not in CRITERIA:
            raise ValueError("Unsupported criterion")
        self.status = CriterionStatus(self.status)
        if self.criterion in {"PP5", "BP6"} and self.status != CriterionStatus.UNKNOWN:
            raise ValueError("PP5 and BP6 are retained as explicit UNKNOWN criteria")
        if self.status != CriterionStatus.MET and self.strength is not None:
            raise ValueError("Only MET can carry applied strength; keep recommendations in metadata")
        if (self.criterion == "PVS1" and self.status == CriterionStatus.MET and
                self.strength not in {"very_strong", "strong", "moderate", "supporting"}):
            raise ValueError("PVS1 MET requires a supported applied strength")

    def to_dict(self):
        value = asdict(self)
        if self.criterion != "PVS1":
            for key in ("evaluation_context", "decision_trace", "rules_used", "warnings"):
                value.pop(key)
        return value
