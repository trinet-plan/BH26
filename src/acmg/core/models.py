"""Small explicit domain contracts; source labels never become rule decisions."""

from dataclasses import asdict, dataclass, field
from enum import StrEnum
import re


CRITERIA = (
    "PVS1", "PS1", "PM1", "PM2", "PM4", "PM5", "PP2", "PP3", "PP5",
    "BA1", "BS1", "BP1", "BP3", "BP4", "BP6", "BP7",
)


class Status(StrEnum):
    MET = "MET"
    NOT_MET = "NOT_MET"
    NOT_EVALUATED = "NOT_EVALUATED"
    NOT_APPLICABLE = "NOT_APPLICABLE"
    MANUAL_REVIEW = "MANUAL_REVIEW"
    DEPRECATED = "DEPRECATED"


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
    status: Status
    variant: dict
    summary: str
    strength: str | None = None
    direction: str | None = None
    evidence_outcome: str | None = None
    evidence: list = field(default_factory=list)
    missing_inputs: list[str] = field(default_factory=list)
    review_points: list[str] = field(default_factory=list)
    conflict_flags: list[str] = field(default_factory=list)
    provenance: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.criterion not in CRITERIA:
            raise ValueError("Unsupported criterion")
        self.status = Status(self.status)
        if self.criterion in {"PP5", "BP6"} and self.status != Status.DEPRECATED:
            raise ValueError("Deprecated criteria cannot be scored")
        if self.criterion == "PVS1" and self.status == Status.MET:
            raise ValueError("PVS1 is provisional")
        if self.status != Status.MET and self.strength is not None:
            raise ValueError("Only MET can carry applied strength; keep recommendations in metadata")

    def to_dict(self):
        return asdict(self)
