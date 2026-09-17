"""Deterministic PP4 point evaluator used by PP4/PP1/BS4 integration.

This module is intentionally independent from JSON files and from any
``reference_data.py`` module.  The caller passes an ordinary Python object
(dataclass or another object) that exposes the fields described by
``PP4ReferenceLike``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional, Protocol


class PP4ReferenceLike(Protocol):
    """Minimal Python-object interface required by ``evaluate_pp4``."""

    reference_id: str
    gene: str
    phenotype_label: str
    locus_model: str
    diagnostic_yield: float
    source_citation: str


# ClinGen 2024 Table 2: diagnostic yield/posterior probability -> Bayesian points.
# For values between rows, use the lower point value (floor to threshold).
DIAGNOSTIC_YIELD_POINT_TABLE: tuple[tuple[float, float], ...] = (
    (0.999, 12.0),
    (0.998, 11.5),
    (0.997, 11.0),
    (0.996, 10.5),
    (0.994, 10.0),
    (0.992, 9.5),
    (0.988, 9.0),
    (0.983, 8.5),
    (0.975, 8.0),
    (0.965, 7.5),
    (0.950, 7.0),
    (0.930, 6.5),
    (0.902, 6.0),
    (0.864, 5.5),
    (0.816, 5.0),
    (0.754, 4.5),
    (0.680, 4.0),
    (0.596, 3.5),
    (0.506, 3.0),
    (0.415, 2.5),
    (0.330, 2.0),
    (0.254, 1.5),
    (0.191, 1.0),
)

PP4_MIN_POSTERIOR = 0.191
LOCUS_EVIDENCE_CAP = 5.0


@dataclass
class PP4EvaluationResult:
    criterion: str = "PP4"
    applicable: bool = False
    reason: str = ""
    gene: Optional[str] = None
    reference_id: Optional[str] = None
    phenotype_label: Optional[str] = None
    phenotype_match: bool = False
    method_comparable: bool = False
    locus_model: Optional[str] = None
    diagnostic_yield: Optional[float] = None
    candidate_variants_on_allele: int = 1
    posterior_probability_for_variant: Optional[float] = None
    points_before_cap: float = 0.0
    points: float = 0.0
    pp1_support_allowed: Optional[bool] = None
    source_citation: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)


def diagnostic_yield_to_points(probability: float) -> float:
    """Convert diagnostic yield/posterior probability to ClinGen Table 2 points."""
    if not 0.0 <= probability <= 1.0:
        raise ValueError("probability must be between 0 and 1")

    for threshold, points in DIAGNOSTIC_YIELD_POINT_TABLE:
        if probability >= threshold:
            return points
    return 0.0


def evaluate_pp4(
    reference: PP4ReferenceLike,
    *,
    phenotype_match: bool,
    method_comparable: bool,
    candidate_variants_on_allele: int = 1,
) -> PP4EvaluationResult:
    """Evaluate PP4 from an ordinary Python reference object.

    No JSON serialization/deserialization is performed here.  ``reference`` may
    be the ``PP4ReferenceRecord`` dataclass defined in
    ``criteria/pp4_pp1_bs4.py``.
    """
    if candidate_variants_on_allele < 1:
        raise ValueError("candidate_variants_on_allele must be >= 1")

    if not 0.0 <= reference.diagnostic_yield <= 1.0:
        raise ValueError("diagnostic_yield must be between 0 and 1")

    result = PP4EvaluationResult(
        gene=reference.gene,
        reference_id=reference.reference_id,
        phenotype_label=reference.phenotype_label,
        phenotype_match=phenotype_match,
        method_comparable=method_comparable,
        locus_model=reference.locus_model,
        diagnostic_yield=reference.diagnostic_yield,
        candidate_variants_on_allele=candidate_variants_on_allele,
        source_citation=reference.source_citation,
    )

    if not phenotype_match:
        result.reason = "phenotype_not_matched_to_reference_cohort"
        return result

    if not method_comparable:
        result.reason = "testing_method_not_comparable_or_unknown"
        return result

    posterior = reference.diagnostic_yield / candidate_variants_on_allele
    result.posterior_probability_for_variant = posterior

    raw_points = diagnostic_yield_to_points(posterior)
    result.points_before_cap = raw_points

    if raw_points < 1.0:
        result.reason = "diagnostic_yield_below_pp4_lower_bound"
        result.pp1_support_allowed = True
        return result

    result.applicable = True
    result.points = min(raw_points, LOCUS_EVIDENCE_CAP)
    result.reason = "pp4_applicable"

    high_yield_homogeneous = (
        reference.locus_model == "homogeneous"
        and reference.diagnostic_yield > 0.90
    )
    result.pp1_support_allowed = not high_yield_homogeneous

    return result
