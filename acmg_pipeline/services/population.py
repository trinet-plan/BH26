"""Resolve population observations without turning an API miss into zero frequency."""

from decimal import Decimal, InvalidOperation


def number(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        result = Decimal(str(value))
        return result if result.is_finite() else None
    except InvalidOperation:
        return None


class PopulationService:
    def __init__(self, providers):
        self.providers = providers

    def get_resolved_evidence(self, variant, context=None):
        observations, failures = [], []
        for provider in self.providers:
            try:
                batch = provider.get_frequency(variant, context)
            except (OSError, ValueError) as exc:
                failures.append({"provider": provider.name, "reason": str(exc)})
                continue
            if batch is None:
                failures.append({"provider": provider.name, "reason": "NO_OBSERVATION"})
            else:
                observations.extend(batch)
        return {"observations": observations, "failures": failures}


def usable_observations(resolved, variant, minimum_an):
    """Keep raw AF, FAF and cohort observations separate. Do not pool overlapping cohorts."""
    valid, rejected = [], []
    for item in resolved["observations"]:
        reasons = []
        if item.get("variant_key") != variant.key:
            reasons.append("VARIANT_MISMATCH")
        if not all(item.get(key) for key in
                   ("source", "source_version", "retrieved_at", "population", "evidence_id")):
            reasons.append("MISSING_PROVENANCE")
        if item.get("quality_status") != "PASS":
            reasons.append("QUALITY_NOT_PASS")
        an, ac, af = (number(item.get(key)) for key in ("AN", "AC", "AF"))
        if an is None or an != an.to_integral_value() or an < minimum_an or an <= 0:
            reasons.append("INSUFFICIENT_AN")
        if ac is None or ac != ac.to_integral_value() or ac < 0 or (an is not None and ac > an):
            reasons.append("INVALID_AC")
        if af is None or not 0 <= af <= 1:
            reasons.append("INVALID_AF")
        if an and ac is not None and af is not None and abs(ac / an - af) > Decimal("0.000001"):
            reasons.append("AC_AN_AF_INCONSISTENT")
        if ac == 0 and item.get("callable") is not True:
            reasons.append("ABSENCE_WITHOUT_CALLABILITY")
        if reasons:
            rejected.append({"evidence_id": item.get("evidence_id"), "reasons": reasons})
        else:
            valid.append(item)
    return valid, rejected
