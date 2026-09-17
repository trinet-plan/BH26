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


# One-sided 95%: 95% of the distribution lies above the bound this returns.
_Z95 = Decimal("1.6448536269514722")
FAF_METHOD = "Wilson score interval, one-sided 95% lower bound, computed from AC/AN"


def faf95(ac, an):
    """Filtering allele frequency: the lower bound of the 95% CI for AC/AN.

    A point-estimate AF says nothing about how well observed it is: 4 alleles in 5,600 and
    400 in 560,000 both read 0.07%, but only the second is evidence that the frequency really
    is that high. BS1 asks whether a variant is *too common* to cause the disease, so the
    honest quantity is the lower bound - the frequency the data supports even in the worst
    case - which is what gnomAD publishes as FAF and what ClinGen SVI asks BS1/BA1 to use.

    gnomAD derives its own faf95 with a Poisson interval and does not return it from the
    GraphQL fields this project requests, so this is computed locally by the Wilson score
    method and is an approximation of that published value, not a copy of it. Every result
    records FAF_METHOD alongside the number.
    """
    ac, an = number(ac), number(an)
    if ac is None or an is None or an <= 0 or ac < 0 or ac > an:
        return None
    if ac == 0:
        return Decimal(0)
    proportion = ac / an
    z_squared = _Z95 * _Z95
    denominator = 1 + z_squared / an
    centre = (proportion + z_squared / (2 * an)) / denominator
    spread = (proportion * (1 - proportion) / an + z_squared / (4 * an * an)).sqrt()
    lower = centre - (_Z95 / denominator) * spread
    return lower if lower > 0 else Decimal(0)


STATISTICS = {"af": "allele frequency", "faf95": "filtering allele frequency"}


def article(measure):
    return ("an " if measure[0] in "aeiou" else "a ") + measure


def observed_frequencies(observations, statistic):
    """Pair each observation with the quantity `statistic` names.

    A frequency threshold and the statistic it was calibrated against are one unit - an AF
    cutoff compared with a FAF is a stricter cutoff than the one that was written down - so
    every criterion that compares against a threshold reads the statistic from its own policy
    and passes it here rather than picking one.
    """
    if statistic == "faf95":
        return [(item, faf95(item.get("AC"), item.get("AN"))) for item in observations]
    return [(item, number(item.get("AF"))) for item in observations]


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
