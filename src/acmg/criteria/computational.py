"""Explicit calibration intervals, never vote counts or cross-predictor thresholds."""

from acmg.core.models import Status
from acmg.criteria.common import annotation_context, get_evidence, result
from acmg.services.population import number


STRENGTH_ORDER = {"supporting": 1, "moderate": 2, "strong": 3, "very_strong": 4}
REQUIRED = ("predictor", "predictor_version", "source", "version", "mechanism", "consequences")


def selected_calibrations(config):
    section = config.get("computational", {})
    names = section.get("selected_calibrations")
    if names is None:
        name = section.get("selected_calibration")
        names = [name] if name else []
    definitions = section.get("calibrations", {})
    return [(name, definitions.get(name, {})) for name in names]


def accepts(calibration, prediction):
    """A calibration applies only to the release it was derived from.

    Ensembl REST serves some predictors without naming the release behind them. Such a score
    is usable only when the calibration itself names that source and states which release is
    being assumed, so the assumption travels with the result instead of disappearing.
    """
    if prediction.get("predictor") != calibration["predictor"]:
        return False, None
    if prediction.get("mechanism") != calibration["mechanism"]:
        return False, None
    if (prediction.get("predictor_version") == calibration["predictor_version"]
            and prediction.get("calibration_eligible") is not False):
        return True, None
    assertion = calibration.get("version_assertion")
    if not isinstance(assertion, dict):
        return False, None
    if any(not assertion.get(field) for field in
           ("source", "unreported_version", "asserted_version", "asserted_by", "justification")):
        return False, None
    if (prediction.get("source") == assertion["source"]
            and prediction.get("predictor_version") == assertion["unreported_version"]):
        return True, assertion
    return False, None


def band_strength(calibration, code, score):
    strengths = []
    for band in calibration.get("bands", {}).get(code) or []:
        low, high = number(band.get("min")), number(band.get("max"))
        strength = band.get("strength")
        if low is None or high is None or low > high or strength not in STRENGTH_ORDER:
            raise ValueError("Invalid computational calibration interval")
        if low <= score <= high:
            strengths.append(strength)
    return max(strengths, key=STRENGTH_ORDER.get) if strengths else None


def evaluate_prediction(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    calibrations = [(name, item) for name, item in selected_calibrations(config)
                    if all(item.get(key) for key in REQUIRED)]
    if not calibrations:
        return result(code, input_data, Status.NOT_EVALUATED, "Calibrated predictor policy unavailable",
                      missing=["computational.selected_calibrations"])
    consequences = set(annotation["consequences"])
    in_scope = [(name, item) for name, item in calibrations
                if consequences & set(item["consequences"])]
    if not in_scope:
        return result(code, input_data, Status.NOT_APPLICABLE, "Outside predictor calibration scope",
                      evidence=[annotation])
    available = get_evidence("computational", input_data, services)
    applied, missing, evidence, assertions = [], [], [annotation], []
    for name, calibration in in_scope:
        matched = []
        for prediction in available:
            accepted, assertion = accepts(calibration, prediction)
            if accepted:
                matched.append(prediction)
                if assertion and assertion not in assertions:
                    assertions.append(assertion)
        if not matched:
            missing.append(calibration["predictor"])
            continue
        scores = {number(item.get("score")) for item in matched}
        if None in scores or len(scores) != 1:
            return result(code, input_data, Status.MANUAL_REVIEW, "Invalid or conflicting prediction scores",
                          evidence=[annotation, *matched], review=["Resolve prediction score provenance"])
        score = next(iter(scores))
        minimum, maximum = number(calibration.get("score_min")), number(calibration.get("score_max"))
        if minimum is None or maximum is None or not minimum <= score <= maximum:
            return result(code, input_data, Status.NOT_EVALUATED, "Score outside calibrated domain",
                          evidence=[annotation, *matched], missing=["valid_score_domain"])
        if not calibration.get("bands", {}).get(code):
            return result(code, input_data, Status.NOT_EVALUATED, "Criterion calibration intervals missing",
                          missing=[f"computational.bands.{code}"], evidence=[annotation, *matched])
        evidence.extend(matched)
        applied.append({
            "calibration": name, "mechanism": calibration["mechanism"],
            "predictor": calibration["predictor"], "score": str(score),
            "strength": band_strength(calibration, code, score),
        })
    if not applied:
        return result(code, input_data, Status.NOT_EVALUATED, "Matching predictor score unavailable",
                      missing=sorted(set(missing)), evidence=[annotation])
    provenance = {"applied_calibrations": [
        {key: item[key] for key in ("calibration", "mechanism", "predictor", "score", "strength")}
        for item in applied]}
    if assertions:
        provenance["version_assertions"] = assertions
    if missing:
        provenance["unavailable_predictors"] = sorted(set(missing))
    supporting = [item for item in applied if item["strength"]]
    # A predicted effect on either protein or splicing can carry PP3, but benign computational
    # evidence requires every applicable mechanism to show no effect: a low protein score says
    # nothing about a disrupted splice site (ClinGen SVI, Walker et al. 2023).
    blocking = [item for item in applied if not item["strength"]] if code == "BP4" else []
    if blocking:
        provenance["benign_blocked_by"] = [
            {"predictor": item["predictor"], "mechanism": item["mechanism"], "score": item["score"]}
            for item in blocking]
        return result(code, input_data, Status.NOT_MET,
                      "Another predicted mechanism does not support benign computational evidence",
                      evidence=evidence, provenance=provenance)
    if not supporting:
        return result(code, input_data, Status.NOT_MET, "Score does not meet criterion calibration",
                      evidence=evidence, provenance=provenance)
    if code == "BP4" and annotation.get("high_confidence_null_or_splice") is True:
        return result(code, input_data, Status.MANUAL_REVIEW, "Benign prediction conflicts with null/splice annotation",
                      evidence=evidence, review=["Reconcile null/splice evidence"], provenance=provenance)
    # Mechanisms are alternatives, never additive: the strongest applicable interval is used.
    strength = max((item["strength"] for item in supporting), key=STRENGTH_ORDER.get)
    return result(code, input_data, Status.MET, "Score meets calibrated interval",
                  strength=strength, evidence=evidence, provenance=provenance)
