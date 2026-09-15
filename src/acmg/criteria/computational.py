"""Explicit calibration intervals, never vote counts or cross-predictor thresholds."""

from acmg.core.models import Status
from acmg.criteria.common import annotation_context, get_evidence, result
from acmg.services.population import number


STRENGTH_ORDER = {"supporting": 1, "moderate": 2, "strong": 3, "very_strong": 4}


def evaluate_prediction(code, input_data, services, config):
    early, annotation = annotation_context(code, input_data, services)
    if early:
        return early
    selected = config.get("computational", {}).get("selected_calibration")
    calibration = config.get("computational", {}).get("calibrations", {}).get(selected, {})
    required = ("predictor", "predictor_version", "source", "version", "mechanism", "consequences")
    if not all(calibration.get(key) for key in required):
        return result(code, input_data, Status.NOT_EVALUATED, "Calibrated predictor policy unavailable",
                      missing=["computational.selected_calibration"])
    if not set(annotation["consequences"]) & set(calibration["consequences"]):
        return result(code, input_data, Status.NOT_APPLICABLE, "Outside predictor calibration scope",
                      evidence=[annotation])
    predictions = [p for p in get_evidence("computational", input_data, services)
                   if p.get("predictor") == calibration["predictor"]
                   and p.get("predictor_version") == calibration["predictor_version"]
                   and p.get("mechanism") == calibration["mechanism"]]
    if not predictions:
        return result(code, input_data, Status.NOT_EVALUATED, "Matching predictor score unavailable",
                      missing=[calibration["predictor"]], evidence=[annotation])
    scores = {number(p.get("score")) for p in predictions}
    if None in scores or len(scores) != 1:
        return result(code, input_data, Status.MANUAL_REVIEW, "Invalid or conflicting prediction scores",
                      evidence=predictions, review=["Resolve prediction score provenance"])
    score = next(iter(scores))
    minimum, maximum = number(calibration.get("score_min")), number(calibration.get("score_max"))
    if minimum is None or maximum is None or not minimum <= score <= maximum:
        return result(code, input_data, Status.NOT_EVALUATED, "Score outside calibrated domain",
                      evidence=predictions, missing=["valid_score_domain"])
    bands = calibration.get("bands", {}).get(code)
    if not bands:
        return result(code, input_data, Status.NOT_EVALUATED, "Criterion calibration intervals missing",
                      missing=[f"computational.bands.{code}"])
    strengths = []
    for band in bands:
        low, high = number(band.get("min")), number(band.get("max"))
        strength = band.get("strength")
        if low is None or high is None or low > high or strength not in STRENGTH_ORDER:
            raise ValueError("Invalid computational calibration interval")
        if low <= score <= high:
            strengths.append(strength)
    evidence = [annotation, *predictions]
    provenance = {"calibration": calibration, "mechanism": calibration["mechanism"]}
    if not strengths:
        return result(code, input_data, Status.NOT_MET, "Score does not meet criterion calibration",
                      evidence=evidence, provenance=provenance)
    strength = max(strengths, key=STRENGTH_ORDER.get)
    if code == "BP4" and annotation.get("high_confidence_null_or_splice") is True:
        return result(code, input_data, Status.MANUAL_REVIEW, "Benign prediction conflicts with null/splice annotation",
                      evidence=evidence, review=["Reconcile null/splice evidence"], provenance=provenance)
    return result(code, input_data, Status.MET, "Score meets calibrated interval",
                  strength=strength, evidence=evidence, provenance=provenance)
