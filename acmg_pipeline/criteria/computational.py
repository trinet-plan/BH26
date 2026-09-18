from acmg_pipeline.constants import CriterionStatus
"""Explicit calibration intervals, never vote counts or cross-predictor thresholds."""

from acmg_pipeline.criteria.common import annotation_context, get_evidence, result
from acmg_pipeline.services.population import number


STRENGTH_ORDER = {"supporting": 1, "moderate": 2, "strong": 3, "very_strong": 4}
# The criterion a calibration would support if this score pointed the other way. A
# calibration states both directions over one score range, so the opposing band is what
# distinguishes a predictor that contradicts this criterion from one that is simply silent.
OPPOSING_CRITERION = {"PP3": "BP4", "BP4": "PP3"}
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

    An unknown predictor version is never silently treated as the calibration version.  A
    deployment may nevertheless elect to use it through an explicit, source-bound policy;
    the UNKNOWN state and the policy then travel with the result for curator review.
    """
    if prediction.get("predictor") != calibration["predictor"]:
        return False, None
    if prediction.get("mechanism") != calibration["mechanism"]:
        return False, None
    if (prediction.get("predictor_version") == calibration["predictor_version"]
            and prediction.get("calibration_eligible") is not False):
        return True, None
    unknown_policy = calibration.get("unknown_version_policy")
    if not isinstance(unknown_policy, dict):
        # Compatibility with pre-existing cached evidence/configuration.  New provider
        # evidence uses predictor_version=UNKNOWN and the policy below.
        unknown_policy = None
    if unknown_policy is not None:
        if any(not unknown_policy.get(field) for field in
               ("source", "accepted_by", "justification")):
            return False, None
        if (prediction.get("source") == unknown_policy["source"]
                and prediction.get("predictor_version") == "UNKNOWN"):
            return True, {
                "version_status": "UNKNOWN",
                "predictor_version": "UNKNOWN",
                "source": prediction["source"],
                "source_version": prediction.get("source_version"),
                "retrieved_at": prediction.get("retrieved_at"),
                "accepted_by": unknown_policy["accepted_by"],
                "justification": unknown_policy["justification"],
            }
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
        return result(code, input_data, CriterionStatus.UNKNOWN, "Calibrated predictor policy unavailable",
                      missing=["computational.selected_calibrations"])
    consequences = set(annotation["consequences"])
    in_scope = [(name, item) for name, item in calibrations
                if consequences & set(item["consequences"])]
    if not in_scope:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Outside predictor calibration scope",
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
            return result(code, input_data, CriterionStatus.UNKNOWN, "Invalid or conflicting prediction scores",
                          evidence=[annotation, *matched], review=["Resolve prediction score provenance"])
        score = next(iter(scores))
        minimum, maximum = number(calibration.get("score_min")), number(calibration.get("score_max"))
        if minimum is None or maximum is None or not minimum <= score <= maximum:
            return result(code, input_data, CriterionStatus.UNKNOWN, "Score outside calibrated domain",
                          evidence=[annotation, *matched], missing=["valid_score_domain"])
        if not calibration.get("bands", {}).get(code):
            return result(code, input_data, CriterionStatus.UNKNOWN, "Criterion calibration intervals missing",
                          missing=[f"computational.bands.{code}"], evidence=[annotation, *matched])
        evidence.extend(matched)
        applied.append({
            "calibration": name, "mechanism": calibration["mechanism"],
            "predictor": calibration["predictor"], "score": str(score),
            "strength": band_strength(calibration, code, score),
            "opposing": band_strength(calibration, OPPOSING_CRITERION[code], score),
        })
    if not applied:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Matching predictor score unavailable",
                      missing=sorted(set(missing)), evidence=[annotation])
    provenance = {"applied_calibrations": [
        {key: item[key] for key in ("calibration", "mechanism", "predictor", "score", "strength",
                                    "opposing")}
        for item in applied]}
    if assertions:
        provenance["version_assertions"] = assertions
    if missing:
        provenance["unavailable_predictors"] = sorted(set(missing))
    supporting = [item for item in applied if item["strength"]]
    # A predicted effect on either protein or splicing can carry PP3, but benign computational
    # evidence requires every applicable mechanism to show no effect: a low protein score says
    # nothing about a disrupted splice site (ClinGen SVI, Walker et al. 2023).
    #
    # "Shows no effect" is not the same as "produced no evidence", and only the first
    # contradicts BP4. Both calibrations in use leave an explicit gap between their two
    # bands - REVEL 0.290-0.644 (Pejaver et al. 2022), SpliceAI 0.1-0.2 (Walker et al. 2023)
    # - and a score inside it is the calibration declining to call the variant in either
    # direction. Reading that silence as a contradiction made every indeterminate mechanism
    # veto BP4, which is how MYH7 c.4472C>G (REVEL 0.398, SpliceAI 0) came back not_met
    # against a curator's BP4: the protein predictor abstained and was counted as objecting.
    #
    # So only a mechanism scoring inside the OPPOSING criterion's own band blocks. An
    # abstaining mechanism cannot support BP4 either, so it is reported rather than
    # absorbed: the benign call then rests on the mechanisms that did answer.
    blocking = [item for item in applied if item["opposing"]] if code == "BP4" else []
    if blocking:
        provenance["benign_blocked_by"] = [
            {"predictor": item["predictor"], "mechanism": item["mechanism"],
             "score": item["score"], "opposing_criterion": OPPOSING_CRITERION[code],
             "opposing_strength": item["opposing"]}
            for item in blocking]
        names = ", ".join(f"{item['predictor']} {item['score']}" for item in blocking)
        return result(code, input_data, CriterionStatus.NOT_MET,
                      f"Another predicted mechanism supports {OPPOSING_CRITERION[code]} rather "
                      f"than a benign effect ({names})",
                      evidence=evidence, provenance=provenance)
    if not supporting:
        return result(code, input_data, CriterionStatus.NOT_MET, "Score does not meet criterion calibration",
                      evidence=evidence, provenance=provenance)
    uninformative = [item for item in applied if not item["strength"] and not item["opposing"]]
    review = []
    if code == "BP4" and uninformative:
        # A benign call carried by the mechanisms that answered, while another abstained, is
        # narrower than "computational evidence suggests no impact" reads. The curator is told
        # which mechanism was left open rather than having to infer it from the scores.
        provenance["uninformative_mechanisms"] = [
            {"predictor": item["predictor"], "mechanism": item["mechanism"],
             "score": item["score"]} for item in uninformative]
        review = [f"{item['predictor']} ({item['score']}) falls between this calibration's "
                  f"{code} and {OPPOSING_CRITERION[code]} intervals, so {item['mechanism']} "
                  f"impact is neither supported nor excluded" for item in uninformative]
    if code == "BP4" and annotation.get("high_confidence_null_or_splice") is True:
        return result(code, input_data, CriterionStatus.UNKNOWN, "Benign prediction conflicts with null/splice annotation",
                      evidence=evidence, review=["Reconcile null/splice evidence"], provenance=provenance)
    # Mechanisms are alternatives, never additive: the strongest applicable interval is used.
    strength = max((item["strength"] for item in supporting), key=STRENGTH_ORDER.get)
    return result(code, input_data, CriterionStatus.MET, "Score meets calibrated interval",
                  strength=strength, evidence=evidence, review=review, provenance=provenance)
