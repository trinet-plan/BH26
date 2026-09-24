"""
ps3_bs3_judgment.py
PS3 / BS3 judgment logic (prompt design + structured output + human-curator hints)

[Role in the pipeline]
  This is the downstream logic that actually reads the paper text and judges
  direction (damaging/normal) for variants that ps3_bs3_ps4_gate.py's variant
  matching gate flagged with should_run_pipeline=True (i.e., not yet curated in
  ClinGen, or curated but this specific criterion code was not evaluated).

[Carried over from the design document]
  - Target direction only; strength (Very Strong .. Supporting) is out of scope
    (design doc section 1, constraint 1)
  - Variant identification is an explicit, mandatory gate (constraint 2). Here
    this means "does the paper text actually test the target variant?" - a
    second, independent gate from ps3_bs3_ps4_gate.match_variant_in_text
    (different purpose, different input: that one matches against ERepo HGVS
    lists, this one matches against free-text paper content)
  - When judgment is not possible, not_clear is a first-class output value
    (constraint 3)

[Added in this round, reflecting findings from design doc sections 6-11]
  - Cross-checking against the approved-assay list (calls
    ps3_bs3_ps4_gate.check_approved_assay as post-processing)
  - Explicit self-report of "is this judgment based on a single study only"
    (section 9: the three PI3K-pathway variants had PS3 itself excluded from
    evaluation because other strong criteria already applied - a recurring
    administrative-exclusion pattern to watch for)
  - An explicit flag for "does this paper contain multiple assay types, with a
    risk of citing the wrong one" (section 6: the MSH2 p.Gly692Val case, where
    an indirect splicing table was mistakenly used instead of the actual
    functional assay)
  - Human-curator hints (curator_hints) are emitted as a field independent of
    the judgment itself

[Still not implemented at this stage]
  - The actual LLM API call (build_prompt() only returns the prompt string;
    which model / which structured-output mechanism - function calling vs JSON
    mode - is used is the caller's responsibility)
  - PS4-specific evidence aggregation (stacking case counts across multiple
    PMIDs) is out of scope here; this module is PS3/BS3 only.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.inputs import variant_identity
from acmg_pipeline.vcf_record import VariantRecord

from acmg_pipeline.gate import check_approved_assay, AssayApplicability
from acmg_pipeline.common import (
    MatchStatus, VariantMatchingResult, CuratorHint, FinalResult,
    PaperContribution, AggregatedJudgment,
    aggregate_multi_paper_results as _generic_aggregate_multi_paper_results,
)


def _clean_enum_token(raw: str) -> str:
    """
    Defensive normalization against the LLM's habit of appending an
    explanation to an enum value (e.g. "PS3(damaging)", "PS3 (damaging)").
    Strips everything from the first parenthesis onward and trims whitespace,
    keeping only the leading token.

    Added after gemma-4 was observed to actually output "PS3(damaging)" in
    practice. The prompt (PROMPT_TEMPLATE) has since been tightened to avoid
    this ambiguity, but since other models could exhibit the same habit, this
    parsing-side safeguard is kept as a second line of defense.
    """
    if not isinstance(raw, str):
        return raw
    cleaned = re.split(r"[\(（]", raw.strip())[0].strip()
    return cleaned


# ============================================================================
# 1. Structured output schema
# ============================================================================

class ResultDirection(str, Enum):
    FUNCTIONALLY_ABNORMAL = "functionally_abnormal"
    FUNCTIONALLY_NORMAL = "functionally_normal"
    INTERMEDIATE = "intermediate"
    MIXED = "mixed"
    UNCLEAR = "unclear"


class OverallDirection(str, Enum):
    PS3 = "PS3"
    BS3 = "BS3"
    NOT_CLEAR = "not_clear"


@dataclass
class ExperimentExtraction:
    assay_type: str
    experimental_system: str
    readout: str
    comparator: str
    result_direction: ResultDirection
    key_findings: list[str] = field(default_factory=list)
    model_system_caveat: Optional[str] = None


@dataclass
class OverallEvidence:
    direction: OverallDirection
    strength_hint: str = "not_clear"  # fixed; this design does not target strength
    rationale: str = ""


@dataclass
class PS3BS3Judgment:
    """Data class mirroring the LLM's structured output directly."""
    variant_matching: VariantMatchingResult
    experiments: list[ExperimentExtraction]
    overall_evidence: OverallEvidence
    # Added in this round: have the LLM self-report a map of the paper's assays
    multiple_assay_types_in_paper: bool = False
    single_study_only: bool = True  # self-assessed: no other study found as evidence, based on citations/mentions

    @staticmethod
    def from_json(data: dict) -> "PS3BS3Judgment":
        vm = data["variant_matching"]
        if isinstance(vm, str):
            # Defensive fallback: the LLM sometimes flattens variant_matching
            # to a bare string (e.g. "matched") instead of the expected
            # object {"match_status": "matched", ...}. Observed in practice
            # (2026-09-15, gemma-4, English prompt) for both MYH7 and PTEN in
            # the same run, causing "string indices must be integers" crashes
            # when the object-style access below was attempted on a string.
            vm = {"match_status": vm}
        return PS3BS3Judgment(
            variant_matching=VariantMatchingResult(
                match_status=MatchStatus(_clean_enum_token(vm["match_status"])),
                match_type=vm.get("match_type"),
                confidence=vm.get("confidence", "low"),
                notes=vm.get("notes", ""),
            ),
            # `or []`, not `.get(key, [])`: dict.get()'s default only fires when
            # the key is absent, and Claude (observed 2026-09-19, switching
            # LLM_PROVIDER away from gemma-4) emits "experiments": null and
            # "key_findings": null outright for the no-functional-assay case,
            # which .get(key, []) still returns as None, not []. `for e in None`
            # then crashes with "'NoneType' object is not iterable" - same class
            # of per-model output-shape quirk as the bare-string variant_matching
            # fallback above, just the opposite direction (present-but-null
            # instead of a flattened type).
            experiments=[
                ExperimentExtraction(
                    assay_type=e["assay_type"],
                    experimental_system=e["experimental_system"],
                    readout=e["readout"],
                    comparator=e["comparator"],
                    result_direction=ResultDirection(_clean_enum_token(e["result_direction"])),
                    key_findings=e.get("key_findings") or [],
                    model_system_caveat=e.get("model_system_caveat"),
                )
                for e in (data.get("experiments") or [])
            ],
            # `data.get("overall_evidence") or {}`, not `data["overall_evidence"]`:
            # same class of per-model output-shape quirk as `experiments` above -
            # Claude (observed 2026-09-19, PS3/BS3 judgments for an
            # match_status=unsuccessful paper) emits "overall_evidence": null
            # outright rather than omitting the key, and `data["overall_evidence"]`
            # then crashes with "'NoneType' object is not subscriptable" on the
            # ["direction"] access that follows.
            overall_evidence=OverallEvidence(
                direction=OverallDirection(_clean_enum_token(
                    (data.get("overall_evidence") or {}).get("direction", "not_clear"))),
                strength_hint=(data.get("overall_evidence") or {}).get("strength_hint", "not_clear"),
                rationale=(data.get("overall_evidence") or {}).get("rationale", ""),
            ),
            multiple_assay_types_in_paper=data.get("multiple_assay_types_in_paper", False),
            single_study_only=data.get("single_study_only", True),
        )

    def structured_evidence_items(self) -> list[dict]:
        """
        One checklist-style item per extracted experiment, for a curator-
        facing UI. Added 2026-09-15 per a real curator UI requirements doc
        (doc/recs for expert board.docx): its PM1 evidence-card mockup shows
        a definition line plus a checklist of individually checked/unchecked
        findings, each with its own supporting detail (e.g. "15 pathogenic
        variants found in a 21bp region..."). Free-text rationale alone
        can't drive that kind of UI, so this exposes the same per-experiment
        data that's already in `experiments` as a stable, generic shape
        (see va_spec_export.py, which reads this via duck typing - it never
        imports PS3BS3Judgment or any other criterion-specific class).

        `checked` here means "this experiment reached a definitive
        functionally_abnormal/functionally_normal call" (as opposed to
        intermediate/mixed/unclear) - it is a fact about that one
        experiment, independent of the paper's overall adopted direction.
        """
        items = []
        for exp in self.experiments:
            checked = exp.result_direction in (
                ResultDirection.FUNCTIONALLY_ABNORMAL, ResultDirection.FUNCTIONALLY_NORMAL,
            )
            detail_parts = [f"{exp.experimental_system}, readout: {exp.readout}, vs {exp.comparator}"]
            if exp.key_findings:
                detail_parts.append("; ".join(exp.key_findings))
            if exp.model_system_caveat:
                detail_parts.append(f"(caveat: {exp.model_system_caveat})")
            items.append({
                "label": f"{exp.assay_type}: {exp.result_direction.value}",
                "checked": checked,
                "detail": " - ".join(detail_parts),
            })
        return items


# ============================================================================
# 2. Prompt template (extended from design doc section 2-4)
# ============================================================================
#
# Switched to English on 2026-09-15: the LLM exchange itself will ultimately
# run in English (the target production model reasons in English), so the
# prompt sent to the model is now in English. Code comments and docstrings in
# this repository are also in English per the same request - this is a
# separate decision from what language the chat with the user is conducted in.

PROMPT_TEMPLATE = """\
You are assisting a clinical genetics curator in evaluating ACMG/AMP PS3/BS3 functional evidence.

Target variant: {gene} {hgvsc} ({hgvsp}, aliases: {equivalents})
Full text of the paper: {full_text}

Follow these steps to reach a judgment.

1. Variant identification: does this paper directly test the target variant?
   match_status: matched / heuristic / single_variant_study / unsuccessful
   * If unsuccessful, set everything below to null and stop here.
   * Even if the paper's primary focus is a different variant, matched is
     acceptable as long as the target variant is explicitly tested as a
     named comparator alongside it.
   * IMPORTANT - large-scale screens: some papers report a saturation
     mutagenesis or massively-parallel functional screen covering thousands
     of variants at once (e.g., deep mutational scanning, MITE libraries,
     VAMP-seq-style abundance screens). For these, the running text usually
     only calls out a handful of representative variants by name/number; the
     specific quantitative result (a score, a Z-score, an enrichment value)
     for most individual variants - including possibly the target variant -
     is often available only in a supplementary data table that is not part
     of the text you were given. If you cannot find the target variant's own
     specific reported value or classification anywhere in the given text
     (only a general description of the screen's methodology, or a mention
     that the dataset "covers all positions" without stating this variant's
     own result), you must set match_status to unsuccessful, even though the
     paper's overall subject matter is clearly relevant. Do not infer this
     variant's result by generalizing from the paper's aggregate findings or
     from other variants' results.

2. Map of the assays in the paper (important - check this before anything else):
   Determine whether this paper uses more than one distinct experimental
   approach for the target variant (e.g., a splicing assay and a separate
   functional assay). If so, set multiple_assay_types_in_paper: true, and
   before proceeding, explicitly identify which experiment should actually be
   cited as the basis for PS3/BS3 (normally the one that most directly
   measures function; do not base the judgment on an indirect result such as
   a splicing-pattern table alone).

3. Experiment extraction (only if match_status is matched / heuristic /
   single_variant_study):
   For each experiment, report: assay_type, experimental_system (cell line /
   mouse / patient-derived, etc.), readout, comparator, and result_direction
   (functionally_abnormal/functionally_normal/intermediate/mixed/unclear).
   IMPORTANT: judge each experiment's result_direction independently, based
   solely on that experiment's own reported numbers or findings compared to
   its own wild-type/control comparator. Do NOT adjust an individual
   experiment's result_direction to match the overall conclusion or the
   direction of other experiments in the same paper. For example, if an
   assay's own reported value is described as wild-type-like, not
   significantly different from control, or within the normal range, that
   experiment's result_direction must be functionally_normal (or
   intermediate/unclear as appropriate) even if other experiments in the same
   paper point the other way, and even if your overall conclusion will end up
   being PS3 or BS3. The overall_evidence.direction field in step 6 is where
   you weigh and combine the experiments; result_direction per experiment is
   not the place to do that weighing.

4. Model system caveat:
   For in vivo animal models (especially mouse knock-ins), note in
   model_system_caveat that this is not direct evidence from human tissue.

5. Self-reported reproducibility:
   Based on citations/mentions in the text, can you identify any other study
   reporting independent functional evidence for the same variant besides
   this paper? If not (i.e., this paper appears to be the only source of
   evidence), set single_study_only: true.

6. Overall judgment:
   direction: output exactly one of the following three strings, with no
   parentheses or added explanation: "PS3" / "BS3" / "not_clear"
   (for reference only, not to be included in the output: PS3 = evidence
   supporting a damaging functional effect, BS3 = evidence supporting a
   normal functional effect, not_clear = evidence is insufficient or
   conflicting)
   strength_hint: fixed to "not_clear" (strength is not targeted by this system)
   rationale: 2-3 sentences of justification

Output a single JSON object with exactly these keys:
variant_matching, experiments, multiple_assay_types_in_paper, single_study_only,
overall_evidence

IMPORTANT: variant_matching, each entry in experiments, and overall_evidence
must each be a JSON object (with the sub-fields described above), never a
bare string. For example, variant_matching must look like
{{"match_status": "matched", "match_type": "protein_notation", "confidence": "high", "notes": "..."}}
- NOT simply "matched" on its own. Below is a minimal example of the overall
expected shape (with placeholder content):

{{
  "variant_matching": {{"match_status": "matched", "match_type": "protein_notation", "confidence": "high", "notes": "..."}},
  "experiments": [
    {{"assay_type": "...", "experimental_system": "...", "readout": "...", "comparator": "...", "result_direction": "functionally_abnormal", "key_findings": ["..."]}}
  ],
  "multiple_assay_types_in_paper": false,
  "single_study_only": true,
  "overall_evidence": {{"direction": "PS3", "strength_hint": "not_clear", "rationale": "..."}}
}}
"""


def build_prompt(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    full_text: str,
) -> str:
    """Build the literature prompt from the shared criterion input types.

    ``clinical_note`` is intentionally accepted at the common boundary even
    though PS3/BS3 is decided from functional literature, not patient data.
    """
    del clinical_note
    gene, hgvsc, hgvsp, equivalents = variant_identity(variant)
    return PROMPT_TEMPLATE.format(
        gene=gene,
        hgvsc=hgvsc,
        hgvsp=hgvsp,
        equivalents=", ".join(equivalents),
        full_text=full_text,
    )


# ============================================================================
# 3. Post-processing: approved-assay cross-check + human-curator hint generation
# ============================================================================

def generate_curator_hints(
    judgment: PS3BS3Judgment,
    gene: str,
    vcep_name: Optional[str],
    criterion: str,  # "PS3" or "BS3"
    pmid: Optional[str] = None,
) -> list[CuratorHint]:
    """
    Does not alter the LLM's judgment itself; generates a separate set of
    hints for a human curator's review. Functions as a checklist that
    encodes the misjudgment patterns discovered in design doc sections 6-11.

    `pmid` is optional (callers working from synthetic/test data may not
    have one) but should be passed whenever available: every hint here talks
    about "this paper"/"this judgment", and without a PMID attached that
    reference is ambiguous the moment the hint is read outside the one
    FinalResult it was generated for (e.g. once it's serialized into a
    VA-Spec EvidenceLine's extensions, or just re-read later in a log).
    """
    hints: list[CuratorHint] = []
    paper_ref = f" (PMID:{pmid})" if pmid else ""

    if judgment.variant_matching.match_status == MatchStatus.UNSUCCESSFUL:
        return hints  # already resolved to not_clear; no hints needed

    # --- Pattern 1: multiple assays coexist in the paper (section 6, the MSH2 p.Gly692Val lesson) ---
    if judgment.multiple_assay_types_in_paper:
        assay_list = "; ".join(
            f"{e.assay_type} ({e.result_direction.value})" for e in judgment.experiments
        ) or "none extracted"
        hints.append(CuratorHint(
            "warning",
            f"This paper{paper_ref} self-reports multiple experimental "
            f"systems: {assay_list}. It is recommended to go back to the "
            "relevant part of the full text and manually confirm which of "
            "these is genuinely the most direct functional evidence for "
            "the overall direction below (there is a past instance of an "
            "indirect results table being cited by mistake).",
        ))

    # --- Pattern 2: judgment based on a single study only (section 9, the PI3K-pathway lesson) ---
    if judgment.single_study_only and judgment.overall_evidence.direction != OverallDirection.NOT_CLEAR:
        hints.append(CuratorHint(
            "caution",
            f"This judgment{paper_ref} is based on a single study only. "
            "For criterion codes on the ClinGen side that require "
            "reproducibility across multiple studies (e.g., PS3_Moderate "
            "or stronger), this evidence alone may not be sufficient in "
            "strength.",
        ))

    # --- Pattern 3: cross-check against the approved-assay list (section 11, the RASopathy lesson) ---
    #
    # NOTE: whether a VCEP has a registered approved-assay spec at all, and
    # whether BS3 itself is inapplicable for that VCEP, are both properties
    # of the VCEP+criterion pair, not of any individual experiment -
    # check_approved_assay() returns NO_SPEC_AVAILABLE / BS3_NOT_APPLICABLE
    # before it ever looks at assay_description (see ps3_bs3_ps4_gate.py),
    # so probing with an empty description here reliably predicts both
    # cases without looking at any real experiment. Checking these two
    # up front - instead of inside the per-experiment loop below - avoids
    # emitting the exact same hint once per extracted experiment (e.g. 4
    # identical "No approved-assay list is registered" copies for a paper
    # with 4 experiments, which was a real bug caught in review).
    if vcep_name is not None:
        vcep_level = check_approved_assay(vcep_name, gene, "", criterion)
        if vcep_level.applicability == AssayApplicability.NO_SPEC_AVAILABLE:
            hints.append(CuratorHint(
                "info",
                f"No approved-assay list is registered for {vcep_name}, so "
                f"institutional applicability could not be checked "
                f"automatically for any of this paper's{paper_ref} "
                "experiments.",
            ))
        elif vcep_level.applicability == AssayApplicability.BS3_NOT_APPLICABLE:
            hints.append(CuratorHint(
                "warning",
                f"BS3 itself is not applicable for {vcep_name} "
                f"({vcep_level.note}). Even if this paper's judgment"
                f"{paper_ref} points toward BS3, it likely cannot be "
                "adopted as such.",
            ))
        else:
            for exp in judgment.experiments:
                check = check_approved_assay(vcep_name, gene, exp.assay_type, criterion)
                if check.applicability == AssayApplicability.APPROVED:
                    hints.append(CuratorHint(
                        "info",
                        f"The assay \"{exp.assay_type}\" matches an approved "
                        f"assay for {vcep_name} ({check.matched_assay.name}). "
                        "No institutional obstacle to applying this criterion "
                        "is apparent on this front.",
                    ))
                elif check.applicability == AssayApplicability.ASSAY_NOT_APPROVED:
                    hints.append(CuratorHint(
                        "warning",
                        f"The assay \"{exp.assay_type}\" did not match any "
                        f"assay approved by {vcep_name}. Even if it is "
                        "biologically compelling, this criterion may not hold "
                        "up institutionally (design doc sections 9-10: this "
                        "pattern was actually confirmed for PIK3R2/PIK3CA).",
                    ))
                elif check.applicability == AssayApplicability.GENE_NOT_LISTED:
                    hints.append(CuratorHint(
                        "caution",
                        f"The assay \"{exp.assay_type}\" is approved for "
                        f"{vcep_name}, but {gene} is not listed among its "
                        f"applicable genes. {check.note}",
                    ))

    # --- Pattern 4: model system limitations (pre-existing design, from Case3-var1) ---
    # De-duplicated (identical caveat text across several experiments in the
    # same paper - e.g. "all assays run in HEK293T cells" - would otherwise
    # repeat verbatim) and attributed to the specific assay it came from,
    # for the same reason as Pattern 3 above: an unattributed "info" hint
    # gives a curator no way to tell which experiment it is about.
    seen_caveats: set[str] = set()
    for exp in judgment.experiments:
        if exp.model_system_caveat and exp.model_system_caveat not in seen_caveats:
            seen_caveats.add(exp.model_system_caveat)
            hints.append(CuratorHint(
                "info",
                f"Model system caveat for \"{exp.assay_type}\": {exp.model_system_caveat}",
            ))

    return hints


# Matches things that look like an actual reported measurement: a decimal
# number (e.g. "1.14"), a percentage (e.g. "86.8%"), a p-value
# (e.g. "p<0.0001", "p=0.05"), or a number immediately followed by a common
# lab unit (uL/mL, min/hr, kDa, uM/mM/nM, or a fold-change "x"). Deliberately
# narrower than "any digit", to avoid matching digits embedded in gene/drug/
# cell-line names (p53, Nutlin-3, A549, HEK293T, CRISPR-Cas9, etc.).
_QUANTITATIVE_VALUE_PATTERN = (
    r"\d+\.\d+"                                  # a decimal number, e.g. 1.14
    r"|\d+(\.\d+)?\s*%"                           # a percentage, e.g. 86.8%
    r"|p\s*[<>=]\s*0?\.\d+"                       # a p-value, e.g. p<0.0001
    r"|\d+(\.\d+)?\s*(u|µ|m)?[lL]\b"              # a volume, e.g. 35 uL, 5 mL
    r"|\d+(\.\d+)?\s*(min|mins|minutes|hr|hrs|hours)\b"  # a duration
    r"|\d+(\.\d+)?\s*(kda|u|µ|m)?[mM]\b"          # a mass/molar concentration, e.g. 10 mM
    r"|\d+(\.\d+)?\s*-?fold\b"                    # a fold-change, e.g. 2-fold
)


def detect_no_quantitative_evidence(judgment: PS3BS3Judgment) -> bool:
    """
    Detects experiments whose key_findings contain no quantitative
    measurement value at all, while a definitive direction was still
    reached. This is a coarse heuristic proxy for "the extraction describes
    the paper's general methodology rather than a specific measurement
    actually reported for the target variant" - genuine primary data (a
    percentage, a p-value, an abundance score, a decimal ratio, a volume, a
    time point, etc.) almost always includes at least one such value
    somewhere across all the experiments; a purely qualitative description
    of what "the screen identifies" or "the assay measures" in the abstract
    does not.

    Matching is restricted to patterns that resemble an actual measured
    value (decimals, percentages, p-values, or a number followed by a unit)
    rather than "any digit anywhere". An earlier version of this check used
    a bare "any digit" test, which turned out to be far too permissive in
    this domain: gene, drug, and cell-line names routinely contain digits
    (e.g., "p53", "Nutlin-3", "A549", "HEK293T"), so almost any biology text
    would satisfy it regardless of whether an actual quantitative
    measurement was present.

    Real-world instance (2026-09-15, gemma-4, English prompt): for TP53
    c.875A>G, both extracted "experiments" described the screen's general
    selection logic ("identifies alleles that lack WT-like activity...")
    with no measured value anywhere in either experiment's key_findings
    (only the digit-bearing identifiers "p53" and "Nutlin-3"), yet the LLM
    still confidently concluded PS3. By contrast, in the same run, MYH7's
    and PTEN's genuinely variant-specific extractions did each include at
    least one concrete measured value (respectively, balloon volumes "31
    and 35 uL" and an abundance score of "1.14").

    This heuristic is intentionally coarse (a single missing value is not
    conclusive proof of over-generalization) and may need refinement as more
    real cases are gathered; for now it is used, like the other detectors in
    this module, as a trigger for a cautious override to not_clear rather
    than as a claim of certainty.
    """
    if not judgment.experiments:
        return False
    all_key_findings_text = " ".join(
        finding for exp in judgment.experiments for finding in exp.key_findings
    )
    return not bool(re.search(_QUANTITATIVE_VALUE_PATTERN, all_key_findings_text, re.IGNORECASE))


def detect_experiment_conflict(judgment: PS3BS3Judgment) -> bool:
    """
    Detects whether both functionally_abnormal and functionally_normal are
    present among the experiment results (i.e., assays within the same paper
    disagree on direction).

    Originally written (2026-09-14) against a PTEN c.112C>T (p.Pro38Ser)
    gemma-4 run where VAMP-seq (protein abundance, functionally_normal) and
    the Akt-activation assay (functionally_abnormal) disagreed. The ground
    truth for PS3 on this exact variant is actually MET (moderate) - see
    test_data.full_criteria_ground_truth - not not_met as this function's
    own history originally assumed, so the case this was built around is
    itself a real, correct MET that a blanket override would suppress. See
    _conflict_reasoned_in_rationale() for the distinction this now draws
    between a paper the LLM read carefully (and reasoned through the
    disagreement) and one it skimmed past.
    """
    directions = {e.result_direction for e in judgment.experiments}
    return (
        ResultDirection.FUNCTIONALLY_ABNORMAL in directions
        and ResultDirection.FUNCTIONALLY_NORMAL in directions
    )


_CONTRASTIVE_MARKERS = (
    "although", "despite", "even though", "whereas", "in contrast",
    "however", "nonetheless", "nevertheless", "but the", "while the",
)


def _conflict_reasoned_in_rationale(judgment: PS3BS3Judgment) -> bool:
    """
    Whether the free-text rationale shows the LLM actually engaged with the
    disagreeing experiments, rather than silently picking one and ignoring
    the other - the real failure mode detect_experiment_conflict was built
    to catch (2026-09-14: the LLM "adopted only" the abnormal result with no
    acknowledgement of the normal one at all).

    A rationale that names a contrastive marker (a paper's own text
    weighing one finding against another - "although VAMP-seq showed
    WT-like abundance... the Akt-activation assay... supports PS3") is
    reasoning through the conflict, not ignoring it. Confirmed against the
    real PTEN c.112C>T case (ground truth PS3=MET) this function was named
    for: its rationale opens with exactly this pattern ("Although VAMP-seq
    showed p.Pro38Ser has WT-like/enhanced protein abundance ... a direct
    functional assay demonstrated ... even though the abundance assay
    alone would suggest a benign/normal result").

    This is a coarse text heuristic, same style/limitations as
    detect_no_quantitative_evidence's pattern matching - it cannot verify
    the reasoning is actually sound, only that the model did not silently
    drop the discordant experiment. detect_experiment_conflict still fires
    and still produces the same warning hint either way; this only decides
    whether the direction is additionally forced to not_clear.
    """
    text = (judgment.overall_evidence.rationale or "").lower()
    if not text:
        return False
    return any(marker in text for marker in _CONTRASTIVE_MARKERS)


# Added 2026-09-15: detect_experiment_conflict above turned out to be
# insufficient on its own. In a second real gemma-4 run (PTEN c.112C>T
# again), the rationale explicitly stated "showing a WT-like abundance score
# (1.14)", yet the result_direction for that same experiment was (along with
# the other experiment) categorized as functionally_abnormal. In other words,
# the LLM retrofits its categorical labels to match its own final
# conclusion, so the conflict never surfaces as a disagreement between
# categorical fields. A separate check is needed to catch the mismatch
# between free-text wording that implies a normal/WT-like finding and a set
# of experiments that are all labeled functionally_abnormal.
_NORMAL_FINDING_PATTERNS = [
    r"wt-?like", r"wild-?type-?like", r"similar to (the )?wild-?type",
    r"not (significantly )?different from (the )?wild-?type",
    r"no significant difference", r"comparable to (the )?wild-?type",
    r"within (the )?normal range",
]

# Words that, if found immediately before a "normal finding" phrase, negate
# its meaning (e.g. "lack WT-like activity" means the OPPOSITE of a normal
# finding - it describes a damaging/loss-of-function allele). Without this
# check, such phrasing is indistinguishable from a genuine normal finding by
# simple substring matching alone.
#
# Real-world instance (2026-09-15, gemma-4, English prompt): for TP53
# c.875A>G, one experiment's key_findings read "The screen identifies
# alleles that lack WT-like activity (loss-of-function)...". This matched
# the "wt-?like" pattern above, causing detect_narrative_categorical_mismatch
# to fire - but for the wrong reason: it happened to catch this specific
# misjudgment, but only by coincidence, not because it correctly understood
# the sentence. The phrase describes damaging alleles, not a normal finding.
_NEGATION_WINDOW_CHARS = 25
_NEGATION_WORDS = ["lack", "lacking", "without", "absence of", "fails to show", "failed to show", "not "]


def _has_uncontradicted_normal_finding(text_lower: str) -> bool:
    """
    True if at least one _NORMAL_FINDING_PATTERNS match exists in text_lower
    that is NOT immediately preceded by a negation word (within
    _NEGATION_WINDOW_CHARS characters). A negated match (e.g., "lack
    WT-like activity") is excluded, since it indicates the opposite of a
    normal finding.
    """
    for pattern in _NORMAL_FINDING_PATTERNS:
        for m in re.finditer(pattern, text_lower):
            window_start = max(0, m.start() - _NEGATION_WINDOW_CHARS)
            preceding = text_lower[window_start:m.start()]
            if not any(neg in preceding for neg in _NEGATION_WORDS):
                return True
    return False


def detect_narrative_categorical_mismatch(judgment: PS3BS3Judgment) -> bool:
    """
    Detects the case where the free-text rationale (and each experiment's
    key_findings) implies a normal/WT-like finding, yet every experiment is
    categorically labeled functionally_abnormal in the structured fields
    (i.e., the conflict never made it into the categorical fields).

    Real-world instance: for PTEN c.112C>T, the rationale explicitly said
    "WT-like abundance score (1.14)", but the VAMP-seq experiment's
    result_direction was functionally_abnormal.
    """
    if not judgment.experiments:
        return False
    text = (judgment.overall_evidence.rationale or "")
    for exp in judgment.experiments:
        text += " " + " ".join(exp.key_findings)
    text_lower = text.lower()
    mentions_normal_finding = _has_uncontradicted_normal_finding(text_lower)
    all_abnormal = all(
        e.result_direction == ResultDirection.FUNCTIONALLY_ABNORMAL for e in judgment.experiments
    )
    return mentions_normal_finding and all_abnormal


def detect_empty_experiments_with_definitive_direction(judgment: PS3BS3Judgment) -> bool:
    """
    Detects the case where match_status is matched/heuristic/single_variant_study
    (i.e., the LLM claims the paper does test the target variant) yet zero
    experiments were extracted, while the overall direction is still a
    definitive PS3 or BS3 (not not_clear). This is internally inconsistent:
    per the prompt's own instructions, a definitive direction should be
    backed by at least one concretely extracted experiment.

    Real-world instance (2026-09-15, gemma-4 run, English prompt): for TP53
    c.875A>G (p.Lys292Arg), the rationale explicitly admitted "the specific
    Z-score for K292R is in the supplementary data rather than the main
    text", yet the LLM still confidently concluded PS3 by reasoning from the
    paper's general methodology rather than variant-specific data. Ground
    truth was BS3 = met. This function catches the empty-experiments symptom
    of that failure mode (it cannot detect the underlying over-generalization
    directly, but an empty experiment list combined with a definitive
    direction is a reliable proxy for it).
    """
    if judgment.variant_matching.match_status == MatchStatus.UNSUCCESSFUL:
        return False  # already resolves to not_clear on its own
    return (
        len(judgment.experiments) == 0
        and judgment.overall_evidence.direction != OverallDirection.NOT_CLEAR
    )


def finalize(
    judgment: PS3BS3Judgment,
    gene: str,
    vcep_name: Optional[str],
    criterion: str,
    pmid: Optional[str] = None,
) -> FinalResult:
    hints = generate_curator_hints(judgment, gene, vcep_name, criterion, pmid=pmid)
    effective_direction = judgment.overall_evidence.direction
    paper_ref = f" (PMID:{pmid})" if pmid else ""

    if detect_experiment_conflict(judgment) and effective_direction != OverallDirection.NOT_CLEAR:
        conflicting = "; ".join(
            f"{e.assay_type} = {e.result_direction.value}" for e in judgment.experiments
        )
        reasoned = _conflict_reasoned_in_rationale(judgment)
        if reasoned:
            hints.append(CuratorHint(
                "caution",
                f"The experiments in this paper{paper_ref} disagree on "
                f"direction ({conflicting}). The rationale shows the LLM "
                f"weighed the disagreement rather than ignoring it, so the "
                f"judgment ({effective_direction.value}) was kept, but a "
                "human curator should still confirm which assay is "
                "genuinely the most direct evidence.",
            ))
        else:
            hints.append(CuratorHint(
                "warning",
                f"The experiments in this paper{paper_ref} disagree on "
                f"direction ({conflicting}), yet the LLM confidently concluded "
                f"{effective_direction.value} without acknowledging the "
                "disagreement in its own rationale. This conflict was "
                "detected automatically and the judgment has been forced "
                "to not_clear.",
            ))
            effective_direction = OverallDirection.NOT_CLEAR
    elif detect_narrative_categorical_mismatch(judgment) and effective_direction != OverallDirection.NOT_CLEAR:
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, the free-text rationale contains "
            f"wording suggesting a normal/WT-like finding, yet every "
            f"experiment is categorized as functionally_abnormal, and the "
            f"LLM confidently concluded {effective_direction.value}. The "
            "categorical labels were likely retrofitted to match the final "
            "conclusion, so the judgment has been forced to not_clear "
            "(real-world instance: PTEN c.112C>T, second run).",
        ))
        effective_direction = OverallDirection.NOT_CLEAR
    elif detect_empty_experiments_with_definitive_direction(judgment) and effective_direction != OverallDirection.NOT_CLEAR:
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, match_status was reported as "
            f"{judgment.variant_matching.match_status.value}, but zero "
            f"experiments were extracted, and the LLM still confidently "
            f"concluded {effective_direction.value}. A definitive "
            "direction with no concrete supporting experiment is "
            "internally inconsistent, so the judgment has been forced to "
            "not_clear (real-world instance: TP53 c.875A>G, where the "
            "rationale admitted the variant-specific data was only in "
            "supplementary material, yet still asserted a direction by "
            "generalizing from the paper's overall methodology).",
        ))
        effective_direction = OverallDirection.NOT_CLEAR
    elif detect_no_quantitative_evidence(judgment) and effective_direction != OverallDirection.NOT_CLEAR:
        hints.append(CuratorHint(
            "warning",
            f"In this paper{paper_ref}, none of the extracted experiments' "
            f"key_findings contain any numeric detail, yet the LLM "
            f"confidently concluded {effective_direction.value}. This "
            "suggests the extraction may describe the paper's general "
            "methodology rather than a measurement specifically reported "
            "for the target variant, so the judgment has been forced to "
            "not_clear (real-world instance: TP53 c.875A>G, where the "
            "paper is a large-scale saturation screen and the target "
            "variant's specific score was only in supplementary data).",
        ))
        effective_direction = OverallDirection.NOT_CLEAR

    return FinalResult(judgment=judgment, curator_hints=hints, effective_direction=effective_direction)


# ============================================================================
# 4. Multi-paper aggregation
# ============================================================================
#
# Added 2026-09-15, in response to a design gap noticed by the user: every
# function above judges ONE paper (one PMID) at a time, but a real ClinGen
# VCEP classification typically cites several papers as evidence for a
# single criterion (a real example seen during design doc section 9: RUNX1
# c.601C>T's ERepo record cited 4 PMIDs in its evidenceLinks). The
# single_study_only field on PS3BS3Judgment only asks the LLM to
# self-report whether it BELIEVES other studies exist, from citations
# mentioned in the one paper it was given - it never actually goes and reads
# those other papers. This section adds a layer that combines the
# per-paper FinalResult objects for multiple actual papers into one
# variant-level conclusion, so "how many papers really support this" is
# computed directly rather than self-reported by the model.

def aggregate_multi_paper_results(contributions: list[PaperContribution]) -> AggregatedJudgment:
    """
    Thin, backward-compatible wrapper around the generic implementation in
    evidence_common.py (extracted there 2026-09-15 when PS4 and PP1/BS4 were
    added, so this logic is not duplicated per criterion module) - just
    supplies this module's own OverallDirection.NOT_CLEAR as the sentinel.
    See evidence_common.aggregate_multi_paper_results for the actual logic
    and docstring.
    """
    return _generic_aggregate_multi_paper_results(contributions, not_clear=OverallDirection.NOT_CLEAR)
