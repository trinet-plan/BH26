"""
va_spec_export.py
Converts this project's judgment output - PS3/BS3 (ps3_bs3_judgment.py),
PS4 (ps4_judgment.py), or PP1/BS4 (segregation_judgment.py) - into a GA4GH
VA-Spec EvidenceLine (as a plain dict, ready for json.dump()).

[v2, 2026-09-15: now built on the real ga4gh.va_spec Pydantic library]
  The first version of this module (see git history / design doc section
  15-2) hand-built plain dicts shaped to match the VA-Spec docs and the
  sample fixture (va_spec_example.yaml), since neither `ga4gh.va_spec` nor
  any other VA-Spec schema library was installed or imported. That worked,
  but was never actually validated against the real schema - it could have
  had undetected drift (wrong field shapes, wrong enum values) anywhere.
  This version imports and constructs the real Pydantic models from
  `ga4gh.va_spec` (`pip install ga4gh.va_spec`) instead, so every
  EvidenceLine this module produces is guaranteed to satisfy the actual
  schema (Pydantic validation raises immediately on construction if not).
  Two real discrepancies were caught by switching to the real models:
    - `Contribution.activityType` is a plain `str`, not a MappableConcept
      (the old code nested it as `{"name": "..."}`, which was wrong).
    - `Contribution.contributor` is an `Agent` model, not an arbitrary dict.
  Everything else (Direction enum values "supports"/"neutral"/"disputes",
  the "ACMG Guidelines, 2015" System string, the Document/Method/
  MappableConcept/Coding/Extension field names) matched what had already
  been hand-built, confirming the original reading of the docs was
  otherwise accurate.

[Criterion-agnostic by design]
  This module only imports the generic shapes from evidence_common.py
  (MatchStatus, CuratorHint, PaperContribution, AggregatedJudgment) - never
  a criterion-specific direction Enum (OverallDirection, or whatever PS4 /
  PP1&BS4 define). "Is this the neutral/inconclusive case" is checked as
  `direction.value == "not_clear"` rather than `direction == SomeEnum.
  NOT_CLEAR`, since every direction Enum in this project uses that literal
  string for its neutral member (see evidence_common.py's own docstring).
  This is what lets build_evidence_line() work unmodified for any of the
  three judgment modules above, and any future one that follows the same
  convention.

[Scope]
  Builds one EvidenceLine per (variant, criterion) - aggregated across every
  cited paper examined for that criterion - not a full VariantPathogenicity
  Statement/Proposition for the whole variant. That wrapping is a separate,
  later concern; the immediate goal agreed with the user is the EvidenceLine
  itself.

[Design decisions, finalized in conversation with the user (2026-09-15),
 superseding the informal v3-era mapping in design doc section 5]
  - Granularity: one top-level EvidenceLine per (variant, criterion),
    corresponding to this pipeline's AggregatedJudgment. Each paper actually
    examined (PaperContribution) becomes a nested EvidenceLine under
    hasEvidenceItems, so the aggregated conclusion and the individual
    per-paper evidence both survive in the same document (the real schema's
    `hasEvidenceItems` explicitly allows nesting another EvidenceLine, not
    just StudyResult/Statement - confirmed via EvidenceLine.model_fields).
  - not_clear -> still emit an EvidenceLine (do not drop it / route it to a
    separate "Evidence Gap" queue as design doc v3 proposed). Use
    directionOfEvidenceProvided = Direction.NEUTRAL, one of exactly three
    real enum members (supports/neutral/disputes - confirmed via
    ga4gh.va_spec.base.enums.Direction) for exactly this case. reportedIn
    (the papers actually examined) and description (why no confident
    direction could be reached) are still populated; strength fields are
    omitted since no strength claim is being made.
  - Curator hints (structured severity+message, no dedicated VA-Spec field)
    go into `extensions`. A short free-text summary goes into the inherited
    `description` field (InformationEntity.description).
  - Paper references (PMID/URL) go into `reportedIn` (Document).
  - Strength: the LLM never actually judges strength (see
    ps3_bs3_judgment.OverallEvidence.strength_hint, which the prompt fixes to
    the literal string "not_clear" - see PROMPT_TEMPLATE - so it carries no
    real signal and is deliberately NOT used here). Instead,
    strengthOfEvidenceProvided/evidenceOutcome are populated from a rough,
    explicitly-flagged heuristic: the count of independently agreeing papers
    that aggregate_multi_paper_results() already computes (1->Supporting,
    2->Moderate, 3+->Strong - all three are real values from
    ga4gh.va_spec.base.enums.STRENGTH_OF_EVIDENCE_PROVIDED_VALUES). This is
    NOT a VCEP-calibrated ACMG strength determination - every EvidenceLine
    that carries a strength value also carries a `strengthEstimationMethod`
    extension saying so, so a human curator does not mistake it for an
    authoritative call.
  - structuredEvidenceItems: a checklist-style extension (see
    _structured_items_extension below) added 2026-09-15 in response to a
    real curator UI requirements doc (doc/recs for expert board.docx).
"""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
from functools import lru_cache
from typing import Any, Optional

from jsonschema import Draft202012Validator, FormatChecker

from acmg_pipeline.criteria.common import DEFAULT_STRENGTH
from acmg_pipeline.automated_va_spec import OUTCOME_PATTERN, extensions_last, output_schema

from ga4gh.core.models import Coding, Extension, MappableConcept
from ga4gh.va_spec.base.core import Agent, Contribution, Direction, Document, EvidenceLine, Method
from ga4gh.va_spec.base.enums import System

from acmg_pipeline.common import (
    AggregatedJudgment, CuratorHint, MatchStatus, PaperContribution,
    is_not_clear, strength_tier_from_paper_count,
)
from acmg_pipeline.constants import IMPLEMENTED_CODES, PATHOGENIC_CODES, CriterionStatus
from acmg_pipeline.criteria import curator_info, reference_links
from acmg_pipeline.vcf_record import VariantRecord

# Codes whose EvidenceLine also gets a `curatorInfo` extension (see
# build_reference_extensions()), and which acmg_pipeline.criteria.curator_info
# function supplies it - the doc's ask for these goes beyond "show a page"
# (reference_links.py) into "check hotspot/nearby benign variants" or
# "check if previously reported", which needs the page's actual content
# fetched, not just a link to it.
#
# PP1 is deliberately absent here: its requested ClinVar page is attached
# as a referenceLink to the real literature EvidenceLine, while curatorInfo
# remains limited to the previously implemented PM1/PM3/PM5 fetchers.
_CODE_CURATOR_INFO_FETCHERS = {
    "PM1": curator_info.uniprot_domain_context,
    "PM5": curator_info.uniprot_domain_context,
    "PM3": curator_info.clinvar_report_context,
}

# The ACMG/AMP 2015 guideline itself, cited as the source of every
# criterion's Method. evidence-cli already attached this to its scored
# lines, and the 1.0.1 output schema requires `specifiedBy.reportedIn`,
# so the literature and workflow lines built here carry it too - without
# it they were the only lines in the document whose Method named no
# source at all.
ACMG_2015_METHOD_DOCUMENT = Document(
    name="Richards et al., 2015, Genet Med.",
    doi="10.1038/gim.2015.30",
    pmid="25741868",
)

PIPELINE_AGENT = Agent(
    id="acmg-literature-llm-pipeline",
    name=(
        "ACMG literature-evidence LLM pipeline (gemma-4 via vLLM) - "
        "AI-generated draft evidence, requires human curator review before use"
    ),
)


def default_strength(code: str) -> str:
    """The ACMG/AMP 2015 default strength of one criterion code.

    evidence-cli's own DEFAULT_STRENGTH map covers only the 16 codes it
    implements, so it cannot answer for PP1/BS4/PS3/BS3/PS4. The ACMG
    prefix rule covers all 28 and agrees with that map everywhere the two
    overlap (asserted in tests), so it is derived rather than duplicated.
    """
    if code == "PVS1":
        return "very_strong"
    if code == "BA1":
        return "stand_alone"
    return {"PS": "strong", "BS": "strong", "PM": "moderate",
            "PP": "supporting", "BP": "supporting"}[code[:2]]


def evidence_line_id(code: str, gene: str, safe_hgvsc: str) -> str:
    """The single EvidenceLine id scheme for the integrated 28-line document.

    Every builder in this module goes through here, including the adapter
    for evidence-cli results (build_automated_evidence_line), which would
    otherwise carry that project's own `urn:bh26:evidence-line:<sha256>`
    form. Two ids for one document was a real problem, not a cosmetic one:

      * the split ran along MET/NOT_MET vs. every other status, not along
        which half produced the line, so the SAME criterion changed id
        FORM depending on its outcome - a consumer could not look a line
        up by (variant, criterion);
      * the urn hashed `evidence_outcome` into the id, so re-running with
        a different threshold (PP3_strong -> PP3_moderate) changed the id
        of the same variant/criterion, breaking run-to-run diffing and
        anything anchored to a line (curator comments, review state).

    `evline:{gene}_{safe_hgvsc}_{code}` is stable across runs and
    reconstructible from (variant, criterion) alone. It stays a valid
    JSON-Schema `iri-reference`, which is all the 1.0.1 output schema asks
    of `id`.

    The urn scheme in acmg_pipeline/automated_va_spec.py (stable_urn(), used
    by to_evidence_line()) is deliberately left alone: the standalone
    `acmg evaluate` CLI has no gene/HGVSC to build this form from (only
    assembly:chrom:pos:ref:alt), and its own output contract and tests are
    pinned to the urn.
    """
    return f"evline:{gene}_{safe_hgvsc}_{code}"


@lru_cache(maxsize=1)
def _integrated_line_schema() -> dict:
    """The 1.0.1 output schema, widened to the fields the integrated document uses.

    `acmg_pipeline/schemas/acmg-evidence-line-1.0.1-output.json` describes
    exactly what evidence-cli emits for a scored criterion, so applying it
    verbatim to all 28 lines is impossible in two ways:

      * it is `additionalProperties: false` and lists neither `reportedIn`
        nor `contributions`, both of which the literature lines carry
        (PMIDs of the papers judged, and the LLM agent's contribution);
      * it requires `evidenceOutcome`, which a workflow line must not have -
        emitting one would assert a judgment that was never made.  Such a
        line is any result that is neither MET nor NOT_MET; see
        build_automated_evidence_line(), which routes those to
        build_workflow_evidence_line() instead.

    [Vocabulary note - corrected 2026-09-17]
      Before the evidence-cli merge (c5b303f) this said the workflow states
      were NOT_EVALUATED / MANUAL_REVIEW / NOT_APPLICABLE, which was the
      six-value `acmg.core.models.Status` enum of the standalone package.
      That enum no longer exists: the merge folded it into the three-value
      `acmg_pipeline.constants.CriterionStatus` (met / not_met / unknown),
      so a workflow line now always carries status "unknown" and those
      three names are never emitted.  VA-Spec documents under
      va_spec_output/ that still show them predate the merge.  A criterion
      that is inapplicable rather than merely undecided is marked inside
      provenance instead - see acmg_pipeline.criteria.common.NOT_APPLICABLE.

    So the schema file itself is NOT edited - its sha256 is recorded as
    provenance by acmg_pipeline.automated_va_spec.output_schema_sha256() and travels
    inside audit envelopes; changing it would silently invalidate those.
    Instead this derives an in-memory superset: the same property schemas,
    plus the two literature-only fields, with `evidenceOutcome` demoted
    from required (it is still cross-checked by _check_acmg_semantics()
    whenever it IS present).
    """
    schema = deepcopy(output_schema())
    schema["required"] = [key for key in schema["required"] if key != "evidenceOutcome"]
    schema["properties"]["reportedIn"] = {"type": "array"}
    schema["properties"]["contributions"] = {"type": "array"}
    # `hasEvidenceItems` in 1.0.1 accepts an IRI string or a cohort allele
    # frequency StudyResult - evidence-cli's only two shapes. The
    # literature lines nest one full EvidenceLine per paper judged (see
    # build_paper_evidence_line), which the real ga4gh.va_spec model
    # explicitly allows and which is where the per-paper rationale and
    # curator hints live. Accept that third shape here rather than
    # flatten it away; the nested objects are already Pydantic-validated
    # as EvidenceLine at construction time.
    schema["properties"]["hasEvidenceItems"]["items"]["oneOf"].append(
        {"type": "object", "properties": {"type": {"const": "EvidenceLine"}},
         "required": ["type"]}
    )
    # `acmgStrength` / `acmgOutcome` are MappableConcepts. evidence-cli
    # hand-builds them without the `type` discriminator; the real
    # ga4gh.core.models.MappableConcept that the literature lines are dumped
    # from emits `"type": "MappableConcept"`. Both are valid VA-Spec, so the
    # discriminator is permitted rather than stripped from the Pydantic
    # output. NOTE: this leaves the two halves emitting slightly different
    # shapes for the same two fields - worth settling one way or the other,
    # but not by deleting a field the upstream model considers correct.
    for name in ("acmgStrength", "acmgOutcome"):
        schema["$defs"][name].setdefault("properties", {})["type"] = {
            "const": "MappableConcept"
        }
    return schema


# PS3/BS3 share one literature judgment (acmg_pipeline.criteria.ps3_bs3):
# the LLM decides the direction actually ESTABLISHED by the evidence (PS3,
# BS3, or not_clear) independent of which of the pair is being tested right
# now - _strength_blocks()'s own docstring is explicit that "testing PS3 but
# finding BS3-direction evidence must be coded BS3_moderate, not
# PS3_moderate". So a "BS3" line whose evidenceOutcome code starts with
# "PS3" (or vice versa) is not a bug, it is the correct way to record "this
# disputes BS3 - the evidence actually shows PS3". Found 2026-09-17 when a
# real demo case (DSG2 c.1592T>G) had clear PS3-direction literature
# evidence and building its BS3 line crashed here.
_ACMG_DISPUTE_SIBLINGS = {"PS3": "BS3", "BS3": "PS3"}


def _check_acmg_semantics(line: dict, criterion: str) -> None:
    """The cross-field ACMG rules from validate_1_0_1(), for any line that scores.

    Same three checks, in the same order, so a line that has an
    evidenceOutcome is held to exactly the standard evidence-cli's scored
    lines already were - previously these applied only to the 16 automated
    codes' MET/NOT_MET lines, leaving the literature codes free to emit
    e.g. direction=supports alongside a *_not_met outcome.
    """
    outcome = line["evidenceOutcome"]["primaryCoding"]["code"]
    direction = line["directionOfEvidenceProvided"]
    if line["specifiedBy"]["methodType"] != criterion or not OUTCOME_PATTERN.fullmatch(outcome):
        raise ValueError(f"VA-Spec ACMG criterion mapping mismatch for {criterion}")
    outcome_code = outcome.split("_", 1)[0]
    sibling = _ACMG_DISPUTE_SIBLINGS.get(criterion)
    if outcome_code != criterion and outcome_code != sibling:
        raise ValueError(f"VA-Spec methodType/evidenceOutcome mismatch for {criterion}")
    if outcome_code == sibling:
        expected = "disputes"
    else:
        expected = "neutral" if outcome.endswith("_not_met") else (
            "disputes" if criterion.startswith("B") else "supports"
        )
    if direction != expected:
        raise ValueError(
            f"VA-Spec direction/evidenceOutcome mismatch for {criterion}: "
            f"{direction!r} with outcome {outcome!r} (expected {expected!r})"
        )


def validate_integrated_line(line: dict, criterion: str) -> dict:
    """Validate one line of the integrated 28-line document. Applied to ALL 28.

    Scored lines (those carrying an evidenceOutcome) get the same schema
    and ACMG semantic checks that acmg_pipeline.automated_va_spec.validate_1_0_1()
    applies to evidence-cli's output; workflow lines get the schema checks
    minus the evidenceOutcome requirement. Raises ValueError on the first
    problem, matching validate_1_0_1()'s behaviour.
    """
    errors = sorted(
        Draft202012Validator(_integrated_line_schema(), format_checker=FormatChecker())
        .iter_errors(line),
        key=lambda item: list(item.path),
    )
    if errors:
        raise ValueError(
            f"VA-Spec 1.0.1 schema validation failed for {criterion}: {errors[0].message}"
        )
    if line.get("evidenceOutcome"):
        _check_acmg_semantics(line, criterion)
    return extensions_last(line)


def _pubmed_url(pmid: str) -> str:
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def _document(pmid: str) -> Document:
    return Document(pmid=pmid, urls=[_pubmed_url(pmid)])


def _direction_of_evidence(direction, criterion: str) -> Direction:
    """direction/disputes here means "for/against pathogenicity", not "for/against
    the criterion under test" - so a BS3 line whose evidence actually establishes
    BS3 (i.e. direction.value == criterion == "BS3") DISPUTES pathogenicity, it
    does not support it. Only a match on a *pathogenic* criterion (PS3, PS4)
    supports pathogenicity; a mismatch (the dispute-sibling case, e.g. PS3
    evidence found while testing BS3) always disputes the criterion under test,
    regardless of that criterion's own P/B prefix. Found 2026-09-18: this always
    returned SUPPORTS on a match, so a real BS3_supporting outcome (MYBPC3-class
    literature, LDLR c.2575G>A) raised acmg_pipeline.export's own
    _check_acmg_semantics() ValueError the first time it was exercised live -
    mirrors the same convention acmg_pipeline.criteria.pp1_bs4_pp4_engine's
    build_evidence_line() already uses via PATHOGENIC_CODES.
    """
    if is_not_clear(direction):
        return Direction.NEUTRAL
    if direction.value != criterion:
        return Direction.DISPUTES
    return Direction.SUPPORTS if criterion in PATHOGENIC_CODES else Direction.DISPUTES


def _hints_extension(hints: list[CuratorHint]) -> Optional[Extension]:
    if not hints:
        return None
    return Extension(
        name="curatorHints",
        value=[{"severity": h.severity, "message": h.message} for h in hints],
    )


def _structured_items_extension(judgment) -> Optional[Extension]:
    """
    Pulls the criterion-specific checklist items off `judgment` via duck
    typing - `structured_evidence_items()` is an OPTIONAL method each
    criterion module's Judgment class may implement (PS3BS3Judgment,
    PS4Judgment, SegregationJudgment all do); this module never imports any
    of those classes, so a getattr probe is how it stays criterion-agnostic
    the same way common.is_not_clear() does for direction.

    Each item is {"label": str, "checked": bool, "detail": str} - see
    design doc section 15-4-1 for the full convention.
    """
    getter = getattr(judgment, "structured_evidence_items", None)
    if not callable(getter):
        return None
    items = getter()
    if not items:
        return None
    return Extension(name="structuredEvidenceItems", value=items)


def _strength_blocks(
    direction, n_agreeing_papers: int,
) -> tuple[Optional[MappableConcept], Optional[MappableConcept], Optional[Extension]]:
    """
    Returns (strengthOfEvidenceProvided, evidenceOutcome, heuristic-disclosure
    extension), any of which may be None.

    The ACMG code is built from the direction actually established
    (`direction.value`, e.g. "PS3" or "BS3"), NOT from whichever criterion
    the caller happened to be testing for. Using the tested criterion
    unconditionally would mislabel a disputing result - e.g. testing PS3 but
    finding BS3-direction evidence must be coded "BS3_moderate: Criterion
    Met", not "PS3_moderate: Criterion Met" (which would claim the opposite
    of what was found).
    """
    if is_not_clear(direction):
        return None, None, None
    tier = strength_tier_from_paper_count(n_agreeing_papers)
    if tier is None:
        return None, None, None
    code_base = direction.value
    strength = MappableConcept(primaryCoding=Coding(code=tier, system=System.ACMG.value))
    # An ACMG outcome code carries a strength suffix only when the strength
    # DIFFERS from the criterion's own default - "PP1_supporting" is not a
    # legal code, because supporting is what plain PP1 already means. This
    # is the same rule evidence-cli applies in
    # acmg_pipeline/criteria/common.py result(); before the merge nothing
    # checked it on this side, and
    # PP1/BS4 at the supporting tier emitted the illegal form.
    outcome_code = code_base if tier == default_strength(code_base) else f"{code_base}_{tier}"
    outcome = MappableConcept(
        primaryCoding=Coding(code=outcome_code, system=System.ACMG.value),
        name=f"ACMG 2015 {code_base} {tier} Criterion Met (heuristic estimate, unconfirmed)",
    )
    disclosure = Extension(
        name="strengthEstimationMethod",
        value=(
            f"heuristic: {n_agreeing_papers} independent paper(s) with agreeing "
            "direction (1=supporting, 2=moderate, 3+=strong). This is a rough "
            "proxy, NOT a VCEP-calibrated ACMG strength determination - the "
            "underlying LLM judgment only assesses direction, never strength "
            "- and requires human curator confirmation."
        ),
    )
    return strength, outcome, disclosure


def build_paper_evidence_line(contribution: PaperContribution, criterion: str, index: int) -> EvidenceLine:
    """Builds the nested, per-paper EvidenceLine for one PaperContribution."""
    result = contribution.result
    hints_ext = _hints_extension(result.curator_hints)
    items_ext = _structured_items_extension(result.judgment)
    extensions = [e for e in (hints_ext, items_ext) if e] or None

    return EvidenceLine(
        id=f"evline:paper-{index}-{contribution.pmid}",
        directionOfEvidenceProvided=_direction_of_evidence(result.effective_direction, criterion),
        reportedIn=[_document(contribution.pmid)],
        description=result.judgment.overall_evidence.rationale or None,
        extensions=extensions,
    )


def build_evidence_line(
    aggregated: AggregatedJudgment,
    gene: str,
    hgvsc: str,
    criterion: str,
    vcep_name: Optional[str] = None,
    variant: Optional[VariantRecord] = None,
) -> dict:
    """
    Builds the top-level, aggregated VA-Spec EvidenceLine for one
    (variant, criterion) pair from this pipeline's AggregatedJudgment, and
    returns it as a JSON-ready dict (`.model_dump(mode="json")` on a real,
    schema-validated ga4gh.va_spec.base.core.EvidenceLine instance - not a
    hand-built dict).
    """
    direction = aggregated.aggregated_direction
    all_pmids = [c.pmid for c in aggregated.contributions]
    relevant_pmids = [
        c.pmid for c in aggregated.contributions
        if c.result.judgment.variant_matching.match_status != MatchStatus.UNSUCCESSFUL
        and not is_not_clear(c.result.effective_direction)
    ]

    strength, outcome, strength_disclosure = _strength_blocks(direction, len(relevant_pmids))
    hints_ext = _hints_extension(aggregated.aggregation_hints)
    status = (
        CriterionStatus.UNKNOWN.value if is_not_clear(direction)
        else CriterionStatus.MET.value if direction.value == criterion
        else CriterionStatus.NOT_MET.value
    )
    assessment_ext = Extension(
        name="bh26AssessmentDetails",
        value={
            "criterion": criterion,
            "status": status,
            "summary": " ".join(h.message for h in aggregated.aggregation_hints)
            or f"Literature evidence was evaluated for {criterion}.",
        },
    )
    extensions = [e for e in (hints_ext, strength_disclosure, assessment_ext) if e]
    if variant is not None:
        extensions.extend(build_reference_extensions(criterion, variant))

    description = " ".join(h.message for h in aggregated.aggregation_hints) or None

    safe_hgvsc = (
        hgvsc.replace(">", "_").replace(".", "_").replace("+", "p").replace("-", "m")
    )

    evidence_line = EvidenceLine(
        id=evidence_line_id(criterion, gene, safe_hgvsc),
        directionOfEvidenceProvided=_direction_of_evidence(direction, criterion),
        reportedIn=[_document(pmid) for pmid in all_pmids] or None,
        description=description,
        hasEvidenceItems=[
            build_paper_evidence_line(c, criterion, i)
            for i, c in enumerate(aggregated.contributions, start=1)
        ] or None,
        strengthOfEvidenceProvided=strength,
        evidenceOutcome=outcome,
        specifiedBy=Method(
            methodType=criterion,
            name=(
                f"{vcep_name or 'Unspecified VCEP'} ACMG/AMP {criterion} "
                "evidence assessment (LLM-assisted draft)"
            ),
            reportedIn=ACMG_2015_METHOD_DOCUMENT,
        ),
        contributions=[Contribution(
            contributor=PIPELINE_AGENT,
            activityType="automated evidence evaluation",
            date=datetime.now(timezone.utc),
        )],
        extensions=extensions or None,
    )
    return validate_integrated_line(
        evidence_line.model_dump(mode="json", exclude_none=True), criterion
    )


def _serialize_curator_info(ctx) -> Optional[dict]:
    """
    Plain-dict serialization for whichever acmg_pipeline.criteria.
    curator_info context type `ctx` is (UniprotDomainContext for PM1/PM5,
    ClinVarReportContext for PM3) - not ctx.__dict__ directly, since both
    contain lists of nested dataclasses that Extension.value (must be
    JSON-serializable) can't hold as-is. Returns None if `ctx` itself is
    None (the fetcher found nothing - e.g. no VariationID resolved, or a
    network error - see each fetcher's own docstring for when that
    happens) - same "no data means no extension" convention as
    reference_url_for_criterion() returning None upstream of this.
    """
    if ctx is None:
        return None
    if isinstance(ctx, curator_info.UniprotDomainContext):
        return {
            "uniprotAccession": ctx.accession,
            "proteinPosition": ctx.position,
            "coveringDomains": [
                {"type": d.feature_type, "start": d.start, "end": d.end, "description": d.description}
                for d in ctx.covering_domains
            ],
            "nearbyVariants": [
                {"position": v.position, "description": v.description}
                for v in ctx.nearby_variants
            ],
        }
    if isinstance(ctx, curator_info.ClinVarReportContext):
        return {
            "clinvarVariationId": ctx.variation_id,
            "segregationMentions": [
                {"text": m.text, "matchedKeyword": m.matched_keyword} for m in ctx.segregation_mentions
            ],
            "transMentions": [
                {"text": m.text, "matchedKeyword": m.matched_keyword} for m in ctx.trans_mentions
            ],
        }
    raise TypeError(f"no curatorInfo serialization defined for {type(ctx).__name__}")


def build_reference_extensions(code: str, variant: VariantRecord) -> list[Extension]:
    """Build auxiliary curator links/facts without making them evidence claims."""
    extensions: list[Extension] = []
    try:
        url = reference_links.reference_url_for_criterion(code, variant)
    except Exception:
        # Reference pages are optional navigation aids. Their network/lookup
        # failure must never erase the criterion assessment itself.
        url = None
    if url is not None:
        extensions.append(Extension(name="referenceLink", value=url))

    fetcher = _CODE_CURATOR_INFO_FETCHERS.get(code)
    if fetcher is not None:
        try:
            info_value = _serialize_curator_info(fetcher(variant))
        except Exception:
            info_value = None
        if info_value is not None:
            extensions.append(Extension(name="curatorInfo", value=info_value))
    return extensions


def _variant_identity(variant: VariantRecord) -> tuple[str, str, str]:
    gene = str(variant.info.get("GENE", ""))
    hgvsc = str(variant.info.get("HGVSC", ""))
    fallback = f"{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"
    safe_hgvsc = (hgvsc or fallback).replace(">", "_").replace(".", "_").replace("+", "p").replace("-", "m")
    return gene or "unknown-gene", hgvsc, safe_hgvsc


def build_workflow_evidence_line(
    code: str,
    variant: VariantRecord,
    *,
    status: str,
    description: str,
    details: Optional[dict[str, Any]] = None,
) -> dict:
    """Emit a neutral VA-Spec line for a non-scoreable workflow state.

    `assessment` carries no `summary`: it would duplicate the EvidenceLine's
    own top-level `description` below verbatim (see automated_va_spec.
    assessment_details()'s docstring for the same call made on the
    scored-line side, 2026-09-18).
    """
    gene, _hgvsc, safe_hgvsc = _variant_identity(variant)
    assessment = {
        "criterion": code,
        "status": status,
    }
    if details:
        assessment.update(details)
    extensions = [Extension(name="bh26AssessmentDetails", value=assessment)]
    extensions.extend(build_reference_extensions(code, variant))
    line = EvidenceLine(
        id=evidence_line_id(code, gene, safe_hgvsc),
        directionOfEvidenceProvided=Direction.NEUTRAL,
        description=description,
        specifiedBy=Method(
            methodType=code,
            name=f"ACMG/AMP {code} assessment ({status.lower().replace('_', ' ')})",
            reportedIn=ACMG_2015_METHOD_DOCUMENT,
        ),
        extensions=extensions,
    )
    return validate_integrated_line(
        line.model_dump(mode="json", exclude_none=True), code
    )


def build_automated_evidence_line(result, variant: VariantRecord) -> dict:
    """Map one evidence-cli result and attach main's curator-reference extensions."""
    from acmg_pipeline.automated_va_spec import (
        _curator_hints_from_result,
        assessment_details,
        to_evidence_line,
        validate_1_0_1,
    )

    status = CriterionStatus(result.status)
    if status in {CriterionStatus.MET, CriterionStatus.NOT_MET}:
        line = to_evidence_line(result)
        additions = [e.model_dump(mode="json", exclude_none=True)
                     for e in build_reference_extensions(result.criterion, variant)]
        if additions:
            line.setdefault("extensions", []).extend(additions)
        # Re-key onto the integrated document's single id scheme. The urn
        # to_evidence_line() produced is correct for the standalone
        # `acmg evaluate` CLI, but inside the 28-line document it would be
        # the only line addressed differently - and it hashes the outcome,
        # so it moves whenever a threshold changes. See evidence_line_id().
        gene, _hgvsc, safe_hgvsc = _variant_identity(variant)
        line["id"] = evidence_line_id(result.criterion, gene, safe_hgvsc)
        # validate_1_0_1() still runs: these lines always score, so they
        # are held to the strict schema (no relaxation for a missing
        # evidenceOutcome) as before the merge.
        return validate_1_0_1(line, result.criterion)

    details = assessment_details(result)
    summary = result.summary or f"{result.criterion} was not scored ({status.value})."
    line = build_workflow_evidence_line(
        result.criterion,
        variant,
        status=status.value,
        description=summary,
        details={key: value for key, value in details.items()
                 if key not in {"criterion", "status", "summary"}},
    )
    curator_hints = _curator_hints_from_result(result)
    if curator_hints:
        line.setdefault("extensions", []).append({
            "name": "curatorHints",
            "value": curator_hints,
        })
    return validate_integrated_line(line, result.criterion)


def build_stub_evidence_line(code: str, variant: VariantRecord) -> dict:
    """
    A minimal EvidenceLine for one of the seven codes this integrated
    pipeline doesn't implement - NOT a judgment. Carries a `referenceLink`
    extension when one is available (see
    acmg_pipeline.criteria.reference_links, built from doc/recs for expert
    board.docx's per-criterion "what to show the curator" asks) so a
    curator or another team's tool has a direct link to check, with no
    evidenceOutcome/strengthOfEvidenceProvided/contributions claiming a
    judgment this pipeline never made. PM1/PM5 additionally get a
    `curatorInfo` extension (see acmg_pipeline.criteria.curator_info) -
    structured UniProt domain/nearby-variant FACTS, still not a hotspot/
    benign-nearby verdict.

    [Why `extensions`, not `reportedIn` - decided 2026-09-16, don't "fix"
     this without re-reading]
      `Document` (what EvidenceLine.reportedIn holds) validates fine with
      only `urls` set, no `pmid` required - confirmed against the real
      ga4gh.va_spec Pydantic model - so reportedIn=[Document(urls=[url])]
      was considered as a more spec-native alternative to a custom
      extension. Deliberately NOT used: `reportedIn` means "the source this
      evidence was reported in" - using it for a gnomAD/ClinVar/UniProt/
      AutoPVS1 page this pipeline never actually read or evaluated would
      overstate this line's own claim (the whole point of this function is
      that NO evaluation happened). `extensions` correctly marks the link
      as auxiliary/navigational, matching this function's own "NOT a
      judgment" framing above.

    `directionOfEvidenceProvided` is a required VA-Spec EvidenceLine field
    (confirmed against the real ga4gh.va_spec Pydantic model - omitting it
    raises a validation error) with no "not evaluated" option in Direction
    (supports/neutral/disputes only - see ga4gh.va_spec.base.core.Direction);
    NEUTRAL is used here, same as this module's own _direction_of_evidence()
    does for a not_clear literature judgment - but unlike that case, nothing
    was actually evaluated here, so `description` says so explicitly rather
    than leaving NEUTRAL to imply "assessed as neutral evidence".

    A line is always returned, even if no reference URL exists, so the
    integrated interface has exactly one explicit state for all 28 codes.

    Raises ValueError for `code` in IMPLEMENTED_CODES (currently just PP1
    among reference_links.CODES_WITH_REFERENCE_URL - the other 4
    implemented codes have no doc/recs-for-expert-board.docx reference URL
    at all) - that code has a REAL judgment via build_evidence_line(), and
    a stub line here would collide on the exact same EvidenceLine id
    (both build f"evline:{gene}_{safe_hgvsc}_{code}"), silently shadowing
    real evidence with a placeholder if a caller ever included both. The
    caller is responsible for only invoking this for codes it has no real
    AggregatedJudgment for - see acmg_pipeline.pipeline.classify_variant_
    from_structured_input() once it's wired to call this function.
    """
    if code in IMPLEMENTED_CODES:
        raise ValueError(
            f"{code!r} is an implemented code (has a real judgment via "
            "build_evidence_line()) - build_stub_evidence_line() is only for "
            "codes this project has no real evidence for."
        )
    return build_workflow_evidence_line(
        code,
        variant,
        status=CriterionStatus.UNKNOWN.value,
        description=(
            f"{code} is not evaluated by this integrated pipeline. No judgment was "
            "made; any referenceLink or curatorInfo extension is an auxiliary "
            "curation aid only."
        ),
        details={"missingInputs": ["criterion implementation"]},
    )
