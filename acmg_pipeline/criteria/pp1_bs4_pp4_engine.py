"""
pp1_bs4_pp4_engine.py

Wires acmg_pipeline.criteria.pp4_pp1_bs4.evaluate_locus_evidence() (pulled
into main 2026-09-17 from r-kobayashi's pp4_pp1_bs4 branch - see that
module's own docstring) into this project's classify()/VA-Spec pipeline,
per the user's explicit direction (2026-09-17): "pipelineと繋げて。テスト
も流れるようにして" (connect it to the pipeline, make sure tests run too).

[Why this is a separate module from pp4_pp1_bs4.py]
  pp4_pp1_bs4.py is the reviewed judgment logic itself (ClinGen 2024
  Bayesian points, phenotype matching, family segregation scoring) and is
  left exactly as pulled in, unmodified. Everything specific to THIS
  project's own integration - the points-to-Strength mapping and the
  CriterionEvidence/EvidenceLine conversion - lives here instead, the same
  separation segregation.py (judgment logic) already keeps from export.py
  (VA-Spec conversion).

[The points-to-Strength mapping follows the paper's own Table 4]
  evaluator.py's DIAGNOSTIC_YIELD_POINT_TABLE and pp4_pp1_bs4.py's
  segregation scoring both produce points on the same Tavtigian-compatible
  scale acmg_pipeline.classification already uses (1/2/4/8 = Supporting/
  Moderate/Strong/Very_Strong). to_criterion_evidence() (below) calls
  acmg_pipeline.criteria.pp1_pp4_strength_table.combined_pp1_pp4_strength()
  to floor PP4's and PP1's points to a Strength TOGETHER (not
  independently, which can report a combined strength the shared
  +5.0-point cap does not actually support) - see that module's own
  docstring for the table itself and its tie-break rule.

[No curated reference registry - PP4/PP1/BS4 are always answered from a
 live literature search (2026-09-17, superseding this module's earlier
 config/pp4_reference_records.json design from the same day)]
  evaluate_locus_evidence() needs a PP4ReferenceRecord (gene-phenotype
  diagnostic-yield data) - a fact no tool can compute from a single
  patient. This module originally read that fact from a hand-curated
  config/pp4_reference_records.json registry (DRAFT/APPROVED/
  AUTO_EXTRACTED), the same "only approved entries load" convention
  config/bs1_thresholds_draft.json still uses for BS1. The user's explicit
  direction (2026-09-17), after several rounds of discussion, was to drop
  that registry entirely:
    - Confirmed there is no "AI drafts, a human approves before use"
      principle in this project - PS3/BS3/PS4 already show the real
      pattern: a live LLM judgment against real literature becomes the
      CriterionEvidence immediately, carrying a disclosure rather than
      being withheld pending a separate approval step.
    - A registry that only a human can populate reintroduces exactly the
      "no data until someone does the curation work first" cold-start
      problem this project is meant to help with in the first place - the
      tool should find candidate evidence itself, not wait for it.
    - An auto-caching layer that persists search results back into the
      registry (an earlier, briefly-tried middle ground) was rejected too:
      it does not match PS3/BS3/PS4's own ephemeral per-call caching, and
      a live PubMed+LLM search is judged rare enough (PP4 does not apply
      to most variants) that repeating it per call is an acceptable cost,
      not a real one.
  evaluate() therefore always calls acmg_pipeline.pp4_literature_search.
  search_diagnostic_yield() (PubMed MCP + this project's own LLM) using
  clinical_note.diagnosis as the phenotype description - the same
  free-text field acmg_pipeline.export previously only fed into VA-Spec's
  objectCondition, now put to a second use. No diagnosis, or nothing
  confirmed found, means an honest UNKNOWN, not a guess.

[Why phenotype_match is assumed True, not computed, for a literature-
 found reference]
  The search query was itself "this gene + this patient's diagnosis," so
  a match is true by construction. (An earlier, now-removed design tried
  two alternative matchers here - an exact match against a curated
  required-HPO-term list, which needs curated data this path has none of,
  and a PubCaseFinder gene-ranking proxy, which real testing against case3
  (MYH7) showed can rank the correct, ClinVar-established gene outside its
  own top 10 for a real but sparse HPO profile. Neither fit a literature-
  derived reference, so this module no longer tries to compute the
  question at all here.)

[method_comparable is assumed True, disclosed as unconfirmed]
  Whether the CURRENT test's methodology is comparable to the cited
  study's is not something a literature search can verify. Rather than
  block PP4 entirely (the same reasoning PS3/BS3/PS4 use for their own
  heuristic, unconfirmed strength estimates), this is assumed True and the
  reference's own source_citation always says so explicitly.

[Why the denominator matters - confirmed with two real experiments]
  Two real searches (case3/MYH7, against both a GeneReviews entry and a
  live PubMed hit) found the same trap twice: a gene's reported percentage
  is very often "share among already gene-positive patients," not the
  true overall diagnostic yield - using the former would badly overstate
  the evidence. acmg_pipeline.pp4_literature_search (not this file)
  discards any percentage not explicitly confirmed as the latter, rather
  than guessing which denominator applies.
"""

from __future__ import annotations

from typing import Optional

from acmg_pipeline.classification import CriterionEvidence, Strength
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.criteria.pp1_pp4_strength_table import combined_pp1_pp4_strength
from acmg_pipeline.criteria.pp4_pp1_bs4 import (
    LocusEvidenceResult,
    PhenotypeMatchResult,
    PP4ReferenceRecord,
    evaluate_locus_evidence,
)
from acmg_pipeline.vcf_record import VariantRecord

async def _build_reference_from_literature(gene: str, diagnosis: str) -> Optional[PP4ReferenceRecord]:
    """Ask acmg_pipeline.pp4_literature_search for a confirmed overall
    diagnostic-yield statistic and, if found, wrap it as a PP4ReferenceRecord.

    Returns None when nothing usable was found - callers must treat that
    as "not evaluable," never as a negative (benign) observation.

    [locus_model is approximated from yield alone - a known, documented
     simplification (2026-09-18)]
      evaluator.py's high_yield_homogeneous check (which suppresses PP1 so
      it isn't double-counted with a locus-homogeneous, high-yield PP4 -
      the paper's own central CTNS/cystinosis argument) requires
      reference.locus_model == "homogeneous". A live literature search has
      no way to confirm the paper's actual test for this ("is the
      diagnostic yield substantially specific to ONE gene, i.e. is there
      literature evidence no other gene is known to cause this phenotype" -
      see Biesecker et al. 2024's own CTNS/FBN1 examples, each backed by an
      explicit "no other gene known" statement, not just a yield number;
      the paper's own step 2.b explicitly allows "locus homogeneity but
      low diagnostic yield" as a distinct case, so yield and homogeneity
      are not the same fact). Confirming genuine locus homogeneity would
      need a further LLM literature judgment (a new hallucination surface)
      or a curated per-disease homogeneity list (the registry design this
      project's PP4 rework deliberately dropped - see this module's own
      docstring). Per the user's explicit direction (2026-09-18), this
      project instead approximates: yield alone already implies the
      homogeneity treatment when it exceeds the same 90% threshold
      evaluator.py's high_yield_homogeneous check already uses.
      This approximation's real-world impact is bounded: every yield that
      can trigger it (>90%) already independently saturates PP4's own
      +5.0 cap on its own (Table 2 reaches +5.0 at just 81.6% yield), so
      getting the homogeneity call "wrong" here never changes the total
      combined PP1+PP4 point value or final classification - it only
      changes whether PP1 is reported as an independently-confirmed MET
      criterion (when it should have been suppressed as redundant with
      PP4) or correctly NOT_MET.
    """
    from acmg_pipeline import pp4_literature_search

    result = await pp4_literature_search.search_diagnostic_yield(gene, diagnosis)
    if not result.found or result.yield_fraction is None:
        return None

    caution = (
        f" [CAUTION: small sample size n={result.sample_size} - interpret with caution]"
        if result.sample_size is not None and result.sample_size < 20 else ""
    )
    return PP4ReferenceRecord(
        reference_id=f"pubmed:{result.pmid}",
        gene=gene,
        phenotype_label=diagnosis,
        phenotype_hpo=(),
        locus_model="homogeneous" if result.yield_fraction > 0.90 else "heterogeneous",
        diagnostic_yield=result.yield_fraction,
        testing_method="unspecified",
        source_citation=(
            f"PMID:{result.pmid} (auto-extracted via literature search, unconfirmed - "
            f"requires human curator review; denominator: {result.denominator_description}; "
            f"sample_size={result.sample_size}){caution}"
        ),
        source_note=result.quote or "",
    )


def to_criterion_evidence(result: LocusEvidenceResult) -> dict[str, CriterionEvidence]:
    """Map one LocusEvidenceResult onto separate PP1/BS4/PP4 CriterionEvidence.

    PP4's and PP1's raw point contributions (0 when not applicable/not
    evaluable) are floored to a Strength TOGETHER, via
    pp1_pp4_strength_table.combined_pp1_pp4_strength() - Table 4 from
    Biesecker et al., 2024 - rather than independently: PP1 is scored from
    result.segregation.pp1_points_used, already zeroed by
    evaluate_locus_evidence() itself when high-yield locus homogeneity
    makes it redundant with PP4, so reading it directly (rather than
    pp1_points_raw) correctly reports PP1 as NOT_MET in that case, not as a
    fabricated MET that double-counts the same evidence. BS4 is scored from
    segregation.bs4_met/bs4_points (always exactly -4.0 in the current
    segregation scoring logic - see pp4_pp1_bs4.py - so no points-to-
    strength floor is needed for it; Table 4 does not cover BS4).
    """
    evidence: dict[str, CriterionEvidence] = {}

    pp4_points = result.pp4.points if result.pp4 is not None and result.pp4.applicable else 0.0
    pp1_points = result.segregation.pp1_points_used if result.segregation.pp1_evaluable else 0.0
    strengths = combined_pp1_pp4_strength(pp1_points=pp1_points, pp4_points=pp4_points)

    if result.pp4 is not None and result.pp4.applicable:
        evidence["PP4"] = (
            CriterionEvidence(code="PP4", status=CriterionStatus.MET, strength=strengths.pp4,
                              source=f"pp1_bs4_pp4_engine: {result.pp4.reference_id} - {result.pp4.source_citation}")
            if strengths.pp4 is not None
            else CriterionEvidence(code="PP4", status=CriterionStatus.NOT_MET,
                                    source=f"pp1_bs4_pp4_engine: {result.pp4.reference_id} (below Supporting threshold)")
        )
    elif result.pp4 is not None:
        evidence["PP4"] = CriterionEvidence(
            code="PP4", status=CriterionStatus.NOT_MET,
            source=f"pp1_bs4_pp4_engine: {result.pp4.reason} ({result.phenotype_match.reason})",
        )
    else:
        evidence["PP4"] = CriterionEvidence(
            code="PP4", status=CriterionStatus.UNKNOWN,
            source=f"pp1_bs4_pp4_engine: phenotype not evaluable ({result.phenotype_match.reason})",
        )

    if not result.segregation.pp1_evaluable:
        evidence["PP1"] = CriterionEvidence(
            code="PP1", status=CriterionStatus.UNKNOWN,
            source="pp1_bs4_pp4_engine: segregation not evaluable (no family data or unknown inheritance mode)",
        )
    else:
        evidence["PP1"] = (
            CriterionEvidence(code="PP1", status=CriterionStatus.MET, strength=strengths.pp1,
                              source="pp1_bs4_pp4_engine: family segregation")
            if strengths.pp1 is not None
            else CriterionEvidence(code="PP1", status=CriterionStatus.NOT_MET,
                                    source="pp1_bs4_pp4_engine: family segregation below Supporting threshold")
        )

    if not result.segregation.bs4_evaluable:
        evidence["BS4"] = CriterionEvidence(
            code="BS4", status=CriterionStatus.UNKNOWN,
            source="pp1_bs4_pp4_engine: segregation not evaluable (no family data or unknown inheritance mode)",
        )
    elif result.segregation.bs4_met:
        evidence["BS4"] = CriterionEvidence(
            code="BS4", status=CriterionStatus.MET, strength=Strength.STRONG,
            source="pp1_bs4_pp4_engine: clear non-segregation",
        )
    else:
        evidence["BS4"] = CriterionEvidence(
            code="BS4", status=CriterionStatus.NOT_MET,
            source="pp1_bs4_pp4_engine: no non-segregation observed",
        )

    return evidence


def build_evidence_line(code: str, evidence: CriterionEvidence, variant: VariantRecord) -> dict:
    """Convert one PP1/BS4/PP4 CriterionEvidence into a real VA-Spec EvidenceLine.

    Mirrors acmg_pipeline.export.build_evidence_line() (the literature
    engine's own builder) closely enough to slot into the same 28-line
    document, but works from a plain CriterionEvidence instead of an
    AggregatedJudgment - there is no per-paper contribution list here, so
    `hasEvidenceItems`/`reportedIn` are omitted rather than fabricated.
    """
    from acmg_pipeline.export import (
        ACMG_2015_METHOD_DOCUMENT, _variant_identity, build_reference_extensions,
        evidence_line_id, validate_integrated_line,
    )
    from acmg_pipeline.automated_va_spec import details_as_extensions
    from acmg_pipeline.classification import PATHOGENIC_CODES
    from ga4gh.core.models import Extension
    from ga4gh.va_spec.base.core import Direction, EvidenceLine, Method

    gene, _hgvsc, safe_hgvsc = _variant_identity(variant)
    status = evidence.status

    if status == CriterionStatus.MET:
        direction = Direction.SUPPORTS if code in PATHOGENIC_CODES else Direction.DISPUTES
        outcome_code = code if evidence.strength and evidence.strength.value == _default_strength(code) else f"{code}_{evidence.strength.value}"
        evidence_outcome = {
            "name": f"ACMG 2015 {code} criterion met",
            "primaryCoding": {"system": "ACMG Guidelines, 2015", "code": outcome_code},
        }
        strength_block = {"primaryCoding": {"system": "ACMG Guidelines, 2015", "code": evidence.strength.value}}
    elif status == CriterionStatus.NOT_MET:
        direction = Direction.NEUTRAL
        evidence_outcome = {
            "name": f"ACMG 2015 {code} criterion not met",
            "primaryCoding": {"system": "ACMG Guidelines, 2015", "code": f"{code}_not_met"},
        }
        strength_block = None
    else:
        direction = Direction.NEUTRAL
        evidence_outcome = None
        strength_block = None

    # No `summary`/`criterion` keys: they would duplicate `description` below
    # (both come from evidence.source) and `specifiedBy.methodType` (both are
    # `code`) verbatim. Unlike automated_va_spec.assessment_details() (see its
    # docstring), nothing here needs `criterion` as a list-disambiguator.
    # An UNKNOWN status (no family/segregation/phenotype data to decide from)
    # is reported as not_met, not unknown - PP1/BS4/PP4 are implemented, this
    # is "ran, couldn't reach a verdict", the same reframing export.
    # build_evidence_line()/build_automated_evidence_line() apply, 2026-09-18.
    # A curatorHint discloses the real reason.
    curator_hints = []
    reported_status = status
    if status == CriterionStatus.UNKNOWN:
        reported_status = CriterionStatus.NOT_MET
        curator_hints.append({
            "severity": "caution", "category": "unevaluated",
            "message": f"{code} could not actually be evaluated from the "
                       "available data (reported as not_met rather than left "
                       "unknown).",
        })
    # `status` sits as its own top-level extension, not grouped under one
    # bh26AssessmentDetails object - see details_as_extensions(), 2026-09-18.
    extensions = [Extension(**d) for d in details_as_extensions({"status": reported_status.value})]
    if curator_hints:
        extensions.append(Extension(name="curatorHints", value=curator_hints))
    ref_extensions, reference_evidence_items = build_reference_extensions(code, variant)
    extensions.extend(ref_extensions)

    line = EvidenceLine(
        id=evidence_line_id(code, gene, safe_hgvsc),
        directionOfEvidenceProvided=direction,
        description=evidence.source or None,
        hasEvidenceItems=reference_evidence_items or None,
        specifiedBy=Method(
            methodType=code,
            name=f"ACMG/AMP {code} assessment (family segregation / phenotype specificity)",
            reportedIn=ACMG_2015_METHOD_DOCUMENT,
        ),
        extensions=extensions,
    )
    payload = line.model_dump(mode="json", exclude_none=True)
    if evidence_outcome is not None:
        payload["evidenceOutcome"] = evidence_outcome
    if strength_block is not None:
        payload["strengthOfEvidenceProvided"] = strength_block
    return validate_integrated_line(payload, code)


def _default_strength(code: str) -> str:
    from acmg_pipeline.export import default_strength
    return default_strength(code)


async def evaluate(
    variant: VariantRecord,
    clinical_note: ClinicalNoteExtraction,
    config: dict,
) -> dict[str, CriterionEvidence]:
    """Return {code: CriterionEvidence} for PP1, BS4, PP4 - the entry point
    acmg_pipeline.pipeline.evaluate_variant_evidence_lines() calls.

    `config` is accepted for signature parity with the automated engine's
    `automated_config` convention but currently unused - PP4's diagnostic-
    yield input always comes from a live literature search now (see this
    module's own docstring), not from a caller-supplied path/override.

    PP1/BS4 (family co-segregation) never need the literature search at
    all - they are scored from clinical_note.family.relatives alone, per
    ClinGen 2024 Table 3 (pp4_pp1_bs4._score_family_segregation()). Only
    PP4 needs a diagnostic-yield statistic. So when there is no diagnosis
    to search for, or the search finds nothing confirmed, this function
    still calls evaluate_locus_evidence() - with phenotype_match_override
    set to matched=None (not True) so PP4 correctly comes back "not
    evaluable" - rather than short-circuiting to UNKNOWN for all three;
    evaluate_pp4() is skipped whenever phenotype.matched is None (see
    pp4_pp1_bs4.evaluate_locus_evidence()), so PP1/BS4's segregation
    scoring below is unaffected either way.

    HPO normalization (acmg_pipeline.hpo_mondo_extraction.normalize_hpo(), a
    TogoMCP + LLM round trip) and the literature search only run once
    clinical_note.diagnosis is non-empty - with nothing to search for,
    PP4 is already decided (not evaluable), so spending real network/LLM
    cost first would be wasted. Both imports are local (not top-level)
    because hpo_mondo_extraction and pp4_literature_search each require
    VLLM_BASE_URL/VLLM_API_KEY at import time (same reason clinical_note.py
    lazily imports clinical_extraction.py instead of importing it at
    module load).
    """
    gene = str(variant.info.get("GENE", ""))
    diagnosis = (clinical_note.diagnosis or "").strip()

    reference = await _build_reference_from_literature(gene, diagnosis) if diagnosis else None

    if reference is not None:
        phenotype_match_override = PhenotypeMatchResult(
            matched=True, reason="phenotype_match_assumed_from_literature_search_query",
        )
    else:
        phenotype_match_override = PhenotypeMatchResult(
            matched=None,
            reason=(
                "no_diagnosis_available_to_search_literature_for" if not diagnosis
                else "literature_search_found_no_confirmed_overall_diagnostic_yield_statistic"
            ),
        )
        # PP1/BS4's segregation scoring below needs SOME PP4ReferenceRecord
        # to call evaluate_locus_evidence() with, but evaluate_pp4() is
        # never invoked while phenotype_match_override.matched is None (see
        # pp4_pp1_bs4.py) - diagnostic_yield=0.0 here is inert, not a guess.
        reference = PP4ReferenceRecord(
            reference_id="none", gene=gene, phenotype_label=diagnosis,
            phenotype_hpo=(), locus_model="heterogeneous", diagnostic_yield=0.0,
            testing_method="unspecified", source_citation="",
        )

    from acmg_pipeline import hpo_mondo_extraction
    clinical_note = await hpo_mondo_extraction.normalize_hpo(clinical_note)

    result = evaluate_locus_evidence(
        clinical_note,
        reference,
        method_comparable=True,  # unconfirmed - disclosed in reference.source_citation
        inheritance_mode=None,  # falls back to clinical_note.family.inheritance_pattern
        fully_penetrant=None,
        low_phenocopy=None,
        ar_case_mode=None,
        phenotype_match_override=phenotype_match_override,
    )
    return to_criterion_evidence(result)
