from acmg_pipeline.constants import CriterionStatus
"""ClinGen general PVS1 decision tree for sequence variants.

Evidence acquisition is intentionally separate. This evaluator consumes versioned, reviewed
transcript, NMD, splice, protein-region, population-LoF, and initiation assessments and never
turns a missing assessment into a negative conclusion.
"""

from copy import deepcopy

from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import (
    annotation_context, condition_scope, normalize_inheritance, ontology_related,
    resolved_condition, result as _base_result, reviewed_or_automated,
)
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.providers.initiation import IC01_METHOD as INITIATION_IC01_METHOD
from acmg_pipeline.providers.protein_region import UNIPROT_METHOD
from acmg_pipeline.providers.splice_default import METHOD as SPLICE_DEFAULT_METHOD
from acmg_pipeline.vcf_record import VariantRecord


def _region_relevance_review(region):
    """The "this is a flagged prediction, not a curator review" caveat for
    NF04/NF06, only when protein_region.py's UniProt-derived first pass (not a
    curator) is what actually answered critical_region_disrupted/
    region_biologically_relevant - see that module's own docstring.
    """
    if not region or region.get("region_relevance_method") != UNIPROT_METHOD:
        return []
    return [
        "NF04/NF06 (critical/biologically relevant region) answered from UniProt feature "
        f"overlap ({region.get('region_relevance_source_version')}), not a curator's own "
        "review - confirm before this MET is used in a classification."
    ]


TRUNCATING = {"stop_gained": "STOP_GAINED", "frameshift_variant": "FRAMESHIFT"}
CANONICAL_SPLICE = {"splice_donor_variant", "splice_acceptor_variant"}
OTHER_LOF = {"transcript_ablation"}


def result(code, input_data, status, summary, **kwargs):
    """Keep PVS1's curator-facing explanation with its decision tree."""
    assert code == "PVS1"
    outcome = {
        CriterionStatus.MET: "Outcome: MET; the available evidence satisfies PVS1.",
        CriterionStatus.NOT_MET: "Outcome: NOT_MET; PVS1 was evaluated but its requirements were not satisfied.",
        CriterionStatus.UNKNOWN: "Outcome: UNKNOWN; no MET or NOT_MET judgment is made from the available information.",
    }[status]
    message = (
        f"{summary.rstrip('.')}. PVS1 evaluates predicted loss-of-function variants only when "
        f"loss of function is an established disease mechanism for the gene. {outcome}"
    )
    return _base_result(code, input_data, status, message, **kwargs)


def _node(node_id, name, node_result, value=None, evidence=(), source="clingen_pvs1_2018",
          next_node=None, note=None):
    return {
        "node_id": node_id,
        "name": name,
        "result": node_result,
        "value": value,
        "evidence_ids": [item["evidence_id"] for item in evidence if item and item.get("evidence_id")],
        "rule_id": node_id,
        "rule_source": source,
        "next_node": next_node,
        "note": note,
    }


def _base_context(input_data, gene=None):
    condition = input_data.get("condition")
    return {
        "gene": gene,
        "condition_id": condition,
        "condition_label": input_data.get("condition_label"),
        "condition_status": "PROVIDED" if condition else "NOT_PROVIDED",
        "condition_specific": False,
        "mechanism_scope": "UNKNOWN",
        # How the curated mechanism was matched to this case, reported rather than assumed.
        # Only EXACT is produced today: normalizing OMIM/Orphanet identifiers and reading
        # ClinGen lumping/splitting decisions are what would produce EQUIVALENT and INCLUDED,
        # and neither is available here, so anything short of an identifier match stays
        # UNKNOWN instead of being upgraded on a resemblance.
        "disease_match": "UNKNOWN",
        "mechanism_source": None,
        "inheritance": normalize_inheritance(input_data.get("inheritance")),
        "moi_match": "UNKNOWN",
        "applicability": "NOT_EVALUATED",
    }


def _rules(config):
    policy = config.get("PVS1")
    if not isinstance(policy, dict):
        return None
    ruleset = policy.get("ruleset")
    rules = policy.get("rules")
    if not isinstance(ruleset, dict) or not isinstance(rules, dict):
        return None
    sources = ruleset.get("sources")
    threshold = rules.get("protein_loss_threshold")
    if not all(isinstance(ruleset.get(key), str) and ruleset[key]
               for key in ("name", "version")):
        return None
    if not isinstance(sources, list) or not sources or any(
            not isinstance(item, dict) or not item.get("id") for item in sources):
        return None
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not 0 < threshold <= 1:
        return None
    precedence = rules.get("mechanism_source_precedence", [])
    if not isinstance(precedence, list) or any(
            not isinstance(item, str) or not item for item in precedence):
        return None
    return policy


def _candidates(category, input_data, services, transcript=None):
    records = services.evidence.get_candidates(category, Variant(**input_data["variant"]))
    records = [item for item in records if reviewed_or_automated(item)]
    if transcript is not None:
        records = [item for item in records
                   if not item.get("transcript") or item.get("transcript") == transcript]
    return records


def _select_context_record(category, input_data, services, transcript):
    records = _candidates(category, input_data, services, transcript)
    condition = input_data.get("condition")
    if condition:
        exact = [item for item in records if item.get("condition") == condition]
        selected = exact or [item for item in records if not item.get("condition")]
    else:
        selected = [item for item in records if not item.get("condition")]
    if not selected:
        return None, "missing"
    if len(selected) != 1:
        return selected, "conflict"
    return selected[0], None


# An unqualified mode is the less specific statement of the same thing, so it is compatible
# with either qualification of it. ACMG asks for compatible inheritance rather than identical
# inheritance, and the two differ exactly here: a source that records "X-linked" without
# saying which zygosity is affected has not contradicted a case assessed as X-linked
# recessive. Autosomal dominant and autosomal recessive are not related this way and never
# stand in for one another.
COMPATIBLE_MODES = {
    "x_linked": {"x_linked", "x_linked_dominant", "x_linked_recessive"},
}


def _moi_compatible(declared, inheritance):
    return inheritance in COMPATIBLE_MODES.get(declared, {declared})


def _moi_applicable(record, inheritance):
    """Whether a curated mechanism record may be read as this case's mechanism.

    The same gene and the same disease can carry different mechanisms under different
    inheritance modes, so a mechanism curated for one mode is not evidence about another.
    A record that names no mode is not scoped to one and stays usable, exactly as a record
    that names no condition stays usable across conditions. A record that does name one is
    usable only when this case names a compatible mode - including when this case names none,
    because an unknown mode cannot be confirmed compatible with anything, and assuming it is
    would borrow the very evidence the scoping exists to keep apart.
    """
    declared = record.get("inheritance")
    if not declared:
        return True
    mode = normalize_inheritance(declared)
    return mode is not None and _moi_compatible(mode, inheritance)


def _mechanism_rank(record, precedence):
    """Where a mechanism record sits in the source order ACMG lays out.

    A reviewed assessment outranks every derived one, whatever produced it. Among derived
    records the configured order decides, because "which source wins" is policy that a result
    has to be able to name, not a fact about this code - it is recorded in the rule set and
    travels into the provenance with everything else. A derived record from a source the
    policy does not rank sits below every source it does, rather than above them by accident.
    """
    if record.get("assessment_method") != "automated":
        return 0
    method = record.get("method")
    return 1 + (precedence.index(method) if method in precedence else len(precedence))


def _resolve_mechanism(input_data, services, annotation, context, precedence=()):
    records = _candidates("gene_disease", input_data, services, annotation["transcript"])
    gene = annotation.get("gene")
    gene_records = [item for item in records if item.get("gene") == gene]
    if records and not gene_records:
        return None, records, "gene_mismatch"
    inheritance = context["inheritance"]
    applicable = [item for item in gene_records if _moi_applicable(item, inheritance)]
    other_mode = [item for item in gene_records if item not in applicable]
    if other_mode and not applicable:
        # Unstated and contradicted are different questions for a curator - one is a field
        # nobody filled in, the other is a judgment about whether a mode's curation carries -
        # so the result says which, and neither is silently treated as the other.
        context["moi_match"] = "UNSTATED" if inheritance is None else "MISMATCH"
        return None, other_mode, "moi_mismatch"
    condition, case_via = resolved_condition(input_data)
    if condition:
        matched = [(item, resolved_condition(item)) for item in applicable]
        # An expert panel that looked at this phenotype and kept it out of the disease it
        # curated has answered the question outright, so its record is set aside before any
        # weaker relation is considered - a resemblance cannot reinstate a decided exclusion.
        ruled_out = {id(item) for item, _ in matched
                     if condition in condition_scope(item, "excluded")}
        considered = [(item, pair) for item, pair in matched if id(item) not in ruled_out]
        exact = [item for item, (value, _) in considered if value == condition]
        # Lumped into the curated disease by the same panel, so the curation is about this
        # phenotype - an explicit decision, unlike an ontology relation.
        included = [item for item, _ in considered
                    if condition in condition_scope(item, "included")
                    and id(item) not in {id(value) for value in exact}]
        if not exact and not included:
            if ruled_out:
                context["disease_match"] = "EXCLUDED"
                return None, [item for item, _ in matched if id(item) in ruled_out], "excluded"
            # A curation for a parent or a child disease is on file. It is not this disease,
            # so it cannot decide the mechanism, but reporting "no mechanism" would send a
            # curator hunting for evidence that is sitting under the neighbouring term.
            related = [item for item, (value, _) in considered
                       if ontology_related(input_data, condition, item, value)]
            if related:
                context["disease_match"] = "PARENT_CHILD"
                return None, related, "parent_child"
            # The gene is curated, but for diseases none of the match levels tie to this one.
            # That is not the same as having nothing: an expert panel has said what loss of
            # function does in this gene, and whether that carries to the disease being
            # assessed is a question about disease entities, which is a curator's to answer.
            # Both the scoped and unscoped records go with it, so the answer is decidable
            # from what the result carries.
            if any(value for _, (value, _) in considered):
                return None, [item for item, _ in considered], "other_disease"
        selected = exact or included or [item for item in applicable
                                         if not item.get("condition")]
        scope = "CONDITION_SPECIFIC" if (exact or included) else "GENE_LEVEL"
        # An identifier match on both sides is EXACT; anything that needed a mapping to line
        # the two up is EQUIVALENT, which is a weaker statement and is reported as one.
        vias = {case_via, *(how for item, (_, how) in considered
                            if id(item) in {id(value) for value in exact})}
        match_level = ("EXACT" if vias == {"identity"} else "EQUIVALENT") if exact else "INCLUDED"
    else:
        selected = [item for item in applicable if not item.get("condition")]
        scope = "GENE_LEVEL"
        match_level = "UNKNOWN"
    if not selected:
        return None, [], "missing"
    # Source precedence: a record that disagrees with a higher-ranked one is superseded by
    # it rather than in conflict with it, so only records at the best rank present decide.
    # Every record stays in the evidence, because a curator should see that they disagreed
    # and which source the answer came from. Two records at the same rank that disagree are a
    # genuine conflict, and still reported as one.
    ranks = {id(item): _mechanism_rank(item, list(precedence)) for item in selected}
    best = min(ranks.values())
    deciding = [item for item in selected if ranks[id(item)] == best]
    values = {item.get("lof_mechanism_established") for item in deciding}
    if len(values) > 1:
        return None, selected, "conflict"
    context["mechanism_scope"] = scope
    context["condition_specific"] = scope == "CONDITION_SPECIFIC"
    context["disease_match"] = match_level if scope == "CONDITION_SPECIFIC" else "UNKNOWN"
    context["mechanism_source"] = deciding[0].get("source")
    context["moi_match"] = ("MATCHED" if any(item.get("inheritance") for item in deciding)
                            else "NOT_SCOPED")
    return deciding[0], selected, None


def _deduplicate(items):
    values = {}
    for item in items:
        values[item.get("evidence_id", id(item))] = item
    return list(values.values())


def _applicability(status, missing, review):
    """Which of ACMG's non-judgments this is, without widening CriterionStatus.

    MET is a judgment. Everything else here is UNKNOWN, and UNKNOWN covers three situations a
    curator acts on differently: the criterion does not apply to this variant or disease, the
    information needed to decide it was not available, or the information conflicts and a
    person has to choose. Each exit already records which of those it is - an unresolved
    requirement, a review point, or neither - so this reads what is there rather than asking
    every exit to restate itself and risking the two drifting apart.
    """
    if status == CriterionStatus.MET:
        return "APPLICABLE"
    if review:
        return "MANUAL_REVIEW"
    if missing:
        return "NOT_EVALUATED"
    return "NOT_APPLICABLE"


def _finish(input_data, state, status, summary, *, strength=None, missing=(), review=(),
            extra_evidence=(), provenance=None):
    unresolved = list(missing)
    state["context"]["applicability"] = _applicability(status, unresolved, review)
    # Kept unwrapped for the preliminary assessment, whose reader must not be handed the
    # criterion-level wording ("Outcome: MET") for a run that reached no verdict.
    state["summary"] = summary
    scope = state["context"]["mechanism_scope"]
    condition_assessment = "MATCHED" if state["context"]["condition_specific"] else "NOT_EVALUATED"
    assessment_scope = {
        "CONDITION_SPECIFIC": "condition_specific",
        "GENE_LEVEL": "gene_level",
        "UNKNOWN": "unknown",
    }[scope]
    return result(
        "PVS1", input_data, status, summary, strength=strength,
        evidence=_deduplicate([*state["evidence"], *extra_evidence]), missing=unresolved,
        review=list(dict.fromkeys([*state.get("pending_review", []), *review])),
        provenance={
            "assessment_scope": assessment_scope,
            "condition_assessment": condition_assessment,
            "pvs1_ruleset": state["ruleset"],
            **(provenance or {}),
        },
        evaluation_context=state["context"], decision_trace=state["trace"],
        rules_used=state["rules_used"], warnings=state["warnings"],
        unresolved_requirements=unresolved,
    )


def _missing_record(input_data, state, category, node_id, name):
    state["trace"].append(_node(node_id, name, "UNKNOWN"))
    return _finish(input_data, state, CriterionStatus.UNKNOWN,
                   f"Reviewed {category} evidence unavailable", missing=[category])


def _conflicting_record(input_data, state, category, records, node_id, name):
    state["trace"].append(_node(node_id, name, "MANUAL_REVIEW", evidence=records))
    return _finish(input_data, state, CriterionStatus.UNKNOWN,
                   f"Conflicting {category} evidence", review=[f"Resolve {category} assessments"],
                   extra_evidence=records)


def _transcript(input_data, services, annotation, state, node_id="NF01"):
    assessment, issue = _select_context_record(
        "transcript_assessment", input_data, services, annotation["transcript"])
    if issue == "missing":
        return None, _missing_record(
            input_data, state, "transcript_assessment", node_id, "transcript_relevance")
    if issue == "conflict":
        return None, _conflicting_record(
            input_data, state, "transcript_assessment", assessment, node_id, "transcript_relevance")
    state["evidence"].append(assessment)
    relevance = assessment.get("relevance")
    if relevance == "NOT_RELEVANT":
        state["trace"].append(_node(
            node_id, "transcript_relevance", "NOT_APPLICABLE", relevance, [assessment]))
        return None, _finish(input_data, state, CriterionStatus.UNKNOWN,
                             "Variant is not on a biologically relevant transcript")
    if relevance != "RELEVANT":
        state["trace"].append(_node(node_id, "transcript_relevance", "UNKNOWN", relevance, [assessment]))
        return None, _finish(input_data, state, CriterionStatus.UNKNOWN,
                             "Transcript relevance unresolved", missing=["transcript_relevance"])
    state["trace"].append(_node(node_id, "transcript_relevance", "PASS", relevance, [assessment]))
    return assessment, None


def _region_path(input_data, services, annotation, state):
    region, issue = _select_context_record("protein_region", input_data, services, annotation["transcript"])
    if issue == "conflict":
        return _conflicting_record(
            input_data, state, "protein_region", region, "NF04", "critical_functional_region")
    if issue == "missing":
        region = None
    else:
        state["evidence"].append(region)

    critical = region.get("critical_region_disrupted") if region else None
    critical_result = "PASS" if critical is True else "FAIL" if critical is False else "UNKNOWN"
    state["trace"].append(_node(
        "NF04", "critical_functional_region", critical_result,
        critical, [region] if region else []))
    if critical is True:
        summary = ("LoF disrupts a UniProt-annotated critical functional region"
                   if region.get("region_relevance_method") == UNIPROT_METHOD
                   else "LoF disrupts a reviewed critical functional region")
        return _finish(input_data, state, CriterionStatus.MET, summary, strength="strong",
                       review=_region_relevance_review(region))

    population, pop_issue = _select_context_record(
        "population_lof", input_data, services, annotation["transcript"])
    if pop_issue == "conflict":
        return _conflicting_record(
            input_data, state, "population_lof", population, "NF05", "population_lof_frequency")
    if pop_issue == "missing":
        population = None
    else:
        state["evidence"].append(population)
    frequent = population.get("lof_variants_frequent") if population else None
    state["trace"].append(_node(
        "NF05", "population_lof_frequency", "FAIL" if frequent is True else "PASS",
        frequent, [population] if population else []))
    if frequent is True:
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "LoF variants are frequent in the affected exon or region")

    relevance = region.get("region_biologically_relevant") if region else None
    if relevance in (False, "NOT_RELEVANT"):
        state["trace"].append(_node(
            "NF06", "region_biological_relevance", "NOT_APPLICABLE", relevance, [region]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Affected protein region is not biologically relevant")
    if relevance not in (True, "RELEVANT"):
        state["trace"].append(_node(
            "NF06", "region_biological_relevance", "UNKNOWN", relevance,
            [region] if region else []))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Affected protein region relevance unresolved",
                       missing=["region_biological_relevance"])
    state["trace"].append(_node(
        "NF06", "region_biological_relevance", "PASS", relevance, [region]))

    lost = region.get("lost_residues")
    total = region.get("total_protein_length")
    if type(lost) is not int or type(total) is not int:
        state["trace"].append(_node(
            "NF07", "protein_loss_fraction", "UNKNOWN", None, [region]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Protein loss fraction unavailable", missing=["protein_loss_fraction"])
    if lost < 0 or total <= 0 or lost > total:
        state["trace"].append(_node(
            "NF07", "protein_loss_fraction", "MANUAL_REVIEW",
            {"lost_residues": lost, "total_protein_length": total}, [region]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Protein loss inputs are internally inconsistent",
                       review=["Correct lost residues and total protein length"])
    fraction = lost / total
    threshold = state["policy"]["rules"]["protein_loss_threshold"]
    strength = "strong" if fraction > threshold else "moderate"
    state["trace"].append(_node(
        "NF07", "protein_loss_fraction", "PASS", fraction, [region],
        note=f"Threshold is > {threshold}"))
    return _finish(input_data, state, CriterionStatus.MET,
                   "NMD-escaping LoF affects a biologically relevant protein region",
                   strength=strength, review=_region_relevance_review(region),
                   provenance={"protein_loss_fraction": fraction})


def _truncating_path(input_data, services, annotation, state):
    transcript, early = _transcript(input_data, services, annotation, state)
    if early:
        return early
    nmd, issue = _select_context_record("nmd_prediction", input_data, services, annotation["transcript"])
    if issue == "missing":
        return _missing_record(input_data, state, "nmd_prediction", "NF02", "nmd_prediction")
    if issue == "conflict":
        return _conflicting_record(
            input_data, state, "nmd_prediction", nmd, "NF02", "nmd_prediction")
    state["evidence"].append(nmd)
    if nmd.get("position_type") == "intron":
        state["pending_review"].append(
            "NF02/NF03 (NMD prediction / exon relevance) answered from the affected "
            "intron's position, a flagged approximation for a canonical splice donor/"
            "acceptor loss (see acmg_pipeline.providers.nmd's own docstring) - not the "
            "exact exon a real skipped-exon or retained-intron outcome would fall in.")
    predicted = nmd.get("predicted")
    if type(predicted) is not bool or not nmd.get("rule_source"):
        state["trace"].append(_node("NF02", "nmd_prediction", "UNKNOWN", predicted, [nmd]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Transcript-specific NMD prediction unresolved", missing=["nmd_prediction"])
    state["trace"].append(_node(
        "NF02", "nmd_prediction", "PASS", predicted, [nmd]))
    if not predicted:
        return _region_path(input_data, services, annotation, state)

    exon_relevance = transcript.get("exon_relevance")
    if exon_relevance in (False, "NOT_RELEVANT"):
        state["trace"].append(_node(
            "NF03", "exon_relevance", "NOT_APPLICABLE", exon_relevance, [transcript]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Affected exon is absent from biologically relevant transcripts")
    if exon_relevance not in (True, "RELEVANT"):
        state["trace"].append(_node(
            "NF03", "exon_relevance", "UNKNOWN", exon_relevance, [transcript]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Affected exon relevance unresolved", missing=["exon_relevance"])
    state["trace"].append(_node(
        "NF03", "exon_relevance", "PASS", exon_relevance, [transcript]))
    return _finish(input_data, state, CriterionStatus.MET,
                   "LoF is predicted to undergo NMD in a relevant transcript and exon",
                   strength="very_strong")


def _splice_path(input_data, services, annotation, state, rna):
    if rna and rna.get("lof_effect_confirmed") is True:
        assessment = deepcopy(rna)
        assessment["used_by"] = list(dict.fromkeys([*assessment.get("used_by", []), "PVS1"]))
        source = "clingen_splicing_2023"
    else:
        assessment, issue = _select_context_record(
            "splice_assessment", input_data, services, annotation["transcript"])
        if issue == "missing":
            return _missing_record(
                input_data, state, "splice_assessment", "SP01", "expected_splice_consequence")
        if issue == "conflict":
            return _conflicting_record(
                input_data, state, "splice_assessment", assessment, "SP01",
                "expected_splice_consequence")
        source = "clingen_pvs1_2018"
    state["evidence"].append(assessment)
    if assessment.get("method") == SPLICE_DEFAULT_METHOD:
        state["pending_review"].append(
            "SP01/SP02 (expected splice consequence/reading frame) answered from a default "
            "policy for canonical splice donor/acceptor variants, not a curator's own review "
            "or a variant-specific RNA/sequence finding - confirm before this result is used "
            "in a classification.")
    rescue = assessment.get("alternative_rescue")
    if rescue is True:
        state["trace"].append(_node(
            "SP01", "expected_splice_consequence", "NOT_APPLICABLE", "RESCUED",
            [assessment], source))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Alternative splicing is expected to rescue the LoF consequence")
    if type(rescue) is not bool:
        state["trace"].append(_node(
            "SP01", "expected_splice_consequence", "UNKNOWN", None, [assessment], source))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Alternative splice rescue unresolved", missing=["alternative_rescue"])
    outcome = assessment.get("splice_outcome")
    if outcome not in {"OUT_OF_FRAME_PTC", "IN_FRAME"}:
        state["trace"].append(_node(
            "SP01", "expected_splice_consequence", "UNKNOWN", outcome, [assessment], source))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Splice consequence unresolved", missing=["splice_outcome"])
    disrupted = assessment.get("reading_frame_disrupted")
    expected = outcome == "OUT_OF_FRAME_PTC"
    if type(disrupted) is not bool or disrupted != expected:
        state["trace"].append(_node(
            "SP02", "reading_frame_disrupted", "MANUAL_REVIEW", disrupted,
            [assessment], source, note=f"Outcome implies {expected}"))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Splice outcome and reading-frame assessment conflict",
                       review=["Reconcile splice consequence and reading frame"])
    state["trace"].append(_node(
        "SP01", "expected_splice_consequence", "PASS", outcome, [assessment], source))
    state["trace"].append(_node(
        "SP02", "reading_frame_disrupted", "PASS", disrupted, [assessment], source))
    if disrupted:
        return _truncating_path(input_data, services, annotation, state)
    _, early = _transcript(input_data, services, annotation, state, "NF01")
    return early or _region_path(input_data, services, annotation, state)


def _start_loss_path(input_data, services, annotation, state):
    _, early = _transcript(input_data, services, annotation, state, "IC00")
    if early:
        return early
    assessment, issue = _select_context_record(
        "initiation_assessment", input_data, services, annotation["transcript"])
    if issue == "missing":
        return _missing_record(
            input_data, state, "initiation_assessment", "IC01", "intact_alternative_transcript")
    if issue == "conflict":
        return _conflicting_record(
            input_data, state, "initiation_assessment", assessment, "IC01",
            "intact_alternative_transcript")
    state["evidence"].append(assessment)
    if assessment.get("alternative_transcript_method") == INITIATION_IC01_METHOD:
        state["pending_review"].append(
            "IC01 (intact alternative transcript) answered from an Ensembl-derived "
            "prediction (does another protein_coding transcript's start codon sit at a "
            "different genomic position), not a curator's own review, and not validated "
            "against a known real case - confirm before this result is used in a "
            "classification.")
    alternative = assessment.get("intact_alternative_transcript")
    if alternative is True:
        state["trace"].append(_node(
            "IC01", "intact_alternative_transcript", "NOT_APPLICABLE", True, [assessment]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "An intact biologically relevant alternative transcript is available")
    if type(alternative) is not bool:
        state["trace"].append(_node(
            "IC01", "intact_alternative_transcript", "UNKNOWN", alternative, [assessment]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Alternative transcript status unresolved",
                       missing=["intact_alternative_transcript"])
    state["trace"].append(_node(
        "IC01", "intact_alternative_transcript", "PASS", False, [assessment]))
    downstream = assessment.get("downstream_in_frame_start")
    if type(downstream) is not bool:
        state["trace"].append(_node(
            "IC02", "downstream_in_frame_start", "UNKNOWN", downstream, [assessment]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Downstream in-frame start status unresolved",
                       missing=["downstream_in_frame_start"])
    state["trace"].append(_node(
        "IC02", "downstream_in_frame_start", "PASS", downstream, [assessment]))
    pathogenic = assessment.get("upstream_pathogenic_evidence")
    if type(pathogenic) is not bool:
        state["trace"].append(_node(
            "IC03", "upstream_pathogenic_evidence", "UNKNOWN", pathogenic, [assessment]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Supporting pathogenic evidence unresolved",
                       missing=["upstream_pathogenic_evidence"])
    strength = "moderate" if pathogenic else "supporting"
    state["trace"].append(_node(
        "IC03", "upstream_pathogenic_evidence", "PASS", pathogenic, [assessment]))
    return _finish(input_data, state, CriterionStatus.MET,
                   "Initiation codon loss evaluated using the general PVS1 framework",
                   strength=strength)


def _preliminary_assessment(input_data, services, annotation, variant_type, state, rna):
    """What the variant-level tree concludes before any disease mechanism is established.

    Transcript relevance, NMD, exon and region significance, splice rescue and alternative
    initiation are decidable from the variant alone. Discarding that work because no disease
    mechanism was established would throw away exactly what a curator needs in order to
    supply the missing context, so it is computed and kept.

    It runs on its own state so nothing it records reaches the criterion's evidence, trace or
    warnings. What comes back is what PVS1 *would* conclude if loss of function were an
    established mechanism for the disease - which is not a claim that it is, and is why the
    strength it carries is named `candidate_strength` and never becomes the result's own.
    """
    scratch = {
        "context": dict(state["context"]),
        "trace": [],
        "evidence": [],
        "warnings": [],
        "rules_used": state["rules_used"],
        "ruleset": state["ruleset"],
        "policy": state["policy"],
        "summary": None,
    }
    if variant_type in {"STOP_GAINED", "FRAMESHIFT"}:
        outcome = _truncating_path(input_data, services, annotation, scratch)
    elif variant_type in {"CANONICAL_SPLICE", "SPLICE_LOF_CONFIRMED"}:
        outcome = _splice_path(input_data, services, annotation, scratch, rna)
    else:
        outcome = _start_loss_path(input_data, services, annotation, scratch)
    return {
        "eligible_lof_variant": True,
        "variant_type": variant_type,
        # The node the variant-level tree stopped at is the decision path, read from the
        # trace rather than tracked separately so it cannot disagree with it.
        "decision_path": outcome.decision_trace[-1]["node_id"] if outcome.decision_trace else None,
        "candidate_strength": outcome.strength,
        "candidate_applicability": outcome.evaluation_context["applicability"],
        "conclusion_if_mechanism_established": scratch["summary"],
        "unresolved_requirements": list(outcome.unresolved_requirements),
        "decision_trace": outcome.decision_trace,
    }


_WORD = __import__("re").compile(r"[a-z0-9]+")


def _phrase(value):
    """A disease name reduced to its words, for comparing two ways of writing one."""
    return " ".join(_WORD.findall(str(value).lower())) if value else ""


def _stated_diagnosis_match(diagnosis, label):
    """How a candidate's name relates to the diagnosis the note states, or None.

    Text, and only text. ACMG's disease-match levels put name resemblance outside what may
    be applied automatically, so this never becomes a condition and never reaches the disease
    gate - it orders and marks the suggestions a curator is choosing between, which is the one
    place a resemblance is useful and harmless.

    Two relations are reported and nothing in between. The names are the same phrase, or the
    stated one appears whole inside the candidate's ("hypertrophic cardiomyopathy" inside
    "MYBPC3-related hypertrophic cardiomyopathy"). Partial word overlap is not reported at
    all: "cardiomyopathy" is shared by diseases that are not the same disease, and a score
    invented here would be a threshold nobody set.
    """
    stated, candidate = _phrase(diagnosis), _phrase(label)
    if not stated or not candidate:
        return None
    if stated == candidate:
        return "LABEL_EQUAL"
    if f" {stated} " in f" {candidate} ":
        return "LABEL_CONTAINS"
    return None


def _candidate_conditions(input_data, services, annotation):
    """The diseases this gene is curated for, as a list to show a curator - never a choice.

    With no condition supplied there is no disease whose mechanism could be checked, and PVS1
    stops. What it can still say is which diseases the gene has been curated for, because
    that is what a curator needs in order to supply the one this case is about.

    Two kinds of record answer that, and they are not interchangeable. A `gene_disease`
    mechanism assessment says what loss of function does in that disease, and is what would
    let PVS1 proceed. A `gene_disease_validity` curation says only that the gene and the
    disease are related, which is a weaker statement and deliberately never a mechanism - it
    is read here, for suggestions, and nowhere else. Validity records are not filtered through
    reviewed_or_automated() for that reason: they are not assessments to be relied on, and
    holding them to the bar for evidence that decides a criterion would only drop context.

    They are offered and never taken. In ClinGen's curation set 85% of curated genes have a
    single disease, which is exactly the shape that invites picking it automatically - and the
    remaining genes are why that would be wrong: among genes curated for more than one
    disease, mechanisms disagree often enough (ABCC9, ACTB, ATP1A2 and CACNA1D each carry both
    gain and loss of function) that choosing for the curator would sometimes choose the
    mechanism too. A single candidate is not evidence that the case is about that disease
    either; it is evidence about what has been curated.
    """
    gene = annotation.get("gene")
    variant = Variant(**input_data["variant"])
    merged = {}

    def entry(condition):
        return merged.setdefault(condition, {
            "condition": condition, "condition_label": None,
            "lof_mechanism_established": None, "inheritance": None,
            "inheritance_as_recorded": None, "gene_disease_validity": None, "sources": [],
        })

    def label(item, found):
        return found["condition_label"] or item.get("condition_label")

    for item in _candidates("gene_disease", input_data, services, annotation["transcript"]):
        if item.get("gene") != gene or not item.get("condition"):
            continue
        found = entry(item["condition"])
        found["condition_label"] = label(item, found)
        found["sources"].append({
            "source": item.get("source"), "source_version": item.get("source_version"),
            "assessment_method": item.get("assessment_method"),
            "lof_mechanism_established": item.get("lof_mechanism_established"),
        })
        if found["inheritance"] is None:
            found["inheritance"] = normalize_inheritance(item.get("inheritance"))
            found["inheritance_as_recorded"] = item.get("inheritance")

    for item in services.evidence.get_candidates("gene_disease_validity", variant):
        if item.get("gene") != gene or not item.get("condition"):
            continue
        found = entry(item["condition"])
        found["condition_label"] = label(item, found)
        found["gene_disease_validity"] = item.get("classification")
        found["sources"].append({
            "source": item.get("source"), "source_version": item.get("source_version"),
            "classification": item.get("classification"),
        })
        if found["inheritance"] is None:
            found["inheritance"] = normalize_inheritance(item.get("moi"))
            found["inheritance_as_recorded"] = item.get("moi")

    diagnosis = (input_data.get("clinical_note") or {}).get("diagnosis")
    for found in merged.values():
        # Only when the sources that spoke about the mechanism agree. Two that disagree are
        # left unstated rather than resolved here, and both stay visible under `sources`.
        stated = {item["lof_mechanism_established"] for item in found["sources"]
                  if "lof_mechanism_established" in item
                  and item["lof_mechanism_established"] is not None}
        found["lof_mechanism_established"] = stated.pop() if len(stated) == 1 else None
        found["stated_diagnosis_match"] = _stated_diagnosis_match(
            diagnosis, found["condition_label"])

    # What the note named comes first, then what would let PVS1 proceed, then the rest. All
    # three are the curator's to choose from; the order only says which to look at first.
    return sorted(merged.values(),
                  key=lambda item: (item["stated_diagnosis_match"] is None,
                                    item["lof_mechanism_established"] is not True,
                                    item["condition"]))


def _not_evaluated(input_data, services, annotation, variant_type, state, rna, node, summary,
                   missing, *, evidence=(), review=(), candidates=()):
    """Stop short of a PVS1 verdict, keeping the variant-level work that is still valid."""
    state["trace"].append(node)
    return _finish(
        input_data, state, CriterionStatus.UNKNOWN, summary, missing=[missing],
        extra_evidence=evidence, review=review,
        provenance={"preliminary_assessment": _preliminary_assessment(
            input_data, services, annotation, variant_type, state, rna),
            **({"candidate_conditions": list(candidates)} if candidates else {})})


def _evaluate_input(input_data, services, config):
    context = _base_context(input_data)
    warnings = []
    if not input_data.get("condition"):
        warnings.extend([
            "Condition was not provided.",
            "Variant-level loss-of-function assessment was retained as a preliminary result; "
            "PVS1 itself was not evaluated.",
        ])
    state = {
        "context": context,
        "trace": [_node("C01", "condition_status",
                        "PASS" if input_data.get("condition") else "UNKNOWN",
                        context["condition_status"], source="acmg_amp_2015")],
        "evidence": [],
        "warnings": warnings,
        "rules_used": [],
        "ruleset": None,
        "policy": {},
        # Caveats raised partway through the tree (e.g. "SP01 answered from a default
        # policy, not a curator") that must survive into whichever _finish() call
        # actually returns, even several path functions later - see
        # _region_relevance_review()'s and the splice-default caveat's own call sites.
        "pending_review": [],
    }

    early, annotation = annotation_context("PVS1", input_data, services)
    if early:
        early.evaluation_context = context
        early.decision_trace = state["trace"]
        early.warnings = warnings
        early.unresolved_requirements = list(early.missing_inputs)
        return early
    context["gene"] = annotation.get("gene")
    state["evidence"].append(annotation)

    consequences = set(annotation["consequences"])
    rna, rna_issue = _select_context_record(
        "rna_assay", input_data, services, annotation["transcript"])
    if rna_issue == "conflict":
        return _conflicting_record(input_data, state, "rna_assay", rna, "V01", "variant_type")
    confirmed_rna = rna if rna and rna.get("lof_effect_confirmed") is True else None
    if consequences & set(TRUNCATING):
        variant_type = next(TRUNCATING[item] for item in TRUNCATING if item in consequences)
    elif consequences & CANONICAL_SPLICE:
        variant_type = "CANONICAL_SPLICE"
    elif "start_lost" in consequences:
        variant_type = "START_LOST"
    elif confirmed_rna:
        variant_type = "SPLICE_LOF_CONFIRMED"
    elif consequences & OTHER_LOF:
        variant_type = "OTHER_LOF"
    else:
        variant_type = "OTHER"

    if variant_type == "OTHER":
        state["trace"].append(_node("V01", "variant_type", "NOT_APPLICABLE", variant_type,
                                    [annotation]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "PVS1 is not applicable: the annotation does not indicate a truncating, "
                       "canonical-splice, start-loss, or RNA-confirmed loss-of-function consequence.")
    if variant_type == "OTHER_LOF":
        state["trace"].append(_node("V01", "variant_type", "UNKNOWN", variant_type, [annotation]))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "No implemented general PVS1 path for this LoF consequence",
                       missing=["supported_pvs1_variant_path"])

    policy = _rules(config)
    if policy is None:
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Versioned PVS1 rule set unavailable", missing=["PVS1.ruleset"])
    state["policy"] = policy
    state["ruleset"] = {
        "name": policy["ruleset"]["name"], "version": policy["ruleset"]["version"]}
    state["rules_used"] = policy["ruleset"]["sources"]

    # A gene can carry loss of function for one disease and gain of function for another, so
    # the gene alone never settles whether PVS1 applies. Without a condition there is no
    # disease whose mechanism could be checked, and a mechanism record curated for no disease
    # in particular does not become disease-specific by being the only one on file. Both stop
    # short of a verdict; neither discards the variant-level work.
    if not input_data.get("condition"):
        candidates = _candidate_conditions(input_data, services, annotation)
        # Named so a curator can supply the disease context, not so the criterion can pick
        # one - see _candidate_conditions() for why a single candidate is still not a choice.
        named = [item for item in candidates if item["stated_diagnosis_match"]]
        diagnosis = (input_data.get("clinical_note") or {}).get("diagnosis")
        if named:
            # The note named a disease and the gene is curated for something that reads like
            # it. Still a resemblance between two strings, so the curator confirms the
            # identifier rather than the pipeline adopting it.
            review = [f"The note states {diagnosis!r}; confirm whether that is "
                      f"{', '.join(item['condition'] for item in named)} and supply it as "
                      f"the disease context"]
        elif candidates:
            review = [f"Supply the disease context; {annotation.get('gene')} has curated "
                      f"evidence for "
                      f"{', '.join(sorted({item['condition'] for item in candidates}))}"]
        else:
            review = []
        return _not_evaluated(
            input_data, services, annotation, variant_type, state, confirmed_rna,
            _node("D01", "disease_match", "UNKNOWN", "NOT_PROVIDED"),
            "Condition was not provided, so no disease-specific loss-of-function mechanism "
            "could be established",
            "condition", review=review, candidates=candidates)

    mechanism, mechanism_records, mechanism_issue = _resolve_mechanism(
        input_data, services, annotation, context,
        policy["rules"].get("mechanism_source_precedence", []))
    if mechanism_issue == "gene_mismatch":
        state["trace"].append(_node(
            "G01", "lof_mechanism_available", "MANUAL_REVIEW", evidence=mechanism_records))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Gene-mechanism mapping conflicts with the annotation",
                       review=["Resolve gene-disease mechanism"], extra_evidence=mechanism_records)
    if mechanism_issue == "conflict":
        state["trace"].append(_node(
            "G01", "lof_mechanism_available", "MANUAL_REVIEW", evidence=mechanism_records))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Conflicting LoF mechanism assessments",
                       review=["Resolve LoF mechanism assessments"], extra_evidence=mechanism_records)
    if mechanism_issue == "excluded":
        # Not a review point: ClinGen already reviewed this phenotype and kept it out, so
        # asking a curator to decide again would be asking a settled question. The node says
        # NOT_APPLICABLE because this curation is not applicable here; the criterion stays
        # NOT_EVALUATED, because a mechanism for this disease could still come from elsewhere
        # and PVS1 has not been ruled out, only this route to it.
        return _not_evaluated(
            input_data, services, annotation, variant_type, state, confirmed_rna,
            _node("D01", "disease_match", "NOT_APPLICABLE", "EXCLUDED", mechanism_records),
            f"The available loss-of-function mechanism evidence is curated for a disease "
            f"whose expert panel explicitly excluded {input_data['condition']!r} from it",
            "disease-specific loss-of-function mechanism",
            evidence=mechanism_records)
    if mechanism_issue == "other_disease":
        diseases = sorted({str(item["condition"]) for item in mechanism_records
                           if item.get("condition")})
        return _not_evaluated(
            input_data, services, annotation, variant_type, state, confirmed_rna,
            _node("D01", "disease_match", "MANUAL_REVIEW", "OTHER_DISEASE_CURATED",
                  mechanism_records),
            f"Loss-of-function mechanism evidence exists for this gene under "
            f"{', '.join(diseases)}, and none of those is {input_data['condition']!r} or "
            f"relates to it in MONDO",
            "disease-specific loss-of-function mechanism",
            evidence=mechanism_records,
            review=[f"Decide whether the mechanism curated for {', '.join(diseases)} "
                    f"applies to {input_data['condition']!r}"])
    if mechanism_issue == "parent_child":
        return _not_evaluated(
            input_data, services, annotation, variant_type, state, confirmed_rna,
            _node("D01", "disease_match", "MANUAL_REVIEW", "PARENT_CHILD", mechanism_records),
            f"The available loss-of-function mechanism evidence is curated for a disease "
            f"related to {input_data['condition']!r} in MONDO, but not for it",
            "disease-specific loss-of-function mechanism",
            evidence=mechanism_records,
            review=["Decide whether the related disease's mechanism applies to this one"])
    if mechanism_issue == "moi_mismatch":
        modes = sorted({str(item.get("inheritance")) for item in mechanism_records})
        raw = input_data.get("inheritance")
        # An unrecognized mode and an absent one both fail to match, but a curator fixes them
        # differently: one is a spelling this pipeline does not know, the other is a field
        # nobody filled in. The message says which.
        declared = (f"{context['inheritance']!r}" if context["inheritance"]
                    else f"unrecognized {raw!r}" if raw else "unstated")
        # Not borrowed, and not buried either. A mechanism curated for another mode is
        # material a curator can act on, and the two ways of failing to match call for
        # different actions: record the mode this case was assessed under, or decide whether
        # the curated mode's mechanism carries to the one it was. Withholding quietly would
        # be the safer-looking choice and the less useful one - the records are attached, the
        # variant-level tree is reported, and the question goes to the person who can answer
        # it. PVS1 is still not applied.
        question = ("Record the inheritance mode this case was assessed under, or decide "
                    f"whether the mechanism curated for {', '.join(modes)} applies to it"
                    if context["inheritance"] is None else
                    f"Decide whether the mechanism curated for {', '.join(modes)} applies to "
                    f"a case assessed as {declared}")
        return _not_evaluated(
            input_data, services, annotation, variant_type, state, confirmed_rna,
            _node("G01", "lof_mechanism_available", "MANUAL_REVIEW", declared,
                  mechanism_records),
            f"The available loss-of-function mechanism evidence is curated for "
            f"{', '.join(modes)} inheritance; this case's inheritance mode is {declared}",
            "loss-of-function disease mechanism for this inheritance mode",
            evidence=mechanism_records, review=[question])
    if mechanism_issue == "missing":
        state["trace"].append(_node("G01", "lof_mechanism_available", "UNKNOWN"))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Loss-of-function disease mechanism unavailable",
                       missing=["loss-of-function disease mechanism"])
    if context["mechanism_scope"] == "GENE_LEVEL":
        return _not_evaluated(
            input_data, services, annotation, variant_type, state, confirmed_rna,
            _node("D01", "disease_match", "UNKNOWN", context["disease_match"],
                  mechanism_records),
            f"The available loss-of-function mechanism evidence is not curated for "
            f"{input_data['condition']!r}, so no disease-specific mechanism could be "
            f"established",
            "disease-specific loss-of-function mechanism",
            evidence=mechanism_records)
    state["trace"].append(_node(
        "D01", "disease_match", "PASS", context["disease_match"], mechanism_records))
    state["evidence"].extend(mechanism_records)
    state["trace"].append(_node(
        "G01", "lof_mechanism_available", "PASS", True, mechanism_records))
    established = mechanism.get("lof_mechanism_established")
    if established is False:
        state["trace"].append(_node(
            "G02", "lof_mechanism_established", "NOT_APPLICABLE", False, mechanism_records))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Loss of function is not an established disease mechanism")
    if established is not True:
        state["trace"].append(_node(
            "G02", "lof_mechanism_established", "UNKNOWN", established, mechanism_records))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Loss-of-function disease mechanism unresolved",
                       missing=["lof_mechanism_established"])
    state["trace"].append(_node(
        "G02", "lof_mechanism_established", "PASS", True, mechanism_records))
    state["trace"].append(_node(
        "V01", "variant_type", "PASS", variant_type, [annotation],
        "clingen_splicing_2023" if variant_type == "SPLICE_LOF_CONFIRMED" else "clingen_pvs1_2018"))

    if variant_type in {"STOP_GAINED", "FRAMESHIFT"}:
        return _truncating_path(input_data, services, annotation, state)
    if variant_type in {"CANONICAL_SPLICE", "SPLICE_LOF_CONFIRMED"}:
        return _splice_path(input_data, services, annotation, state, confirmed_rna)
    if variant_type == "START_LOST":
        return _start_loss_path(input_data, services, annotation, state)
    raise AssertionError(f"Unhandled PVS1 variant type: {variant_type}")


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    return _evaluate_input(criterion_input(variant, clinical_note), services, config)
