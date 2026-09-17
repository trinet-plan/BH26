from acmg_pipeline.constants import CriterionStatus
"""ClinGen general PVS1 decision tree for sequence variants.

Evidence acquisition is intentionally separate. This evaluator consumes versioned, reviewed
transcript, NMD, splice, protein-region, population-LoF, and initiation assessments and never
turns a missing assessment into a negative conclusion.
"""

from copy import deepcopy

from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.criteria.common import annotation_context, result as _base_result, reviewed_or_automated
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


def _resolve_mechanism(input_data, services, annotation, context):
    records = _candidates("gene_disease", input_data, services, annotation["transcript"])
    gene = annotation.get("gene")
    gene_records = [item for item in records if item.get("gene") == gene]
    if records and not gene_records:
        return None, records, "gene_mismatch"
    condition = input_data.get("condition")
    if condition:
        exact = [item for item in gene_records if item.get("condition") == condition]
        selected = exact or [item for item in gene_records if not item.get("condition")]
        scope = "CONDITION_SPECIFIC" if exact else "GENE_LEVEL"
    else:
        selected = [item for item in gene_records if not item.get("condition")]
        scope = "GENE_LEVEL"
    if not selected:
        return None, [], "missing"
    values = {item.get("lof_mechanism_established") for item in selected}
    if len(values) > 1:
        return None, selected, "conflict"
    context["mechanism_scope"] = scope
    context["condition_specific"] = scope == "CONDITION_SPECIFIC"
    return selected[0], selected, None


def _deduplicate(items):
    values = {}
    for item in items:
        values[item.get("evidence_id", id(item))] = item
    return list(values.values())


def _finish(input_data, state, status, summary, *, strength=None, missing=(), review=(),
            extra_evidence=(), provenance=None):
    unresolved = list(missing)
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


def _evaluate_input(input_data, services, config):
    context = _base_context(input_data)
    warnings = []
    if not input_data.get("condition"):
        warnings.extend([
            "Condition was not provided.",
            "PVS1 was evaluated using gene-level loss-of-function mechanism evidence.",
        ])
    state = {
        "context": context,
        "trace": [_node("C01", "condition_status", "PASS", context["condition_status"],
                        source="acmg_amp_2015")],
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

    mechanism, mechanism_records, mechanism_issue = _resolve_mechanism(
        input_data, services, annotation, context)
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
    if mechanism_issue == "missing":
        state["trace"].append(_node("G01", "lof_mechanism_available", "UNKNOWN"))
        return _finish(input_data, state, CriterionStatus.UNKNOWN,
                       "Loss-of-function disease mechanism unavailable",
                       missing=["loss-of-function disease mechanism"])
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
