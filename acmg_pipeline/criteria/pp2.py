from dataclasses import replace

from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.mechanism import evaluate_mechanism
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def _summary(result):
    outcome = {"met": "the available evidence satisfies PP2.", "not_met": "PP2 was evaluated but its requirements were not satisfied."}.get(result.status.value, "no MET or NOT_MET judgment is made from the available information.")
    return replace(result, summary=(f"{result.summary.rstrip('.')}. PP2 evaluates missense variants in genes where missense is an established disease mechanism and benign missense variation is constrained. Outcome: {result.status.value.upper()}; {outcome}"))


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    return _summary(evaluate_mechanism("PP2", input_data, services, config))
