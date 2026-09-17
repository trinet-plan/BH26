from dataclasses import replace

from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.comparator import evaluate_comparator
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def _summary(result):
    outcome = {"met": "the available evidence satisfies PM5.", "not_met": "PM5 was evaluated but its requirements were not satisfied."}.get(result.status.value, "no MET or NOT_MET judgment is made from the available information.")
    return replace(result, summary=(f"{result.summary.rstrip('.')}. PM5 compares a missense variant with an independently established pathogenic variant causing a different amino-acid substitution at the same residue. Outcome: {result.status.value.upper()}; {outcome}"))


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    return _summary(evaluate_comparator("PM5", input_data, services, config))
