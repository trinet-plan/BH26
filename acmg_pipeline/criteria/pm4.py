from dataclasses import replace

from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.regions import evaluate_region
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def _summary(result):
    outcome = {"met": "the available evidence satisfies PM4.", "not_met": "PM4 was evaluated but its requirements were not satisfied."}.get(result.status.value, "no MET or NOT_MET judgment is made from the available information.")
    return replace(result, summary=(f"{result.summary.rstrip('.')}. PM4 evaluates protein-length changes caused by an in-frame insertion, in-frame deletion, or stop-loss variant. Outcome: {result.status.value.upper()}; {outcome}"))


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    return _summary(evaluate_region("PM4", input_data, services, config))
