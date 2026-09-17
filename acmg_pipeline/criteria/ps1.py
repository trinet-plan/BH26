from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.comparator import evaluate_comparator
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    return evaluate_comparator("PS1", input_data, services, config)
