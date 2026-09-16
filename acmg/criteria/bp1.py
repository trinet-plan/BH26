from acmg.core.interface import criterion_input
from acmg.criteria.mechanism import evaluate_mechanism
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    return evaluate_mechanism("BP1", input_data, services, config)
