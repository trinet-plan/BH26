from acmg.core.interface import criterion_input
from acmg.criteria.regions import evaluate_region
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    return evaluate_region("PM4", input_data, services, config)
