from acmg.core.interface import criterion_input
from acmg.core.models import Status
from acmg.criteria.common import result
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    return result("BP6", input_data, Status.DEPRECATED,
                  "ClinGen General policy: retrieve primary evidence instead of scoring external assertions")
