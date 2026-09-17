from acmg_pipeline.constants import CriterionStatus
from acmg_pipeline.automated_core.interface import criterion_input
from acmg_pipeline.criteria.common import result as _base_result
from acmg_pipeline.clinical_note import ClinicalNoteExtraction
from acmg_pipeline.vcf_record import VariantRecord


def result(code, input_data, status, summary, **kwargs):
    assert code == "PP5"
    message = (f"{summary.rstrip('.')}. PP5 is not scored from an external assertion alone; "
               f"the primary evidence must be reviewed. Outcome: UNKNOWN; no MET or NOT_MET judgment is made.")
    return _base_result(code, input_data, status, message, **kwargs)


def evaluate(variant: VariantRecord, clinical_note: ClinicalNoteExtraction, services, config):
    input_data = criterion_input(variant, clinical_note)
    return result("PP5", input_data, CriterionStatus.UNKNOWN,
                  "ClinGen General policy: retrieve primary evidence instead of scoring external assertions")
