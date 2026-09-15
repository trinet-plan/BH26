from acmg.core.models import Status
from acmg.criteria.common import result


def evaluate(input_data, services, config):
    return result("PP5", input_data, Status.DEPRECATED,
                  "ClinGen General policy: retrieve primary evidence instead of scoring external assertions")
