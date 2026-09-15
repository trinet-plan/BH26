from acmg.criteria.regions import evaluate_region


def evaluate(input_data, services, config):
    return evaluate_region("BP3", input_data, services, config)
