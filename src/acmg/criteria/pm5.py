from acmg.criteria.comparator import evaluate_comparator


def evaluate(input_data, services, config):
    return evaluate_comparator("PM5", input_data, services, config)
