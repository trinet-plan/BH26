from acmg.criteria.computational import evaluate_prediction


def evaluate(input_data, services, config):
    return evaluate_prediction("BP4", input_data, services, config)
