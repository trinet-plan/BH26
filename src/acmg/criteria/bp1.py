from acmg.criteria.mechanism import evaluate_mechanism


def evaluate(input_data, services, config):
    return evaluate_mechanism("BP1", input_data, services, config)
