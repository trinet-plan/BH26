from acmg.criteria.mechanism import evaluate_mechanism


def evaluate(input_data, services, config):
    return evaluate_mechanism("PP2", input_data, services, config)
