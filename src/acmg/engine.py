"""Independent criteria plus non-scoring conflict annotations."""

from importlib import import_module
from types import SimpleNamespace

from acmg.core.models import CRITERIA, Status, Variant
from acmg.providers.local import LocalPopulationProvider
from acmg.services.evidence import EvidenceService
from acmg.services.population import PopulationService


CONFLICTS = (
    ("PVS1", "PM4"), ("PVS1", "PP3"), ("PS1", "PP3"), ("PS1", "PM5"),
    ("PM4", "BP3"), ("PP3", "BP4"), ("PP3", "BP7"), ("BP4", "BP7"),
    ("PP2", "BP1"), ("PM2", "BA1"), ("PM2", "BS1"), ("BA1", "BS1"),
)


def make_services(evidence, population_providers=None):
    observations = [item for item in evidence if item.get("category") == "population"]
    names = population_providers if population_providers is not None else sorted({
        item.get("source", "unknown") for item in observations})
    providers = [LocalPopulationProvider(name, [o for o in observations if o.get("source") == name])
                 for name in names]
    return SimpleNamespace(population=PopulationService(providers), evidence=EvidenceService(evidence))


def evaluate_record(input_data, services, config, criteria=CRITERIA):
    Variant(**input_data["variant"])
    if len(criteria) != len(set(criteria)) or any(code not in CRITERIA for code in criteria):
        raise ValueError("Criteria must be unique supported codes")
    results = [import_module(f"acmg.criteria.{code.lower()}").evaluate(input_data, services, config)
               for code in criteria]
    by_code = {r.criterion: r for r in results}
    for left, right in CONFLICTS:
        pair = [by_code.get(left), by_code.get(right)]
        if all(r is not None and r.status == Status.MET for r in pair):
            flag = f"REVIEW_OVERLAP:{left}:{right}"
            for value in pair:
                value.conflict_flags.append(flag)
    return results
