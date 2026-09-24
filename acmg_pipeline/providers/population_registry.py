"""Configuration-driven registry for independent population-frequency databases."""

from acmg_pipeline.providers.gnomad import GnomadProvider
from acmg_pipeline.providers.togovar import API_VERSION, TogoVarProvider


DEFAULT_POPULATION_SOURCES = [{
    "provider": "togovar",
    "api_version": API_VERSION,
    "frequency_sources": ["gnomad", "tommo"],
}]


def _togovar(client, definition):
    allowed = {"provider", "api_version", "frequency_sources"}
    unknown = set(definition) - allowed
    if unknown:
        raise ValueError(f"Unknown TogoVar provider settings: {sorted(unknown)}")
    return TogoVarProvider(
        client,
        api_version=definition.get("api_version", API_VERSION),
        frequency_sources=definition.get("frequency_sources"),
    )


def _gnomad(client, definition):
    """The only source in this registry that reports per-subpopulation observations
    (GnomadProvider.PRIMARY_POPULATIONS), not just one "global" figure per dataset -
    TogoVar's own API has no subpopulation breakdown at all (confirmed against real cached
    responses, 2026-09-22: every TogoVar frequency entry is tagged "<dataset>:global").
    Most real ClinGen VCEP BA1/BS1 specifications are calibrated against gnomAD's popmax
    (the highest observed subpopulation frequency), which TogoVar-only sourcing can never
    supply - see criteria/ba1.py's/bs1.py's own gene-specific thresholds. observed_
    frequencies()'s own max() over every resolved observation already does the "popmax"
    selection once gnomAD's per-population observations are in the pool; no criterion-side
    change is needed, only this registry entry and a population_sources list that includes
    both "togovar" and "gnomad".
    """
    allowed = {"provider", "dataset", "release"}
    unknown = set(definition) - allowed
    if unknown:
        raise ValueError(f"Unknown gnomAD provider settings: {sorted(unknown)}")
    kwargs = {}
    if "dataset" in definition:
        kwargs["dataset"] = definition["dataset"]
    if "release" in definition:
        kwargs["release"] = definition["release"]
    return GnomadProvider(client, **kwargs)


POPULATION_PROVIDER_FACTORIES = {
    "togovar": _togovar,
    "gnomad": _gnomad,
}


def register_population_provider(name, factory):
    """Register one reviewed adapter factory during application startup."""
    if not isinstance(name, str) or not name or not callable(factory):
        raise ValueError("A provider registration requires a name and callable factory")
    if name in POPULATION_PROVIDER_FACTORIES:
        raise ValueError(f"Population provider already registered: {name}")
    POPULATION_PROVIDER_FACTORIES[name] = factory


def build_population_providers(client, definitions=None):
    """Build enabled providers from JSON-compatible configuration.

    A separate frequency database needs an adapter implementing ``name`` and
    ``get_frequency(variant)`` plus an explicit registry entry. Dynamic imports
    are deliberately avoided so configuration cannot execute arbitrary code.
    """
    definitions = DEFAULT_POPULATION_SOURCES if definitions is None else definitions
    if not isinstance(definitions, list) or not definitions:
        raise ValueError("population_sources must be a nonempty list")
    providers = []
    configured = set()
    for definition in definitions:
        if not isinstance(definition, dict):
            raise ValueError("Each population source must be an object")
        provider_name = definition.get("provider")
        factory = POPULATION_PROVIDER_FACTORIES.get(provider_name)
        if factory is None:
            raise ValueError(f"Unknown population provider: {provider_name!r}")
        if provider_name in configured:
            raise ValueError(f"Duplicate population provider: {provider_name}")
        configured.add(provider_name)
        provider = factory(client, definition)
        if (not isinstance(getattr(provider, "name", None), str)
                or not callable(getattr(provider, "get_frequency", None))):
            raise TypeError(f"Invalid population provider adapter: {provider_name}")
        providers.append(provider)
    return providers
