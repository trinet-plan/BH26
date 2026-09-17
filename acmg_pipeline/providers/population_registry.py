"""Configuration-driven registry for independent population-frequency databases."""

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


POPULATION_PROVIDER_FACTORIES = {
    "togovar": _togovar,
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
