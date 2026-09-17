import pytest

from acmg_pipeline.providers.population_registry import (
    POPULATION_PROVIDER_FACTORIES,
    build_population_providers,
)
from acmg_pipeline.providers.togovar import TogoVarProvider


def test_default_registry_builds_togovar_with_gnomad_and_tommo_groups():
    providers = build_population_providers(object())
    assert len(providers) == 1
    assert isinstance(providers[0], TogoVarProvider)
    assert providers[0].frequency_sources == ("gnomad", "tommo")


def test_configuration_builds_a_separate_registered_frequency_database(monkeypatch):
    class OtherFrequencyDb:
        name = "OtherFrequencyDB"

        def __init__(self, client, release):
            self.client = client
            self.release = release

        def get_frequency(self, variant):
            return None

    def factory(client, definition):
        if set(definition) != {"provider", "release"}:
            raise ValueError("Unexpected OtherFrequencyDB settings")
        return OtherFrequencyDb(client, definition["release"])

    monkeypatch.setitem(POPULATION_PROVIDER_FACTORIES, "other_frequency_db", factory)
    client = object()
    providers = build_population_providers(client, [
        {
            "provider": "togovar",
            "api_version": "0.9.1",
            "frequency_sources": ["gnomad", "tommo"],
        },
        {"provider": "other_frequency_db", "release": "2026.1"},
    ])

    assert [provider.name for provider in providers] == ["TogoVar", "OtherFrequencyDB"]
    assert providers[1].client is client
    assert providers[1].release == "2026.1"


@pytest.mark.parametrize("definitions, message", [
    ([], "nonempty"),
    ([{"provider": "unknown"}], "Unknown population provider"),
    ([{"provider": "togovar"}, {"provider": "togovar"}], "Duplicate"),
    (["togovar"], "must be an object"),
])
def test_invalid_population_source_configuration_is_rejected(definitions, message):
    with pytest.raises(ValueError, match=message):
        build_population_providers(object(), definitions)
