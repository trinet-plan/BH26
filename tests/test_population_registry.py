import pytest

from acmg_pipeline.providers.gnomad import GnomadProvider
from acmg_pipeline.providers.population_registry import (
    POPULATION_PROVIDER_FACTORIES,
    build_population_providers,
)
from acmg_pipeline.providers.togovar import TogoVarProvider


def test_gnomad_is_registered_alongside_togovar():
    """Added 2026-09-22: TogoVar's own API never reports a per-subpopulation
    frequency (only one "global" figure per dataset - confirmed against real cached
    responses), so a population_sources list combining "togovar" and "gnomad" is what
    lets BA1/BS1's gene-specific VCEP thresholds (most of which are popmax-calibrated)
    see a real popmax candidate - see population_registry._gnomad()'s own docstring."""
    providers = build_population_providers(object(), [
        {"provider": "togovar", "api_version": "0.9.1", "frequency_sources": ["gnomad"]},
        {"provider": "gnomad", "release": "4.1.1"},
    ])
    assert [provider.name for provider in providers] == ["TogoVar", "gnomAD"]
    assert isinstance(providers[1], GnomadProvider)
    assert providers[1].release == "4.1.1"


def test_gnomad_factory_rejects_an_unknown_setting():
    with pytest.raises(ValueError, match="Unknown gnomAD provider settings"):
        build_population_providers(object(), [{"provider": "gnomad", "bogus": True}])


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
