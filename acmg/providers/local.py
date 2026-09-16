"""Normalized local evidence contracts used by offline runs and provider tests."""


class LocalPopulationProvider:
    def __init__(self, name, observations):
        self.name = name
        self.observations = observations

    def get_frequency(self, variant, context=None):
        matches = [item for item in self.observations if item.get("variant_key") == variant.key]
        return matches or None


class ProviderRegistry:
    def __init__(self):
        self.providers = {}

    def register(self, category, name, provider):
        if (category, name) in self.providers:
            raise ValueError(f"Duplicate provider: {category}/{name}")
        self.providers[category, name] = provider

    def get_enabled(self, category):
        return [provider for (kind, _), provider in self.providers.items() if kind == category]
