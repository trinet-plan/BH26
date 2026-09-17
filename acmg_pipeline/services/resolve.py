"""Per-variant evidence resolution, shared by the CLI and the integrated pipeline.

Before this module, `acmg/cli.py`'s `prepare-demo-online` was the only
place that knew how to drive the providers, and the integrated pipeline
took evidence as an argument instead. Both now go through
`VariantProviderSuite` for everything that is decided one variant at a
time (dbNSFP scores, the ClinVar VCV record, the PS1/PM5 comparator
searches, the PM1 hotspot search), so that logic exists once.

[What the CLI still owns, deliberately]
  * Identity resolution. `prepare-demo-online` calls Ensembl to RESOLVE
    each record's coordinates before any Variant exists, and needs the
    returned candidate for `reconcile()`. That is an identity step, not an
    evidence step.
  * The gnomAD batch. `GnomadProvider.get_frequencies` issues ONE GraphQL
    request covering every variant in the run. Splitting it into
    per-variant `get_frequency` calls would change the request body, and
    the HTTP cache is keyed on that body - the committed offline fixture
    cache holds the batched response, and `cache_set_sha256` in the
    identity manifest pins the exact cache-entry set. Per-variant calls
    would miss the cache entirely and move that hash, so the batch stays.
  `ProviderEvidenceResolver` therefore does those two itself, per variant,
  for the single-variant case the integrated pipeline has.

[Why the caller no longer supplies evidence]
  An evidence record is only usable if it carries `evidence_id` (a
  resolvable IRI), `source`, `source_version`, `retrieved_at` and
  `quality_status` - `acmg_pipeline.services.evidence.EvidenceService` filters on
  all five, and `acmg_pipeline.automated_va_spec.evidence_reference()` refuses to
  export an item without an IRI. Only a provider can supply those. In
  particular they cannot be read off the VCF INFO column: the demo VCFs
  carry identity/context (GENE/TRANSCRIPT/HGVSC/HGVSP/ZYGOSITY/
  DISEASE_ASSOCIATION) plus already-decided answers (CLNSIG, ACMG_CODES)
  and uncalibrated scores (AM_*/AG_*), but no population frequency and no
  source version for anything. INFO stays identity and context only.

[Absent is not zero]
  A variant missing from gnomAD yields no population record at all, so
  PM2/BA1/BS1 report NOT_EVALUATED ("no reliable population observation")
  rather than being handed AF=0.

[Failures are reported, not swallowed]
  One provider being unreachable must not look like "this variant has no
  such evidence". Every call returns its error alongside its records, so
  the CLI can put it in that record's `issues` and the pipeline can keep
  it in `ResolvedEvidence.failures`.
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone

from acmg_pipeline.automated_core.models import Variant
from acmg_pipeline.providers.clinvar import (
    VCV, ClinVarComparatorProvider, ClinVarHotspotProvider, ClinVarProvider,
)
from acmg_pipeline.providers.dbnsfp import DbnsfpProvider
from acmg_pipeline.providers.ensembl import EnsemblIdentityProvider
from acmg_pipeline.providers.gnomad import GnomadProvider
from acmg_pipeline.providers.http import CachedHttpClient, FetchError

PROVIDER_ERRORS = (ValueError, FetchError)


@dataclass
class ProviderOutcome:
    """One provider call: what it produced, and why it produced nothing."""

    records: list[dict] = field(default_factory=list)
    error: str | None = None
    # search_ps1/search_pm5/search_hotspot each return a query record that
    # documents the search itself plus zero or more matches. The counts the
    # CLI manifest reports distinguish the two, so they are kept apart.
    searched: bool = False
    matches: int = 0


@dataclass
class ResolvedEvidence:
    """Normalized evidence for one variant, plus whatever could not be fetched."""

    records: list[dict] = field(default_factory=list)
    failures: list[dict] = field(default_factory=list)

    def absorb(self, provider: str, outcome: ProviderOutcome) -> ProviderOutcome:
        self.records.extend(outcome.records)
        if outcome.error is not None:
            self.failures.append({"provider": provider, "error": outcome.error})
        return outcome


def splice_score_for(predictions, transcript):
    """The SpliceAI score the PS1/PM5 comparators use to gate a splicing effect."""
    return next(
        (item.get("score") for item in predictions
         if item.get("predictor") == "SpliceAI" and item.get("transcript") == transcript),
        None,
    )


class VariantProviderSuite:
    """The provider calls that are made once per variant, in one place.

    Constructed with an already-configured HTTP client so the caller keeps
    control of caching and offline mode, and with the Ensembl provider the
    comparator needs for its own lookups.
    """

    def __init__(self, client, *, clinvar_release, ensembl_provider, hotspot_policy=None):
        self._client = client
        self._clinvar_release = clinvar_release
        self._ensembl = ensembl_provider
        self._hotspot_policy = hotspot_policy
        # One instance for the whole run: DbnsfpProvider caches the release
        # it discovered in `.metadata`, and the manifest reports that as
        # `provider_version`. A fresh instance per variant would throw it
        # away and leave the release unrecorded - which matters, because
        # PP3/BP4 calibration is only valid for the dbNSFP release that
        # produced the scores.
        self.dbnsfp = DbnsfpProvider(client)

    def dbnsfp_version(self):
        return (self.dbnsfp.metadata or {}).get("dbnsfp_version")

    def predictions(self, variant, transcript=None) -> ProviderOutcome:
        try:
            return ProviderOutcome(records=list(self.dbnsfp.get_predictions(variant, transcript)))
        except PROVIDER_ERRORS as exc:
            return ProviderOutcome(error=str(exc))

    def clinvar_record(self, accession, variant) -> tuple[ProviderOutcome, dict | None]:
        """Returns (outcome, identity evidence) - the latter belongs to the record, not the evidence set."""
        if not accession or not VCV.fullmatch(str(accession)):
            return ProviderOutcome(), None
        provider = ClinVarProvider(self._client, self._clinvar_release)
        try:
            item, identity = provider.get_record(str(accession), variant)
        except PROVIDER_ERRORS as exc:
            return ProviderOutcome(error=str(exc)), None
        return ProviderOutcome(records=[item]), identity

    def comparator(self, criterion, annotation, variant, splice_score) -> ProviderOutcome:
        """PS1 (same amino acid change) or PM5 (same residue, different change)."""
        provider = ClinVarComparatorProvider(
            self._client, self._clinvar_release, self._ensembl
        )
        search = {"PS1": provider.search_ps1, "PM5": provider.search_pm5}[criterion]
        try:
            query, matches = search(annotation, variant, splice_score)
        except PROVIDER_ERRORS as exc:
            return ProviderOutcome(error=str(exc))
        return ProviderOutcome(
            records=[query, *matches], searched=True, matches=len(matches)
        )

    def hotspot(self, annotation, variant) -> ProviderOutcome:
        provider = ClinVarHotspotProvider(
            self._client, self._clinvar_release, self._hotspot_policy
        )
        try:
            query, region = provider.search_hotspot(annotation, variant)
        except PROVIDER_ERRORS as exc:
            return ProviderOutcome(error=str(exc))
        records = [query] if region is None else [query, region]
        return ProviderOutcome(records=records, searched=True, matches=int(region is not None))


class StaticEvidenceResolver:
    """Return a fixed evidence list. For tests and for replaying a cached run.

    Keeps the integrated entry point hermetic without reintroducing an
    evidence argument on it: the seam moves from "every caller passes
    evidence" to "a test injects a resolver".
    """

    def __init__(self, records: list[dict]):
        if not isinstance(records, list) or not all(isinstance(r, dict) for r in records):
            raise TypeError("records must be a list of dictionaries")
        self._records = records

    def resolve(self, identity: dict, variant: Variant) -> ResolvedEvidence:
        return ResolvedEvidence(records=list(self._records))


class ProviderEvidenceResolver:
    """Query Ensembl, gnomAD, dbNSFP and ClinVar for one variant.

    Every request goes through `CachedHttpClient`, so a run against a
    populated cache directory is reproducible offline and byte-identical -
    pass `offline=True` to make a cache miss an error instead of a network
    call.
    """

    def __init__(
        self,
        cache_dir,
        *,
        offline: bool = False,
        ensembl_release: str | None = None,
        gnomad_release: str = "4.1.1",
        clinvar_release: str | None = None,
        evidence_cache_dir=None,
        hotspot_policy: dict | None = None,
    ):
        self._client = CachedHttpClient(cache_dir, offline=offline)
        self._external = (
            CachedHttpClient(evidence_cache_dir, offline=offline)
            if evidence_cache_dir else self._client
        )
        self._ensembl_release = ensembl_release
        self._gnomad_release = gnomad_release
        self._clinvar_release = (
            clinvar_release or datetime.now(timezone.utc).date().isoformat()
        )
        self._hotspot_policy = hotspot_policy
        self._ensembl = None

    def _ensembl_provider(self) -> EnsemblIdentityProvider:
        if self._ensembl is None:
            release = self._ensembl_release or EnsemblIdentityProvider.current_release(self._client)
            self._ensembl = EnsemblIdentityProvider(self._client, release)
        return self._ensembl

    def _suite(self) -> VariantProviderSuite:
        return VariantProviderSuite(
            self._external,
            clinvar_release=self._clinvar_release,
            ensembl_provider=self._ensembl_provider(),
            hotspot_policy=self._hotspot_policy,
        )

    def resolve(self, identity: dict, variant: Variant) -> ResolvedEvidence:
        """`identity` is the VCF INFO view: GENE/TRANSCRIPT/HGVSC/CLNVARIATIONID."""
        resolved = ResolvedEvidence()
        annotation, predictions = self._annotate(identity, variant, resolved)
        if annotation is not None:
            resolved.records.append(annotation)
            resolved.records.extend(predictions)
        self._add_population(variant, resolved)

        suite = self._suite()
        transcript = annotation.get("transcript") if annotation else None
        resolved.absorb(DbnsfpProvider.name, suite.predictions(variant, transcript))
        outcome, _identity_evidence = suite.clinvar_record(
            identity.get("CLNVARIATIONID"), variant
        )
        resolved.absorb(ClinVarProvider.name, outcome)

        if annotation is None or "missense_variant" not in (annotation.get("consequences") or []):
            return resolved
        # PS1/PM5/PM1 all compare protein-level changes, so they are only
        # meaningful for a missense annotation.
        splice_score = splice_score_for(predictions, annotation.get("transcript"))
        for criterion in ("PS1", "PM5"):
            resolved.absorb(
                f"{ClinVarComparatorProvider.name}:{criterion}",
                suite.comparator(criterion, annotation, variant, splice_score),
            )
        if self._hotspot_policy and annotation.get("protein_start"):
            resolved.absorb(ClinVarHotspotProvider.name, suite.hotspot(annotation, variant))
        return resolved

    def _annotate(self, identity, variant, resolved):
        record = {"identity": identity}
        try:
            _candidate, annotation, predictions = (
                self._ensembl_provider().map_record_with_evidence(record)
            )
        except PROVIDER_ERRORS as exc:
            resolved.failures.append({"provider": "Ensembl", "error": str(exc)})
            return None, []
        # The transcript HGVS the VCF names and the coordinates the caller
        # is evaluating must describe the same variant. A silent mismatch
        # here would attach another variant's annotation and predictions to
        # this one, so it is a failure, never a warning.
        if annotation.get("variant_key") != variant.key:
            resolved.failures.append({
                "provider": "Ensembl",
                "error": (
                    f"HGVS {identity.get('TRANSCRIPT')}:{identity.get('HGVSC')} resolves to "
                    f"{annotation.get('variant_key')}, not {variant.key}"
                ),
            })
            return None, []
        return annotation, predictions

    def _add_population(self, variant, resolved):
        provider = GnomadProvider(self._external, release=self._gnomad_release)
        try:
            observations = provider.get_frequency(variant)
        except PROVIDER_ERRORS as exc:
            resolved.failures.append({"provider": provider.name, "error": str(exc)})
            return
        if observations:  # None = not in gnomAD; contribute nothing, never AF=0
            resolved.records.extend(observations)
