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
  * Batch preparation. The CLI owns batching across demo variants; the
    integrated pipeline resolves one variant at a time through TogoVar.

[Why the caller no longer supplies evidence]
  An evidence record is only usable if it carries `evidence_id` (a
  resolvable IRI), `source`, `source_version`, `retrieved_at` and
  `quality_status` - `acmg_pipeline.services.evidence.EvidenceService` filters on
  all five, and `acmg_pipeline.automated_va_spec.evidence_reference()` refuses to
  export an item without an IRI. Only a provider can supply those. In
  particular they cannot be read off the demo VCFs' INFO column: those
  files carry variant identity plus fields such as ZYGOSITY/DISEASE_ASSOCIATION,
  but disease context is now accepted only from the clinical-note parser; they also carry
  already-decided answers (CLNSIG, ACMG_CODES)
  and uncalibrated scores (AM_*/AG_*), but no population frequency and no
  source version for anything.

  That was about those files, not about VCF as a format, and it still
  holds for them. A file written by VEP or bcftools annotate carries the
  same frequencies and consequences the providers fetch, and names the
  tool and release that produced them in its headers. Where it does,
  providers/vcf_annotation.py reads those fields under this same contract
  and drops whatever cannot meet it - a field with no pinnable release is
  skipped and fetched instead, a curator's conclusion cannot be mapped at
  all, and a predictor score keeps the version that decides whether a
  calibration accepts it. It is opt-in, and every record it produces
  carries origin="vcf_input", because a file handed to us with the variant
  is not the same evidence as a source this pipeline queried.

[Absent is not zero]
  A variant missing from TogoVar yields no population record at all, so
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
from acmg_pipeline.providers.clingen_dosage import (
    METHOD as DOSAGE_METHOD, ClinGenDosageProvider,
)
from acmg_pipeline.providers.clingen_gene_validity import ClinGenGeneValidityProvider
from acmg_pipeline.providers.clingen_lumping import ClinGenLumpingProvider
from acmg_pipeline.providers.cspec_applicability import (
    METHOD as CSPEC_METHOD, CSpecApplicabilityProvider,
)
from acmg_pipeline.providers.gene2phenotype import (
    METHOD as G2P_METHOD, Gene2PhenotypeProvider,
)
from acmg_pipeline.providers.mondo import MondoMappingProvider
from acmg_pipeline.providers.mondo_hierarchy import MondoHierarchyProvider
from acmg_pipeline.providers.clinvar_spectrum import ClinvarSpectrumProvider
from acmg_pipeline.providers.gene_disease_draft import GeneDiseaseDraftProvider
from acmg_pipeline.providers.gnomad_constraint import GnomadConstraintProvider
from acmg_pipeline.providers.initiation import InitiationProvider
from acmg_pipeline.providers.mane import ManeTranscriptProvider
from acmg_pipeline.providers.nmd import TRUNCATING as NMD_TRUNCATING, NmdPredictionProvider
from acmg_pipeline.providers.protein_region import ProteinRegionProvider
from acmg_pipeline.providers.splice_default import SpliceDefaultProvider
from acmg_pipeline.providers.upstream_pathogenic import UpstreamPathogenicProvider
from acmg_pipeline.providers.dbnsfp import DbnsfpProvider
from acmg_pipeline.providers.ensembl import EnsemblIdentityProvider
from acmg_pipeline.providers.http import CachedHttpClient, FetchError
from acmg_pipeline.providers.population_registry import build_population_providers

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
    """Normalized evidence for one variant, what could not be fetched, and from where.

    `manifest` is one entry per provider that ran. It is not logging: `use_restriction`
    records what the evidence may be used for (a dosage score is PVS1's mechanism gate and
    nothing else; a gene-disease association is a suggestion and never a mechanism), and
    `unresolved` records what was asked for and came back empty, which a caller cannot
    otherwise tell from what was never asked. Both were recorded on the CLI's own path and
    would have been lost by moving to this one.
    """

    records: list[dict] = field(default_factory=list)
    failures: list[dict] = field(default_factory=list)
    manifest: list[dict] = field(default_factory=list)
    # Provenance for the identity itself (which ClinVar record this variant was matched to),
    # which belongs on the audit record rather than in the evidence.
    identity_evidence: list[dict] = field(default_factory=list)

    def absorb(self, provider: str, outcome: ProviderOutcome) -> ProviderOutcome:
        self.records.extend(outcome.records)
        if outcome.error is not None:
            self.failures.append({"provider": provider, "error": outcome.error})
        return outcome

    def note(self, provider, *, use_restriction, records=(), error=None, **details):
        """Record that `provider` ran, what it produced, and what that may be used for."""
        entry = {"provider": provider, "evidence": len(records),
                 "use_restriction": use_restriction,
                 "errors": [error] if error else []}
        entry.update({key: value for key, value in details.items() if value is not None})
        existing = next((item for item in self.manifest
                         if item["provider"] == provider), None)
        if existing is None:
            self.manifest.append(entry)
            return
        # One resolver instance answers many variants on the batch path, and a manifest that
        # reported only the last of them would understate what ran.
        existing["evidence"] += entry["evidence"]
        existing["errors"].extend(entry["errors"])
        for key, value in details.items():
            if isinstance(value, int) and isinstance(existing.get(key), int):
                existing[key] += value
            elif isinstance(value, list):
                existing[key] = _extend_unique(existing.get(key) or [], value)
            elif value is not None:
                existing.setdefault(key, value)


def _extend_unique(existing, addition):
    """Append what is not already there, for lists whose items need not be hashable."""
    combined = list(existing)
    for item in addition:
        if item not in combined:
            combined.append(item)
    return combined


def merge_manifests(manifests):
    """Combine per-variant provider manifests into one run-level manifest.

    A resolver answers one variant at a time, so a batch produces one manifest per variant
    and reporting only the last would understate what ran. Counts add up, lists are unioned,
    and the constants each provider states about itself - its version, its method, what the
    evidence may be used for - are kept as they are.
    """
    combined = {}
    for manifest in manifests:
        for entry in manifest:
            existing = combined.get(entry["provider"])
            if existing is None:
                combined[entry["provider"]] = {**entry,
                                               "errors": list(entry.get("errors") or [])}
                continue
            for key, value in entry.items():
                if key == "provider":
                    continue
                if isinstance(value, bool) or not isinstance(value, (int, list)):
                    existing.setdefault(key, value)
                elif isinstance(value, list):
                    existing[key] = _extend_unique(existing.get(key) or [], value)
                else:
                    existing[key] = (existing.get(key) or 0) + value
    return list(combined.values())


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
    """Query Ensembl, TogoVar, dbNSFP and ClinVar for one variant.

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
        population_sources: list[dict] | None = None,
        clinvar_release: str | None = None,
        evidence_cache_dir=None,
        hotspot_policy: dict | None = None,
        with_clingen_dosage: bool = False,
        with_gene_disease_draft: bool = False,
        gene_disease_draft_policy: dict | None = None,
        with_clinvar_spectrum: bool = False,
        with_splice_default: bool = False,
        splice_default_policy_version: str | None = None,
        with_initiation_assessment: bool = False,
        with_pvs1_transcript_gates: bool = False,
        with_gene2phenotype: bool = False,
        cspec_applicability_path=None,
        with_disease_matching: bool = False,
        with_gene_disease_associations: bool = False,
    ):
        self._client = CachedHttpClient(cache_dir, offline=offline)
        self._external = (
            CachedHttpClient(evidence_cache_dir, offline=offline)
            if evidence_cache_dir else self._client
        )
        self._ensembl_release = ensembl_release
        self._population_sources = population_sources
        self._population = None
        self._clinvar_release = (
            clinvar_release or datetime.now(timezone.utc).date().isoformat()
        )
        self._hotspot_policy = hotspot_policy
        self._with_clingen_dosage = with_clingen_dosage
        self._with_gene_disease_draft = with_gene_disease_draft
        self._gene_disease_draft_policy = gene_disease_draft_policy
        self._with_clinvar_spectrum = with_clinvar_spectrum
        self._with_splice_default = with_splice_default
        self._splice_default_policy_version = splice_default_policy_version
        self._with_initiation_assessment = with_initiation_assessment
        self._with_pvs1_transcript_gates = with_pvs1_transcript_gates
        self._with_gene2phenotype = with_gene2phenotype
        # A path rather than a boolean: this producer reads a collected snapshot the project
        # commits, so a run that wants it has to say which snapshot it means.
        self._cspec_applicability_path = cspec_applicability_path
        self._cspec_applicability = None
        self._with_gene_disease_associations = with_gene_disease_associations
        # One switch for the three sources that only ever refine the disease match: mapping
        # OMIM/Orphanet to MONDO, MONDO ancestry, and ClinGen's lumping decisions. They
        # answer one question between them and are useless apart, so they are asked for once.
        self._with_disease_matching = with_disease_matching
        self._mondo = None
        self._mane = None
        self._clingen_gene_validity = None
        self._gnomad_constraint = None
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

    def _population_providers(self):
        if self._population is None:
            # An empty list is "no population sources", distinct from None, which means the
            # registry's own defaults. A caller that supplies its own frequencies, or wants
            # none, has to be able to say so without the defaults firing behind it.
            self._population = (
                [] if self._population_sources == []
                else build_population_providers(self._external, self._population_sources)
            )
        return self._population

    @property
    def identity_provider(self) -> EnsemblIdentityProvider:
        """The Ensembl provider this resolver annotates with.

        Shared rather than rebuilt so a caller that has to resolve a record's identity before
        it has a variant to ask about - which is what turns a VCF row into a prepared record -
        does it against the same release and the same cache the evidence comes from.
        """
        return self._ensembl_provider()

    @property
    def cache_clients(self):
        """The HTTP caches this resolver replays from, as {label: client}.

        Exposed so a caller can record what a run actually fetched - which cache entries, and
        whether the network was touched at all. An offline run that silently reached the
        network would otherwise look identical to one that did not.
        """
        return {"identity": self._client, "external": self._external}

    def resolve(self, identity: dict, variant: Variant) -> ResolvedEvidence:
        """Resolve variant identity plus the parser-owned CONDITION supplied by the caller."""
        resolved = ResolvedEvidence()
        annotation, predictions = self._annotate(identity, variant, resolved)
        if annotation is not None:
            resolved.records.append(annotation)
            resolved.records.extend(predictions)
        self._add_population(variant, resolved)
        self._add_lof_mechanism(annotation, variant, resolved)
        self._add_gene2phenotype_mechanism(annotation, variant, resolved)
        self._add_cspec_mechanism(annotation, variant, resolved)
        self._add_gene_disease_associations(annotation, variant, resolved)
        self._add_disease_matching(annotation, identity, resolved)
        self._add_gene_disease_draft(annotation, variant, identity, resolved)
        self._add_splice_default(annotation, variant, resolved)
        self._add_initiation_assessment(annotation, variant, identity, resolved)
        self._add_pvs1_transcript_gates(annotation, variant, resolved)

        try:
            suite = self._suite()
        except PROVIDER_ERRORS as exc:
            # Building the suite asks Ensembl which release is current, and that lookup can
            # fail on its own. Outside a guard it escaped resolve() entirely, taking the
            # population evidence already gathered with it and surfacing to the API as an
            # ExceptionGroup with the cause buried. It is a provider failure like any other.
            resolved.failures.append({"provider": "Ensembl", "error": str(exc)})
            return resolved
        transcript = annotation.get("transcript") if annotation else None
        dbnsfp = resolved.absorb(DbnsfpProvider.name, suite.predictions(variant, transcript))
        resolved.note(DbnsfpProvider.name, records=dbnsfp.records, error=dbnsfp.error,
                      provider_version=suite.dbnsfp_version(), queried_variants=1,
                      calibration_use="PP3_BP4_WITH_CONFIGURED_CALIBRATION",
                      use_restriction="PP3_BP4_WITH_CONFIGURED_CALIBRATION")
        accession = identity.get("CLNVARIATIONID")
        outcome, identity_evidence = suite.clinvar_record(accession, variant)
        resolved.absorb(ClinVarProvider.name, outcome)
        if identity_evidence:
            resolved.identity_evidence.append(identity_evidence)
        # Only a VCV-shaped accession was ever a query; anything else was skipped without
        # asking, and counting it as a match would report lookups that never happened.
        if accession and VCV.fullmatch(str(accession)):
            resolved.note(ClinVarProvider.name, records=outcome.records,
                          provider_version=self._clinvar_release,
                          queried_accessions=1,
                          matched_records=1 if outcome.error is None else 0,
                          error=outcome.error,
                          classification_use="NOT_PP5_BP6",
                          ps1_comparator_errors=[],
                          use_restriction="NOT_PP5_BP6")

        if annotation is None or "missense_variant" not in (annotation.get("consequences") or []):
            return resolved
        # PS1/PM5/PM1 all compare protein-level changes, so they are only
        # meaningful for a missense annotation.
        splice_score = splice_score_for(predictions, annotation.get("transcript"))
        counts = {"PS1": ("ps1_comparator_searches", "ps1_comparator_evidence"),
                  "PM5": ("pm5_residue_searches", "pm5_comparator_evidence")}
        for criterion in ("PS1", "PM5"):
            comparison = resolved.absorb(
                f"{ClinVarComparatorProvider.name}:{criterion}",
                suite.comparator(criterion, annotation, variant, splice_score),
            )
            searches, matches = counts[criterion]
            if comparison.error is not None:
                resolved.note(ClinVarProvider.name, use_restriction="NOT_PP5_BP6",
                              ps1_comparator_errors=[{"variant_key": variant.key,
                                                      "criterion": criterion,
                                                      "error": comparison.error}])
                continue
            resolved.note(ClinVarProvider.name, use_restriction="NOT_PP5_BP6",
                          **{searches: 1 if comparison.searched else 0,
                             matches: comparison.matches})
        if self._hotspot_policy and annotation.get("protein_start"):
            hotspot = resolved.absorb(
                ClinVarHotspotProvider.name, suite.hotspot(annotation, variant))
            policy = self._hotspot_policy or {}
            resolved.note(ClinVarHotspotProvider.name, records=hotspot.records,
                          provider_version=self._clinvar_release,
                          error=hotspot.error, region_evidence=hotspot.matches,
                          # The policy the density was counted under travels with the count:
                          # the same window and thresholds are what make it reproducible.
                          policy_version=policy.get("policy_version"),
                          policy_source=policy.get("policy_source"),
                          window_aa=policy.get("window_aa"),
                          min_pathogenic=policy.get("min_pathogenic"),
                          max_benign=policy.get("max_benign"),
                          use_restriction="PM1_HOTSPOT_ROUTE_ONLY")
        return resolved

    def _add_lof_mechanism(self, annotation, variant, resolved):
        """PVS1's G01 gate, from ClinGen dosage rather than a per-gene review.

        Off unless asked for: it supplies an automated stand-in for a curated judgment, so a
        run has to opt into it. A curated gene_disease record for the same gene still wins -
        _resolve_mechanism() prefers a condition-specific one, and these are gene-level.
        """
        if not self._with_clingen_dosage or annotation is None:
            return
        provider = ClinGenDosageProvider(self._external)
        try:
            records = provider.get_mechanism(variant, annotation.get("gene"))
        except PROVIDER_ERRORS as exc:
            resolved.failures.append({"provider": ClinGenDosageProvider.name, "error": str(exc)})
            resolved.note(ClinGenDosageProvider.name, error=str(exc), genes_queried=1,
                          use_restriction="PVS1_LOF_MECHANISM_GATE_ONLY")
            return
        resolved.records.extend(records)
        resolved.note(ClinGenDosageProvider.name, records=records, genes_queried=1,
                      method=DOSAGE_METHOD,
                      provider_version=records[0]["source_version"] if records else None,
                      use_restriction="PVS1_LOF_MECHANISM_GATE_ONLY")

    def _add_gene2phenotype_mechanism(self, annotation, variant, resolved):
        """PVS1's mechanism gate from G2P, which curates per gene-disease pair.

        Off unless asked for, like the dosage stand-in beside it. Where dosage scores one
        haploinsufficiency judgment per gene, G2P states the molecular mechanism for each
        disease and names its MONDO term, so it answers for a gene that loses function in one
        disease and gains it in another - the case a per-gene score cannot separate.
        """
        if not self._with_gene2phenotype or annotation is None:
            return
        gene = annotation.get("gene")
        if not gene:
            return
        try:
            records = Gene2PhenotypeProvider(self._external).get_mechanism(variant, gene)
        except PROVIDER_ERRORS as exc:
            resolved.failures.append(
                {"provider": Gene2PhenotypeProvider.name, "error": str(exc)})
            resolved.note(Gene2PhenotypeProvider.name, error=str(exc), genes_queried=1,
                          use_restriction="PVS1_LOF_MECHANISM_GATE_ONLY")
            return
        resolved.records.extend(records)
        resolved.note(Gene2PhenotypeProvider.name, records=records, genes_queried=1,
                      method=G2P_METHOD,
                      use_restriction="PVS1_LOF_MECHANISM_GATE_ONLY")

    def _add_cspec_mechanism(self, annotation, variant, resolved):
        """PVS1's mechanism gate from a VCEP's own criteria specification.

        Off unless a snapshot path is given, like the two producers beside it. Where dosage
        scores one gene and G2P curates gene-disease pairs it has reached, this answers for
        any gene a VCEP wrote a specification for - which is how a recessive, non-dosage-
        scored gene gets an answer at all. Reading PVS1 applicability as the mechanism is a
        project decision the provider records on every record; see its module docstring.
        """
        if not self._cspec_applicability_path or annotation is None:
            return
        gene = annotation.get("gene")
        if not gene:
            return
        if self._cspec_applicability is None:
            self._cspec_applicability = CSpecApplicabilityProvider(
                self._cspec_applicability_path)
        try:
            records = self._cspec_applicability.get_mechanism(variant, gene)
        except (OSError, ValueError) as exc:
            resolved.failures.append(
                {"provider": CSpecApplicabilityProvider.name, "error": str(exc)})
            resolved.note(CSpecApplicabilityProvider.name, error=str(exc), genes_queried=1,
                          use_restriction="PVS1_LOF_MECHANISM_GATE_ONLY")
            return
        resolved.records.extend(records)
        resolved.note(CSpecApplicabilityProvider.name, records=records, genes_queried=1,
                      method=CSPEC_METHOD,
                      provider_version=records[0]["source_version"] if records else None,
                      use_restriction="PVS1_LOF_MECHANISM_GATE_ONLY")

    def _add_gene_disease_associations(self, annotation, variant, resolved):
        """The diseases a gene is curated for, as candidates a curator can be offered.

        Off unless asked for. These are ClinGen gene-disease validity rows, and they carry
        their own `gene_disease_validity` category so a mechanism lookup cannot reach them: a
        validity classification says the gene and the disease are related and never that loss
        of function is why, which is the conflation gene_disease_draft.py exists to prevent.

        The rows are per gene, so an id of their own would collapse every variant in one gene
        into the first one the run saw - the same defect the dosage provider's ids were fixed
        for - and the variant key is appended to keep them apart.
        """
        if not self._with_gene_disease_associations or annotation is None:
            return
        gene = annotation.get("gene")
        if not gene:
            return
        try:
            rows = self._clingen_gene_validity_provider().get_validity(gene)
        except PROVIDER_ERRORS as exc:
            resolved.failures.append(
                {"provider": ClinGenGeneValidityProvider.name, "error": str(exc)})
            resolved.note(ClinGenGeneValidityProvider.name, error=str(exc), genes_queried=1,
                          use_restriction="CANDIDATE_CONDITION_SUGGESTION_ONLY")
            return
        records = [{**row, "variant_key": variant.key,
                    "evidence_id": f"{row['evidence_id']}:{variant.key}"}
                   for row in rows
                   if str(row.get("condition") or "").startswith("MONDO:")]
        resolved.records.extend(records)
        resolved.note(ClinGenGeneValidityProvider.name, records=records, genes_queried=1,
                      provider_version=records[0]["source_version"] if records else None,
                      use_restriction="CANDIDATE_CONDITION_SUGGESTION_ONLY")

    def _mondo_provider(self) -> MondoMappingProvider:
        if self._mondo is None:
            self._mondo = MondoMappingProvider(self._external)
        return self._mondo

    def normalize_condition(self, condition):
        """The case's condition resolved to MONDO, or None when it is already one or cannot be.

        Called by the caller that owns the case's context rather than returned as evidence:
        this says which disease the case is about, not something retrieved about the variant.
        """
        if not self._with_disease_matching or not condition:
            return None
        if str(condition).startswith("MONDO:"):
            return None
        try:
            return self._mondo_provider().normalize(condition)
        except PROVIDER_ERRORS:
            return None

    def _add_disease_matching(self, annotation, identity, resolved):
        """What PVS1's disease gate needs to compare two diseases that are written differently.

        Two things here, and the third (the case's own condition, resolved to MONDO) through
        normalize_condition() because it belongs to the case rather than to the variant. All
        of them refinements of one comparison and none of them a mechanism: the
        case's condition resolved to MONDO when it was recorded in OMIM or Orphanet, the MONDO
        ancestry of every disease named on a mechanism record, and ClinGen's lumping decisions
        for each curated gene-disease pair. Without them the gate can only answer on identical
        identifiers, and withholds PVS1 for evidence that is on file under a neighbouring term.

        They travel beside the records rather than replacing anything: the original identifier
        is kept, and an exclusion or an ancestry relation is attached to the curation that
        made it.
        """
        if not self._with_disease_matching or annotation is None:
            return
        try:
            lumping = ClinGenLumpingProvider(self._external, self._mondo_provider())
            tree = MondoHierarchyProvider(self._external)
            asked, scoped, related, unresolved = 0, 0, 0, []
            for record in resolved.records:
                if record.get("category") != "gene_disease":
                    continue
                curated = record.get("condition")
                if not str(curated or "").startswith("MONDO:"):
                    continue
                asked += 1
                scope = lumping.get_scope(record.get("gene"), curated)
                if scope:
                    record["condition_scope"] = scope
                    scoped += 1
                ancestry = tree.ancestors(curated)
                if ancestry:
                    record["condition_ancestors"] = ancestry["ancestors"]
                    related += 1
                if not scope and not ancestry:
                    unresolved.append(curated)
        except PROVIDER_ERRORS as exc:
            resolved.failures.append(
                {"provider": "PVS1 disease matching", "error": str(exc)})
            resolved.note("PVS1 disease matching", error=str(exc),
                          use_restriction="DISEASE_MATCH_ONLY")
            return
        resolved.note("PVS1 disease matching", conditions_queried=asked,
                      lumping_scopes=scoped, ancestries=related,
                      unresolved=sorted(set(unresolved)),
                      use_restriction="DISEASE_MATCH_ONLY")

    def _clingen_gene_validity_provider(self) -> ClinGenGeneValidityProvider:
        if self._clingen_gene_validity is None:
            self._clingen_gene_validity = ClinGenGeneValidityProvider(self._external)
        return self._clingen_gene_validity

    def _gnomad_constraint_provider(self) -> GnomadConstraintProvider:
        if self._gnomad_constraint is None:
            self._gnomad_constraint = GnomadConstraintProvider(self._external)
        return self._gnomad_constraint

    def _add_gene_disease_draft(self, annotation, variant, identity, resolved):
        """PP2/BP1/PVS1's gene-disease mechanism, as a statistical suggestion.

        Off unless asked for, same convention as _add_lof_mechanism(): a
        ClinGen Gene-Disease Validity + gnomAD constraint suggestion is a
        lower bar than a curator's own reviewed gene_disease record, so a
        run has to opt in. Always produced under category
        "gene_disease_draft", never "gene_disease" -
        GeneDiseaseDraftProvider's own docstring is explicit that renaming
        it is a reviewer's decision, not this resolver's - a criterion that
        wants it (acmg_pipeline.criteria.mechanism's evaluate_mechanism())
        reads that category as an explicitly flagged fallback.

        Needs a disease context to pick the right one of a gene's several
        ClinGen-curated diseases (e.g. MYH7 has separate curations for
        cardiomyopathy and skeletal myopathy) - `identity["CONDITION"]`
        (a MONDO ID, from _IDENTITY_INFO_KEYS) supplies it. Without one,
        GeneDiseaseDraftProvider.build() still runs, but its own matching
        (no condition on the request side) cannot select a single disease
        out of several, so it correctly comes back as a low-confidence
        suggestion rather than guessing which disease was meant.
        """
        if not self._with_gene_disease_draft or annotation is None:
            return
        gene = annotation.get("gene")
        if not gene:
            return
        condition = identity.get("CONDITION") if isinstance(identity, dict) else None
        try:
            validity = self._clingen_gene_validity_provider().get_validity(gene)
        except PROVIDER_ERRORS as exc:
            resolved.failures.append(
                {"provider": ClinGenGeneValidityProvider.name, "error": str(exc)})
            return
        if not validity:
            return
        # A case's condition and ClinGen's own curated condition for the same gene are often
        # written at different MONDO depths (e.g. a case recorded under "familial adenomatous
        # polyposis 1", MONDO:0021056, against ClinGen's curated "classic or attenuated
        # familial adenomatous polyposis", MONDO:0021057, its direct parent) - an exact-string
        # match against `condition` below would then find nothing, the same gap
        # _add_disease_matching() exists to close for PVS1's own mechanism gate. Reusing
        # MondoHierarchyProvider.related() here (not duplicating the comparison) remaps
        # `condition` to whichever curated validity row it is the same-or-related to, so
        # GeneDiseaseDraftProvider's own exact-match selection still works unmodified.
        if condition and self._with_disease_matching:
            exact = any(row.get("condition") == condition for row in validity)
            if not exact:
                try:
                    tree = MondoHierarchyProvider(self._external)
                    matched = next(
                        (row["condition"] for row in validity
                         if row.get("condition")
                         and tree.related(condition, row["condition"])),
                        None,
                    )
                except PROVIDER_ERRORS as exc:
                    resolved.failures.append(
                        {"provider": MondoHierarchyProvider.name, "error": str(exc)})
                    matched = None
                if matched:
                    condition = matched
        try:
            constraint = self._gnomad_constraint_provider().get_constraint(gene)
        except PROVIDER_ERRORS as exc:
            resolved.failures.append(
                {"provider": GnomadConstraintProvider.name, "error": str(exc)})
            constraint = None
        clinvar_spectrum = None
        if self._with_clinvar_spectrum:
            try:
                clinvar_spectrum = ClinvarSpectrumProvider(
                    self._external, self._clinvar_release).get_spectrum(gene)
            except PROVIDER_ERRORS as exc:
                resolved.failures.append(
                    {"provider": ClinvarSpectrumProvider.name, "error": str(exc)})
        try:
            provider = GeneDiseaseDraftProvider(self._gene_disease_draft_policy or {})
            drafts = provider.build(
                variant_keys=[variant.key], gene=gene,
                transcript=annotation.get("transcript"),
                generated_at=datetime.now(timezone.utc).isoformat(),
                validity=validity, condition=condition, constraint=constraint,
                clinvar_spectrum=clinvar_spectrum,
            )
        except ValueError as exc:
            resolved.failures.append({"provider": GeneDiseaseDraftProvider.name, "error": str(exc)})
            return
        resolved.records.extend(drafts)

    # Same set acmg_pipeline.criteria.pvs1.CANONICAL_SPLICE names - duplicated rather than
    # imported to keep this module decoupled from criteria internals, the same way it never
    # imports from pvs1.py for anything else.
    _CANONICAL_SPLICE_CONSEQUENCES = frozenset({"splice_donor_variant", "splice_acceptor_variant"})

    def _add_splice_default(self, annotation, variant, resolved):
        """PVS1's SP01/SP02 default first-pass answer for canonical splice variants.

        Off unless asked for, same convention as the other optional steps above. Only
        for splice_donor_variant/splice_acceptor_variant - see
        acmg_pipeline.providers.splice_default's own docstring for why this is a
        literature base-rate default (checked against MYBPC3 c.2905+1G>A), not a
        per-variant sequence computation.
        """
        if not self._with_splice_default or annotation is None:
            return
        if not (set(annotation.get("consequences") or []) & self._CANONICAL_SPLICE_CONSEQUENCES):
            return
        try:
            provider = SpliceDefaultProvider(self._splice_default_policy_version)
            records = provider.get_splice_assessment(
                variant, annotation.get("transcript"),
                generated_at=datetime.now(timezone.utc).isoformat(),
            )
        except ValueError as exc:
            resolved.failures.append({"provider": SpliceDefaultProvider.name, "error": str(exc)})
            return
        resolved.records.extend(records)

    def _add_initiation_assessment(self, annotation, variant, identity, resolved):
        """PVS1's start-loss path (IC01/IC02/IC03) - wired for the first time here.

        InitiationProvider and UpstreamPathogenicProvider have existed since before this
        method, but only automated_cli.py's separate batch path ever called them - the
        integrated pipeline never did, so a live start_lost variant always came back
        UNKNOWN here regardless of what those two providers could answer. Same "off
        unless asked for" convention as the other optional steps above.
        """
        if not self._with_initiation_assessment or annotation is None:
            return
        if "start_lost" not in (annotation.get("consequences") or []):
            return
        gene, transcript = annotation.get("gene"), annotation.get("transcript")
        hgvsc = identity.get("HGVSC") if isinstance(identity, dict) else None
        if not (gene and transcript and hgvsc):
            return
        release = self._ensembl_release or self._ensembl_provider().release
        try:
            records = InitiationProvider(self._external, release).get_initiation_assessment(
                variant, gene, transcript, hgvsc)
        except PROVIDER_ERRORS as exc:
            resolved.failures.append({"provider": InitiationProvider.name, "error": str(exc)})
            return
        if not records:
            return
        record = records[0]
        downstream_start_codon = record.get("downstream_start_codon")
        try:
            upstream = UpstreamPathogenicProvider(self._external, self._clinvar_release) \
                .get_upstream_evidence(variant, gene, transcript, downstream_start_codon)
        except PROVIDER_ERRORS as exc:
            resolved.failures.append({"provider": UpstreamPathogenicProvider.name, "error": str(exc)})
            upstream = {}
        resolved.records.append({**record, **upstream})

    def _mane_provider(self) -> ManeTranscriptProvider:
        if self._mane is None:
            self._mane = ManeTranscriptProvider.from_directory(self._external)
        return self._mane

    def _add_pvs1_transcript_gates(self, annotation, variant, resolved):
        """PVS1's NF01/NF02/NF03/NF04/NF06/NF07 - MANE Select relevance, NMD
        prediction, and the protein-loss measurement (with its own UniProt-
        derived NF04/NF06 first pass - see providers/protein_region.py).

        Off unless asked for, same convention as the other optional steps
        above. Mirrors automated_cli.py's --with-mane-transcript/--with-nmd-
        prediction batch path (mane.get_transcript_assessment(),
        nmd.get_nmd_prediction(), region.get_protein_region()), which existed
        long before this method but was never called from here - a live
        PVS1 evaluation through the integrated pipeline always stopped at
        NF01 (transcript relevance unresolved) regardless of what those
        three already-implemented providers could answer.
        """
        if not self._with_pvs1_transcript_gates or annotation is None:
            return
        gene, transcript = annotation.get("gene"), annotation.get("transcript")
        hgvsc = annotation.get("hgvsc")
        if not (gene and transcript):
            return
        # NF02/NF04/NF06/NF07 also apply to a canonical splice variant that SP02 (see
        # _add_splice_default()) determines disrupts the reading frame - _truncating_path()
        # is reached from both V01=TRUNCATING and SP02=disrupted, and NF02 asks the same
        # NMD question either way (it is about position, not the original consequence type).
        consequences = set(annotation.get("consequences") or [])
        truncating = bool(consequences & (NMD_TRUNCATING | self._CANONICAL_SPLICE_CONSEQUENCES))
        release = self._ensembl_release or self._ensembl_provider().release

        exon = None
        if truncating and hgvsc:
            try:
                exon = NmdPredictionProvider(self._external, release).exon_on_transcript(
                    gene, transcript, hgvsc)
            except PROVIDER_ERRORS:
                exon = None
        try:
            resolved.records.extend(self._mane_provider().get_transcript_assessment(
                variant, gene, transcript, exon=exon))
        except PROVIDER_ERRORS as exc:
            resolved.failures.append({"provider": ManeTranscriptProvider.name, "error": str(exc)})

        if not (truncating and hgvsc):
            return
        try:
            nmd_provider = NmdPredictionProvider(self._external, release)
            resolved.records.extend(
                nmd_provider.get_nmd_prediction(variant, gene, transcript, hgvsc))
        except PROVIDER_ERRORS as exc:
            resolved.failures.append({"provider": NmdPredictionProvider.name, "error": str(exc)})
            return
        try:
            resolved.records.extend(ProteinRegionProvider(self._external, release)
                                    .get_protein_region(variant, gene, transcript, hgvsc,
                                                        annotation.get("protein_start")))
        except PROVIDER_ERRORS as exc:
            resolved.failures.append({"provider": ProteinRegionProvider.name, "error": str(exc)})

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
        for provider in self._population_providers():
            try:
                observations = provider.get_frequency(variant)
            except PROVIDER_ERRORS as exc:
                resolved.failures.append({"provider": provider.name, "error": str(exc)})
                continue
            if observations:  # Provider absence contributes nothing, never AF=0.
                resolved.records.extend(observations)
