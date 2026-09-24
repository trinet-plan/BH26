"""CLI entry point. Audit is usable before external dependencies are provisioned."""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from acmg_pipeline.automated_core.input import audit_demo
from acmg_pipeline.automated_core.identity import evaluation_inputs, reconcile
from acmg_pipeline.automated_core.reference import FastaReference
from acmg_pipeline.automated_core.models import CRITERIA, Variant
from acmg_pipeline.gene_disease import build_assessment_document, build_draft_document
from acmg_pipeline.automated_output import run_internal
from acmg_pipeline.criteria.pvs1 import CANONICAL_SPLICE
from acmg_pipeline.providers.clingen_dosage import METHOD as DOSAGE_METHOD, ClinGenDosageProvider
from acmg_pipeline.providers.gene2phenotype import (
    METHOD as G2P_METHOD, Gene2PhenotypeProvider,
)
from acmg_pipeline.providers.clingen_gene_validity import ClinGenGeneValidityProvider
from acmg_pipeline.providers.clingen_lumping import (
    METHOD as LUMPING_METHOD, ClinGenLumpingProvider,
)
from acmg_pipeline.providers.mondo import METHOD as MONDO_METHOD, MondoMappingProvider
from acmg_pipeline.providers.splice_default import (
    METHOD as SPLICE_DEFAULT_METHOD, SpliceDefaultProvider,
)
from acmg_pipeline.providers.mondo_hierarchy import (
    METHOD as MONDO_TREE_METHOD, MondoHierarchyProvider,
)
from acmg_pipeline.providers.clinvar import (
    VCV, ClinVarComparatorProvider, ClinVarHotspotProvider, ClinVarProvider,
)
from acmg_pipeline.providers.ensembl import EnsemblIdentityProvider
from acmg_pipeline.providers.gnomad import GnomadProvider
from acmg_pipeline.providers.mane import METHOD as MANE_METHOD, ManeTranscriptProvider
from acmg_pipeline.providers.initiation import (
    METHOD as INITIATION_METHOD, START_LOST, InitiationProvider,
)
from acmg_pipeline.providers.upstream_pathogenic import (
    METHOD as UPSTREAM_METHOD, UpstreamPathogenicProvider,
)
from acmg_pipeline.providers.protein_region import (
    METHOD as REGION_METHOD, ProteinRegionProvider,
)
from acmg_pipeline.providers.population_registry import (
    DEFAULT_POPULATION_SOURCES,
)
from acmg_pipeline.providers.nmd import (
    METHOD as NMD_METHOD, RULE_SOURCE as NMD_RULE_SOURCE, TRUNCATING as NMD_TRUNCATING,
    NmdPredictionProvider,
)
from acmg_pipeline.providers.http import CachedHttpClient
from acmg_pipeline.providers.togovar import API_VERSION as TOGOVAR_API_VERSION, TogoVarProvider
from acmg_pipeline.services.resolve import (
    ProviderEvidenceResolver, VariantProviderSuite, merge_manifests, splice_score_for,
)


def _population_sources(args):
    """The population sources this run reads frequencies from, or [] for none.

    The rule set is where they belong: which frequency databases a run consults decides what
    PM2 and BS1 see, so it is a recorded policy rather than a command-line default. This
    command used to build a TogoVar provider directly and leave frequency_sources unset,
    which takes every group TogoVar offers - a wider set than config/demo-rules.json asks
    for, and a different one from what the same variant gets through the API.

    --togovar-api-version still overrides the version, because pinning the API a run replayed
    is the flag's job; which databases to consult is not.
    """
    if not args.with_togovar:
        return []
    configured = (json.loads(args.rules.read_text(encoding="utf-8")).get("population_sources")
                  if args.rules else None)
    if not configured:
        return DEFAULT_POPULATION_SOURCES
    return [{**source, "api_version": args.togovar_api_version}
            if source.get("provider") == "togovar" else source
            for source in configured]


def _hotspot_policy(args):
    if not args.rules:
        raise ValueError("--with-pm1-hotspot requires --rules with PM1.hotspot")
    return json.loads(args.rules.read_text(encoding="utf-8")).get("PM1", {}).get("hotspot", {})


def _splice_default_version(args):
    if not args.with_splice_default:
        return None
    if not args.rules:
        raise ValueError(
            "--with-splice-default requires --rules with PVS1.splice_default_policy_version")
    return (json.loads(args.rules.read_text(encoding="utf-8"))
            .get("PVS1", {}).get("splice_default_policy_version"))


def _gene_disease_draft_policy(args):
    if not args.with_gene_disease_draft:
        return None
    if not args.rules:
        raise ValueError("--with-gene-disease-draft requires --rules with a gene_disease_draft policy")
    policy = json.loads(args.rules.read_text(encoding="utf-8")).get("gene_disease_draft")
    if not policy:
        raise ValueError("--rules has no gene_disease_draft policy")
    return policy


def _apply_context_conditions(records, context_path):
    """Set identity['CONDITION'] on every record from --context's record_contexts,
    BEFORE evidence resolution - see --context's own help text for why this has
    to happen here rather than at evaluate time.

    Mutates each record's identity dict in place and returns records unchanged
    otherwise; a record_id --context has nothing for is left alone (no
    condition, same as not passing --context at all for that record).
    """
    if context_path is None:
        return records
    document = json.loads(context_path.read_text(encoding="utf-8"))
    record_contexts = document.get("record_contexts") or {}
    for record in records:
        entry = record_contexts.get(record["record_id"])
        condition = entry.get("condition") if entry else None
        if condition:
            record["identity"]["CONDITION"] = condition
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(prog="acmg")
    sub = parser.add_subparsers(dest="command", required=True)
    audit = sub.add_parser("audit-demo", help="Audit original VCFs without claiming verified identity")
    audit.add_argument("--input-dir", type=Path, default=Path("democase"))
    audit.add_argument("--output-dir", type=Path, required=True)
    prepare = sub.add_parser("prepare-demo", help="Resolve demo identities against indexed GRCh38 FASTA")
    prepare.add_argument("--input-dir", type=Path, default=Path("democase"))
    prepare.add_argument("--reference", type=Path, required=True)
    prepare.add_argument("--identity-evidence", type=Path, required=True,
                         help="JSON mapping record_id to independently retrieved identity candidates")
    prepare.add_argument("--output-dir", type=Path, required=True)
    online = sub.add_parser("prepare-demo-online",
                            help="Resolve demo HGVS identities with cached Ensembl GRCh38 data")
    online.add_argument("--input-dir", type=Path, default=Path("democase"))
    online.add_argument("--cache-dir", type=Path, required=True)
    online.add_argument("--output-dir", type=Path, required=True)
    online.add_argument("--ensembl-release")
    online.add_argument("--evidence-cache-dir", type=Path)
    population = online.add_mutually_exclusive_group()
    population.add_argument(
        "--with-togovar", action="store_true",
        help="Fetch population frequencies through the TogoVar GRCh38 API",
    )
    population.add_argument(
        "--with-gnomad", action="store_true",
        help="Legacy direct gnomAD fetch retained for replaying existing cached runs",
    )
    online.add_argument("--togovar-api-version", default=TOGOVAR_API_VERSION)
    online.add_argument("--gnomad-release", default="4.1.1")
    online.add_argument("--with-clinvar", action="store_true")
    online.add_argument("--with-dbnsfp", action="store_true",
                        help="Fetch dbNSFP meta-predictor scores pinned to their dbNSFP release")
    online.add_argument("--with-pm1-hotspot", action="store_true",
                        help="Count ClinVar missense density around each residue as PM1 hotspot proxy")
    online.add_argument("--with-clingen-dosage", action="store_true",
                        help="Derive PVS1's LoF-mechanism gate from ClinGen haploinsufficiency "
                             "scores (automated stand-in for a curated gene_disease record)")
    online.add_argument("--with-gene2phenotype", action="store_true",
                        help="Take PVS1's LoF-mechanism gate from G2P curation, which is "
                             "scoped to one gene-disease pair and names its MONDO disease")
    online.add_argument("--cspec-applicability", default=None, metavar="PATH",
                        help="Take PVS1's LoF-mechanism gate from a collected snapshot of "
                             "ClinGen VCEP criteria specifications "
                             "(config/cspec_applicability.json), which answers for any "
                             "gene a VCEP wrote a specification for - including recessive "
                             "genes the dosage score and G2P leave unanswered")
    online.add_argument("--with-mondo-mapping", action="store_true",
                        help="Resolve each record's OMIM/Orphanet condition to MONDO so "
                             "PVS1's disease gate can compare it with curated evidence")
    online.add_argument("--with-clingen-gene-validity", action="store_true",
                        help="Attach ClinGen's curated gene-disease associations, so a case "
                             "with no condition can be offered the diseases the gene is "
                             "curated for. Never a mechanism - see doc and PVS1's own gate")
    online.add_argument("--with-clingen-lumping", action="store_true",
                        help="Attach the phenotypes each curated disease lumps in or keeps "
                             "out, so an included case is matched and an excluded one is "
                             "answered instead of being sent to review")
    online.add_argument("--with-mondo-hierarchy", action="store_true",
                        help="Look up MONDO ancestry so a mechanism curated for a parent or "
                             "child disease reaches a curator instead of being reported as "
                             "no mechanism at all")
    online.add_argument("--with-mane-transcript", action="store_true",
                        help="Assert PVS1's transcript-relevance gate when the evaluated "
                             "transcript is the gene's MANE Select (automated stand-in)")
    online.add_argument("--with-splice-default", action="store_true",
                        help="Answer PVS1's SP01/SP02 for canonical splice donor/acceptor "
                             "variants from the configured default policy, flagged as a "
                             "prediction rather than a curator's review")
    online.add_argument("--with-nmd-prediction", action="store_true",
                        help="Predict NMD from VEP exon numbering for PVS1's NF02 gate "
                             "(no record for the last two exons, where the rule needs a distance)")
    online.add_argument("--with-region-assessment", action="store_true",
                        help="Answer PM4/BP3's repeat-region question from UniProt Repeat/"
                             "Domain feature overlap over the altered protein interval, "
                             "flagged as an automated first pass rather than a curator's "
                             "own review")
    online.add_argument("--with-bp7-splice-assessment", action="store_true",
                        help="Answer BP7 from VEP's splice_region_variant consequence plus "
                             "the SpliceAI score already fetched for every variant (ClinGen "
                             "SVI 2023, PMID:37352859) - no new HTTP request of its own")
    online.add_argument("--with-clinvar-spectrum", action="store_true",
                        help="Attach each gene's ClinVar missense/truncating P/LP counts, an "
                             "input to --with-gene-disease-draft's suggestion")
    online.add_argument("--with-gene-disease-draft", action="store_true",
                        help="Derive PP2/BP1/PVS1's gene-disease mechanism as a flagged "
                             "statistical suggestion (ClinGen Gene-Disease Validity + gnomAD "
                             "constraint, thresholds from --rules' gene_disease_draft policy) "
                             "when no curator-reviewed gene_disease record exists - see "
                             "acmg_pipeline.criteria.mechanism's evaluate_mechanism(). Needs "
                             "--context too: ClinGenGeneValidityProvider can only pick the "
                             "right one of a gene's several curated diseases with a condition "
                             "to match against.")
    online.add_argument("--rules", type=Path,
                        help="Rules JSON supplying PM1.hotspot / gene_disease_draft thresholds")
    online.add_argument("--context", type=Path,
                        help="Context JSON (same record_contexts shape as evaluate's --context) "
                             "supplying each record's disease condition BEFORE evidence "
                             "resolution runs. Without this, identity['CONDITION'] is never set "
                             "during prepare-demo-online (only evaluate's own --context sees a "
                             "condition, too late for --with-gene-disease-draft's generation "
                             "step) - see acmg_pipeline.services.resolve.ProviderEvidenceResolver."
                             "_add_gene_disease_draft()'s own docstring")
    online.add_argument("--clinvar-release", default=datetime.now(timezone.utc).date().isoformat())
    online.add_argument("--offline", action="store_true")
    evaluate = sub.add_parser("evaluate", help="Evaluate independently sourced evidence for prepared variants")
    evaluate.add_argument("--input", type=Path, required=True)
    evaluate.add_argument("--evidence", type=Path)
    evaluate.add_argument("--config", type=Path)
    evaluate.add_argument("--context", type=Path,
                          help="Curated clinical context: condition, disease thresholds, BA1 exceptions")
    evaluate.add_argument("--criteria", default="all")
    evaluate.add_argument("--output-dir", type=Path, required=True)
    evaluate.add_argument("--offline", action="store_true")
    evaluate.add_argument("--internal-only", action="store_true",
                          help="Explicitly omit VA-Spec export while its dependency gate is pending")
    drafts = sub.add_parser(
        "build-gene-disease-drafts",
        help="Build review-only PP2/BP1/PVS1 mechanism candidates from versioned sources",
    )
    drafts.add_argument("--input", type=Path, required=True)
    drafts.add_argument("--output", type=Path, required=True)
    assessments = sub.add_parser(
        "build-gene-disease-evidence",
        help="Expand source-anchored mechanism decisions and optionally merge base evidence",
    )
    assessments.add_argument("--input", type=Path, required=True)
    assessments.add_argument("--base-evidence", type=Path)
    assessments.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "build-gene-disease-evidence":
            if args.output.exists():
                raise ValueError(f"Output already exists: {args.output}")
            source = json.loads(args.input.read_text(encoding="utf-8-sig"))
            base = (json.loads(args.base_evidence.read_text(encoding="utf-8-sig"))
                    if args.base_evidence else None)
            payload = build_assessment_document(source, base)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            print(f"Generated combined evidence with {len(payload['evidence'])} records")
            return 0
        if args.command == "build-gene-disease-drafts":
            if args.output.exists():
                raise ValueError(f"Output already exists: {args.output}")
            payload = build_draft_document(json.loads(args.input.read_text(encoding="utf-8-sig")))
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            print(f"Generated {len(payload['evidence'])} review-only gene-disease drafts")
            return 0
        if args.command == "evaluate":
            criteria = CRITERIA if args.criteria == "all" else tuple(args.criteria.split(","))
            if not criteria or any(code not in CRITERIA for code in criteria) or len(criteria) != len(set(criteria)):
                raise ValueError("Unknown or duplicate criterion")
            payload = run_internal(args.input, args.evidence, args.config, args.output_dir, criteria,
                                   va_spec=not args.internal_only, context_path=args.context)
            print(f"Evaluated {len(payload['records'])} records; {len(payload['input_errors'])} errors; VA-Spec {payload['va_spec_export']}")
            return 2 if payload["input_errors"] else 0
        if args.command in {"audit-demo", "prepare-demo", "prepare-demo-online"}:
            records = audit_demo(args.input_dir)
            # Keyed by the identifier as written, because that is what a record carries and
            # what a curator will look for in the manifest. Only the online command fills it.
            condition_mappings = {}
            condition_ancestry = {}
            if args.command == "prepare-demo":
                reference = FastaReference(args.reference)
                candidates = json.loads(args.identity_evidence.read_text(encoding="utf-8-sig"))
                if not isinstance(candidates, dict):
                    raise ValueError("Identity evidence must be an object keyed by record_id")
                unknown = set(candidates) - {r["record_id"] for r in records}
                if unknown:
                    raise ValueError(f"Unknown identity record IDs: {sorted(unknown)}")
                records = [reconcile(r, candidates.get(r["record_id"], []), reference) for r in records]
            if args.command == "prepare-demo-online":
                # One evidence implementation, shared with the API path. The CLI keeps only
                # what the resolver has no business in: resolving each record's identity, and
                # the gnomAD batch, whose single multi-variant request is the shape the
                # committed offline cache holds - a per-variant lookup would replay none of it.
                resolver = ProviderEvidenceResolver(
                    args.cache_dir, offline=args.offline,
                    evidence_cache_dir=args.evidence_cache_dir or args.cache_dir,
                    ensembl_release=args.ensembl_release,
                    clinvar_release=args.clinvar_release,
                    population_sources=_population_sources(args),
                    hotspot_policy=_hotspot_policy(args) if args.with_pm1_hotspot else None,
                    with_clingen_dosage=args.with_clingen_dosage,
                    with_gene2phenotype=args.with_gene2phenotype,
                    cspec_applicability_path=args.cspec_applicability,
                    with_disease_matching=(args.with_mondo_mapping or args.with_clingen_lumping
                                           or args.with_mondo_hierarchy),
                    with_gene_disease_associations=args.with_clingen_gene_validity,
                    with_pvs1_transcript_gates=(args.with_mane_transcript
                                                or args.with_nmd_prediction),
                    with_initiation_assessment=args.with_nmd_prediction,
                    with_splice_default=args.with_splice_default,
                    splice_default_policy_version=_splice_default_version(args),
                    with_clinvar_spectrum=args.with_clinvar_spectrum,
                    with_gene_disease_draft=args.with_gene_disease_draft,
                    gene_disease_draft_policy=_gene_disease_draft_policy(args),
                    with_region_assessment=args.with_region_assessment,
                    with_bp7_splice_assessment=args.with_bp7_splice_assessment,
                )
                provider = resolver.identity_provider
                mapped = {}
                for record in records:
                    try:
                        candidate, _annotation, _predictions = (
                            provider.map_record_with_evidence(record))
                        mapped[record["record_id"]] = [candidate]
                    except ValueError as exc:
                        record["issues"].append(f"IDENTITY_PROVIDER_ERROR: {exc}")
                        mapped[record["record_id"]] = []
                records = [reconcile(r, mapped[r["record_id"]], provider.reference)
                           for r in records]
                records = _apply_context_conditions(records, args.context)

                evidence = []
                external_manifest = []
                variants = {}
                by_variant_key = {}
                for row in records:
                    if not row["resolution"]:
                        continue
                    variant = Variant(**row["resolution"]["variant"])
                    variants[variant.key] = variant
                    by_variant_key.setdefault(variant.key, []).append(row)
                # Once per variant, not once per row. Several ALT records can resolve to the
                # same variant, and asking for the same evidence again would replay from
                # cache but count twice in the manifest.
                for key, rows in by_variant_key.items():
                    resolved_evidence = resolver.resolve(rows[0]["identity"], variants[key])
                    evidence.extend(resolved_evidence.records)
                    external_manifest = merge_manifests(
                        [external_manifest, resolved_evidence.manifest])
                    for row in rows:
                        # Which ClinVar record this variant was matched to is provenance for
                        # the identity, not evidence about the variant, so it goes on the row.
                        row["resolution"]["evidence"].extend(
                            resolved_evidence.identity_evidence)
                        for failure in resolved_evidence.failures:
                            row["issues"].append(
                                f"{failure['provider'].upper()}_PROVIDER_ERROR: "
                                f"{failure['error']}")
                if args.with_gnomad:
                    gnomad = GnomadProvider(
                        CachedHttpClient(args.evidence_cache_dir or args.cache_dir,
                                         offline=args.offline),
                        release=args.gnomad_release)
                    try:
                        batches = gnomad.get_frequencies(variants.values())
                    except ValueError as exc:
                        for row in records:
                            if row["resolution"]:
                                row["issues"].append(f"GNOMAD_PROVIDER_ERROR: {exc}")
                        batches = {}
                    observations = [item for batch in batches.values() if batch
                                    for item in batch]
                    evidence.extend(observations)
                    external_manifest.append({
                        "provider": gnomad.name, "provider_version": gnomad.release,
                        "queried_variants": len(variants), "observed_variants": sum(
                            batch is not None for batch in batches.values()),
                        "evidence": len(observations),
                    })
            args.output_dir.mkdir(parents=True, exist_ok=False)
            output = args.output_dir / "audit.json"
            output.write_text(json.dumps({"schema_version": "1.0", "records": records},
                                         ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            resolved = evaluation_inputs(records)
            for record in resolved:
                mapping = condition_mappings.get(record.get("condition"))
                if mapping:
                    record["condition_mapping"] = mapping
                # Ancestry is keyed by the MONDO term, which is what the mapping resolved to.
                term = (mapping or {}).get("normalized_condition") or record.get("condition")
                ancestors = condition_ancestry.get(term)
                if ancestors:
                    record["condition_ancestors"] = ancestors
            if args.command in {"prepare-demo", "prepare-demo-online"}:
                (args.output_dir / "variants.json").write_text(
                    json.dumps({"schema_version": "1.0", "records": resolved},
                               ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            if args.command == "prepare-demo-online":
                resolved_keys = {record["resolution"]["variant"]["assembly"] + ":" +
                                 record["resolution"]["variant"]["chrom"] + ":" +
                                 str(record["resolution"]["variant"]["pos"]) + ":" +
                                 record["resolution"]["variant"]["ref"] + ":" +
                                 record["resolution"]["variant"]["alt"]
                                 for record in records if record["resolution"]}
                evidence = [item for item in evidence if item["variant_key"] in resolved_keys]
                evidence = list({item["evidence_id"]: item for item in evidence}.values())
                annotation_count = sum(item.get("category") == "annotation" for item in evidence)
                prediction_count = sum(item.get("category") == "computational" for item in evidence)
                (args.output_dir / "evidence.json").write_text(
                    json.dumps({"schema_version": "1.0", "evidence": evidence},
                               ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
                # The resolver owns the caches now, so the reproducibility record is read
                # off it: which entries a run replayed, and whether it touched the network
                # at all. An offline run that silently did is the thing this catches.
                caches = resolver.cache_clients
                client, external_client = caches["identity"], caches["external"]
                manifest = {
                    "schema_version": "1.0", "provider": provider.name,
                    "provider_version": provider.release,
                    "retrieved_at": datetime.now(timezone.utc).isoformat(),
                    "offline": args.offline,
                    "network_used": client.network_used or external_client.network_used,
                    "cache_entries": client.used,
                    "cache_set_sha256": hashlib.sha256(
                        json.dumps(client.used, sort_keys=True).encode()).hexdigest(),
                    "records": len(records), "resolved": len(resolved),
                    "annotation_evidence": annotation_count,
                    "computational_evidence": prediction_count,
                    "external_providers": external_manifest,
                    "external_network_used": external_client.network_used,
                    "external_cache_entries": external_client.used,
                    "external_cache_set_sha256": hashlib.sha256(
                        json.dumps(external_client.used, sort_keys=True).encode()).hexdigest(),
                    "total_evidence": len(evidence),
                }
                (args.output_dir / "identity-manifest.json").write_text(
                    json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            print(f"Audited {len(records)} ALT records; {len(resolved)} identities resolved: {output}")
        return 0
    except (ValueError, OSError) as exc:
        parser.exit(2, f"acmg: {exc}\n")


# Without this, `python -m acmg_pipeline.automated_cli ...` imports the module, defines
# main(), runs nothing and exits 0 - a silent success that looks like a completed run and
# leaves no output directory. Every documented invocation goes through -m, so the entry
# point has to exist here.
if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
