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
from acmg_pipeline.providers.clinvar import (
    VCV, ClinVarComparatorProvider, ClinVarHotspotProvider, ClinVarProvider,
)
from acmg_pipeline.providers.ensembl import EnsemblIdentityProvider
from acmg_pipeline.providers.gnomad import GnomadProvider
from acmg_pipeline.providers.http import CachedHttpClient
from acmg_pipeline.providers.togovar import API_VERSION as TOGOVAR_API_VERSION, TogoVarProvider
from acmg_pipeline.services.resolve import VariantProviderSuite, splice_score_for


def main(argv=None):
    parser = argparse.ArgumentParser(prog="acmg")
    sub = parser.add_subparsers(dest="command", required=True)
    audit = sub.add_parser("audit-demo", help="Audit original VCFs without claiming verified identity")
    audit.add_argument("--input-dir", type=Path, default=Path("demo-data"))
    audit.add_argument("--output-dir", type=Path, required=True)
    prepare = sub.add_parser("prepare-demo", help="Resolve demo identities against indexed GRCh38 FASTA")
    prepare.add_argument("--input-dir", type=Path, default=Path("demo-data"))
    prepare.add_argument("--reference", type=Path, required=True)
    prepare.add_argument("--identity-evidence", type=Path, required=True,
                         help="JSON mapping record_id to independently retrieved identity candidates")
    prepare.add_argument("--output-dir", type=Path, required=True)
    online = sub.add_parser("prepare-demo-online",
                            help="Resolve demo HGVS identities with cached Ensembl GRCh38 data")
    online.add_argument("--input-dir", type=Path, default=Path("demo-data"))
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
    online.add_argument("--rules", type=Path,
                        help="Rules JSON supplying PM1.hotspot thresholds for --with-pm1-hotspot")
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
                client = CachedHttpClient(args.cache_dir, offline=args.offline)
                release = args.ensembl_release or EnsemblIdentityProvider.current_release(client)
                provider = EnsemblIdentityProvider(client, release)
                mapped = {}
                annotations = []
                predictions = []
                for record in records:
                    try:
                        candidate, annotation, record_predictions = provider.map_record_with_evidence(record)
                        mapped[record["record_id"]] = [candidate]
                        annotations.append(annotation)
                        predictions.extend(record_predictions)
                    except ValueError as exc:
                        record["issues"].append(f"IDENTITY_PROVIDER_ERROR: {exc}")
                        mapped[record["record_id"]] = []
                records = [reconcile(r, mapped[r["record_id"]], provider.reference) for r in records]
                evidence = [*annotations, *predictions]
                external_client = CachedHttpClient(
                    args.evidence_cache_dir or args.cache_dir, offline=args.offline
                )
                # Everything decided one variant at a time (dbNSFP, the
                # ClinVar VCV record, the PS1/PM5 comparators) goes through
                # the same suite the integrated pipeline uses; only the
                # gnomAD batch below stays here, because its single
                # multi-variant request is what the committed offline cache
                # holds. See acmg/services/resolve.py.
                suite = VariantProviderSuite(
                    external_client, clinvar_release=args.clinvar_release,
                    ensembl_provider=provider,
                )
                external_manifest = []
                variants = {
                    row["resolution"]["variant"]["assembly"] + ":" +
                    row["resolution"]["variant"]["chrom"] + ":" +
                    str(row["resolution"]["variant"]["pos"]) + ":" +
                    row["resolution"]["variant"]["ref"] + ":" +
                    row["resolution"]["variant"]["alt"]: Variant(**row["resolution"]["variant"])
                    for row in records if row["resolution"]
                }
                if args.with_togovar:
                    togovar = TogoVarProvider(
                        external_client, api_version=args.togovar_api_version
                    )
                    observations = []
                    observed_variants = 0
                    errors = []
                    for key, variant in sorted(variants.items()):
                        try:
                            batch = togovar.get_frequency(variant)
                        except ValueError as exc:
                            errors.append({"variant_key": key, "error": str(exc)})
                            continue
                        if batch:
                            observed_variants += 1
                            observations.extend(batch)
                    evidence.extend(observations)
                    external_manifest.append({
                        "provider": togovar.name,
                        "provider_version": f"API {togovar.api_version}",
                        "queried_variants": len(variants),
                        "observed_variants": observed_variants,
                        "evidence": len(observations),
                        "errors": errors,
                    })
                if args.with_gnomad:
                    gnomad = GnomadProvider(external_client, release=args.gnomad_release)
                    try:
                        batches = gnomad.get_frequencies(variants.values())
                    except ValueError as exc:
                        for row in records:
                            if row["resolution"]:
                                row["issues"].append(f"GNOMAD_PROVIDER_ERROR: {exc}")
                        batches = {}
                    observations = [item for batch in batches.values() if batch for item in batch]
                    evidence.extend(observations)
                    external_manifest.append({
                        "provider": gnomad.name, "provider_version": gnomad.release,
                        "queried_variants": len(variants), "observed_variants": sum(
                            batch is not None for batch in batches.values()),
                        "evidence": len(observations),
                    })
                if args.with_dbnsfp:
                    dbnsfp = suite.dbnsfp
                    transcripts = {item["variant_key"]: item["transcript"] for item in annotations}
                    scores = 0
                    dbnsfp_errors = []
                    for key, variant in sorted(variants.items()):
                        outcome = suite.predictions(variant, transcripts.get(key))
                        if outcome.error is not None:
                            dbnsfp_errors.append({"variant_key": key, "error": outcome.error})
                            continue
                        evidence.extend(outcome.records)
                        scores += len(outcome.records)
                    external_manifest.append({
                        "provider": dbnsfp.name,
                        "provider_version": suite.dbnsfp_version(),
                        "queried_variants": len(variants), "evidence": scores,
                        "errors": dbnsfp_errors,
                        "calibration_use": "PP3_BP4_WITH_CONFIGURED_CALIBRATION",
                    })
                if args.with_clinvar:
                    clinvar = ClinVarProvider(external_client, args.clinvar_release)
                    seen = set()
                    clinvar_results = {}
                    clinvar_count = 0
                    for row in records:
                        accession = row["identity"].get("CLNVARIATIONID")
                        if not accession or not VCV.fullmatch(accession) or not row["resolution"]:
                            continue
                        variant = Variant(**row["resolution"]["variant"])
                        lookup = (accession, variant.key)
                        if lookup not in seen:
                            seen.add(lookup)
                            outcome, identity = suite.clinvar_record(accession, variant)
                            clinvar_results[lookup] = (outcome, identity)
                            if outcome.error is None:
                                clinvar_count += 1
                        outcome, identity = clinvar_results[lookup]
                        if outcome.error is None:
                            evidence.extend(outcome.records)
                            row["resolution"]["evidence"].append(identity)
                        else:
                            row["issues"].append(f"CLINVAR_PROVIDER_ERROR: {outcome.error}")
                    clinvar_manifest = {
                        "provider": clinvar.name, "provider_version": clinvar.release,
                        "queried_accessions": len(seen), "matched_records": clinvar_count,
                        "evidence": clinvar_count,
                        "classification_use": "NOT_PP5_BP6",
                    }
                    comparator = ClinVarComparatorProvider(
                        external_client, args.clinvar_release, provider
                    )
                    unique_annotations = {
                        (item["variant_key"], item["transcript"]): item
                        for item in annotations if "missense_variant" in item["consequences"]
                    }
                    comparator_searches = 0
                    comparator_evidence = 0
                    comparator_errors = []
                    pm5_searches = 0
                    pm5_evidence = 0
                    tally = {"PS1": [0, 0], "PM5": [0, 0]}
                    for annotation in unique_annotations.values():
                        variant = variants[annotation["variant_key"]]
                        splice_score = splice_score_for(
                            [item for item in predictions
                             if item.get("variant_key") == variant.key],
                            annotation["transcript"],
                        )
                        for criterion in ("PS1", "PM5"):
                            outcome = suite.comparator(
                                criterion, annotation, variant, splice_score
                            )
                            if outcome.error is not None:
                                comparator_errors.append({
                                    "variant_key": variant.key, "error": outcome.error,
                                    "criterion": criterion,
                                })
                                continue
                            evidence.extend(outcome.records)
                            tally[criterion][0] += 1
                            tally[criterion][1] += outcome.matches
                    comparator_searches, comparator_evidence = tally["PS1"]
                    pm5_searches, pm5_evidence = tally["PM5"]
                    clinvar_manifest["pm5_residue_searches"] = pm5_searches
                    clinvar_manifest["pm5_comparator_evidence"] = pm5_evidence
                    clinvar_manifest["ps1_comparator_searches"] = comparator_searches
                    clinvar_manifest["ps1_comparator_evidence"] = comparator_evidence
                    clinvar_manifest["ps1_comparator_errors"] = comparator_errors
                    external_manifest.append(clinvar_manifest)
                if args.with_pm1_hotspot:
                    if not args.rules:
                        raise ValueError("--with-pm1-hotspot requires --rules with PM1.hotspot")
                    policy = json.loads(args.rules.read_text(encoding="utf-8"))                         .get("PM1", {}).get("hotspot", {})
                    hotspot = ClinVarHotspotProvider(external_client, args.clinvar_release, policy)
                    hotspot_suite = VariantProviderSuite(
                        external_client, clinvar_release=args.clinvar_release,
                        ensembl_provider=provider, hotspot_policy=policy,
                    )
                    hotspot_regions = 0
                    hotspot_errors = []
                    for annotation in {
                        (item["variant_key"], item["transcript"]): item
                        for item in annotations
                        if "missense_variant" in item["consequences"] and item.get("protein_start")
                    }.values():
                        outcome = hotspot_suite.hotspot(
                            annotation, variants[annotation["variant_key"]]
                        )
                        if outcome.error is not None:
                            hotspot_errors.append({"variant_key": annotation["variant_key"],
                                                   "error": outcome.error})
                            continue
                        evidence.extend(outcome.records)
                        hotspot_regions += outcome.matches
                    external_manifest.append({
                        "provider": hotspot.name, "provider_version": args.clinvar_release,
                        "policy_version": policy.get("policy_version"),
                        "policy_source": policy.get("policy_source"),
                        "window_aa": policy.get("window_aa"),
                        "min_pathogenic": policy.get("min_pathogenic"),
                        "max_benign": policy.get("max_benign"),
                        "region_evidence": hotspot_regions, "errors": hotspot_errors,
                        "use_restriction": "PM1_HOTSPOT_ROUTE_ONLY",
                    })
            args.output_dir.mkdir(parents=True, exist_ok=False)
            output = args.output_dir / "audit.json"
            output.write_text(json.dumps({"schema_version": "1.0", "records": records},
                                         ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            resolved = evaluation_inputs(records)
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
