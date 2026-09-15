"""CLI entry point. Audit is usable before external dependencies are provisioned."""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from acmg.core.input import audit_demo
from acmg.core.identity import evaluation_inputs, reconcile
from acmg.core.reference import FastaReference
from acmg.core.models import CRITERIA, Variant
from acmg.output import run_internal
from acmg.providers.clinvar import VCV, ClinVarProvider
from acmg.providers.ensembl import EnsemblIdentityProvider
from acmg.providers.gnomad import GnomadProvider
from acmg.providers.http import CachedHttpClient


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
    online.add_argument("--with-gnomad", action="store_true")
    online.add_argument("--gnomad-release", default="4.1.1")
    online.add_argument("--with-clinvar", action="store_true")
    online.add_argument("--clinvar-release", default=datetime.now(timezone.utc).date().isoformat())
    online.add_argument("--offline", action="store_true")
    evaluate = sub.add_parser("evaluate", help="Evaluate independently sourced evidence for prepared variants")
    evaluate.add_argument("--input", type=Path, required=True)
    evaluate.add_argument("--evidence", type=Path)
    evaluate.add_argument("--config", type=Path)
    evaluate.add_argument("--criteria", default="all")
    evaluate.add_argument("--output-dir", type=Path, required=True)
    evaluate.add_argument("--offline", action="store_true")
    evaluate.add_argument("--internal-only", action="store_true",
                          help="Explicitly omit VA-Spec export while its dependency gate is pending")
    args = parser.parse_args(argv)
    try:
        if args.command == "evaluate":
            criteria = CRITERIA if args.criteria == "all" else tuple(args.criteria.split(","))
            if not criteria or any(code not in CRITERIA for code in criteria) or len(criteria) != len(set(criteria)):
                raise ValueError("Unknown or duplicate criterion")
            payload = run_internal(args.input, args.evidence, args.config, args.output_dir, criteria,
                                   va_spec=not args.internal_only)
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
                for record in records:
                    try:
                        candidate, annotation = provider.map_record_with_annotation(record)
                        mapped[record["record_id"]] = [candidate]
                        annotations.append(annotation)
                    except ValueError as exc:
                        record["issues"].append(f"IDENTITY_PROVIDER_ERROR: {exc}")
                        mapped[record["record_id"]] = []
                records = [reconcile(r, mapped[r["record_id"]], provider.reference) for r in records]
                evidence = list(annotations)
                external_client = CachedHttpClient(
                    args.evidence_cache_dir or args.cache_dir, offline=args.offline
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
                            try:
                                clinvar_results[lookup] = clinvar.get_record(accession, variant)
                                clinvar_count += 1
                            except ValueError as exc:
                                clinvar_results[lookup] = exc
                        result = clinvar_results[lookup]
                        if isinstance(result, tuple):
                            item, identity = result
                            evidence.append(item)
                            row["resolution"]["evidence"].append(identity)
                        else:
                            row["issues"].append(f"CLINVAR_PROVIDER_ERROR: {result}")
                    external_manifest.append({
                        "provider": clinvar.name, "provider_version": clinvar.release,
                        "queried_accessions": len(seen), "matched_records": clinvar_count,
                        "evidence": clinvar_count,
                        "classification_use": "NOT_PP5_BP6",
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
