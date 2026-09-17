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
from acmg_pipeline.providers.clingen_dosage import METHOD as DOSAGE_METHOD, ClinGenDosageProvider
from acmg_pipeline.providers.gene2phenotype import (
    METHOD as G2P_METHOD, Gene2PhenotypeProvider,
)
from acmg_pipeline.providers.mondo import METHOD as MONDO_METHOD, MondoMappingProvider
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
from acmg_pipeline.providers.nmd import (
    METHOD as NMD_METHOD, RULE_SOURCE as NMD_RULE_SOURCE, TRUNCATING as NMD_TRUNCATING,
    NmdPredictionProvider,
)
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
    online.add_argument("--with-clingen-dosage", action="store_true",
                        help="Derive PVS1's LoF-mechanism gate from ClinGen haploinsufficiency "
                             "scores (automated stand-in for a curated gene_disease record)")
    online.add_argument("--with-gene2phenotype", action="store_true",
                        help="Take PVS1's LoF-mechanism gate from G2P curation, which is "
                             "scoped to one gene-disease pair and names its MONDO disease")
    online.add_argument("--with-mondo-mapping", action="store_true",
                        help="Resolve each record's OMIM/Orphanet condition to MONDO so "
                             "PVS1's disease gate can compare it with curated evidence")
    online.add_argument("--with-mane-transcript", action="store_true",
                        help="Assert PVS1's transcript-relevance gate when the evaluated "
                             "transcript is the gene's MANE Select (automated stand-in)")
    online.add_argument("--with-nmd-prediction", action="store_true",
                        help="Predict NMD from VEP exon numbering for PVS1's NF02 gate "
                             "(no record for the last two exons, where the rule needs a distance)")
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
            # Keyed by the identifier as written, because that is what a record carries and
            # what a curator will look for in the manifest. Only the online command fills it.
            condition_mappings = {}
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
                if args.with_clingen_dosage:
                    # PVS1 stops at G01 without a gene_disease record. This supplies one
                    # from published dosage curation, marked automated so it can never be
                    # mistaken for the per-gene review it stands in for.
                    dosage = ClinGenDosageProvider(external_client)
                    dosage_records, dosage_errors, seen_genes = [], [], set()
                    for annotation in annotations:
                        gene = annotation.get("gene")
                        key = (annotation["variant_key"], gene)
                        if not gene or key in seen_genes:
                            continue
                        seen_genes.add(key)
                        try:
                            dosage_records.extend(
                                dosage.get_mechanism(variants[annotation["variant_key"]], gene))
                        except (FetchError, ValueError) as exc:
                            dosage_errors.append(f"{gene}: {exc}")
                    evidence.extend(dosage_records)
                    external_manifest.append({
                        "provider": dosage.name,
                        "provider_version": dosage_records[0]["source_version"] if dosage_records else None,
                        "method": DOSAGE_METHOD,
                        "genes_queried": len(seen_genes), "evidence": len(dosage_records),
                        "errors": dosage_errors,
                        "use_restriction": "PVS1_LOF_MECHANISM_GATE_ONLY",
                    })
                if args.with_gene2phenotype:
                    # Disease-scoped where ClinGen dosage is gene-scoped, so it separates a
                    # gene that loses function in one disease from one that gains it in
                    # another - the case a single per-gene score cannot.
                    g2p = Gene2PhenotypeProvider(external_client)
                    g2p_records, g2p_errors, g2p_seen = [], [], set()
                    for annotation in annotations:
                        gene = annotation.get("gene")
                        key = (annotation["variant_key"], gene)
                        if not gene or key in g2p_seen:
                            continue
                        g2p_seen.add(key)
                        try:
                            g2p_records.extend(
                                g2p.get_mechanism(variants[annotation["variant_key"]], gene))
                        except (FetchError, ValueError, KeyError) as exc:
                            g2p_errors.append(f"{gene}: {exc}")
                    evidence.extend(g2p_records)
                    external_manifest.append({
                        "provider": g2p.name,
                        "provider_version": sorted(
                            {item["source_version"] for item in g2p_records}) or None,
                        "method": G2P_METHOD,
                        "genes_queried": len(g2p_seen), "evidence": len(g2p_records),
                        "errors": g2p_errors,
                        "use_restriction": "PVS1_LOF_MECHANISM_GATE_ONLY",
                    })
                if args.with_mondo_mapping:
                    # PVS1's disease gate compares identifiers, so a case in OMIM and a
                    # curation in MONDO have to be resolved to one vocabulary first. The
                    # original identifier is kept; the mapping travels beside it.
                    mondo = MondoMappingProvider(external_client)
                    mondo_errors = []
                    for record in records:
                        condition = record.get("condition")
                        if not condition or condition in condition_mappings:
                            continue
                        try:
                            mapping = mondo.normalize(condition)
                        except (FetchError, ValueError) as exc:
                            mondo_errors.append(f"{condition}: {exc}")
                            continue
                        if mapping:
                            condition_mappings[condition] = mapping
                    resolvable = {record.get("condition") for record in records
                                  if record.get("condition")}
                    external_manifest.append({
                        "provider": mondo.name,
                        "provider_version": next(
                            (item["source_version"] for item in condition_mappings.values()
                             if item.get("source_version")), None),
                        "method": MONDO_METHOD,
                        "conditions_queried": len(resolvable),
                        "evidence": len(condition_mappings),
                        "unresolved": sorted(resolvable - set(condition_mappings)),
                        "errors": mondo_errors,
                        "use_restriction": "DISEASE_MATCH_EQUIVALENCE_ONLY",
                    })
                if args.with_mane_transcript:
                    # PVS1's NF01 gate. Only a MANE Select match produces a record; see
                    # providers/mane.py for why a non-match is not NOT_RELEVANT.
                    mane = ManeTranscriptProvider.from_directory(external_client)
                    mane_records, mane_errors, seen = [], [], set()
                    for annotation in annotations:
                        key = (annotation["variant_key"], annotation.get("transcript"))
                        if key in seen:
                            continue
                        seen.add(key)
                        try:
                            # The exon also answers NF03 - see providers/mane.py. It is
                            # looked up only when this run is deriving NMD anyway, so no
                            # extra VEP request is made for a run that is not.
                            exon = None
                            if args.with_nmd_prediction:
                                exon = NmdPredictionProvider(
                                    external_client,
                                    args.ensembl_release or provider.release,
                                ).exon_on_transcript(
                                    annotation.get("gene"), annotation.get("transcript"),
                                    annotation.get("hgvsc"))
                            mane_records.extend(mane.get_transcript_assessment(
                                variants[annotation["variant_key"]],
                                annotation.get("gene"), annotation.get("transcript"),
                                exon=exon))
                        except (FetchError, ValueError) as exc:
                            mane_errors.append(f"{annotation.get('gene')}: {exc}")
                    evidence.extend(mane_records)
                    external_manifest.append({
                        "provider": mane.name, "provider_version": f"MANE v{mane.release}",
                        "method": MANE_METHOD,
                        "transcripts_queried": len(seen), "evidence": len(mane_records),
                        "errors": mane_errors,
                        "use_restriction": "PVS1_TRANSCRIPT_RELEVANCE_GATE_ONLY",
                    })
                if args.with_nmd_prediction:
                    # PVS1's NF02 gate. Its own VEP request, so the committed offline
                    # annotation cache keeps replaying unchanged - see providers/nmd.py.
                    nmd = NmdPredictionProvider(external_client, args.ensembl_release
                                                or provider.release)
                    nmd_records, nmd_errors, seen = [], [], set()
                    for annotation in annotations:
                        key = (annotation["variant_key"], annotation.get("transcript"))
                        if key in seen or not set(annotation.get("consequences") or []) & NMD_TRUNCATING:
                            continue
                        seen.add(key)
                        try:
                            nmd_records.extend(nmd.get_nmd_prediction(
                                variants[annotation["variant_key"]], annotation.get("gene"),
                                annotation.get("transcript"), annotation.get("hgvsc")))
                        except (FetchError, ValueError) as exc:
                            nmd_errors.append(f"{annotation.get('hgvsc')}: {exc}")
                    evidence.extend(nmd_records)
                    external_manifest.append({
                        "provider": nmd.name, "provider_version": nmd.release,
                        "method": NMD_METHOD, "rule_source": NMD_RULE_SOURCE,
                        "variants_queried": len(seen), "evidence": len(nmd_records),
                        "errors": nmd_errors,
                        "use_restriction": "PVS1_NMD_GATE_ONLY",
                    })
                    # PVS1's NF07 measurement, for the variants the NMD rule says escape
                    # decay. The region gates NF04/NF06 stay unanswered - see
                    # providers/protein_region.py.
                    region = ProteinRegionProvider(external_client, nmd.release)
                    region_records, region_errors = [], []
                    for annotation in annotations:
                        if not set(annotation.get("consequences") or []) & NMD_TRUNCATING:
                            continue
                        try:
                            region_records.extend(region.get_protein_region(
                                variants[annotation["variant_key"]], annotation.get("gene"),
                                annotation.get("transcript"), annotation.get("hgvsc"),
                                annotation.get("protein_start")))
                        except (FetchError, ValueError) as exc:
                            region_errors.append(f"{annotation.get('hgvsc')}: {exc}")
                    evidence.extend(region_records)
                    external_manifest.append({
                        "provider": region.name, "provider_version": region.release,
                        "method": REGION_METHOD,
                        "evidence": len(region_records), "errors": region_errors,
                        "use_restriction": "PVS1_PROTEIN_LOSS_MEASUREMENT_ONLY",
                    })
                    # PVS1's IC02 gate. IC01 and IC03 stay unanswered - see
                    # providers/initiation.py.
                    initiation = InitiationProvider(external_client, nmd.release)
                    init_records, init_errors = [], []
                    for annotation in annotations:
                        if START_LOST not in (annotation.get("consequences") or []):
                            continue
                        try:
                            init_records.extend(initiation.get_initiation_assessment(
                                variants[annotation["variant_key"]], annotation.get("gene"),
                                annotation.get("transcript"), annotation.get("hgvsc")))
                        except (FetchError, ValueError) as exc:
                            init_errors.append(f"{annotation.get('hgvsc')}: {exc}")
                    evidence.extend(init_records)
                    external_manifest.append({
                        "provider": initiation.name, "provider_version": initiation.release,
                        "method": INITIATION_METHOD,
                        "evidence": len(init_records), "errors": init_errors,
                        "use_restriction": "PVS1_DOWNSTREAM_START_ONLY",
                    })
                    # PVS1's IC03 gate, bounded by the codon IC02 found.
                    upstream = UpstreamPathogenicProvider(external_client, args.clinvar_release)
                    # Merged into the record IC02 built: PVS1 reads two
                    # initiation_assessment records as a conflict.
                    answered, up_errors = 0, []
                    for record in init_records:
                        codon = record.get("downstream_start_codon")
                        if not codon:
                            continue
                        try:
                            fields = upstream.get_upstream_evidence(
                                variants[record["variant_key"]], record["gene"],
                                record["transcript"], codon)
                        except (FetchError, ValueError) as exc:
                            up_errors.append(f"{record['gene']}: {exc}")
                            continue
                        if fields:
                            record.update(fields)
                            answered += 1
                    external_manifest.append({
                        "provider": upstream.name, "provider_version": args.clinvar_release,
                        "method": UPSTREAM_METHOD,
                        "evidence": answered, "errors": up_errors,
                        "use_restriction": "PVS1_UPSTREAM_PATHOGENIC_ONLY",
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
