"""Auditable internal outputs. VA-Spec export is a separate validation gate."""

import csv
import hashlib
import json
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

from acmg import __version__
from acmg.core.context import apply_context, context_summary, load_context
from acmg.core.input import sha256_file
from acmg.core.models import CRITERIA, Variant
from acmg.engine import evaluate_record, make_services


def load_object(path):
    def reject_constant(value):
        raise ValueError(f"Invalid JSON numeric constant: {value}")
    def reject_duplicates(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    value = json.loads(Path(path).read_text(encoding="utf-8-sig"), parse_constant=reject_constant,
                       object_pairs_hook=reject_duplicates)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def run_internal(input_path, evidence_path, config_path, output_dir, criteria=CRITERIA,
                 *, va_spec=False, context_path=None):
    """Evaluate prepared JSON. No source annotations are promoted into independent evidence."""
    document = load_object(input_path)
    evidence_document = load_object(evidence_path) if evidence_path else {"evidence": []}
    config = load_object(config_path) if config_path else {}
    context = load_context(load_object(context_path)) if context_path else None
    records = document.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("Prepared input must contain a nonempty records list")
    evidence = evidence_document.get("evidence")
    if not isinstance(evidence, list) or any(not isinstance(item, dict) for item in evidence):
        raise ValueError("Evidence must be a list of normalized evidence objects")
    services = make_services(evidence, config.get("population_providers"))
    outputs, errors, seen = [], [], set()
    for index, record in enumerate(records):
        record_id = record.get("record_id") if isinstance(record, dict) else None
        try:
            if not isinstance(record, dict) or not isinstance(record_id, str) or not record_id:
                raise ValueError("Record requires a string record_id")
            if record_id in seen:
                raise ValueError("Duplicate record_id")
            seen.add(record_id)
            Variant(**record["variant"])
            identity = record.get("identity_provenance")
            if not isinstance(identity, list) or not identity:
                raise ValueError("Prepared input requires identity_provenance from verified mapping")
            results = evaluate_record(apply_context(record, context), services, config, criteria)
            outputs.append({"record_id": record_id, "source": record.get("source", {}),
                            "variant": record["variant"], "identity_provenance": identity,
                            "results": [result.to_dict() for result in results]})
        except (ValueError, KeyError, TypeError) as exc:
            errors.append({"record_index": index, "record_id": record_id, "error": str(exc)})
    va_document = None
    if va_spec:
        from acmg.va_spec.mapper import export_document

        # Treat model validation as a gate: do not create a partially written run directory.
        va_document = export_document(outputs)
    output_dir = Path(output_dir)
    # Reuse is explicit at the caller: never silently overwrite evidence from an older run.
    output_dir.mkdir(parents=True, exist_ok=False)
    export_status = "VALIDATED" if va_spec else "NOT_PERFORMED"
    payload = {"schema_version": "1.1", "records": outputs, "input_errors": errors,
               "va_spec_export": export_status}
    result_path = output_dir / "results.json"
    result_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (output_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["record_id", "variant", "criterion", "status", "strength", "summary", "conflict_flags"])
        for record in outputs:
            for value in record["results"]:
                writer.writerow([record["record_id"], Variant(**record["variant"]).key, value["criterion"],
                                 value["status"], value["strength"] or "", value["summary"],
                                 ";".join(value["conflict_flags"])])
    va_path = None
    va_instances = []
    if va_spec:
        va_path = output_dir / "evidence-lines.json"
        va_path.write_text(json.dumps(va_document, ensure_ascii=False, indent=2,
                                      allow_nan=False) + "\n", encoding="utf-8")
        va_dir = output_dir / "va-spec-1.0.1"
        va_dir.mkdir()
        for record in va_document["records"]:
            safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", record["record_id"])
            for wrapped in record["evidence_lines"]:
                path = va_dir / f"{safe_id}--{wrapped['criterion']}.json"
                path.write_text(json.dumps(wrapped["evidence_line"], ensure_ascii=False,
                                           indent=2, allow_nan=False) + "\n", encoding="utf-8")
                va_instances.append({"path": str(path), "sha256": sha256_file(path)})
    paths = {"input": input_path, "evidence": evidence_path, "config": config_path,
             "context": context_path}
    manifest = {
        "schema_version": "1.0", "tool_version": __version__, "python_version": platform.python_version(),
        "evaluated_at": datetime.now(timezone.utc).isoformat(), "network_used": False,
        "inputs": {name: {"path": str(path), "sha256": sha256_file(path)} for name, path in paths.items() if path},
        "criteria": list(criteria), "evaluated_records": len(outputs), "input_errors": len(errors),
        "va_spec_export": export_status, "final_classification": "OUT_OF_SCOPE",
        "result_sha256": sha256_file(result_path), "curated_context": context_summary(context),
        "evidence_sha256": hashlib.sha256(json.dumps(evidence, sort_keys=True, allow_nan=False).encode()).hexdigest(),
    }
    if va_path:
        manifest["va_spec"] = {
            **va_document["validated_by"], "envelope_path": str(va_path),
            "envelope_sha256": sha256_file(va_path), "instance_count": len(va_instances),
            "criterion_assessment_count": sum(
                len(record["criterion_assessments"]) for record in va_document["records"]
            ),
            "referenced_evidence_count": sum(
                len(record["referenced_evidence"]) for record in va_document["records"]
            ),
            "instances": va_instances,
        }
    (output_dir / "run-manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return payload
