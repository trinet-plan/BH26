"""Curated clinical context attached to prepared records at evaluation time.

Disease context, disease-specific thresholds and the BA1 exception list are human decisions
rather than retrieved data, so they arrive as a separate, versioned document instead of being
mixed into the evidence a provider produced. A list that has not been fully transcribed says
so, and an incomplete list resolves nothing rather than asserting absence.
"""

from acmg.core.models import Variant


CONTEXT_FIELDS = ("condition", "condition_label", "inheritance", "disease_frequency_threshold")


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def load_context(document):
    """Validate a curated-context document and return the parts records can be given."""
    _require(isinstance(document, dict), "Curated context must be a JSON object")
    version = document.get("context_version")
    _require(isinstance(version, str) and version, "Curated context requires a context_version")
    records = document.get("records", {})
    _require(isinstance(records, dict), "Curated context records must be an object")
    parsed = {}
    for key, value in records.items():
        _require(isinstance(value, dict), f"Curated context for {key} must be an object")
        _require(len(key.split(":")) == 5, f"Curated context key must be a variant key: {key}")
        unknown = set(value) - set(CONTEXT_FIELDS)
        _require(not unknown, f"Unsupported curated context fields for {key}: {sorted(unknown)}")
        parsed[key] = dict(value)
    record_contexts = document.get("record_contexts", {})
    _require(isinstance(record_contexts, dict),
             "Curated context record_contexts must be an object")
    parsed_record_contexts = {}
    for record_id, value in record_contexts.items():
        _require(isinstance(record_id, str) and record_id,
                 "Curated record context requires a nonempty record_id")
        _require(isinstance(value, dict),
                 f"Curated record context for {record_id} must be an object")
        unknown = set(value) - set(CONTEXT_FIELDS)
        _require(not unknown,
                 f"Unsupported curated record context fields for {record_id}: {sorted(unknown)}")
        parsed_record_contexts[record_id] = dict(value)
    exceptions = document.get("ba1_exceptions")
    if exceptions is not None:
        _require(isinstance(exceptions, dict), "ba1_exceptions must be an object")
        _require(all(exceptions.get(field) for field in ("source", "source_version", "reviewed_at")),
                 "ba1_exceptions requires source, source_version and reviewed_at")
        # Where the list came from matters as much as its content: a hand-transcribed table
        # carries a different risk than a retrieved file, so every run records which it was.
        _require(exceptions.get("entry_method") in {"manual_transcription", "retrieved"},
                 "ba1_exceptions must state entry_method (manual_transcription or retrieved)")
        _require(isinstance(exceptions.get("complete"), bool),
                 "ba1_exceptions must state whether the list is complete")
        listed = exceptions.get("variants", [])
        _require(isinstance(listed, list), "ba1_exceptions variants must be a list")
        for item in listed:
            _require(isinstance(item, dict), "Each BA1 exception must be an object")
            key = item.get("variant_key")
            _require(isinstance(key, str) and len(key.split(":")) == 5,
                     f"BA1 exception requires a variant_key: {item}")
            # The published list is a table of transcript HGVS, so record how each entry was
            # resolved to a GRCh38 allele; a wrong key would silently exempt the wrong variant.
            _require(all(item.get(field) for field in ("gene", "hgvs_c", "caid", "resolved_by")),
                     f"BA1 exception {key} requires gene, hgvs_c, caid and resolved_by")
    return {"context_version": version, "records": parsed,
            "record_contexts": parsed_record_contexts, "ba1_exceptions": exceptions,
            "source": document.get("source")}


def apply_context(record, context):
    """Attach curated context to one prepared record without touching retrieved evidence."""
    if not context:
        return record
    key = Variant(**record["variant"]).key
    updated = {
        **record,
        **context["records"].get(key, {}),
        **context.get("record_contexts", {}).get(record.get("record_id"), {}),
    }
    exceptions = context.get("ba1_exceptions")
    # An incomplete list cannot say a variant is absent from it, so it resolves nothing.
    if exceptions and exceptions["complete"]:
        updated["ba1_exception_assessment"] = {
            "source": exceptions["source"], "source_version": exceptions["source_version"],
            "reviewed_at": exceptions["reviewed_at"],
            "is_exception": key in {item["variant_key"] for item in exceptions.get("variants", [])},
            "list_size": len(exceptions.get("variants", [])),
        }
    return updated


def context_summary(context):
    if not context:
        return None
    exceptions = context.get("ba1_exceptions")
    summary = {"context_version": context["context_version"], "source": context.get("source"),
               "records": len(context["records"]),
               "record_contexts": len(context.get("record_contexts", {}))}
    if exceptions:
        summary["ba1_exceptions"] = {
            "source": exceptions["source"], "source_version": exceptions["source_version"],
            "complete": exceptions["complete"], "variants": len(exceptions.get("variants", [])),
            "entry_method": exceptions["entry_method"],
            "transcription": exceptions.get("transcription"),
        }
    return summary
