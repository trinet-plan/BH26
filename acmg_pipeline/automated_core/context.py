"""Curated clinical context attached to prepared records at evaluation time.

Disease context, disease-specific thresholds and the BA1 exception list are human decisions
rather than retrieved data, so they arrive as a separate, versioned document instead of being
mixed into the evidence a provider produced. A list that has not been fully transcribed says
so, and an incomplete list resolves nothing rather than asserting absence.
"""

from acmg_pipeline.automated_core.models import Variant


CONTEXT_FIELDS = ("condition", "condition_label", "inheritance", "disease_frequency_threshold",
                  "ba1_threshold_override")
# The two fields a gene_frequency_thresholds entry may name - each mirrors one criterion's
# own threshold shape (BA1's is checked by criteria/ba1.py's threshold_policy(), BS1's IS
# disease_frequency_threshold, the same field a per-variant/per-record curated context
# already supplies - see apply_context()'s GENE_THRESHOLD_FIELDS mapping below).
_GENE_THRESHOLD_KEYS = ("ba1", "bs1")


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
    gene_thresholds = document.get("gene_frequency_thresholds")
    if gene_thresholds is not None:
        _require(isinstance(gene_thresholds, dict), "gene_frequency_thresholds must be an object")
        for gene, entry in gene_thresholds.items():
            _require(isinstance(gene, str) and gene, "gene_frequency_thresholds key must be a gene symbol")
            _require(isinstance(entry, dict), f"gene_frequency_thresholds[{gene!r}] must be an object")
            unknown = set(entry) - set(_GENE_THRESHOLD_KEYS)
            _require(not unknown,
                     f"Unsupported gene_frequency_thresholds fields for {gene!r}: {sorted(unknown)}")
            for criterion_key, threshold in entry.items():
                _require(isinstance(threshold, dict),
                         f"gene_frequency_thresholds[{gene!r}][{criterion_key!r}] must be an object")
                required = ("source", "source_version", "reviewed_at", "frequency_statistic")
                _require(all(threshold.get(field) for field in required),
                         f"gene_frequency_thresholds[{gene!r}][{criterion_key!r}] requires "
                         f"{', '.join(required)}")
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
            "gene_frequency_thresholds": gene_thresholds, "source": document.get("source")}


# Which context field each gene_frequency_thresholds sub-entry feeds - "bs1" feeds the same
# disease_frequency_threshold field a per-variant/per-record curated context already supplies
# (criteria/bs1.py's own disease_specific_threshold() reads it, condition-matching included -
# a gene_frequency_thresholds entry names its own real curated condition, same as a per-
# variant one would), "ba1" feeds a new field criteria/ba1.py's threshold_policy() checks
# before falling back to config["BA1"]'s global default.
GENE_THRESHOLD_TARGET_FIELD = {"bs1": "disease_frequency_threshold", "ba1": "ba1_threshold_override"}


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
    gene_entry = (context.get("gene_frequency_thresholds") or {}).get(record.get("gene"))
    if gene_entry:
        # Lowest priority: a per-variant or per-record entry above already named the same
        # field is a more specific curation and is never overwritten by a gene-wide default.
        for source_key, target_field in GENE_THRESHOLD_TARGET_FIELD.items():
            if target_field not in updated and gene_entry.get(source_key):
                updated[target_field] = gene_entry[source_key]
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
               "record_contexts": len(context.get("record_contexts", {})),
               "gene_frequency_thresholds": len(context.get("gene_frequency_thresholds") or {})}
    if exceptions:
        summary["ba1_exceptions"] = {
            "source": exceptions["source"], "source_version": exceptions["source_version"],
            "complete": exceptions["complete"], "variants": len(exceptions.get("variants", [])),
            "entry_method": exceptions["entry_method"],
            "transcription": exceptions.get("transcription"),
        }
    return summary
