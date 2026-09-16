"""Validated local evidence transport for annotation, prediction and curation."""


class EvidenceService:
    def __init__(self, records):
        self.records = records

    def get(self, category, variant, context):
        selected = []
        for item in self.get_candidates(category, variant):
            if item.get("transcript") and item["transcript"] != context.get("transcript"):
                continue
            if item.get("condition") and item["condition"] != context.get("condition"):
                continue
            selected.append(item)
        return selected

    def get_candidates(self, category, variant):
        """Return provenance-complete records before transcript/condition resolution.

        Most criteria can use the default context filter in ``get``. PVS1 must explicitly
        prefer condition-specific mechanism evidence and then fall back to gene-level evidence,
        so it needs to inspect both scopes without accidentally borrowing another condition.
        """
        return [item for item in self.records
                if item.get("category") == category
                and item.get("variant_key") == variant.key
                and item.get("quality_status") == "PASS"
                and all(item.get(key) for key in
                        ("evidence_id", "source", "source_version", "retrieved_at"))]
