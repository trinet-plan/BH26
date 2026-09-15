"""Validated local evidence transport for annotation, prediction and curation."""


class EvidenceService:
    def __init__(self, records):
        self.records = records

    def get(self, category, variant, context):
        selected = []
        for item in self.records:
            if item.get("category") != category or item.get("variant_key") != variant.key:
                continue
            if item.get("quality_status") != "PASS":
                continue
            if not all(item.get(key) for key in
                       ("evidence_id", "source", "source_version", "retrieved_at")):
                continue
            if item.get("transcript") and item["transcript"] != context.get("transcript"):
                continue
            if item.get("condition") and item["condition"] != context.get("condition"):
                continue
            selected.append(item)
        return selected
