"""Map internal workflow results to validated GA4GH VA-Spec Evidence Lines."""

from importlib.metadata import version

from ga4gh.va_spec.acmg_2015 import VariantPathogenicityEvidenceLine

from acmg.core.models import Status


SYSTEM = "ACMG Guidelines, 2015"
METHOD_TYPES = {
    "PVS1": "Null variant assessment",
    "PS1": "Same amino acid change assessment",
    "PM1": "Mutational hot spot and functional domain assessment",
    "PM2": "Population Data Assessment",
    "PM4": "Protein length change assessment",
    "PM5": "Novel missense position assessment",
    "PP2": "Variant spectrum assessment",
    "PP3": "In silico functional impact assessment",
    "PP5": "Reputable Source Assessment",
    "BA1": "Population Data Assessment",
    "BS1": "Population Data Assessment",
    "BP1": "Variant spectrum assessment",
    "BP3": "Protein length change assessment",
    "BP4": "In silico functional impact assessment",
    "BP6": "Reputable Source Assessment",
    "BP7": "Predicted silent variant assessment",
}
STRENGTHS = {
    "stand_alone": "standalone",
    "standalone": "standalone",
    "very_strong": "very strong",
    "strong": "strong",
    "moderate": "moderate",
    "supporting": "supporting",
}


def concept(code):
    return {"primaryCoding": {"system": SYSTEM, "code": code}}


def evidence_reference(item):
    identifier = item.get("evidence_id") if isinstance(item, dict) else None
    if not isinstance(identifier, str) or ":" not in identifier:
        raise ValueError("Every exported evidence item requires an IRI evidence_id")
    return identifier


def to_evidence_line(result):
    """Return None for workflow-only states; validate every exported line with VA-SPEC-Python."""
    status = Status(result.status)
    if status not in {Status.MET, Status.NOT_MET}:
        return None
    direction = result.direction
    if status == Status.NOT_MET:
        direction = "neutral"  # VA-Spec machine-readable enum; internal status remains NOT_MET.
    if direction not in {"supports", "disputes", "neutral"}:
        raise ValueError(f"Invalid VA-Spec direction for {result.criterion}")
    method_type = METHOD_TYPES[result.criterion]
    payload = {
        "type": "EvidenceLine",
        "name": f"{result.criterion} assessment for {result.variant['assembly']}:{result.variant['chrom']}:{result.variant['pos']}:{result.variant['ref']}:{result.variant['alt']}",
        "description": result.summary,
        "specifiedBy": {
            "type": "Method",
            "name": "ACMG/AMP 2015 with ClinGen General Guidance",
            "methodType": method_type,
            "reportedIn": {
                "type": "Document",
                "name": "Richards et al., 2015, Genet Med.",
                "doi": "10.1038/gim.2015.30",
                "pmid": "25741868",
                "urls": ["https://pubmed.ncbi.nlm.nih.gov/25741868/"],
            },
        },
        "directionOfEvidenceProvided": direction,
        "evidenceOutcome": concept(result.evidence_outcome),
    }
    if status == Status.MET:
        if result.strength not in STRENGTHS:
            raise ValueError(f"Unknown evidence strength for {result.criterion}")
        payload["strengthOfEvidenceProvided"] = concept(STRENGTHS[result.strength])
    references = list(dict.fromkeys(evidence_reference(item) for item in result.evidence))
    if references:
        payload["hasEvidenceItems"] = references
    model = VariantPathogenicityEvidenceLine.model_validate(payload)
    return model.model_dump(mode="json", exclude_none=True)


def export_record(record):
    lines = []
    evidence = {}
    for result in record["results"]:
        # Rehydrate only the fields needed by the mapper while retaining the validated internal result.
        from acmg.core.models import CriterionResult

        value = CriterionResult(**result)
        line = to_evidence_line(value)
        if line is None:
            continue
        lines.append({"criterion": value.criterion, "evidence_line": line})
        for item in value.evidence:
            identifier = evidence_reference(item)
            if identifier in evidence and evidence[identifier] != item:
                raise ValueError(f"Conflicting evidence object: {identifier}")
            evidence[identifier] = item
    return {"record_id": record["record_id"], "variant": record["variant"],
            "evidence_lines": lines, "referenced_evidence": evidence}


def export_document(records):
    return {
        "envelope_schema_version": "1.0",
        "profile": "Variant Pathogenicity Evidence Line (ACMG 2015)",
        "validated_by": {
            "package": "ga4gh.va-spec",
            "package_version": version("ga4gh.va-spec"),
            "model": "VariantPathogenicityEvidenceLine",
            "model_schema_id": VariantPathogenicityEvidenceLine.schema_id(),
        },
        "records": [export_record(record) for record in records],
    }
