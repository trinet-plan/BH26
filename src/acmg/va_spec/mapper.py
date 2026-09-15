"""Map internal workflow results to validated GA4GH VA-Spec Evidence Lines."""

import hashlib
import json
import re
from functools import lru_cache
from importlib.metadata import version
from importlib.resources import files

from ga4gh.va_spec.acmg_2015 import VariantPathogenicityEvidenceLine
from jsonschema import Draft202012Validator, FormatChecker

from acmg.core.models import Status


SYSTEM = "ACMG Guidelines, 2015"
SCHEMA_VERSION = "1.0.1"
SCHEMA_ID = ("https://w3id.org/ga4gh/schema/va-spec/1.0.1/acmg-2015/json/"
             "VariantPathogenicityEvidenceLine")
OUTCOME_PATTERN = re.compile(
    r"^((?:PVS1)(?:_(?:not_met|(?:strong|moderate|supporting)))?|"
    r"(?:PS[1-4]|BS[1-4])(?:_(?:not_met|(?:very_strong|moderate|supporting)))?|"
    r"BA1(?:_not_met)?|(?:PM[1-6])(?:_(?:not_met|(?:very_strong|strong|supporting)))?|"
    r"(PP[1-5]|BP[1-7])(?:_(?:not_met|very_strong|strong|moderate))?)$"
)
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


@lru_cache(maxsize=1)
def output_schema():
    resource = files("acmg.va_spec").joinpath("schemas/acmg-evidence-line-1.0.1-output.json")
    schema = json.loads(resource.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return schema


def output_schema_sha256():
    resource = files("acmg.va_spec").joinpath("schemas/acmg-evidence-line-1.0.1-output.json")
    return hashlib.sha256(resource.read_bytes()).hexdigest()


def validate_1_0_1(line, criterion):
    """Validate the emitted subset and ACMG cross-field semantics for VA-Spec 1.0.1."""
    errors = sorted(Draft202012Validator(
        output_schema(), format_checker=FormatChecker()
    ).iter_errors(line), key=lambda item: list(item.path))
    if errors:
        raise ValueError(f"VA-Spec 1.0.1 schema validation failed: {errors[0].message}")
    method_type = line["specifiedBy"]["methodType"]
    outcome = line["evidenceOutcome"]["primaryCoding"]["code"]
    direction = line["directionOfEvidenceProvided"]
    if method_type != criterion or not OUTCOME_PATTERN.fullmatch(outcome):
        raise ValueError("VA-Spec 1.0.1 ACMG criterion mapping mismatch")
    if outcome.split("_", 1)[0] != criterion:
        raise ValueError("VA-Spec 1.0.1 methodType/evidenceOutcome mismatch")
    expected = "neutral" if outcome.endswith("_not_met") else (
        "disputes" if criterion.startswith("B") else "supports"
    )
    if direction != expected:
        raise ValueError("VA-Spec 1.0.1 direction/evidenceOutcome mismatch")
    return line


def to_evidence_line(result):
    """Return None for workflow states; normalize then validate the 1.0.1 output."""
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
    # The installed reference model normalizes core 1.0-compatible structures, but its
    # newer snapshot uses semantic method categories. VA-Spec 1.0.1 uses ACMG codes here.
    model = VariantPathogenicityEvidenceLine.model_validate(payload)
    line = model.model_dump(mode="json", exclude_none=True)
    line["specifiedBy"]["methodType"] = result.criterion
    # GKS-Core 1.0.0 (referenced by VA-Spec 1.0.1) predates the required
    # MappableConcept type discriminator added in newer GKS-Core snapshots.
    line["evidenceOutcome"].pop("type", None)
    if "strengthOfEvidenceProvided" in line:
        line["strengthOfEvidenceProvided"].pop("type", None)
    return validate_1_0_1(line, result.criterion)


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
            "schema_version": SCHEMA_VERSION,
            "schema_id": SCHEMA_ID,
            "validator": "jsonschema Draft202012Validator plus ACMG cross-field checks",
            "validation_scope": "Each evidence_line; enclosing records document is a BH26 envelope",
            "output_profile_schema_sha256": output_schema_sha256(),
            "structural_normalizer": "ga4gh.va-spec",
            "structural_normalizer_version": version("ga4gh.va-spec"),
            "model": "VariantPathogenicityEvidenceLine",
            "structural_normalizer_model_schema_id": VariantPathogenicityEvidenceLine.schema_id(),
        },
        "records": [export_record(record) for record in records],
    }
