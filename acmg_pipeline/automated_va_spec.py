from acmg_pipeline.constants import CriterionStatus
"""Map internal workflow results to validated GA4GH VA-Spec Evidence Lines."""

import hashlib
import json
import re
from decimal import Decimal, InvalidOperation
from functools import lru_cache
from importlib.metadata import version
from importlib.resources import files

from ga4gh.va_spec.acmg_2015 import VariantPathogenicityEvidenceLine
from jsonschema import Draft202012Validator, FormatChecker

from acmg_pipeline.automated_core.models import CRITERIA


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


def concept(code, name=None):
    value = {"primaryCoding": {"system": SYSTEM, "code": code}}
    if name:
        value["name"] = name
    return value


def evidence_reference(item):
    identifier = item.get("evidence_id") if isinstance(item, dict) else None
    if not isinstance(identifier, str) or ":" not in identifier:
        raise ValueError("Every exported evidence item requires an IRI evidence_id")
    return identifier


def stable_urn(kind, value):
    digest = hashlib.sha256(value.encode()).hexdigest()
    return f"urn:bh26:{kind}:{digest}"


def population_study_result(item, variant):
    """Represent a normalized population observation like the official gnomAD example."""
    try:
        ac = int(item["AC"])
        an = int(item["AN"])
        af = float(Decimal(str(item["AF"])))
    except (KeyError, TypeError, ValueError, InvalidOperation) as exc:
        raise ValueError("Population Evidence cannot form a VA-Spec StudyResult") from exc
    if ac < 0 or an <= 0 or not 0 <= af <= 1:
        raise ValueError("Invalid population values for VA-Spec StudyResult")
    source = item.get("source")
    source_version = item.get("source_version")
    population = item.get("population")
    if not all(isinstance(value, str) and value for value in
               (source, source_version, population)):
        raise ValueError("Population Evidence provenance is incomplete")
    variant_key = f"{variant['assembly']}:{variant['chrom']}:{variant['pos']}:{variant['ref']}:{variant['alt']}"
    dataset_iri = ("https://gnomad.broadinstitute.org/"
                   f"?dataset=gnomad_r4&version={source_version}")
    return {
        "id": evidence_reference(item),
        "type": "CohortAlleleFrequencyStudyResult",
        "name": f"{source} {population} allele frequency for {variant_key}",
        "focusAllele": stable_urn("variant", variant_key),
        "focusAlleleFrequency": af,
        "focusAlleleCount": ac,
        "locusAlleleCount": an,
        "sourceDataSet": {
            "id": dataset_iri, "type": "DataSet",
            "name": f"{source} v{source_version}", "version": source_version,
        },
        "cohort": {
            "id": stable_urn("cohort", f"{source}:{source_version}:{population}"),
            "type": "StudyGroup", "name": population,
        },
        "specifiedBy": {
            "type": "Method", "name": "gnomAD browser allele frequency calculation",
            "reportedIn": {
                "type": "Document", "name": "gnomAD browser help",
                "urls": ["https://gnomad.broadinstitute.org/help"],
            },
        },
        "qualityMeasures": {
            "qualityStatus": item.get("quality_status"),
            "callable": item.get("callable"),
            "filters": item.get("filters", []),
            "variantFlags": item.get("variant_flags", []),
            "retrievedAt": item.get("retrieved_at"),
        },
    }


def va_evidence_item(item, variant):
    if isinstance(item, dict) and item.get("category") == "population":
        return population_study_result(item, variant)
    return evidence_reference(item)


_EVIDENCE_NAMES = {
    "annotation": "Transcript variant annotation",
    "comparator": "Pathogenic comparator assessment",
    "comparator_search": "ClinVar comparator search result",
    "computational": "Computational prediction result",
    "gene_disease": "Gene-disease mechanism assessment",
    "hotspot_search": "ClinVar regional search result",
    "nmd_prediction": "Nonsense-mediated decay prediction",
    "population_lof": "Regional loss-of-function population assessment",
    "protein_region": "Protein region impact assessment",
    "region": "Protein region assessment",
    "rna_assay": "RNA assay result",
    "splice_assessment": "Splice consequence assessment",
    "synonymous_assessment": "Synonymous variant assessment",
    "transcript_assessment": "Transcript and exon relevance assessment",
}


def evidence_catalog_item(item):
    """Build a resolvable audit object for an EvidenceLine IRI reference.

    Only profiles standardized by VA-Spec are embedded in the validated EvidenceLine itself.
    Other normalized evidence records are represented in the BH26 envelope as a generic
    StudyResult-shaped object. Domain facts live in an Extension so they are not presented as
    fields from an official VA-Spec StudyResult profile.
    """
    identifier = evidence_reference(item)
    category = item.get("category", "evidence")
    source = item.get("source")
    version_value = item.get("source_version")
    facts = {
        key: value for key, value in item.items()
        if key not in {
            "evidence_id", "category", "source", "source_version", "retrieved_at",
            "quality_status", "curator", "reviewed_at", "assessment_method", "method",
            "policy_version",
        }
    }
    value = {
        "id": identifier,
        "type": "StudyResult",
        "name": _EVIDENCE_NAMES.get(category, f"{category.replace('_', ' ').title()} result"),
        "description": (
            f"Normalized {category} evidence used by the criterion evaluator; "
            "domain-specific values are carried in the observations extension."
        ),
        "extensions": [
            {"name": "bh26EvidenceCategory", "value": category},
            {"name": "observations", "value": facts},
            {"name": "normalizedEvidenceSha256", "value": hashlib.sha256(
                json.dumps(item, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=False, allow_nan=False).encode()
            ).hexdigest()},
        ],
        "qualityMeasures": {
            "qualityStatus": item.get("quality_status"),
            "retrievedAt": item.get("retrieved_at"),
        },
    }
    if isinstance(source, str) and source:
        dataset_key = f"{source}:{version_value or 'unversioned'}"
        value["sourceDataSet"] = {
            "id": stable_urn("dataset", dataset_key),
            "type": "DataSet",
            "name": source if not version_value else f"{source} {version_value}",
        }
        if isinstance(version_value, str) and version_value:
            value["sourceDataSet"]["version"] = version_value
    method = item.get("method") or item.get("assessment_method")
    if method:
        value["specifiedBy"] = {
            "type": "Method",
            "name": str(method),
            "extensions": ([{"name": "policyVersion", "value": item["policy_version"]}]
                           if item.get("policy_version") else []),
        }
    if item.get("curator"):
        contribution = {
            "type": "Contribution",
            "contributor": {"type": "Agent", "name": item["curator"]},
            "activityType": "evidence evaluation",
        }
        if item.get("reviewed_at"):
            contribution["date"] = item["reviewed_at"]
        value["contributions"] = [contribution]
    reported = item.get("primary_evidence")
    if isinstance(reported, list) and reported:
        value["reportedIn"] = list(dict.fromkeys(
            entry for entry in reported if isinstance(entry, str) and entry
        ))
    return value


def assessment_details(result):
    """Return the complete workflow explanation shared by all criterion outputs."""
    value = {
        "criterion": result.criterion,
        "status": result.status.value,
        "summary": result.summary,
        "evidenceItemIds": list(dict.fromkeys(
            evidence_reference(item) for item in result.evidence
        )),
        "provenance": result.provenance,
    }
    optional = {
        "strength": result.strength,
        "direction": result.direction,
        "evidenceOutcome": result.evidence_outcome,
        "missingInputs": result.missing_inputs,
        "reviewPoints": result.review_points,
        "conflictFlags": result.conflict_flags,
        "evaluationContext": result.evaluation_context,
        "decisionTrace": result.decision_trace,
        "rulesUsed": result.rules_used,
        "warnings": result.warnings,
        "unresolvedRequirements": result.unresolved_requirements,
    }
    value.update({key: field_value for key, field_value in optional.items() if field_value})
    return value


@lru_cache(maxsize=1)
def output_schema():
    resource = files("acmg_pipeline").joinpath("schemas/acmg-evidence-line-1.0.1-output.json")
    schema = json.loads(resource.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return schema


def output_schema_sha256():
    resource = files("acmg_pipeline").joinpath("schemas/acmg-evidence-line-1.0.1-output.json")
    return hashlib.sha256(resource.read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def envelope_schema():
    resource = files("acmg_pipeline").joinpath("schemas/bh26-audit-envelope-1.1.json")
    schema = json.loads(resource.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return schema


def envelope_schema_sha256():
    resource = files("acmg_pipeline").joinpath("schemas/bh26-audit-envelope-1.1.json")
    return hashlib.sha256(resource.read_bytes()).hexdigest()


def validate_envelope(document):
    errors = sorted(Draft202012Validator(
        envelope_schema(), format_checker=FormatChecker()
    ).iter_errors(document), key=lambda item: list(item.path))
    if errors:
        raise ValueError(f"BH26 audit envelope validation failed: {errors[0].message}")
    for record in document["records"]:
        assessments = record["criterion_assessments"]
        codes = [item["criterion"] for item in assessments]
        if len(codes) != len(set(codes)) or any(code not in CRITERIA for code in codes):
            raise ValueError("BH26 audit envelope has unknown or duplicate criterion assessment")
        by_code = {item["criterion"]: item for item in assessments}
        catalog = record["referenced_evidence"]
        for identifier, item in catalog.items():
            if item["id"] != identifier:
                raise ValueError("BH26 audit Evidence Item key/id mismatch")
        for assessment in assessments:
            missing = set(assessment["evidenceItemIds"]) - set(catalog)
            if missing:
                raise ValueError(f"Unresolved Evidence Item references: {sorted(missing)}")
        line_ids = set()
        for wrapped in record["evidence_lines"]:
            criterion = wrapped["criterion"]
            if criterion not in by_code or wrapped["assessment_details"] != by_code[criterion]:
                raise ValueError("EvidenceLine assessment details disagree with criterion audit")
            line = wrapped["evidence_line"]
            if line["id"] in line_ids:
                raise ValueError("Duplicate EvidenceLine id in one record")
            line_ids.add(line["id"])
            details = next((item["value"] for item in line.get("extensions", [])
                            if item.get("name") == "bh26AssessmentDetails"), None)
            if details != by_code[criterion]:
                raise ValueError("EvidenceLine extension disagrees with criterion audit")
    return document


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
    status = CriterionStatus(result.status)
    if status not in {CriterionStatus.MET, CriterionStatus.NOT_MET}:
        return None
    direction = result.direction
    if status == CriterionStatus.NOT_MET:
        direction = "neutral"  # VA-Spec machine-readable enum; internal status remains NOT_MET.
    if direction not in {"supports", "disputes", "neutral"}:
        raise ValueError(f"Invalid VA-Spec direction for {result.criterion}")
    method_type = METHOD_TYPES[result.criterion]
    payload = {
        "type": "EvidenceLine",
        "id": stable_urn(
            "evidence-line",
            f"{result.criterion}:{result.evidence_outcome}:{json.dumps(result.variant, sort_keys=True)}",
        ),
        "name": f"{result.criterion} assessment for {result.variant['assembly']}:{result.variant['chrom']}:{result.variant['pos']}:{result.variant['ref']}:{result.variant['alt']}",
        "description": result.summary,
        "extensions": [{"name": "bh26AssessmentDetails",
                        "value": assessment_details(result)}],
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
        "evidenceOutcome": concept(
            result.evidence_outcome,
            f"ACMG 2015 {result.criterion} criterion " +
            ("met" if status == CriterionStatus.MET else "not met"),
        ),
    }
    if status == CriterionStatus.MET:
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
    if result.evidence:
        line["hasEvidenceItems"] = [va_evidence_item(item, result.variant)
                                    for item in result.evidence]
    # GKS-Core 1.0.0 (referenced by VA-Spec 1.0.1) predates the required
    # MappableConcept type discriminator added in newer GKS-Core snapshots.
    line["evidenceOutcome"].pop("type", None)
    if "strengthOfEvidenceProvided" in line:
        line["strengthOfEvidenceProvided"].pop("type", None)
    return validate_1_0_1(line, result.criterion)


def export_record(record):
    lines = []
    evidence = {}
    assessments = []
    for result in record["results"]:
        # Rehydrate only the fields needed by the mapper while retaining the validated internal result.
        from acmg_pipeline.automated_core.models import CriterionResult

        value = CriterionResult(**result)
        details = assessment_details(value)
        assessments.append(details)
        for item in value.evidence:
            identifier = evidence_reference(item)
            if identifier in evidence and evidence[identifier] != item:
                raise ValueError(f"Conflicting evidence object: {identifier}")
            evidence[identifier] = item
        line = to_evidence_line(value)
        if line is None:
            continue
        lines.append({"criterion": value.criterion, "evidence_line": line,
                      "assessment_details": details})
    return {"record_id": record["record_id"], "variant": record["variant"],
            "criterion_assessments": assessments, "evidence_lines": lines,
            "referenced_evidence": {
                identifier: evidence_catalog_item(item) for identifier, item in evidence.items()
            }}


def export_document(records):
    document = {
        "envelope_schema_version": "1.1",
        "profile": "Variant Pathogenicity Evidence Line (ACMG 2015)",
        "validated_by": {
            "schema_version": SCHEMA_VERSION,
            "schema_id": SCHEMA_ID,
            "validator": "jsonschema Draft202012Validator plus ACMG cross-field checks",
            "validation_scope": "Each evidence_line; enclosing records document is a BH26 envelope",
            "output_profile_schema_sha256": output_schema_sha256(),
            "audit_envelope_schema_version": "1.1",
            "audit_envelope_schema_id": envelope_schema()["$id"],
            "audit_envelope_schema_sha256": envelope_schema_sha256(),
            "structural_normalizer": "ga4gh.va-spec",
            "structural_normalizer_version": version("ga4gh.va-spec"),
            "model": "VariantPathogenicityEvidenceLine",
            "structural_normalizer_model_schema_id": VariantPathogenicityEvidenceLine.schema_id(),
        },
        "records": [export_record(record) for record in records],
    }
    return validate_envelope(document)
