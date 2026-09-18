"""
clinical_extraction.py

Extract structured clinical-note facts into the shared ClinicalNoteExtraction
model defined in acmg_pipeline.clinical_note.

Primary API:
    extract_clinical_note(...)-> ClinicalNoteExtraction

The LLM returns JSON text only as an interchange format. The public function
parses that JSON into ClinicalNoteExtraction and returns the class instance.

Important:
- HPO IDs are NOT generated here. hpo_id must be null at this stage.
- ACMG criteria are NOT assigned here.
- affected/family-history/variant status are never guessed from missing data.
- If target_variant_label is omitted, all target-specific variant_status and
  zygosity fields are forced to None.

Because ApiCaseInput carries exactly one VCF variant per case, a caller may
supply a human-readable target_variant_label for that single variant. Then
variant_status/zygosity refer only to that target variant.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import asdict
from pathlib import Path
from typing import Any, Optional

from openai import OpenAI

from acmg_pipeline.api_input import ApiCaseInput
from acmg_pipeline.clinical_note import ClinicalNoteExtraction


# parents[1]: this file lives at <repo root>/acmg_pipeline/clinical_extraction.py,
# so one parent up is <repo root> (where .env lives). The pp4_pp1_bs4 branch had
# parents[3] here, which only resolves correctly under a different, deeper
# checkout layout than this repo's own - fixed when pulling this file into main
# (2026-09-17).
ROOT_DIR = Path(__file__).resolve().parents[1]


def _load_env_file(path: Path = ROOT_DIR / ".env") -> None:
    if not path.exists():
        return

    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip())


_load_env_file()

VLLM_BASE_URL = os.environ.get("VLLM_BASE_URL", "")
VLLM_API_KEY = os.environ.get("VLLM_API_KEY", "")
VLLM_MODEL = os.environ.get("VLLM_MODEL", "google/gemma-4-26B-A4B-it")

if not VLLM_BASE_URL or not VLLM_API_KEY:
    raise RuntimeError(
        "VLLM_BASE_URL / VLLM_API_KEY are not configured in the repository-root .env"
    )


SYSTEM_PROMPT = """\
You are a clinical genetics information-extraction assistant.

Extract only facts explicitly supported by the supplied clinical note.
Return data that matches the ClinicalNoteExtraction schema exactly.

Rules:
- Do not invent or repair missing information.
- Do not generate HPO IDs. Every hpo_id must be null.
- Do not generate a MONDO identifier. condition_id must always be null; it is
  resolved separately, from the diagnosis text you extract, against a real
  ontology - never guess or recall one from memory.
- Do not assign ACMG/AMP criteria or variant pathogenicity.
- Do not infer inheritance pattern from the pedigree. Populate inheritance_pattern
  only when the note explicitly states an inheritance pattern.
- Do not infer biological relationship confirmation. Parent/child wording alone is
  not confirmation; biological_relationship_confirmed must remain null unless the
  note explicitly confirms or refutes it.
- Do not infer paternity or maternity confirmation.
- affected_status is true only when affected status is explicitly supported.
- affected_status is false only when the note explicitly says unaffected or gives a
  direct clinical assessment showing absence of the relevant disease phenotype.
- If affected status is unclear, use null.
- "no reported history" alone is not sufficient to call a person unaffected.
- family_history_status is true only for an explicit relevant positive family history,
  false only for an explicit relevant negative family history, otherwise null.
- clinical_features must contain positive clinical findings only. The current
  ClinicalFeature model has no present/absent field, so do not put a negated finding
  such as "no myocardial hypertrophy" into clinical_features.
- age and age_of_onset should be short strings close to the wording in the note,
  e.g. "63 years" or "childhood". Use null if not stated.
- variant_status and zygosity refer only to the target variant supplied by the user.
- If no target variant is supplied, variant_status and zygosity must be null for the
  proband, relatives, father, and mother.
- If a target variant is supplied but the note does not explicitly say whether that
  person carries it, variant_status must be null, not false.
- father_phenotype and mother_phenotype are concise source-grounded descriptions;
  use null if unavailable.
- When the note describes a parent-related fact (e.g. age at death, cause of death) without
  stating which parent it applies to, do not guess or assign it to a specific parent. Record
  it as a family.relatives entry with relationship "parent (unspecified)" instead of "father"
  or "mother".
- diagnosis is the proband's overall clinical diagnosis / condition name, using the
  note's own wording (e.g. "hypertrophic cardiomyopathy"), not an ACMG/variant
  classification. Extract it whenever the note states one, including a hedged
  impression such as "likely familial HCM" - still extract the named condition
  ("hypertrophic cardiomyopathy") in that case, since the hedge is about certainty,
  not about which disease is being discussed. Use null only when no diagnosis or
  named condition is discussed at all.

Return ONLY valid JSON. Do not wrap it in Markdown.
"""


def build_user_prompt(
    clinical_note: str,
    target_variant_label: Optional[str] = None,
) -> str:
    target_text = (
        target_variant_label.strip()
        if target_variant_label and target_variant_label.strip()
        else "NONE PROVIDED"
    )

    return f"""\
Target variant for target-specific genotype fields:
{target_text}

Use exactly this JSON shape:
{{
  "proband": {{
    "phenotype": {{
      "affected_status": null,
      "clinical_features": [
        {{
          "label": "...",
          "hpo_id": null
        }}
      ],
      "age": null,
      "age_of_onset": null
    }},
    "genotype": {{
      "variant_status": null,
      "zygosity": null
    }}
  }},
  "family": {{
    "inheritance_pattern": null,
    "family_history_status": null,
    "relatives": [
      {{
        "relationship": "...",
        "biological_relationship_confirmed": null,
        "affected_status": null,
        "clinical_features": [
          {{
            "label": "...",
            "hpo_id": null
          }}
        ],
        "age": null,
        "age_of_onset": null,
        "variant_status": null,
        "zygosity": null
      }}
    ]
  }},
  "de_novo": {{
    "father_variant_status": null,
    "mother_variant_status": null,
    "father_phenotype": null,
    "mother_phenotype": null,
    "paternity_confirmed": null,
    "maternity_confirmed": null
  }},
  "diagnosis": null,
  "condition_id": null
}}

Boolean semantics:
- affected_status: true=affected, false=explicitly unaffected, null=unknown.
- variant_status: true=target variant present, false=target variant explicitly absent,
  null=unknown/not assessed.
- family_history_status: true=explicit relevant positive family history,
  false=explicit relevant negative family history, null=unknown.

----- BEGIN CLINICAL NOTE -----
{clinical_note}
----- END CLINICAL NOTE -----
"""


def make_client() -> OpenAI:
    return OpenAI(base_url=VLLM_BASE_URL, api_key=VLLM_API_KEY)


def _extract_json_object(text: str) -> dict[str, Any]:
    text = (text or "").strip()
    if not text:
        raise ValueError("LLM returned an empty response")

    text = re.sub(r"^```(?:json)?\\s*", "", text, flags=re.IGNORECASE)
    text = re.sub(r"\\s*```$", "", text)

    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        start = text.find("{")
        if start < 0:
            raise ValueError("LLM output did not contain a JSON object")
        value, _ = json.JSONDecoder().raw_decode(text[start:])

    if not isinstance(value, dict):
        raise ValueError("LLM output JSON is not an object")
    return value


def _force_no_target_variant(data: dict[str, Any]) -> None:
    """Clear target-specific fields when no VCF target was supplied."""
    proband = data.get("proband")
    if isinstance(proband, dict):
        genotype = proband.get("genotype")
        if isinstance(genotype, dict):
            genotype["variant_status"] = None
            genotype["zygosity"] = None

    family = data.get("family")
    if isinstance(family, dict):
        relatives = family.get("relatives")
        if isinstance(relatives, list):
            for relative in relatives:
                if isinstance(relative, dict):
                    relative["variant_status"] = None
                    relative["zygosity"] = None

    de_novo = data.get("de_novo")
    if isinstance(de_novo, dict):
        de_novo["father_variant_status"] = None
        de_novo["mother_variant_status"] = None


def _force_hpo_ids_null(data: dict[str, Any]) -> None:
    """Never trust HPO IDs emitted during the LLM extraction step."""
    proband = data.get("proband")
    if isinstance(proband, dict):
        phenotype = proband.get("phenotype")
        if isinstance(phenotype, dict):
            for feature in phenotype.get("clinical_features") or []:
                if isinstance(feature, dict):
                    feature["hpo_id"] = None

    family = data.get("family")
    if isinstance(family, dict):
        for relative in family.get("relatives") or []:
            if not isinstance(relative, dict):
                continue
            for feature in relative.get("clinical_features") or []:
                if isinstance(feature, dict):
                    feature["hpo_id"] = None


def _force_condition_id_null(data: dict[str, Any]) -> None:
    """Never trust a MONDO ID emitted during the LLM extraction step.

    Same reasoning as _force_hpo_ids_null(): this step has no ontology to
    check itself against, so any condition_id it produced would be a raw
    LLM guess - confirmed in real testing to be wrong for at least one demo
    case (see acmg_pipeline.hpo_mondo_extraction.resolve_diagnosis_mondo()'s
    own docstring). condition_id is resolved separately, against real EBI
    OLS4 candidates, from the diagnosis text this step does extract.
    """
    data["condition_id"] = None


def extract_clinical_note(
    clinical_note: str,
    *,
    target_variant_label: Optional[str] = None,
    client: Optional[OpenAI] = None,
    model: str = VLLM_MODEL,
) -> ClinicalNoteExtraction:
    """Extract one clinical note into the shared ClinicalNoteExtraction model."""
    if not clinical_note.strip():
        return ClinicalNoteExtraction()

    llm_client = client or make_client()
    response = llm_client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {
                "role": "user",
                "content": build_user_prompt(
                    clinical_note,
                    target_variant_label=target_variant_label,
                ),
            },
        ],
    )

    raw_text = (response.choices[0].message.content or "").strip()
    data = _extract_json_object(raw_text)

    _force_hpo_ids_null(data)
    _force_condition_id_null(data)
    if not target_variant_label:
        _force_no_target_variant(data)

    return ClinicalNoteExtraction.from_json(data)


def extract_api_case(
    case_input: ApiCaseInput,
    *,
    target_variant_label: Optional[str] = None,
    client: Optional[OpenAI] = None,
    model: str = VLLM_MODEL,
) -> ClinicalNoteExtraction:
    """Extract the clinical_note field from one ApiCaseInput fixture."""
    return extract_clinical_note(
        case_input.clinical_note,
        target_variant_label=target_variant_label,
        client=client,
        model=model,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract ClinicalNoteExtraction from one case api_input JSON fixture."
    )
    parser.add_argument(
        "api_input",
        type=Path,
        help="Path to democase/case*_api_input_*.json",
    )
    parser.add_argument(
        "--target-variant",
        default=None,
        help="Human-readable description of the single VCF variant for target-specific status fields.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional debug JSON output path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    case_input = ApiCaseInput.from_json_file(args.api_input)
    result = extract_api_case(
        case_input,
        target_variant_label=args.target_variant,
    )

    # JSON is only for CLI/debug display. The Python API returns the class object.
    payload = json.dumps(asdict(result), ensure_ascii=False, indent=2)
    print(payload)

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
        print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
