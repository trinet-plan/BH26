"""
clinical_note.py
Classes for the structured data extracted from the free-text
"clinical_note" field of a democase/case*_api_input_*.json file: proband
phenotype/genotype, family history (relatives), and de novo (trio) data.

Extraction itself is done elsewhere via an LLM - this module only defines
the shape that LLM output is parsed into (ClinicalNoteExtraction.from_json).
Every field not directly stated in the note should come back None/empty
rather than guessed.

[Caveat: hpo_id, variant_status and zygosity are unverified]
  hpo_id (e.g. "HP:0001639") and zygosity/variant_status are whatever the
  LLM reports, taken at face value from the note's own wording - none of
  these are cross-checked against the real HPO ontology or any variant
  database; `label` and the free-text fields are the only ones directly
  grounded in the source text.

[Where this fits]
  The `de_novo` section (father/mother variant_status + phenotype +
  paternity/maternity confirmation) is exactly the data PS2 ("de novo,
  both parents confirmed, no family history of disease") and PM6 ("assumed
  de novo, parentage not confirmed") need - both are currently Layer-3
  "Evidence Gap / Manual Review Required" codes for this project (see
  CLAUDE.md's table), not yet automated; this module only extracts the raw
  fields, it doesn't itself decide PS2/PM6.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ClinicalFeature:
    label: str
    hpo_id: Optional[str] = None  # e.g. "HP:0001639"; None if not confidently known

    @staticmethod
    def from_json(data: dict) -> "ClinicalFeature":
        return ClinicalFeature(label=data["label"], hpo_id=data.get("hpo_id"))


@dataclass
class ProbandPhenotype:
    affected_status: Optional[bool] = None
    clinical_features: list[ClinicalFeature] = field(default_factory=list)
    age: Optional[str] = None
    age_of_onset: Optional[str] = None

    @staticmethod
    def from_json(data: dict) -> "ProbandPhenotype":
        return ProbandPhenotype(
            affected_status=data.get("affected_status"),
            clinical_features=[ClinicalFeature.from_json(f) for f in (data.get("clinical_features") or [])],
            age=data.get("age"),
            age_of_onset=data.get("age_of_onset"),
        )


@dataclass
class ProbandGenotype:
    variant_status: Optional[bool] = None
    zygosity: Optional[str] = None  # e.g. "heterozygous", "homozygous", "compound_heterozygous"

    @staticmethod
    def from_json(data: dict) -> "ProbandGenotype":
        return ProbandGenotype(variant_status=data.get("variant_status"), zygosity=data.get("zygosity"))


@dataclass
class Proband:
    phenotype: ProbandPhenotype = field(default_factory=ProbandPhenotype)
    genotype: ProbandGenotype = field(default_factory=ProbandGenotype)

    @staticmethod
    def from_json(data: dict) -> "Proband":
        return Proband(
            phenotype=ProbandPhenotype.from_json(data.get("phenotype") or {}),
            genotype=ProbandGenotype.from_json(data.get("genotype") or {}),
        )


@dataclass
class Relative:
    relationship: str  # e.g. "father", "mother", "son", "sibling", as stated in the note
    biological_relationship_confirmed: Optional[bool] = None
    affected_status: Optional[bool] = None
    clinical_features: list[ClinicalFeature] = field(default_factory=list)
    age: Optional[str] = None
    age_of_onset: Optional[str] = None
    variant_status: Optional[bool] = None
    zygosity: Optional[str] = None

    @staticmethod
    def from_json(data: dict) -> "Relative":
        return Relative(
            relationship=data.get("relationship", ""),
            biological_relationship_confirmed=data.get("biological_relationship_confirmed"),
            affected_status=data.get("affected_status"),
            clinical_features=[ClinicalFeature.from_json(f) for f in (data.get("clinical_features") or [])],
            age=data.get("age"),
            age_of_onset=data.get("age_of_onset"),
            variant_status=data.get("variant_status"),
            zygosity=data.get("zygosity"),
        )


@dataclass
class Family:
    inheritance_pattern: Optional[str] = None  # e.g. "autosomal dominant", "de novo", as stated/inferred in the note
    family_history_status: Optional[bool] = None
    relatives: list[Relative] = field(default_factory=list)

    @staticmethod
    def from_json(data: dict) -> "Family":
        return Family(
            inheritance_pattern=data.get("inheritance_pattern"),
            family_history_status=data.get("family_history_status"),
            relatives=[Relative.from_json(r) for r in (data.get("relatives") or [])],
        )


@dataclass
class DeNovo:
    father_variant_status: Optional[bool] = None
    mother_variant_status: Optional[bool] = None
    father_phenotype: Optional[str] = None
    mother_phenotype: Optional[str] = None
    paternity_confirmed: Optional[bool] = None
    maternity_confirmed: Optional[bool] = None

    @staticmethod
    def from_json(data: dict) -> "DeNovo":
        return DeNovo(
            father_variant_status=data.get("father_variant_status"),
            mother_variant_status=data.get("mother_variant_status"),
            father_phenotype=data.get("father_phenotype"),
            mother_phenotype=data.get("mother_phenotype"),
            paternity_confirmed=data.get("paternity_confirmed"),
            maternity_confirmed=data.get("maternity_confirmed"),
        )


@dataclass
class ClinicalNoteExtraction:
    proband: Proband = field(default_factory=Proband)
    family: Family = field(default_factory=Family)
    de_novo: DeNovo = field(default_factory=DeNovo)
    # 自由文の診断名(例: "hypertrophic cardiomyopathy")。VA-Spec Statementの
    # objectCondition(対象疾患)に使う - see z_tmp_va_spec_statement_decisions.md。
    # コード体系(MedGen/OMIM等)への変換は行わず、この自由文をそのまま
    # MappableConcept.name に載せる想定。
    diagnosis: Optional[str] = None

    @staticmethod
    def from_json(data: dict) -> "ClinicalNoteExtraction":
        return ClinicalNoteExtraction(
            proband=Proband.from_json(data.get("proband") or {}),
            family=Family.from_json(data.get("family") or {}),
            de_novo=DeNovo.from_json(data.get("de_novo") or {}),
            diagnosis=data.get("diagnosis"),
        )


def extract_clinical_note(text: str) -> ClinicalNoteExtraction:
    """Boundary for the clinical-note extraction component.

    The public API intentionally accepts the main project's original free-text
    ``clinical_note`` field.  Until the dedicated extractor supplied by the
    clinical-note workstream is connected here, do not infer clinical facts
    from prose: return an empty extraction and let dependent criteria be
    reported as UNKNOWN.
    """
    if not isinstance(text, str):
        raise TypeError("clinical_note must be a string")
    return ClinicalNoteExtraction()
