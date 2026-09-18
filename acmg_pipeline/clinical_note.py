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

import re
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
    # 表示・文献検索には自由文を、疾患特異的な自動判定にはMONDO IDを使う。
    # 疾患情報はVariantRecord.infoに書き戻さない。
    diagnosis: Optional[str] = None
    # diagnosisをMONDOに正規化したID(例: "MONDO:0007739")。
    #
    # 生産者は2つある。clinical-note parserが抽出時に埋める経路と、
    # acmg_pipeline.hpo_mondo_extraction.resolve_diagnosis_mondo()がTogoMCP経由で
    # 後から埋める経路。どちらも同じ「この症例が対象とする疾患」を指すので、
    # フィールドは1つに保つ(2026-09-18のマージで mondo_id と condition_id の2案が
    # 並行実装されたが、読む側 - automated_core.interface.criterion_input() が
    # data["condition"]へ渡し、PP2/BP1/PM1/PVS1の疾患照合を駆動する - は1つしか
    # 無いため、下流の`condition`と名前を揃えたこちらに寄せた)。
    # 解決できない場合はNoneのまま。diagnosisは決して上書きしない。
    condition_id: Optional[str] = None

    def __post_init__(self) -> None:
        # どちらの生産者が埋めても、ここを通る。OMIM/Orphanet/MedGen/自由文が
        # 紛れ込むと疾患照合が静かに外れるので、形式は型の側で強制する。
        if self.condition_id is not None and not re.fullmatch(r"MONDO:\d+", self.condition_id):
            raise ValueError(
                "clinical-note condition_id must be a MONDO identifier such as MONDO:0005045"
            )

    @staticmethod
    def from_json(data: dict) -> "ClinicalNoteExtraction":
        return ClinicalNoteExtraction(
            proband=Proband.from_json(data.get("proband") or {}),
            family=Family.from_json(data.get("family") or {}),
            de_novo=DeNovo.from_json(data.get("de_novo") or {}),
            diagnosis=data.get("diagnosis"),
            condition_id=data.get("condition_id"),
        )


def extract_clinical_note(text: str) -> ClinicalNoteExtraction:
    """Boundary for the clinical-note extraction component.

    The public API intentionally accepts the main project's original free-text
    ``clinical_note`` field. Delegates to acmg_pipeline.clinical_extraction's
    real LLM-based extractor (pulled into main 2026-09-17, from the
    clinical-note workstream's pp4_pp1_bs4 branch) - imported lazily, inside
    this function, rather than at module level: clinical_extraction.py raises
    RuntimeError at import time if VLLM_BASE_URL/VLLM_API_KEY aren't
    configured, and this module (clinical_note.py) is imported broadly enough
    (by segregation.py, pipeline.py, every criterion module, ...) that forcing
    that requirement onto every one of those imports would be a real
    regression for any caller that never actually needs to extract a note.
    """
    if not isinstance(text, str):
        raise TypeError("clinical_note must be a string")
    from acmg_pipeline.clinical_extraction import extract_clinical_note as _extract
    return _extract(text)
