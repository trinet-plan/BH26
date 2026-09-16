"""Shared structured clinical-note types used at criterion boundaries."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ClinicalFeature:
    label: str
    hpo_id: Optional[str] = None

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
    zygosity: Optional[str] = None

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
    relationship: str
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
    inheritance_pattern: Optional[str] = None
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

    @staticmethod
    def from_json(data: dict) -> "ClinicalNoteExtraction":
        return ClinicalNoteExtraction(
            proband=Proband.from_json(data.get("proband") or {}),
            family=Family.from_json(data.get("family") or {}),
            de_novo=DeNovo.from_json(data.get("de_novo") or {}),
        )
