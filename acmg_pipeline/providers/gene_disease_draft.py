"""Create review-only gene-disease mechanism drafts from aggregate data.

The output is deliberately not ``gene_disease`` evidence.  Constraint metrics,
database coverage and ClinVar counts can prioritize review, but they cannot by
themselves establish a disease mechanism for PP2, BP1 or PVS1.
"""

from __future__ import annotations

import hashlib
import json
import math


POSITIVE_VALIDITY = {"definitive", "strong", "moderate"}
LIMITED_VALIDITY = {"limited"}
NEGATIVE_VALIDITY = {"disputed", "refuted", "no_known_disease_relationship"}
PROVENANCE_FIELDS = ("source", "source_version", "retrieved_at")


def _canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _require_provenance(record, label):
    if not isinstance(record, dict) or not all(record.get(key) for key in PROVENANCE_FIELDS):
        raise ValueError(f"{label} requires source, source_version and retrieved_at")


def _number(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number")
    return float(value)


def _count(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{label} must be a non-negative integer")
    return value


def _normalized_validity(value):
    return str(value).strip().lower().replace(" ", "_").replace("-", "_")


class GeneDiseaseDraftProvider:
    """Generate transparent triage suggestions that always require human review."""

    name = "BH26 gene-disease draft generator"

    def __init__(self, policy):
        if not isinstance(policy, dict):
            raise ValueError("Draft policy must be an object")
        required = (
            "policy_version",
            "policy_source",
            "pp2_min_mis_z",
            "bp1_max_missense_fraction",
            "bp1_min_truncating_count",
            "pvs1_min_p_li",
            "pvs1_max_loeuf",
        )
        if not all(policy.get(key) is not None for key in required):
            raise ValueError("Draft policy is incomplete")
        self.policy = dict(policy)
        for key in (
            "pp2_min_mis_z", "bp1_max_missense_fraction",
            "pvs1_min_p_li", "pvs1_max_loeuf",
        ):
            self.policy[key] = _number(self.policy[key], key)
        self.policy["bp1_min_truncating_count"] = _count(
            self.policy["bp1_min_truncating_count"],
            "bp1_min_truncating_count",
        )

    @staticmethod
    def _validity_state(records, gene, condition):
        if not isinstance(records, list):
            raise ValueError("validity must be a list")
        matched = []
        for record in records:
            _require_provenance(record, "validity record")
            if record.get("gene") != gene:
                raise ValueError("Validity gene does not match the requested gene")
            if not record.get("classification"):
                raise ValueError("Validity record requires classification")
            if condition and record.get("condition") == condition:
                matched.append(record)
            elif not record.get("condition"):
                matched.append(record)
        states = {_normalized_validity(item["classification"]) for item in matched}
        positive = bool(states & POSITIVE_VALIDITY)
        limited = bool(states & LIMITED_VALIDITY)
        negative = bool(states & NEGATIVE_VALIDITY)
        if sum((positive, limited, negative)) > 1:
            state = "CONFLICTING"
        elif positive:
            state = "SUPPORTED"
        elif limited and states <= LIMITED_VALIDITY:
            state = "LIMITED"
        elif negative and states <= NEGATIVE_VALIDITY:
            state = "CURATED_NEGATIVE"
        else:
            state = "UNKNOWN"
        return state, matched

    @staticmethod
    def _constraint(record, gene):
        if record is None:
            return None
        _require_provenance(record, "constraint")
        if record.get("gene") != gene:
            raise ValueError("Constraint gene does not match the requested gene")
        result = dict(record)
        for key in ("mis_z", "p_li", "loeuf"):
            if result.get(key) is not None:
                result[key] = _number(result[key], f"constraint.{key}")
        return result

    @staticmethod
    def _spectrum(record, gene):
        if record is None:
            return None
        _require_provenance(record, "ClinVar spectrum")
        if record.get("gene") != gene:
            raise ValueError("ClinVar spectrum gene does not match the requested gene")
        if not isinstance(record.get("complete"), bool):
            raise ValueError("ClinVar spectrum must declare complete")
        result = dict(record)
        for key in (
            "pathogenic_missense_count",
            "pathogenic_truncating_count",
            "benign_missense_count",
        ):
            result[key] = _count(result.get(key), f"ClinVar spectrum.{key}")
        return result

    @staticmethod
    def _suggestion(status, reasons, requires):
        return {
            "status": status,
            "reasons": reasons,
            "requires_review_of": requires,
            "criterion_met": None,
        }

    def _suggestions(self, validity_state, constraint, spectrum):
        if validity_state == "CURATED_NEGATIVE":
            blocked = self._suggestion(
                "CURATED_NEGATIVE", ["Gene-disease validity is curated negative"],
                ["Resolve the gene-disease relationship before criterion assessment"],
            )
            return {code: dict(blocked) for code in ("PP2", "BP1", "PVS1")}
        if validity_state != "SUPPORTED":
            unresolved = self._suggestion(
                "INSUFFICIENT", ["Supported gene-disease validity is unavailable"],
                ["Confirm the exact gene-disease relationship"],
            )
            return {code: dict(unresolved) for code in ("PP2", "BP1", "PVS1")}

        mis_z = constraint.get("mis_z") if constraint else None
        p_li = constraint.get("p_li") if constraint else None
        loeuf = constraint.get("loeuf") if constraint else None

        if mis_z is None:
            pp2 = self._suggestion(
                "INSUFFICIENT", ["Missense constraint is unavailable"],
                ["missense mechanism", "benign missense spectrum"],
            )
        else:
            candidate = mis_z >= self.policy["pp2_min_mis_z"]
            pp2 = self._suggestion(
                "CANDIDATE" if candidate else "NOT_SUGGESTED",
                [f"misZ={mis_z:g} compared with {self.policy['pp2_min_mis_z']:g}"],
                ["missense mechanism", "low benign missense variation", "spectrum completeness"],
            )

        # Real ClinGen VCEPs determine BP1 from the empirical ratio of curated
        # pathogenic missense vs. truncating variants for the gene (e.g. the
        # LDLR VCEP: BP1 does not apply because most FH-pathogenic LDLR
        # variants are missense) - not from gnomAD population constraint.
        # pLI/mis_z were used here before (2026-09-21, found via the
        # 679-variant ground-truth run) and silently broke for reduced-
        # penetrance, adult-onset genes (APC/BRCA1/BRCA2/PALB2): carriers of
        # a truncating variant in those genes are usually asymptomatic
        # through reproductive age, so gnomAD - a mostly-healthy population
        # sample - does not select against them and pLI comes back near 0,
        # even though truncating variants are the established mechanism.
        # PALB2 confirms the fix: pLI~0 (would have failed the old gate),
        # but 9 pathogenic missense vs. 1092 pathogenic truncating
        # (0.8% missense) is unambiguous by the real method.
        if spectrum is None or not spectrum["complete"]:
            bp1 = self._suggestion(
                "INSUFFICIENT", ["Complete ClinVar spectrum input is required"],
                ["predominantly truncating mechanism", "absence of an established missense mechanism"],
            )
        else:
            pathogenic_missense = spectrum["pathogenic_missense_count"]
            pathogenic_truncating = spectrum["pathogenic_truncating_count"]
            total = pathogenic_missense + pathogenic_truncating
            if pathogenic_truncating < self.policy["bp1_min_truncating_count"]:
                bp1 = self._suggestion(
                    "INSUFFICIENT",
                    [f"ClinVar pathogenic truncating count={pathogenic_truncating} below the "
                     f"{self.policy['bp1_min_truncating_count']} needed to conclude predominance"],
                    ["predominantly truncating mechanism", "ClinVar assertion review"],
                )
            else:
                missense_fraction = pathogenic_missense / total
                candidate = missense_fraction <= self.policy["bp1_max_missense_fraction"]
                bp1 = self._suggestion(
                    "CANDIDATE" if candidate else "NOT_SUGGESTED",
                    [(f"ClinVar pathogenic missense fraction={missense_fraction:.1%} "
                      f"({pathogenic_missense}/{total}) compared with "
                      f"{self.policy['bp1_max_missense_fraction']:.1%}")],
                    ["predominantly truncating mechanism", "ClinVar assertion review", "spectrum completeness"],
                )

        if p_li is None and loeuf is None:
            pvs1 = self._suggestion(
                "INSUFFICIENT", ["LoF constraint is unavailable"],
                ["loss-of-function disease mechanism"],
            )
        else:
            candidate = (
                (p_li is not None and p_li >= self.policy["pvs1_min_p_li"])
                or (loeuf is not None and loeuf <= self.policy["pvs1_max_loeuf"])
            )
            reasons = []
            if p_li is not None:
                reasons.append(f"pLI={p_li:g} compared with {self.policy['pvs1_min_p_li']:g}")
            if loeuf is not None:
                reasons.append(f"LOEUF={loeuf:g} compared with {self.policy['pvs1_max_loeuf']:g}")
            pvs1 = self._suggestion(
                "CANDIDATE" if candidate else "NOT_SUGGESTED", reasons,
                ["loss-of-function disease mechanism", "transcript and variant-level PVS1 decision tree"],
            )
        return {"PP2": pp2, "BP1": bp1, "PVS1": pvs1}

    def build(
        self,
        *,
        variant_keys,
        gene,
        transcript,
        generated_at,
        validity,
        condition=None,
        constraint=None,
        clinvar_spectrum=None,
    ):
        """Return one audit draft per variant key; never return criterion-ready evidence."""
        if not isinstance(variant_keys, list) or not variant_keys:
            raise ValueError("variant_keys must be a nonempty list")
        if any(not isinstance(key, str) or len(key.split(":")) != 5 for key in variant_keys):
            raise ValueError("Every variant key must be assembly:chrom:pos:ref:alt")
        if not all(isinstance(value, str) and value for value in (gene, transcript, generated_at)):
            raise ValueError("gene, transcript and generated_at are required")
        validity_state, matched_validity = self._validity_state(validity, gene, condition)
        constraint = self._constraint(constraint, gene)
        spectrum = self._spectrum(clinvar_spectrum, gene)
        suggestions = self._suggestions(validity_state, constraint, spectrum)
        source_inputs = {
            "validity": matched_validity,
            "constraint": constraint,
            "clinvar_spectrum": spectrum,
        }
        fingerprint = hashlib.sha256(_canonical({
            "gene": gene,
            "condition": condition,
            "transcript": transcript,
            "policy": self.policy,
            "inputs": source_inputs,
        }).encode()).hexdigest()
        common = {
            "category": "gene_disease_draft",
            "quality_status": "DRAFT",
            "assessment_status": "DRAFT",
            "criterion_eligible": False,
            "evidence_id": f"urn:bh26:gene-disease-draft:{fingerprint}",
            "source": self.name,
            "source_version": self.policy["policy_version"],
            "retrieved_at": generated_at,
            "gene": gene,
            "condition": condition,
            "transcript": transcript,
            "validity_state": validity_state,
            "suggestions": suggestions,
            "source_inputs": source_inputs,
            "policy": self.policy,
            "review_instruction": (
                "Do not rename this record to gene_disease until a reviewer confirms the "
                "mechanism fields and records curator/reviewed_at provenance."
            ),
        }
        return [
            {
                **common,
                "variant_key": key,
                "evidence_id": (
                    f"urn:bh26:gene-disease-draft:{fingerprint}:"
                    f"{hashlib.sha256(key.encode()).hexdigest()[:16]}"
                ),
            }
            for key in variant_keys
        ]
