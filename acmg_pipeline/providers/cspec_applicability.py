"""ClinGen VCEP criteria specifications read as PVS1's disease-scoped mechanism gate.

[Why this exists]
  PVS1 stops at G01 unless a `gene_disease` record states `lof_mechanism_established`. The
  two automated producers beside this one both go silent for genes they have not covered:
  ClinGen dosage emits nothing for a haploinsufficiency score of 30 (the autosomal-recessive
  value - see clingen_dosage.py for why that is deliberate), and Gene2Phenotype answers only
  for gene-disease pairs it has curated. DYSF c.3498_3499delinsAA - Pathogenic with PVS1 at
  the ClinGen LGMD VCEP - reaches neither, and comes out of the pipeline as a VUS.

  A VCEP that wrote a criteria specification has already ruled on whether PVS1 applies to
  the gene, and ClinGen publishes those specifications machine-readably. This reads the
  collected snapshot of them.

[The mapping, and that it is a decision rather than a reading]
  A VCEP marking PVS1 "Applicable" is a statement about its own specification. Treating it
  as `lof_mechanism_established` is an interpretation, because the two are not the same
  sentence: PVS1 is the criterion for null variants in a gene where loss of function is a
  known mechanism, so a VCEP declaring it applicable has decided that premise holds - but it
  said so about the criterion, not about the mechanism.

  The project made that mapping deliberately (curator decision, 2026-09-18) on the evidence
  that it reproduces judgments already made by hand: on the three genes where
  config/gene-disease-review-decisions.json and CSpec both have a disease-scoped answer, the
  two agree 3/3 (MYBPC3 true/Applicable, MYH7 false/Not applicable, TNNI3 false/Not
  applicable). Every record says so in `policy_note` and carries POLICY_VERSION, so a result
  can be traced back to the decision that produced it rather than looking like a reading.

  "Not applicable" becomes False, not a missing value. A VCEP that struck PVS1 out of its
  own specification has ruled on the premise, which is the MYH7 case - and reporting that as
  "unknown" would send a curator looking for a mechanism the panel already declined.

[What is deliberately NOT read]
  Only PVS1. The snapshot carries all 28 criteria, but `pp2_applicable` and `bp1_applicable`
  in a reviewed assessment mean the criterion is in scope and leave met/not-met to the
  mechanism fields beside them, while CSpec's "Not applicable" means the VCEP struck the
  criterion out. They are different axes, so carrying those two over would change what the
  fields say. PP2/BP1 keep reporting their inputs as missing.

[Why these records do not outrank a reviewed one]
  They are `assessment_method: "automated"`. The snapshot is machine-collected and no
  curator confirmed that a given gene's record is the right one for a given case - only that
  reading VCEP applicability this way is sound. A reviewed assessment for the same gene and
  disease still wins; PVS1's own `_resolve_mechanism()` decides on the reviewed record and
  keeps this one visible as evidence beside it.

[Scope]
  Records are emitted per gene-disease pair, because the snapshot carries the MONDO term the
  specification was written against. A gene entry naming no MONDO term produces nothing: a
  mechanism has to be about a disease, and matching on a gene alone is what the disease gate
  exists to prevent.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

METHOD = "clingen_cspec_pvs1_applicability"
POLICY_VERSION = "BH26-cspec-pvs1-applicability-v1"
POLICY_NOTE = (
    "Derived from the ClinGen VCEP criteria specification's PVS1 applicability for this "
    "gene-disease pair. Applicable is read as loss of function being an established "
    "mechanism and Not applicable as it not being one - a project decision (2026-09-18), "
    "not a restatement of what the VCEP wrote."
)

# What the snapshot has to carry before anything here will read it.
_REQUIRED_KEYS = ("registry_status", "source", "source_url", "retrieved_at", "entries")

APPLICABILITY_MEANING = {
    "Applicable": True,
    "Not applicable": False,
}


class CSpecApplicabilityProvider:
    """Builds automated, disease-scoped `gene_disease` evidence from collected VCEP specs."""

    name = "ClinGen Criteria Specification Registry (CSpec)"

    def __init__(self, path):
        self.path = Path(path)
        self._document = None

    def _load(self) -> dict:
        """The snapshot, read once per provider instance.

        A malformed or truncated file raises rather than resolving to "no records": silently
        answering nothing here is indistinguishable from a gene the registry does not cover,
        and PVS1 would report a missing mechanism for every variant in the run.
        """
        if self._document is None:
            document = json.loads(self.path.read_text(encoding="utf-8-sig"))
            missing = [key for key in _REQUIRED_KEYS if not document.get(key)]
            if missing:
                raise ValueError(
                    f"CSpec applicability snapshot is missing {', '.join(missing)}")
            if not isinstance(document["entries"], dict):
                raise ValueError("CSpec applicability snapshot's entries must be an object")
            self._document = document
        return self._document

    def get_mechanism(self, variant, gene):
        """Every gene-disease pair this gene's specification settles a PVS1 premise for."""
        if not gene:
            return []
        document = self._load()
        entry = document["entries"].get(gene)
        if not entry:
            return []
        stated = (entry.get("criteria") or {}).get("PVS1")
        established = APPLICABILITY_MEANING.get(stated)
        if established is None:
            return []

        source_version = document.get("source_version") or document["retrieved_at"]
        inheritance = _inheritance(entry.get("inheritance") or [])
        records = []
        for condition in entry.get("mondo") or []:
            records.append({
                "category": "gene_disease",
                "variant_key": variant.key,
                "evidence_id": (
                    f"cspec:{gene}:{condition}:{source_version}:"
                    f"{hashlib.sha256(variant.key.encode()).hexdigest()[:16]}"
                ),
                "source": self.name,
                "source_version": source_version,
                "retrieved_at": document["retrieved_at"],
                "quality_status": "PASS",
                "gene": gene,
                "condition": condition,
                "lof_mechanism_established": established,
                "inheritance": inheritance,
                # Curated by a panel, but read here without a human confirming that this
                # record is the right one for this case, so it says so.
                "assessment_method": "automated",
                "method": METHOD,
                "policy_version": POLICY_VERSION,
                "cspec_pvs1_applicability": stated,
                "cspec_svis": sorted(entry.get("svis") or []),
                "registry_status": document["registry_status"],
                "policy_note": POLICY_NOTE,
                "response_sha256": hashlib.sha256(
                    f"{gene}|{condition}|{stated}|{source_version}".encode()).hexdigest(),
            })
        return records


def _inheritance(modes) -> str | None:
    """The one inheritance mode the specification names, in this project's vocabulary.

    More than one is left unstated rather than picked between: which of them the case is
    under is the disease gate's question, and a guess here would answer it silently.
    """
    normalized = {mode for mode in (_MODES.get(str(item).strip()) for item in modes) if mode}
    return normalized.pop() if len(normalized) == 1 else None


# The specifications record inheritance as the HPO mode-of-inheritance label.
_MODES = {
    "Autosomal recessive inheritance": "autosomal_recessive",
    "Autosomal dominant inheritance": "autosomal_dominant",
    "X-linked inheritance": "x_linked",
    "X-linked recessive inheritance": "x_linked",
    "X-linked dominant inheritance": "x_linked",
    "Mitochondrial inheritance": "mitochondrial",
    "Semidominant inheritance": "semidominant",
}
