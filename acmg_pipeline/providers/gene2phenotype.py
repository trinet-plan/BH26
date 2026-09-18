"""Gene2Phenotype curation read as PVS1's disease-scoped loss-of-function mechanism.

[Why this exists]
  PVS1's disease gate needs a mechanism curated for the disease being assessed, named by an
  identifier it can compare. ClinGen dosage answers that for the genes it has scored, but it
  scores one haploinsufficiency judgment per gene, so a gene that loses function in one
  disease and gains it in another is exactly the case it cannot separate. G2P curates per
  gene-disease pair, states the molecular mechanism directly, and publishes the disease's
  MONDO accession alongside it - the three things the gate asks for.

[What is read, and what is deliberately not]
  `molecular_mechanism` decides the mechanism. `genotype` is G2P's allelic requirement and
  becomes the record's inheritance mode. `disease.ontology_terms` supplies the MONDO
  accession, and a record that names none is not emitted: a disease name is not an identifier,
  and matching on names is the route ACMG's match levels exist to rule out.

  `confidence` is a gene-disease validity statement, not a mechanism one, so it never sets
  `lof_mechanism_established`. It decides only whether a record is emitted at all: a limited,
  disputed or refuted curation is not a foundation to apply PVS1 on, and reading one as
  "loss of function is not the mechanism" would turn validity into mechanism, which is the
  conflation gene_disease_draft.py exists to prevent.

[Why "undetermined" emits nothing]
  G2P says "undetermined" where it has not settled the mechanism. That is not a statement
  that loss of function is ruled out, so it produces no record and PVS1 keeps reporting the
  mechanism as unavailable - the same reading the dosage provider gives a score of 30.
"""

from __future__ import annotations

import hashlib

API = "https://www.ebi.ac.uk/gene2phenotype/api"
METHOD = "gene2phenotype_molecular_mechanism"

# G2P's molecular mechanism vocabulary. None means "emit nothing" - see the module docstring.
MECHANISM_MEANING = {
    "loss of function": True,
    "gain of function": False,
    "dominant negative": False,
    "undetermined non-loss-of-function": False,
    "undetermined": None,
}
# Validity strong enough to act on. Everything else emits nothing rather than a negative.
USABLE_CONFIDENCE = {"definitive", "strong", "moderate"}

# G2P records an allelic requirement where ACMG records an inheritance mode. The autosomal
# terms carry over exactly. The X-linked terms name which zygosity is affected, which is more
# than "X-linked" and less than a claim about dominance in the ACMG sense, so they map to the
# unqualified mode and the comparison treats that as compatible with either.
GENOTYPE_MODES = {
    "monoallelic_autosomal": "autosomal_dominant",
    "biallelic_autosomal": "autosomal_recessive",
    "monoallelic_X_hemizygous": "x_linked",
    "monoallelic_X_heterozygous": "x_linked",
    "biallelic_X": "x_linked",
    "mitochondrial": "mitochondrial",
}


def mondo_accession(disease):
    """The one MONDO term this disease is recorded under, or None.

    More than one is not a tie to break here: it would pick which disease the curation was
    about, which is the curator's call.
    """
    if not isinstance(disease, dict):
        return None
    terms = {term.get("accession") for term in disease.get("ontology_terms") or []
             if isinstance(term, dict) and str(term.get("accession", "")).startswith("MONDO:")}
    return next(iter(terms)) if len(terms) == 1 else None


class Gene2PhenotypeProvider:
    """Builds automated, disease-scoped `gene_disease` evidence from G2P curation."""

    name = "Gene2Phenotype (EBI)"

    def __init__(self, client, *, api=API):
        self.client = client
        self.api = api

    def get_mechanism(self, variant, gene):
        """Every G2P record for `gene` that settles a mechanism for a MONDO-named disease."""
        if not gene:
            return []
        summary = self.client.fetch(f"{self.api}/gene/{gene}/summary/")
        records = []
        for entry in summary["body"].get("records_summary") or []:
            if entry.get("confidence") not in USABLE_CONFIDENCE:
                continue
            established = MECHANISM_MEANING.get(entry.get("molecular_mechanism"))
            if established is None:
                continue
            stable_id = entry.get("stable_id")
            if not stable_id:
                continue
            detail = self.client.fetch(f"{self.api}/lgd/{stable_id}/")
            body = detail["body"]
            condition = mondo_accession(body.get("disease"))
            release = body.get("last_updated")
            if not condition or not release:
                continue
            mechanism = body.get("molecular_mechanism") or {}
            records.append({
                "category": "gene_disease",
                "variant_key": variant.key,
                "evidence_id": f"g2p:{stable_id}:{release}:{variant.key}",
                "source": self.name,
                "source_version": release,
                "retrieved_at": detail["retrieved_at"],
                "quality_status": "PASS",
                "gene": gene,
                "condition": condition,
                # The name G2P curated the disease under, so a candidate list can be read
                # without resolving every identifier by hand.
                "condition_label": (body.get("disease") or {}).get("name"),
                "lof_mechanism_established": established,
                "inheritance": GENOTYPE_MODES.get(entry.get("genotype")),
                # Curated by a panel, but read here without a human confirming that this
                # record is the right one for this case, so it says so.
                "assessment_method": "automated",
                "method": METHOD,
                "policy_version": release,
                "g2p_stable_id": stable_id,
                "g2p_confidence": entry.get("confidence"),
                "g2p_genotype": entry.get("genotype"),
                "g2p_molecular_mechanism": mechanism.get("mechanism"),
                # "inferred" or "evidence": how G2P itself arrived at the mechanism.
                "g2p_mechanism_support": mechanism.get("mechanism_support"),
                "g2p_panels": sorted(entry.get("panels") or []),
                "policy_note": (
                    "Derived from the G2P molecular mechanism for this gene-disease record. "
                    "Confidence decides whether the record is used, never whether loss of "
                    "function is the mechanism."
                ),
                "response_sha256": hashlib.sha256(
                    str(detail.get("body_sha256") or stable_id).encode()).hexdigest(),
            })
        return records
