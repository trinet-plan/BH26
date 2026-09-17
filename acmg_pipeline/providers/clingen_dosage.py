"""ClinGen Dosage Sensitivity haploinsufficiency scores as an automated LoF-mechanism signal.

[Why this exists]
  PVS1 stops at G01 unless a `gene_disease` record states
  `lof_mechanism_established`. That is a curated judgment, and until one exists for a gene
  the whole decision tree is unreachable - in a run over the ground-truth set, all five
  variants whose expert-panel call is PVS1=met stopped at G01 without even reaching the
  variant-type node. This supplies the gate from a published, machine-readable source so
  the rest of the tree can run, following AutoPVS1's approach of substituting ClinGen's
  dosage curation for a per-gene mechanism review.

[Why it is not a curated assessment, and says so]
  gene_disease_draft.py's rule holds: an aggregate metric cannot by itself establish a
  disease mechanism. So every record here carries assessment_method="automated" with the
  method and policy_version that produced it - the shape
  acmg_pipeline.criteria.common.reviewed_or_automated() accepts for exactly this case - and
  never a curator/reviewed_at. A curated record for the same gene wins on its own merits:
  _resolve_mechanism() prefers a condition-specific record, and these are gene-level.

[Why a score of 30 is unresolved rather than false]
  The haploinsufficiency score answers "does losing one copy cause disease", which is not
  the question PVS1 asks. 30 means "gene associated with autosomal recessive phenotype",
  and a recessive disease can have loss of function as its mechanism perfectly well - the
  ClinGen Phenylketonuria VCEP applied PVS1 to PAH c.806delT, and PAH scores 30. Reading 30
  as "not a LoF mechanism" would strike out recessive PVS1 wholesale, so no record is
  emitted and PVS1 keeps reporting the mechanism as unavailable for a curator to resolve.

[Scores that do become false]
  0 (no evidence), 1 (little), 2 (emerging) and 40 (evidence *against* dosage sensitivity)
  are all statements about the dominant mechanism this gate is asking after. They produce
  lof_mechanism_established=False, which PVS1 reports as "loss of function is not an
  established disease mechanism" rather than as missing input.
"""

from __future__ import annotations

import hashlib
import re

CURATION_LIST_URL = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"

# The mapping is policy, not code: it is recorded on every record so a result says which
# reading of the score produced it. None means "emit nothing" - see the module docstring.
HAPLOINSUFFICIENCY_MEANING = {
    "3": True,    # sufficient evidence for dosage pathogenicity
    "2": False,   # emerging evidence
    "1": False,   # little evidence
    "0": False,   # no evidence
    "40": False,  # evidence suggests the gene is NOT dosage sensitive
    "30": None,   # autosomal recessive phenotype - not what this score decides
}
METHOD = "clingen_dosage_haploinsufficiency"


# The file dates itself in a leading comment of its own, e.g. "#16 Sep,2026".
_RELEASE_LINE = re.compile(r"^#\s*(\d{1,2}\s+[A-Za-z]{3,9},?\s*\d{4})\s*$")


def _release(lines):
    """The file's own release date; without it there is no version to pin a policy to, and
    an unversioned policy is not usable evidence here."""
    for line in lines[:20]:
        if not line.startswith("#"):
            break
        match = _RELEASE_LINE.match(line)
        if match:
            return match.group(1).strip()
    return None


def parse(text):
    """Return {gene symbol: {"score": str, "gene_id": str}} plus the file's release."""
    lines = text.splitlines()
    header = next((line for line in lines if line.startswith("#Gene Symbol")), None)
    if header is None:
        raise ValueError("ClinGen dosage list is missing its #Gene Symbol header")
    columns = header.lstrip("#").split("\t")
    try:
        symbol_at = columns.index("Gene Symbol")
        score_at = columns.index("Haploinsufficiency Score")
    except ValueError as error:
        raise ValueError("ClinGen dosage list is missing a required column") from error
    def optional(name):
        return columns.index(name) if name in columns else None

    gene_id_at = optional("Gene ID")
    disease_at = optional("Haploinsufficiency Disease ID")
    evaluated_at = optional("Date Last Evaluated")
    description_at = optional("Haploinsufficiency Description")

    genes = {}
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) <= score_at:
            continue
        symbol = fields[symbol_at].strip()
        if symbol:
            def at(index):
                if index is None or len(fields) <= index:
                    return None
                return fields[index].strip() or None

            genes[symbol] = {
                "score": fields[score_at].strip(),
                "gene_id": at(gene_id_at),
                "disease_id": at(disease_at),
                "last_evaluated": at(evaluated_at),
                "description": at(description_at),
            }
    return genes, _release(lines)


class ClinGenDosageProvider:
    """Builds gene-level `gene_disease` evidence from the ClinGen dosage curation list."""

    name = "ClinGen Dosage Sensitivity Map"

    def __init__(self, client, *, url=CURATION_LIST_URL):
        self.client = client
        self.url = url

    def get_mechanism(self, variant, gene):
        """One automated, gene-level record for `gene`, or none.

        None covers three different situations, and PVS1 reports all of them as an
        unavailable mechanism: the gene is not in the list at all, its score is 30, or the
        file did not state a release to pin the policy to.
        """
        if not gene:
            return []
        response = self.client.fetch(self.url, response_format="text")
        text = response["body"]
        genes, release = parse(text)
        if not release:
            return []
        entry = genes.get(gene)
        if entry is None:
            return []
        established = HAPLOINSUFFICIENCY_MEANING.get(entry["score"])
        if established is None:
            return []
        digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        return [{
            "category": "gene_disease",
            "variant_key": variant.key,
            # The variant key belongs in the id: these records are per variant, and the
            # CLI deduplicates evidence by evidence_id, so a gene-only id collapsed every
            # variant in one gene down to whichever was seen first.
            "evidence_id": f"urn:sha256:{digest}:clingen-dosage:{gene}:{variant.key}",
            "source": self.name,
            "source_version": release,
            "retrieved_at": response["retrieved_at"],
            "quality_status": "PASS",
            "gene": gene,
            "lof_mechanism_established": established,
            # Not a curated judgment, and every consumer can see that it is not.
            "assessment_method": "automated",
            "method": METHOD,
            "policy_version": release,
            "haploinsufficiency_score": entry["score"],
            "haploinsufficiency_description": entry["description"],
            # Recorded as provenance, deliberately NOT as `condition`: a condition here
            # would make _resolve_mechanism() treat this as a condition-specific
            # assessment, which an automated gene-level signal is not.
            "haploinsufficiency_disease_id": entry["disease_id"],
            "gene_last_evaluated": entry["last_evaluated"],
            "clingen_gene_id": entry["gene_id"],
            "policy_note": (
                "Derived from the ClinGen haploinsufficiency score, not from a review of "
                "this gene-disease pair. A score of 30 (autosomal recessive) yields no "
                "record, because it does not answer whether loss of function is the "
                "mechanism."
            ),
            "response_sha256": digest,
        }]
