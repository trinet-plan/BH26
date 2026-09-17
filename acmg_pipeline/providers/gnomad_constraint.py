"""gnomAD GraphQL gene-constraint adapter (mis_z / pLI / LOEUF).

Separate from providers/gnomad.py's per-variant frequency lookups: this is a
per-gene, not per-allele, query - gnomAD reports constraint once per gene per
release, not per variant.
"""

from __future__ import annotations

QUERY = """
query GeneConstraint($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
  gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
    gene_id
    gnomad_constraint { mis_z pLI oe_lof_upper }
  }
}
""".strip()


class GnomadConstraintProvider:
    name = "gnomAD gene constraint"

    def __init__(self, client, *, release="4.1.1", reference_genome="GRCh38"):
        if not release:
            raise ValueError("gnomAD constraint release must be recorded")
        self.client = client
        self.release = release
        self.reference_genome = reference_genome

    def get_constraint(self, gene: str) -> dict | None:
        response = self.client.fetch(
            "https://gnomad.broadinstitute.org/api",
            data={"query": QUERY,
                  "variables": {"geneSymbol": gene, "referenceGenome": self.reference_genome}},
            dataset_version=self.release, allow_application_errors=True,
        )
        body = response["body"]
        if not isinstance(body, dict):
            raise ValueError("Unexpected gnomAD constraint response")
        data = body.get("data")
        gene_data = data.get("gene") if isinstance(data, dict) else None
        errors = [item for item in body.get("errors", []) if isinstance(item, dict)]
        if gene_data is None:
            messages = [item.get("message", "") for item in errors]
            if messages and all("not found" in message.lower() for message in messages):
                return None
            raise ValueError("gnomAD gene constraint query failed")
        constraint = gene_data.get("gnomad_constraint")
        if not isinstance(constraint, dict):
            return None
        for key in ("mis_z", "pLI", "oe_lof_upper"):
            value = constraint.get(key)
            if value is not None and not isinstance(value, (int, float)):
                raise ValueError(f"Invalid gnomAD constraint value for {key}")
        return {
            "category": "gene_constraint", "gene": gene,
            "source": self.name, "source_version": self.release,
            "retrieved_at": response["retrieved_at"],
            "evidence_id": f"https://gnomad.broadinstitute.org/gene/{gene}?dataset=gnomad_r4",
            "quality_status": "PASS",
            "mis_z": constraint.get("mis_z"),
            "p_li": constraint.get("pLI"),
            "loeuf": constraint.get("oe_lof_upper"),
        }
