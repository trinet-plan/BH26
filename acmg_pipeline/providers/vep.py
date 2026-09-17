"""Ensembl VEP region annotation adapter (GRCh38 only)."""

from urllib.parse import quote


class VepProvider:
    name = "vep"

    def __init__(self, client, release):
        if not release:
            raise ValueError("Ensembl release must be recorded")
        self.client = client
        self.release = release

    def annotate(self, variant, context=None):
        # VCF representation is accepted by the POST region endpoint, including anchored indels.
        body = {"variants": [f"{variant.chrom} {variant.pos} . {variant.ref} {variant.alt} . . ."]}
        response = self.client.fetch("https://rest.ensembl.org/vep/human/region?hgvs=1&mane=1&protein=1",
                                     data=body, dataset_version=self.release)
        if not isinstance(response["body"], list) or len(response["body"]) != 1:
            raise ValueError("Unexpected VEP response")
        output = []
        for row in response["body"][0].get("transcript_consequences", []):
            transcript = row.get("transcript_id")
            if not transcript:
                continue
            amino_acids = row.get("amino_acids", "").split("/")
            output.append({"category": "annotation", "variant_key": variant.key,
                           "evidence_id": f"vep:{variant.key}:{transcript}",
                           "source": "Ensembl VEP", "source_version": self.release,
                           "retrieved_at": response["retrieved_at"], "quality_status": "PASS",
                           "transcript": transcript, "gene": row.get("gene_symbol"),
                           "consequences": row.get("consequence_terms", []),
                           "protein_id": row.get("protein_id"), "protein_start": row.get("protein_start"),
                           "protein_end": row.get("protein_end"), "hgvsc": row.get("hgvsc"),
                           "hgvsp": row.get("hgvsp"), "mane_select": row.get("mane_select"),
                           "ref_aa": amino_acids[0] if amino_acids else None,
                           "alt_aa": amino_acids[1] if len(amino_acids) == 2 else None})
        return output

    def map_hgvs(self, transcript, hgvsc):
        """Retrieve mapping response; reference validation is a separate mandatory step."""
        hgvs = quote(f"{transcript}:{hgvsc}", safe="")
        return self.client.fetch(f"https://rest.ensembl.org/vep/human/hgvs/{hgvs}?hgvs=1",
                                 dataset_version=self.release)
