"""Ensembl-backed HGVS identity mapping with reference-validated VCF alleles."""

import hashlib
import json

from acmg.core.models import Variant
from acmg.core.reference import normalize
from acmg.providers.http import canonical_json


def reverse_complement(allele):
    if allele == "-":
        return allele
    return allele.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class EnsemblReference:
    def __init__(self, client, release):
        self.client = client
        self.release = str(release)

    def sequence(self, chrom, start, end):
        if start < 1 or end < start:
            raise ValueError("Reference interval out of bounds")
        region = f"{str(chrom).removeprefix('chr')}:{start}..{end}:1"
        response = self.client.fetch(
            f"https://rest.ensembl.org/sequence/region/human/{region}",
            response_format="text", dataset_version=self.release,
        )
        body = response["body"]
        # The sequence endpoint may serialize JSON despite a text/plain Accept header.
        try:
            decoded = json.loads(body)
        except (json.JSONDecodeError, TypeError):
            decoded = None
        value = decoded.get("seq") if isinstance(decoded, dict) else body
        if not isinstance(value, str):
            raise ValueError("Unexpected Ensembl sequence response")
        value = value.strip().upper()
        if len(value) != end - start + 1 or any(base not in "ACGT" for base in value):
            raise ValueError("Reference interval is unavailable or contains ambiguous bases")
        return value


class EnsemblIdentityProvider:
    name = "Ensembl VEP HGVS"

    def __init__(self, client, release):
        if not release:
            raise ValueError("Ensembl release must be recorded")
        self.client = client
        self.release = str(release)
        self.reference = EnsemblReference(client, release)

    @staticmethod
    def current_release(client):
        response = client.fetch("https://rest.ensembl.org/info/data", dataset_version="current")
        releases = response["body"].get("releases") if isinstance(response["body"], dict) else None
        if not isinstance(releases, list) or not releases or not all(isinstance(v, int) for v in releases):
            raise ValueError("Unexpected Ensembl release response")
        return str(max(releases))

    def map_record(self, record):
        transcript = record["identity"].get("TRANSCRIPT")
        hgvsc = record["identity"].get("HGVSC")
        if not transcript or not hgvsc:
            raise ValueError("Record has no versioned transcript HGVS identity")
        response = self.client.fetch(
            f"https://rest.ensembl.org/vep/human/hgvs/{transcript}%3A{hgvsc}?hgvs=1",
            dataset_version=self.release,
        )
        rows = response["body"]
        if not isinstance(rows, list) or len(rows) != 1 or not isinstance(rows[0], dict):
            raise ValueError("Unexpected Ensembl HGVS response")
        row = rows[0]
        expected = f"{transcript}:{hgvsc}"
        if row.get("input") != expected or row.get("assembly_name") != "GRCh38":
            raise ValueError("Ensembl response does not match requested HGVS/build")
        chrom, start, end = row.get("seq_region_name"), row.get("start"), row.get("end")
        alleles = row.get("allele_string", "").upper().split("/")
        strand = row.get("strand")
        if not isinstance(chrom, str) or not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("Ensembl response lacks genomic coordinates")
        if len(alleles) != 2 or any(not allele or any(base not in "ACGT-" for base in allele)
                                    for allele in alleles):
            raise ValueError("Unsupported Ensembl allele representation")
        if strand not in {-1, 1}:
            raise ValueError("Ensembl response lacks a valid strand")
        if strand == -1:
            alleles = [reverse_complement(allele) for allele in alleles]
        ref, alt = alleles
        if alt == "-":
            deleted = self.reference.sequence(chrom, start, end)
            if deleted != ref:
                raise ValueError("Ensembl deleted allele disagrees with reference")
            anchor = self.reference.sequence(chrom, start - 1, start - 1)
            variant = Variant("GRCh38", chrom, start - 1, anchor + ref, anchor)
        elif ref == "-":
            anchor = self.reference.sequence(chrom, start - 1, start - 1)
            variant = Variant("GRCh38", chrom, start - 1, anchor, anchor + alt)
        else:
            variant = Variant("GRCh38", chrom, start, ref, alt)
        variant = normalize(variant, self.reference)
        return {
            "variant": variant.to_dict(),
            "source": self.name,
            "source_version": self.release,
            "retrieved_at": response["retrieved_at"],
            "response_sha256": hashlib.sha256(canonical_json(rows).encode()).hexdigest(),
            "matched_identifiers": {"TRANSCRIPT": transcript, "HGVSC": hgvsc},
        }
