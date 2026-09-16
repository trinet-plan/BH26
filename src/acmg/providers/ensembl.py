"""Ensembl-backed HGVS identity mapping with reference-validated VCF alleles."""

import hashlib
import json
from urllib.parse import quote

from acmg.core.models import Variant
from acmg.core.reference import normalize
from acmg.providers.http import canonical_json


VEP_OPTIONS = (
    "AlphaMissense=1;Conservation=1;REVEL=1;SpliceAI=2;hgvs=1;protein=1;"
    "refseq=1;transcript_version=1"
)


def _number(value):
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return value
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _amino_acids(consequence):
    value = consequence.get("amino_acids")
    if not isinstance(value, str) or not value:
        return None, None
    if "/" not in value:
        return value, value
    if value.count("/") != 1:
        return None, None
    ref, alt = value.split("/")
    return (ref or None), (alt or None)


def _protein_length_change(consequence, ref_aa, alt_aa):
    terms = set(consequence.get("consequence_terms", []))
    if terms & {"missense_variant", "synonymous_variant"}:
        return 0
    if terms & {"inframe_insertion", "inframe_deletion"} and ref_aa and alt_aa:
        return len(alt_aa.replace("-", "")) - len(ref_aa.replace("-", ""))
    # VEP's changed amino-acid token cannot establish the complete extension or
    # truncation length for stop-loss, frameshift, start-loss, or stop-gain calls.
    return None


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
        candidate, _ = self.map_record_with_annotation(record)
        return candidate

    def map_record_with_annotation(self, record):
        candidate, annotation, _ = self.map_record_with_evidence(record)
        return candidate, annotation

    def map_record_with_evidence(self, record):
        transcript = record["identity"].get("TRANSCRIPT")
        hgvsc = record["identity"].get("HGVSC")
        if not transcript or not hgvsc:
            raise ValueError("Record has no versioned transcript HGVS identity")
        encoded_hgvs = quote(f"{transcript}:{hgvsc}", safe="")
        response = self.client.fetch(
            f"https://rest.ensembl.org/vep/human/hgvs/{encoded_hgvs}?{VEP_OPTIONS}",
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
        response_sha256 = hashlib.sha256(canonical_json(rows).encode()).hexdigest()
        candidate = {
            "variant": variant.to_dict(),
            "source": self.name,
            "source_version": self.release,
            "retrieved_at": response["retrieved_at"],
            "response_sha256": response_sha256,
            "matched_identifiers": {"TRANSCRIPT": transcript, "HGVSC": hgvsc},
        }
        gene = record["identity"].get("GENE")
        transcript_rows = [item for item in row.get("transcript_consequences", [])
                           if isinstance(item, dict) and item.get("transcript_id") == transcript]
        if len(transcript_rows) > 1:
            raise ValueError("Ensembl returned multiple consequences for the requested transcript")
        consequence = transcript_rows[0] if transcript_rows else {}
        if consequence.get("gene_symbol") and gene and consequence["gene_symbol"] != gene:
            raise ValueError("Ensembl transcript consequence disagrees with requested gene")
        consequences = sorted(set(consequence.get("consequence_terms", [])))
        ref_aa, alt_aa = _amino_acids(consequence)
        annotation = {
            "category": "annotation", "variant_key": variant.key,
            "evidence_id": f"urn:sha256:{response_sha256}:annotation:{transcript}",
            "source": self.name, "source_version": self.release,
            "retrieved_at": response["retrieved_at"], "quality_status": "PASS",
            "transcript": transcript, "gene": gene, "hgvsc": expected,
            "consequences": consequences,
            "protein_id": consequence.get("protein_id"),
            "hgvsp": consequence.get("hgvsp"),
            "protein_start": consequence.get("protein_start"),
            "protein_end": consequence.get("protein_end"),
            "exon": consequence.get("exon"), "intron": consequence.get("intron"),
            "cds_start": consequence.get("cds_start"), "cds_end": consequence.get("cds_end"),
            "ref_aa": ref_aa, "alt_aa": alt_aa,
            "protein_length_change": _protein_length_change(consequence, ref_aa, alt_aa),
            "high_confidence_null_or_splice": bool(
                set(consequences) & {"stop_gained", "frameshift_variant", "splice_donor_variant",
                                     "splice_acceptor_variant", "start_lost"}),
        }
        predictions = self._prediction_evidence(
            variant.key, transcript, consequence, response["retrieved_at"], response_sha256
        )
        return candidate, annotation, predictions

    def _prediction_evidence(self, variant_key, transcript, consequence, retrieved_at,
                             response_sha256):
        common = {
            "category": "computational", "variant_key": variant_key,
            "source_version": self.release, "retrieved_at": retrieved_at,
            "quality_status": "PASS", "transcript": transcript,
        }
        records = []

        alpha = consequence.get("alphamissense")
        if isinstance(alpha, dict) and _number(alpha.get("am_pathogenicity")) is not None:
            records.append({
                **common,
                "evidence_id": f"urn:sha256:{response_sha256}:prediction:alphamissense:{transcript}",
                "source": "Ensembl VEP AlphaMissense",
                "predictor": "AlphaMissense",
                "predictor_version": "2023",
                "mechanism": "protein",
                "score": _number(alpha["am_pathogenicity"]),
                "classification": alpha.get("am_class"),
                "calibration_eligible": False,
                "calibration_note": "Model class thresholds are not ACMG PP3/BP4 evidence calibration",
                "version_provenance": "AlphaMissense Database Copyright 2023; Ensembl VEP release " + self.release,
            })

        revel = _number(consequence.get("revel", consequence.get("REVEL")))
        if revel is not None:
            records.append({
                **common,
                "evidence_id": f"urn:sha256:{response_sha256}:prediction:revel:{transcript}",
                "source": "Ensembl VEP REVEL",
                "predictor": "REVEL",
                # Ensembl REST does not expose the backing REVEL file version.
                "predictor_version": f"unreported-Ensembl-{self.release}",
                "mechanism": "protein", "score": revel,
                "calibration_eligible": False,
                "version_note": "Backing REVEL data version is not exposed by Ensembl REST",
            })

        splice = consequence.get("spliceai")
        if isinstance(splice, dict):
            deltas = {key: _number(splice.get(key)) for key in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")}
            valid = [value for value in deltas.values() if value is not None]
            if valid:
                records.append({
                    **common,
                    "evidence_id": f"urn:sha256:{response_sha256}:prediction:spliceai:{transcript}",
                    "source": "Ensembl VEP SpliceAI",
                    "predictor": "SpliceAI",
                    "predictor_version": f"unreported-Ensembl-{self.release}",
                    "mechanism": "splicing", "score": max(valid),
                    "delta_scores": deltas,
                    "dataset": "Ensembl/GENCODE v37 MANE raw scores (REST SpliceAI=2)",
                    "calibration_eligible": False,
                    "version_note": "SpliceAI model version is not exposed by Ensembl REST",
                })

        conservation = _number(consequence.get("conservation"))
        if conservation is not None:
            records.append({
                **common,
                "evidence_id": f"urn:sha256:{response_sha256}:prediction:conservation:{transcript}",
                "source": "Ensembl VEP Conservation",
                "predictor": "Ensembl Compara conservation",
                "predictor_version": self.release,
                "mechanism": "conservation", "score": conservation,
                "calibration_eligible": False,
                "version_note": "REST response does not identify the conservation method/track",
            })
        return records
