"""ClinVar VCV retrieval for identity corroboration; classifications are never PP5/BP6."""

import hashlib
import re
from urllib.parse import urlencode
from xml.etree import ElementTree

VCV = re.compile(r"^VCV(\d{9})(?:\.(\d+))?$")


def local_name(tag):
    return tag.rsplit("}", 1)[-1]


class ClinVarProvider:
    name = "ClinVar"

    def __init__(self, client, release):
        if not release:
            raise ValueError("ClinVar release date/version must be recorded")
        self.client = client
        self.release = release

    def get_record(self, accession, expected_variant):
        match = VCV.fullmatch(accession)
        if not match:
            raise ValueError("Unsupported ClinVar VCV accession")
        variation_id = str(int(match.group(1)))
        query = urlencode({"db": "clinvar", "id": variation_id, "rettype": "vcv",
                           "is_variationid": "true"})
        response = self.client.fetch(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}",
            response_format="text", dataset_version=self.release,
        )
        try:
            root = ElementTree.fromstring(response["body"])
        except ElementTree.ParseError as exc:
            raise ValueError("Invalid ClinVar XML") from exc
        archives = [node for node in root.iter() if local_name(node.tag) == "VariationArchive"]
        if len(archives) != 1:
            raise ValueError("Unexpected ClinVar VCV response")
        archive = archives[0]
        returned_accession = archive.attrib.get("Accession")
        returned_version = archive.attrib.get("Version")
        if returned_accession != accession.split(".")[0]:
            raise ValueError("ClinVar accession mismatch")
        locations = []
        for node in archive.iter():
            attrs = node.attrib
            assembly = attrs.get("Assembly") or attrs.get("assembly")
            chrom = attrs.get("Chr") or attrs.get("chr")
            pos = attrs.get("start") or attrs.get("Start")
            ref = attrs.get("referenceAlleleVCF") or attrs.get("ReferenceAlleleVCF")
            alt = attrs.get("alternateAlleleVCF") or attrs.get("AlternateAlleleVCF")
            if assembly and chrom and pos and ref and alt:
                try:
                    locations.append((assembly, str(chrom).removeprefix("chr"), int(pos), ref, alt))
                except ValueError:
                    continue
        expected = (expected_variant.assembly, expected_variant.chrom, expected_variant.pos,
                    expected_variant.ref, expected_variant.alt)
        if expected not in locations:
            raise ValueError("ClinVar record does not contain expected GRCh38 variant")
        classifications = []
        review_statuses = []
        for node in archive.iter():
            name = local_name(node.tag)
            if name in {"Description", "Classification"} and node.text and node.text.strip():
                classifications.append(node.text.strip())
            review = node.attrib.get("ReviewStatus") or node.attrib.get("reviewStatus")
            if review:
                review_statuses.append(review)
            if name == "ReviewStatus" and node.text and node.text.strip():
                review_statuses.append(node.text.strip())
        digest = hashlib.sha256(response["body"].encode()).hexdigest()
        versioned = returned_accession + (f".{returned_version}" if returned_version else "")
        evidence = {
            "category": "clinvar_record", "variant_key": expected_variant.key,
            "evidence_id": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/",
            "source": self.name, "source_version": f"{self.release}:{versioned}",
            "retrieved_at": response["retrieved_at"], "quality_status": "PASS",
            "accession": versioned, "variation_id": variation_id,
            "requested_accession": accession,
            "requested_version": match.group(2),
            "version_advanced": bool(match.group(2) and returned_version != match.group(2)),
            "classifications": sorted(set(classifications)),
            "review_statuses": sorted(set(review_statuses)), "response_sha256": digest,
            "use_restriction": "IDENTITY_AND_COMPARATOR_CANDIDATE_ONLY_NOT_PP5_BP6",
        }
        identity = {
            "variant": expected_variant.to_dict(), "source": self.name,
            "source_version": f"{self.release}:{versioned}",
            "retrieved_at": response["retrieved_at"], "response_sha256": digest,
            "matched_identifiers": {"CLNVARIATIONID": accession},
        }
        return evidence, identity
