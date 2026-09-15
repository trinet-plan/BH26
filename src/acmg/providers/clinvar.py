"""ClinVar VCV retrieval for identity corroboration; classifications are never PP5/BP6."""

import hashlib
import re
from urllib.parse import urlencode
from xml.etree import ElementTree

from acmg.core.models import Variant
from acmg.providers.http import canonical_json

VCV = re.compile(r"^VCV(\d{9})(?:\.(\d+))?$")


def local_name(tag):
    return tag.rsplit("}", 1)[-1]


def direct_child(node, name):
    return next((child for child in node if local_name(child.tag) == name), None)


def child_path(node, *names):
    for name in names:
        node = direct_child(node, name)
        if node is None:
            return None
    return node


def node_text(node):
    return node.text.strip() if node is not None and node.text and node.text.strip() else None


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


class ClinVarComparatorProvider:
    """Discover PS1 protein matches; final disease relevance remains human-reviewed."""

    name = "ClinVar protein comparator search"
    spliceai_no_impact_threshold = 0.1

    def __init__(self, client, release, ensembl):
        if not release:
            raise ValueError("ClinVar release date/version must be recorded")
        self.client = client
        self.release = release
        self.ensembl = ensembl

    def search_ps1(self, annotation, query_variant, query_splice_score=None):
        required = ("transcript", "protein_id", "protein_start", "ref_aa", "alt_aa", "hgvsp")
        if not all(annotation.get(key) for key in required):
            raise ValueError("PS1 search requires complete transcript/protein annotation")
        term = f'"{annotation["hgvsp"]}"[All Fields]'
        query = urlencode({"db": "clinvar", "term": term, "retmode": "json", "retmax": 100})
        response = self.client.fetch(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{query}",
            dataset_version=self.release,
        )
        body = response["body"]
        result = body.get("esearchresult") if isinstance(body, dict) else None
        ids = result.get("idlist") if isinstance(result, dict) else None
        try:
            count = int(result["count"])
            retmax = int(result["retmax"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Unexpected ClinVar comparator search response") from exc
        if not isinstance(ids, list) or any(not str(value).isdigit() for value in ids):
            raise ValueError("Unexpected ClinVar comparator ID list")
        complete = count <= retmax and len(ids) == count
        search_digest = hashlib.sha256(canonical_json(body).encode()).hexdigest()
        search_evidence = {
            "category": "comparator_search", "variant_key": query_variant.key,
            "evidence_id": (f"urn:sha256:{search_digest}:clinvar-ps1-search:"
                            f'{annotation["transcript"]}'),
            "source": self.name, "source_version": self.release,
            "retrieved_at": response["retrieved_at"], "quality_status": "PASS",
            "transcript": annotation["transcript"], "protein_id": annotation["protein_id"],
            "protein_start": annotation["protein_start"], "ref_aa": annotation["ref_aa"],
            "alt_aa": annotation["alt_aa"], "hgvsp": annotation["hgvsp"],
            "query": term, "returned_count": count, "candidate_ids": [str(value) for value in ids],
            "complete": complete,
        }
        comparators = []
        for variation_id in ids:
            candidate = self._candidate(str(variation_id), annotation, query_variant,
                                        query_splice_score)
            if candidate is not None:
                comparators.append(candidate)
        return search_evidence, comparators

    def _candidate(self, variation_id, annotation, query_variant, query_splice_score):
        query = urlencode({"db": "clinvar", "id": variation_id, "rettype": "vcv",
                           "is_variationid": "true"})
        response = self.client.fetch(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}",
            response_format="text", dataset_version=self.release,
        )
        try:
            root = ElementTree.fromstring(response["body"])
        except ElementTree.ParseError as exc:
            raise ValueError("Invalid ClinVar comparator XML") from exc
        archives = [node for node in root.iter() if local_name(node.tag) == "VariationArchive"]
        if len(archives) != 1:
            raise ValueError("Unexpected ClinVar comparator response")
        archive = archives[0]
        expressions = [node_text(node) for node in archive.iter()
                       if local_name(node.tag) == "Expression"]
        if annotation["hgvsp"] not in expressions:
            return None
        prefix = annotation["transcript"] + ":"
        hgvsc_values = sorted(set(value for value in expressions if value and value.startswith(prefix)))
        if not hgvsc_values:
            return None
        mapped_values = {}
        errors = []
        for expression in hgvsc_values:
            hgvsc = expression.split(":", 1)[1]
            record = {"identity": {"TRANSCRIPT": annotation["transcript"], "HGVSC": hgvsc,
                                   "GENE": annotation.get("gene")}}
            try:
                mapped, candidate_annotation, predictions = \
                    self.ensembl.map_record_with_evidence(record)
            except ValueError as exc:
                errors.append(f"{expression}: {exc}")
                continue
            mapped_values.setdefault(
                Variant(**mapped["variant"]).key,
                (expression, mapped, candidate_annotation, predictions),
            )
        if not mapped_values:
            raise ValueError("ClinVar comparator HGVS could not be mapped: " + "; ".join(errors))
        if len(mapped_values) != 1:
            raise ValueError("ClinVar comparator HGVS aliases map to conflicting variants")
        selected_hgvsc, mapped, candidate_annotation, predictions = next(iter(mapped_values.values()))
        comparator_variant = Variant(**mapped["variant"])
        if comparator_variant == query_variant:
            return None
        if candidate_annotation.get("protein_id") != annotation["protein_id"]:
            return None
        changes = self._changed_residues(candidate_annotation)
        expected_change = (annotation["protein_start"], annotation["ref_aa"], annotation["alt_aa"])
        if changes != [expected_change] or candidate_annotation.get("hgvsp") != annotation["hgvsp"]:
            return None
        classified = direct_child(archive, "ClassifiedRecord")
        germline = child_path(classified, "Classifications", "GermlineClassification") \
            if classified is not None else None
        classification = node_text(direct_child(germline, "Description")) \
            if germline is not None else None
        review_status = node_text(direct_child(germline, "ReviewStatus")) \
            if germline is not None else None
        review_eligible = bool(review_status and (
            "criteria provided" in review_status.lower()
            or "expert panel" in review_status.lower()
            or "practice guideline" in review_status.lower()
        ) and "conflicting" not in review_status.lower())
        candidate_splice = next((item.get("score") for item in predictions
                                 if item.get("predictor") == "SpliceAI"), None)
        splice_checked = isinstance(query_splice_score, (int, float)) \
            and isinstance(candidate_splice, (int, float))
        splice_conflict = None if not splice_checked else (
            query_splice_score > self.spliceai_no_impact_threshold
            or candidate_splice > self.spliceai_no_impact_threshold
        )
        conditions = self._conditions(archive)
        accession = archive.attrib.get("Accession")
        version = archive.attrib.get("Version")
        versioned = accession + (f".{version}" if accession and version else "") \
            if accession else f"VariationID:{variation_id}"
        digest = hashlib.sha256(response["body"].encode()).hexdigest()
        return {
            "category": "comparator", "variant_key": query_variant.key,
            "evidence_id": (
                f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/#ps1-"
                f'{hashlib.sha256((query_variant.key + ":" + annotation["transcript"]).encode()).hexdigest()[:16]}'
            ),
            "source": "ClinVar", "source_version": f"{self.release}:{versioned}",
            "retrieved_at": response["retrieved_at"], "quality_status": "PASS",
            "response_sha256": digest, "accession": versioned, "variation_id": variation_id,
            "transcript": annotation["transcript"], "protein_id": annotation["protein_id"],
            "protein_start": annotation["protein_start"], "ref_aa": annotation["ref_aa"],
            "alt_aa": annotation["alt_aa"], "protein_change": annotation["hgvsp"],
            "comparator_protein_interval": {
                "start": candidate_annotation.get("protein_start"),
                "end": candidate_annotation.get("protein_end"),
                "ref_aa": candidate_annotation.get("ref_aa"),
                "alt_aa": candidate_annotation.get("alt_aa"),
            },
            "comparator_hgvsc": selected_hgvsc, "comparator_hgvsc_aliases": hgvsc_values,
            "comparator_variant": comparator_variant.to_dict(),
            "classification": classification, "review_status": review_status,
            "review_status_eligible": review_eligible, "conditions": conditions,
            "exact_protein_match": True, "different_nucleotide_variant": True,
            "query_spliceai": query_splice_score, "comparator_spliceai": candidate_splice,
            "splice_effect_checked": splice_checked, "splice_conflict": splice_conflict,
            "splice_policy": "ClinGen SVI Walker 2023; no-impact delta score <=0.1",
            "assessment_scope": "protein_level",
        }

    @staticmethod
    def _changed_residues(annotation):
        start = annotation.get("protein_start")
        ref = annotation.get("ref_aa")
        alt = annotation.get("alt_aa")
        if type(start) is not int or not isinstance(ref, str) or not isinstance(alt, str):
            return []
        if "-" in {ref, alt} or len(ref) != len(alt):
            return []
        return [(start + offset, ref_aa, alt_aa)
                for offset, (ref_aa, alt_aa) in enumerate(zip(ref, alt))
                if ref_aa != alt_aa]

    @staticmethod
    def _conditions(archive):
        values = set()
        classified = direct_child(archive, "ClassifiedRecord")
        conditions = child_path(classified, "Classifications", "GermlineClassification",
                                "ConditionList") if classified is not None else None
        if conditions is None:
            return []
        for node in conditions.iter():
            name = local_name(node.tag)
            if name == "XRef" and node.attrib.get("DB") and node.attrib.get("ID"):
                database, identifier = node.attrib["DB"], node.attrib["ID"]
                if identifier.startswith((database + ":", "MONDO:", "HP:")):
                    values.add(identifier)
                else:
                    values.add(f"{database}:{identifier}")
            elif name == "ElementValue" and node.attrib.get("Type") == "Preferred":
                value = node_text(node)
                if value:
                    values.add(value)
        return sorted(values)
