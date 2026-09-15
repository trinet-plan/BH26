"""ClinVar VCV retrieval for identity corroboration; classifications are never PP5/BP6."""

import hashlib
import re
from urllib.parse import urlencode
from xml.etree import ElementTree

from acmg.core.models import Variant
from acmg.providers.http import canonical_json

VCV = re.compile(r"^VCV(\d{9})(?:\.(\d+))?$")
PROTEIN_CHANGE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
PROTEIN_HGVS = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


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


GENE_MISSENSE_RETMAX = 5000
SUMMARY_BATCH = 200


def gene_missense_search(client, release, gene, retmax=GENE_MISSENSE_RETMAX):
    """One gene-wide missense search, shared by the PM1 density and PM5 residue lookups."""
    term = f'{gene}[gene] AND "missense variant"[molecular consequence]'
    query = urlencode({"db": "clinvar", "term": term, "retmode": "json", "retmax": retmax})
    response = client.fetch(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{query}",
        dataset_version=release,
    )
    body = response["body"]
    found = body.get("esearchresult") if isinstance(body, dict) else None
    ids = found.get("idlist") if isinstance(found, dict) else None
    try:
        count = int(found["count"])
        limit = int(found["retmax"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Unexpected ClinVar gene search response") from exc
    if not isinstance(ids, list) or any(not str(value).isdigit() for value in ids):
        raise ValueError("Unexpected ClinVar gene search ID list")
    return {
        "term": term, "ids": [str(value) for value in ids], "count": count,
        "complete": count <= limit and len(ids) == count,
        "retrieved_at": response["retrieved_at"],
        "digest": hashlib.sha256(canonical_json(body).encode()).hexdigest(),
    }


def summary_documents(client, release, ids, batch=SUMMARY_BATCH):
    documents = {}
    retrieved_at = None
    for offset in range(0, len(ids), batch):
        chunk = ids[offset:offset + batch]
        query = urlencode({"db": "clinvar", "id": ",".join(chunk), "retmode": "json"})
        response = client.fetch(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{query}",
            dataset_version=release,
        )
        retrieved_at = response["retrieved_at"] if retrieved_at is None \
            else max(retrieved_at, response["retrieved_at"])
        body = response["body"]
        result = body.get("result") if isinstance(body, dict) else None
        if not isinstance(result, dict):
            raise ValueError("Unexpected ClinVar summary response")
        for uid in chunk:
            document = result.get(uid)
            if not isinstance(document, dict):
                raise ValueError(f"Missing ClinVar summary for {uid}")
            documents[uid] = document
    return documents, retrieved_at


def protein_change_positions(protein_change):
    """Only plain single-residue substitutions such as R248W are positioned."""
    if not isinstance(protein_change, str):
        return []
    positions = []
    for token in protein_change.split(","):
        match = PROTEIN_CHANGE.fullmatch(token.strip())
        if match:
            positions.append(int(match.group(2)))
    return sorted(set(positions))


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
            "complete": complete, "search_scope": "exact_protein_change",
        }
        comparators = []
        for variation_id in ids:
            candidate = self._candidate(str(variation_id), annotation, query_variant,
                                        query_splice_score)
            if candidate is not None:
                comparators.append(candidate)
        return search_evidence, comparators

    def search_pm5(self, annotation, query_variant, query_splice_score=None):
        """Residue-scoped search: PS1's exact-change query cannot show what else is reported."""
        required = ("transcript", "gene", "protein_id", "protein_start", "ref_aa", "alt_aa")
        if not all(annotation.get(key) for key in required):
            raise ValueError("PM5 search requires complete transcript/protein annotation")
        found = gene_missense_search(self.client, self.release, annotation["gene"])
        search = {
            "category": "comparator_search", "variant_key": query_variant.key,
            "evidence_id": (f'urn:sha256:{found["digest"]}:clinvar-pm5-search:'
                            f'{annotation["protein_id"]}:{annotation["protein_start"]}'),
            "source": self.name, "source_version": self.release,
            "retrieved_at": found["retrieved_at"], "quality_status": "PASS",
            "transcript": annotation["transcript"], "protein_id": annotation["protein_id"],
            "protein_start": annotation["protein_start"], "ref_aa": annotation["ref_aa"],
            "gene": annotation["gene"], "query": found["term"],
            "returned_count": found["count"], "search_scope": "residue",
            "complete": found["complete"],
            "protein_change_source": "ClinVar esummary gene-level protein_change",
            "position_confirmation": f'ClinVar VCV protein expression on {annotation["protein_id"]}',
        }
        if not found["complete"]:
            search["quality_status"] = "INCOMPLETE_SEARCH"
            return search, []
        documents, retrieved_at = summary_documents(self.client, self.release, found["ids"])
        residue = annotation["protein_start"]
        candidate_ids = [uid for uid, document in documents.items()
                         if residue in protein_change_positions(document.get("protein_change"))]
        search["retrieved_at"] = max(search["retrieved_at"], retrieved_at or search["retrieved_at"])
        search["candidate_ids"] = candidate_ids
        comparators = []
        for uid in candidate_ids:
            # The summary positions are transcript-agnostic, so each candidate is confirmed
            # against this protein reference before it can support PM5.
            candidate = self._candidate(uid, annotation, query_variant, query_splice_score,
                                        scope="residue")
            if candidate is not None:
                comparators.append(candidate)
        return search, comparators

    def _candidate(self, variation_id, annotation, query_variant, query_splice_score,
                   scope="exact_protein_change"):
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
        if scope == "exact_protein_change" and annotation["hgvsp"] not in expressions:
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
        if scope == "exact_protein_change":
            expected = (annotation["protein_start"], annotation["ref_aa"], annotation["alt_aa"])
            if changes != [expected] or candidate_annotation.get("hgvsp") != annotation["hgvsp"]:
                return None
            comparator_alt = annotation["alt_aa"]
            comparator_change = annotation["hgvsp"]
        else:
            # PM5: the same residue and reference amino acid, a different substitution.
            if len(changes) != 1:
                return None
            position, reference, substitution = changes[0]
            if (position != annotation["protein_start"] or reference != annotation["ref_aa"]
                    or substitution == annotation["alt_aa"]):
                return None
            comparator_alt = substitution
            comparator_change = candidate_annotation.get("hgvsp")
            if not comparator_change:
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
                f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/#{scope}-"
                f'{hashlib.sha256((query_variant.key + ":" + annotation["transcript"]).encode()).hexdigest()[:16]}'
            ),
            "source": "ClinVar", "source_version": f"{self.release}:{versioned}",
            "retrieved_at": response["retrieved_at"], "quality_status": "PASS",
            "response_sha256": digest, "accession": versioned, "variation_id": variation_id,
            "transcript": annotation["transcript"], "protein_id": annotation["protein_id"],
            "protein_start": annotation["protein_start"], "ref_aa": annotation["ref_aa"],
            "alt_aa": comparator_alt, "protein_change": comparator_change,
            "query_alt_aa": annotation["alt_aa"], "query_protein_change": annotation["hgvsp"],
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
            "exact_protein_match": scope == "exact_protein_change",
            "residue_match": scope == "residue", "search_scope": scope,
            "different_nucleotide_variant": True,
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


class ClinVarHotspotProvider:
    """Local missense density around a residue, used only as a hotspot proxy.

    ClinVar density says where pathogenic variants have been *reported*, which tracks how
    often a gene is tested as much as biology. It can therefore support the PM1 hotspot
    route, but never the critical-domain route, and the counted accessions are kept so the
    density can be re-checked against a later ClinVar release.
    """

    name = "ClinVar protein hotspot density"
    method = "clinvar_local_density"
    summary_batch = 200
    search_retmax = 5000
    PATHOGENIC = ("pathogenic", "likely pathogenic")
    BENIGN = ("benign", "likely benign")

    def __init__(self, client, release, policy):
        if not release:
            raise ValueError("ClinVar release date/version must be recorded")
        missing = [field for field in ("window_aa", "min_pathogenic", "max_benign",
                                       "method", "policy_source", "policy_version")
                   if policy.get(field) is None or policy[field] == ""]
        if missing:
            raise ValueError("PM1 hotspot policy is incomplete: " + ",".join(missing))
        if policy["method"] != self.method:
            raise ValueError(f"PM1 hotspot policy method must be {self.method}")
        if type(policy["window_aa"]) is not int or policy["window_aa"] < 0:
            raise ValueError("PM1 hotspot window must be a non-negative integer")
        self.client = client
        self.release = release
        self.policy = policy

    def search_hotspot(self, annotation, query_variant):
        required = ("transcript", "gene", "protein_id", "protein_start")
        if not all(annotation.get(key) for key in required):
            raise ValueError("PM1 hotspot search requires complete transcript/protein annotation")
        window = self.policy["window_aa"]
        first = annotation["protein_start"]
        last = annotation.get("protein_end") or first
        start, end = max(1, first - window), last + window
        found = gene_missense_search(self.client, self.release, annotation["gene"],
                                     self.search_retmax)
        term, count, complete = found["term"], found["count"], found["complete"]
        digest, response = found["digest"], {"retrieved_at": found["retrieved_at"]}
        search = {
            "category": "region_search", "variant_key": query_variant.key,
            "evidence_id": f'urn:sha256:{digest}:clinvar-pm1-search:{annotation["transcript"]}',
            "source": self.name, "source_version": self.release,
            "retrieved_at": response["retrieved_at"], "quality_status": "PASS",
            "transcript": annotation["transcript"], "protein_id": annotation["protein_id"],
            "gene": annotation["gene"], "query": term, "returned_count": count,
            "protein_interval": {"start": start, "end": end},
            "window_aa": window, "method": self.method,
            "policy_version": self.policy["policy_version"], "complete": complete,
        }
        if not complete:
            # A truncated gene search cannot bound the counts, so emit no density evidence.
            search["quality_status"] = "INCOMPLETE_SEARCH"
            return search, None
        counted, retrieved_at = self._counts(found["ids"], annotation, start, end,
                                             response["retrieved_at"], query_variant)
        pathogenic = [item for item in counted if item["bucket"] == "pathogenic"]
        benign = [item for item in counted if item["bucket"] == "benign"]
        search["counted_variants"] = counted
        region = {
            "category": "region", "variant_key": query_variant.key,
            "evidence_id": (f'urn:sha256:{digest}:clinvar-pm1-hotspot:'
                            f'{annotation["protein_id"]}:{start}-{end}'),
            "source": self.name, "source_version": f"{self.release}:{self.policy['policy_version']}",
            "retrieved_at": retrieved_at, "quality_status": "PASS",
            "transcript": annotation["transcript"], "protein_id": annotation["protein_id"],
            "gene": annotation["gene"], "start": start, "end": end,
            "region_type": "mutational_hotspot",
            # Derived from reported density; a human curator never reviewed this region.
            "assessment_method": "automated",
            "method": self.method, "policy_version": self.policy["policy_version"],
            "policy_source": self.policy["policy_source"], "window_aa": window,
            "min_pathogenic": self.policy["min_pathogenic"],
            "max_benign": self.policy["max_benign"],
            "pathogenic_count": len(pathogenic), "benign_count": len(benign),
            "self_excluded": sum(item["bucket"] == "query_variant" for item in counted),
            "unplaced_excluded": sum(item["bucket"] == "unplaced_on_protein" for item in counted),
            "outside_window_excluded": sum(item["bucket"] == "outside_window_on_protein"
                                           for item in counted),
            "counted_variants": counted,
            "counted_conditions": sorted({condition for item in pathogenic
                                          for condition in item["conditions"]}),
            "consequence_scope": "missense variant",
            "protein_change_source": "ClinVar esummary gene-level protein_change",
            "use_restriction": "HOTSPOT_PROXY_ONLY_NOT_CRITICAL_DOMAIN",
        }
        return search, region

    def _counts(self, ids, annotation, start, end, retrieved_at, query_variant):
        documents, summarised_at = summary_documents(self.client, self.release, ids,
                                                     self.summary_batch)
        counted = [item for item in
                   (self._summary(uid, document, annotation, start, end, query_variant)
                    for uid, document in documents.items()) if item is not None]
        return counted, max(retrieved_at, summarised_at or retrieved_at)

    def _summary(self, uid, document, annotation, start, end, query_variant):
        gene = document.get("gene_sort") or ""
        genes = {value.get("symbol") for value in document.get("genes", [])
                 if isinstance(value, dict)}
        if gene != annotation["gene"] and annotation["gene"] not in genes:
            return None
        classification = document.get("germline_classification") \
            or document.get("clinical_significance") or {}
        description = (classification.get("description") or "").strip().lower()
        if "conflicting" in description:
            bucket = "conflicting"
        elif description in self.PATHOGENIC:
            bucket = "pathogenic"
        elif description in self.BENIGN:
            bucket = "benign"
        else:
            return None
        reported = protein_change_positions(document.get("protein_change"))
        if not any(start <= position <= end for position in reported):
            return None
        # The summary lists every transcript's numbering, and isoforms of the same gene can
        # differ by tens of residues, so the window hit is only a candidate until the change
        # is read off this protein reference.
        confirmed = self._protein_positions(uid, annotation["protein_id"])
        in_window = [position for position in confirmed if start <= position <= end]
        # The variant under evaluation must not support its own hotspot density.
        if self._is_query_variant(document, query_variant):
            bucket = "query_variant"
        elif not confirmed:
            bucket = "unplaced_on_protein"
        elif not in_window:
            bucket = "outside_window_on_protein"
        trait = classification.get("trait_set") or document.get("trait_set") or []
        return {
            "variation_id": uid, "accession": document.get("accession"),
            "protein_change": document.get("protein_change"),
            "reported_positions": reported, "protein_positions": in_window or confirmed,
            "position_source": f'ClinVar VCV protein expression on {annotation["protein_id"]}',
            "classification": classification.get("description"),
            "review_status": classification.get("review_status"), "bucket": bucket,
            "conditions": sorted({value.get("trait_name") for value in trait
                                  if isinstance(value, dict) and value.get("trait_name")}),
        }

    def _protein_positions(self, variation_id, protein_id):
        """Residues of single-substitution changes ClinVar reports on this protein."""
        query = urlencode({"db": "clinvar", "id": variation_id, "rettype": "vcv",
                           "is_variationid": "true"})
        response = self.client.fetch(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}",
            response_format="text", dataset_version=self.release,
        )
        try:
            root = ElementTree.fromstring(response["body"])
        except ElementTree.ParseError as exc:
            raise ValueError("Invalid ClinVar hotspot candidate XML") from exc
        prefix = protein_id + ":"
        positions = set()
        for node in root.iter():
            if local_name(node.tag) != "Expression":
                continue
            value = node_text(node)
            if not value or not value.startswith(prefix):
                continue
            match = PROTEIN_HGVS.fullmatch(value[len(prefix):])
            if match:
                positions.add(int(match.group(2)))
        return sorted(positions)

    @staticmethod
    def _is_query_variant(document, query_variant):
        """Allele-level identity only; a different ALT at the same position is real evidence."""
        for entry in document.get("variation_set", []):
            if not isinstance(entry, dict):
                continue
            spdi = entry.get("canonical_spdi")
            if isinstance(spdi, str):
                parts = spdi.rsplit(":", 3)
                if len(parts) == 4 and parts[2] == query_variant.ref and parts[3] == query_variant.alt:
                    try:
                        interbase = int(parts[1])
                    except ValueError:
                        interbase = None
                    if interbase is not None and interbase + 1 == query_variant.pos:
                        return True
            for location in entry.get("variation_loc", []):
                if not isinstance(location, dict) or location.get("assembly_name") != "GRCh38":
                    continue
                if (location.get("chr") == query_variant.chrom
                        and location.get("start") == str(query_variant.pos)
                        and location.get("ref") == query_variant.ref
                        and location.get("alt") == query_variant.alt):
                    return True
        return False

