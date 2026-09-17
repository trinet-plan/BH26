"""
acmg_pipeline/criteria/curator_info.py

Informational-only helpers for the "use AI to check ..." asks in doc/recs
for expert board.docx - NOT ACMG judgment logic. Per the user's explicit
decision (2026-09-16): these functions return structured FACTS for a
curator (or an LLM prompt built elsewhere) to read, never a MET/NOT_MET
verdict or a strength - PM1/PM5/PM3/PP1's real judgment logic stays out of
scope here exactly as it does for acmg_pipeline/criteria/stubs.py's other
Layer-1/clinical-record codes. This module exists because the doc's ask
for these four codes is richer than a single reference URL (see
acmg_pipeline.criteria.reference_links, which covers the "show a page"
half of the same doc) - it specifically asks to summarize/check something
ABOUT that page's content, which needs the content fetched and structured,
not just linked to.

[Source: doc/recs for expert board.docx, verbatim]
  PM1, PM5: "Can show variant location/domain, and can use AI to check if
            hotspot/any nearby benign variants."
  PM3:      "Use AI to check if in trans previously reported (clinvar)"
  PP1:      "Use AI to check if segregation previously reported (clinvar)"

[This module's scope, 2026-09-16]
  uniprot_domain_context() (PM1/PM5): implemented - UniProt's own REST API
    (rest.uniprot.org) carries curated "Natural variant" feature
    annotations for well-studied disease genes, each with a clinical
    description (e.g. "in CMH1"), which already IS most of "hotspot/nearby
    benign variants" without needing a separate ClinVar call or an LLM
    summarization step - see that function's own docstring for the real
    MYH7 p.Arg719Trp validation.

  clinvar_report_context() (PM3/PP1): implemented - ClinVar's own VCV
    record XML (efetch, db=clinvar) carries free-text submitter
    interpretation <Comment> elements that, for a well-studied variant,
    directly narrate exactly what PM3/PP1 ask about in plain English (e.g.
    "segregated with disease in >20 affected relatives", "compound
    heterozygosity for the R719W and M349T mutations") - no separate LLM
    summarization needed to surface it, only a keyword filter over
    comments ClinVar submitters already wrote. See that function's own
    docstring for the real MYH7 p.Arg719Trp validation.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import requests

from acmg_pipeline.criteria.reference_links import gene_to_uniprot_accession
from acmg_pipeline.vcf_record import VariantRecord

_UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
_UNIPROT_TIMEOUT_SEC = 15

# Feature types worth reporting as "location/domain" context - the doc's
# own wording ("variant location/domain") maps to UniProt's structural/
# functional annotation types, not every feature type UniProt has (e.g.
# "Sequence conflict" or "Modified residue" aren't domain/location facts).
_DOMAIN_FEATURE_TYPES = {"Domain", "Region", "Chain", "Coiled coil", "Binding site", "Active site", "Motif"}

_HGVSP_POSITION_RE = re.compile(r"p\.\(?[A-Za-z]{3}(\d+)[A-Za-z*]{1,3}\)?")


@dataclass
class DomainFeature:
    feature_type: str
    start: int
    end: int
    description: str = ""


@dataclass
class NearbyVariant:
    position: int
    description: str  # UniProt's own free-text clinical annotation, e.g. "in CMH1; dbSNP:rs121913641"


@dataclass
class UniprotDomainContext:
    accession: str
    position: int
    covering_domains: list[DomainFeature] = field(default_factory=list)
    nearby_variants: list[NearbyVariant] = field(default_factory=list)  # within `window` residues, position-sorted


def _extract_protein_position(hgvsp: str) -> int | None:
    """
    Pulls the numeric codon position out of an HGVSP string, e.g.
    "p.(Arg719Trp)" -> 719, "p.Arg719Trp" -> 719. Works for missense
    notation specifically (PM1/PM5 are both missense-focused criteria per
    their ACMG definitions); frameshift/nonsense notation like
    "p.(Lys93ArgfsTer3)" also matches (the leading digits before "fs"/"*"),
    but the returned position is less meaningful there since PM1/PM5's own
    "hotspot"/"same-residue-different-change" concepts are missense-
    specific - callers should treat a non-missense HGVSP's result with
    that caveat in mind rather than this function refusing to return one.
    """
    if not hgvsp or hgvsp == "N/A":
        return None
    m = _HGVSP_POSITION_RE.search(hgvsp)
    return int(m.group(1)) if m else None


def uniprot_domain_context(variant: VariantRecord, window: int = 10) -> UniprotDomainContext | None:
    """
    Criteria: PM1, PM5 (doc: "show variant location/domain, use AI to
      check hotspot/any nearby benign variants"). Informational only - no
      hotspot/benign-nearby VERDICT is computed here; this returns the raw
      facts a curator (or an LLM prompt built elsewhere) needs to make
      that call themselves.
    Site: UniProt's REST API (rest.uniprot.org/uniprotkb/{accession}.json) -
      the SAME service uniprot_page_url() in reference_links.py links a
      curator to manually; this function fetches its actual content
      instead of just linking to it.
    Logic:
      1. variant.info["GENE"] -> UniProt accession, via reference_links.
         gene_to_uniprot_accession() (shared, not duplicated).
      2. variant.info["HGVSP"] -> protein position, via
         _extract_protein_position() above.
      3. Fetch the UniProt entry's `features` list. Keep every Domain/
         Region/Chain/Coiled coil/Binding site/Active site/Motif feature
         whose [start, end] range covers the target position
         (covering_domains) - this is the "location/domain" half.
      4. Keep every "Natural variant" feature within `window` residues of
         the target position (nearby_variants), each carrying UniProt's
         own free-text clinical description (e.g. "in CMH1; dbSNP:...") -
         this is the "hotspot/nearby benign variants" half: for a
         well-studied disease gene, UniProt's curators have already
         tagged nearby variants with their disease association, which is
         most of what a hotspot check needs without a separate ClinVar
         fetch or an LLM summarization pass.
      Returns None if GENE/HGVSP/accession/position can't be resolved, or
      the UniProt fetch fails - never raises for a missing/network-error
      case, matching gene_to_uniprot_accession()'s own convention.

    [Validated against real data, 2026-09-16]
      MYH7 p.(Arg719Trp) (accession P12883, position 719): covering_domains
      = [Domain "Myosin motor" 85-778] (PM1's "critical/well-established
      functional domain" - directly confirmed, not inferred). nearby_
      variants (window=10) = 7 entries at positions 712-728, EVERY ONE
      annotated "in CMH1" (hypertrophic cardiomyopathy) with no benign-
      annotated variant among them - including TWO other entries at
      position 719 ITSELF (PM5's own "different missense change at the
      same residue" signal, found directly rather than inferred).
    """
    gene = variant.info.get("GENE")
    hgvsp = variant.info.get("HGVSP")
    if not gene or not hgvsp:
        return None
    position = _extract_protein_position(hgvsp)
    if position is None:
        return None
    accession = gene_to_uniprot_accession(gene)
    if not accession:
        return None

    try:
        resp = requests.get(_UNIPROT_ENTRY_URL.format(accession=accession), timeout=_UNIPROT_TIMEOUT_SEC)
        resp.raise_for_status()
        features = resp.json().get("features", [])
    except (requests.RequestException, ValueError):
        return None

    covering_domains: list[DomainFeature] = []
    nearby_variants: list[NearbyVariant] = []
    for f in features:
        try:
            start = f["location"]["start"]["value"]
            end = f["location"]["end"]["value"]
        except (KeyError, TypeError):
            continue
        if f.get("type") in _DOMAIN_FEATURE_TYPES and start <= position <= end:
            covering_domains.append(DomainFeature(
                feature_type=f["type"], start=start, end=end, description=f.get("description", ""),
            ))
        elif f.get("type") == "Natural variant" and abs(start - position) <= window:
            nearby_variants.append(NearbyVariant(position=start, description=f.get("description", "")))

    nearby_variants.sort(key=lambda v: v.position)
    return UniprotDomainContext(
        accession=accession, position=position,
        covering_domains=covering_domains, nearby_variants=nearby_variants,
    )


# ============================================================================
# PM3, PP1 - "check if previously reported in ClinVar"
# ============================================================================

_CLINVAR_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_CLINVAR_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_CLINVAR_TIMEOUT_SEC = 15

# Keyword filters over ClinVar submitters' own free-text <Comment> elements -
# NOT a judgment, just picking out which of a variant's (often 10+) comments
# are even relevant to segregation (PP1) vs. trans-configuration (PM3), since
# most comments discuss neither and would just be noise if returned wholesale.
_SEGREGATION_KEYWORDS = ("segregat",)  # matches segregate/segregated/segregation/co-segregat(e/ed/ion)
_TRANS_KEYWORDS = ("in trans", "compound heterozygo", "co-occurring")


@dataclass
class ClinVarComment:
    text: str
    matched_keyword: str


@dataclass
class ClinVarReportContext:
    variation_id: str
    segregation_mentions: list[ClinVarComment] = field(default_factory=list)  # for PP1
    trans_mentions: list[ClinVarComment] = field(default_factory=list)        # for PM3


def _resolve_clinvar_variation_id(variant: VariantRecord) -> str | None:
    """
    Prefers variant.info["CLNVARIATIONID"] when it's already a real ID (this
    project's demo VCFs carry it as e.g. "VCV000014104", or the placeholder
    "not_registered_novel" for variants ClinVar has never seen - the latter
    is deliberately rejected here rather than passed through as if it were
    a real ID). Falls back to a ClinVar ESearch (same gene+HGVSc query
    reference_links.clinvar_search_url() builds a link for, but resolved to
    an actual ID here) ONLY when it returns exactly one hit - confirmed
    empirically that a missense variant's query can return >1 hit (MYH7
    c.2155C>T's own query also matches the neighboring p.Arg719Gln variant,
    VariationID 14107, at the same codon); guessing which of several
    candidates is the right one would risk attributing another variant's
    ClinVar comments to this one, so this function returns None rather than
    guess when the search isn't unambiguous.
    """
    raw = variant.info.get("CLNVARIATIONID", "")
    m = re.search(r"(\d+)", raw) if raw and raw != "not_registered_novel" else None
    if m:
        # ClinVar VariationIDs are zero-padded in the "VCV000014104" form
        # this project's demo VCFs use, but NOT zero-padded as bare IDs (the
        # EFetch id= param, and ClinVar's own URLs) - int() then str() strips
        # the padding ("000014104" -> "14104") rather than passing the
        # padded string through, which the bare esearch/efetch ID space
        # doesn't use (confirmed empirically: NCBI's endpoints happened to
        # tolerate the padded form when tested, but normalizing here avoids
        # depending on that leniency).
        return str(int(m.group(1)))

    gene = variant.info.get("GENE")
    hgvsc = variant.info.get("HGVSC")
    if not gene or not hgvsc:
        return None
    try:
        resp = requests.get(
            _CLINVAR_ESEARCH_URL,
            params={"db": "clinvar", "term": f"{gene}[gene] AND {hgvsc}", "retmode": "json"},
            timeout=_CLINVAR_TIMEOUT_SEC,
        )
        resp.raise_for_status()
        id_list = resp.json().get("esearchresult", {}).get("idlist", [])
    except (requests.RequestException, ValueError):
        return None
    return id_list[0] if len(id_list) == 1 else None


def clinvar_report_context(variant: VariantRecord) -> ClinVarReportContext | None:
    """
    Criteria: PM3 (doc: "use AI to check if in trans previously reported
      (clinvar)"), PP1 (doc: "use AI to check if segregation previously
      reported (clinvar)"). Informational only - no PM3/PP1 MET/NOT_MET
      verdict is computed here (in particular, this says nothing about
      whether a trans-observation was with a PATHOGENIC variant, which is
      what PM3 itself actually requires - only that ClinVar submitters
      mention a trans/compound-het observation at all).
    Site: ClinVar (ncbi.nlm.nih.gov/clinvar), via NCBI's public EFetch API
      (VCV record XML, db=clinvar) - the same record reference_links.
      clinvar_search_url() links a curator to search for manually; this
      function fetches its actual submitter-comment content instead.
    Logic:
      1. Resolve a ClinVar VariationID via _resolve_clinvar_variation_id()
         above.
      2. EFetch that VCV record's full XML and collect the text of every
         <Comment> element (ClinVar submitters' own free-text
         interpretation summaries, e.g. "This variant ... segregated with
         disease in >20 affected relatives") AND every
         <Attribute Type="Description"> element (curated case-report-style
         free text ClinVar embeds from sources like OMIM - confirmed
         empirically to be where some trans/compound-het narratives
         actually live, NOT inside <Comment>; see the validation note
         below for why both were needed).
      3. Split into segregation_mentions (text contains a _SEGREGATION_
         KEYWORDS match) and trans_mentions (a _TRANS_KEYWORDS match) - a
         text block can appear in both, or neither and be dropped entirely
         (most of a variant's comments/attributes discuss neither topic).
      Returns None if no VariationID could be resolved or the EFetch fails.

    [Validated against real data, 2026-09-16 - including a self-caught bug]
      First pass only scanned <Comment> elements: MYH7 c.2155C>T
      (VariationID 14104) returned 8 segregation_mentions but 0
      trans_mentions, even though a real compound-heterozygosity narrative
      ("...two missense mutations: one was the R719W mutation and the
      other was an M349T mutation... The authors hypothesized that
      compound heterozygosity for the R719W and M349T mutations resulted
      in the particularly severe phenotype") is genuinely present in this
      record - it lives inside an <Attribute Type="Description"> element,
      not a <Comment>. Expanded to scan both element types, confirmed
      against the same record: now correctly returns that trans_mentions
      hit. (A separate bug in _resolve_clinvar_variation_id() - not
      stripping "VCV000014104"'s zero-padding - was also caught and fixed
      during this same validation pass; see that function's docstring.)
    Field: EvidenceLine.extensions (Extension(name="curatorInfo")) once
      wired through export.build_stub_evidence_line() for PM3/PP1, the
      same convention as uniprot_domain_context() above for PM1/PM5.
    """
    variation_id = _resolve_clinvar_variation_id(variant)
    if variation_id is None:
        return None
    try:
        resp = requests.get(
            _CLINVAR_EFETCH_URL,
            params={"db": "clinvar", "id": variation_id, "is_variationid": "true",
                    "rettype": "vcv", "retmode": "xml"},
            timeout=_CLINVAR_TIMEOUT_SEC,
        )
        resp.raise_for_status()
        root = ET.fromstring(resp.content)
    except (requests.RequestException, ET.ParseError):
        return None

    text_elements = root.findall(".//Comment") + root.findall(".//Attribute[@Type='Description']")

    segregation_mentions: list[ClinVarComment] = []
    trans_mentions: list[ClinVarComment] = []
    for comment_el in text_elements:
        text = "".join(comment_el.itertext()).strip()
        if not text:
            continue
        lower = text.lower()
        for kw in _SEGREGATION_KEYWORDS:
            if kw in lower:
                segregation_mentions.append(ClinVarComment(text=text, matched_keyword=kw))
                break
        for kw in _TRANS_KEYWORDS:
            if kw in lower:
                trans_mentions.append(ClinVarComment(text=text, matched_keyword=kw))
                break

    return ClinVarReportContext(
        variation_id=variation_id,
        segregation_mentions=segregation_mentions,
        trans_mentions=trans_mentions,
    )
