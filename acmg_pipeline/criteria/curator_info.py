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

  PM3/PP1's "check if previously reported in ClinVar" - NOT implemented
    here yet. That needs ClinVar's own variant-level submission data (not
    just a search-results URL, which reference_links.clinvar_search_url()
    already provides), which this module doesn't fetch. Left as an open
    follow-up rather than built speculatively without real data to
    validate a query strategy against, the same discipline used elsewhere
    in this project (e.g. search_candidate_pmids() in pipeline.py was only
    built after confirming real query patterns against real demo data).
"""

from __future__ import annotations

import re
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
