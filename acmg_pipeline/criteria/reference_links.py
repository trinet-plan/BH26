"""
acmg_pipeline/criteria/reference_links.py

URL generation helpers for the curator-facing reference pages named in
doc/recs for expert board.docx (per-criterion "what to show the curator"
requirements) - NOT judgment logic. This project does not implement PVS1/
PS1/PM1/PM2/PM3/PM5/PP1/PP2/BA1/BS1/BS2/BP4's own decision logic (see
acmg_pipeline/criteria/stubs.py - Layer-1 codes are another team's
responsibility, PP1's literature/clinical-note judgment is this project's
own but is separate from this doc's "show the ClinVar page" UI ask). These
functions only build a reference URL a curator (or an AI assistant acting
on the curator's behalf) can open - they never fetch, parse, or interpret
the destination page's content.

Takes VariantRecord (acmg_pipeline.vcf_record) directly as the parameter,
per the same "argument is a pulled-in entity class, not a bespoke wrapper"
convention as acmg_pipeline.criteria.segregation.from_clinical_note_family() -
a later VCF INFO field addition needs no signature change here.

[Source: doc/recs for expert board.docx, verbatim per-criterion asks]
  PVS1:            "Add autoPVS1 chart data/result" -> https://autopvs1.bgi.com
  PS1:              "show the clinvar page/summary for that codon"
  PM1, PM5:         "show variant location/domain, use AI to check hotspot/
                     nearby benign variants" -> the domain/hotspot judgment
                     itself is still not URL-able, but the UniProt entry
                     page needed to show it is, via a gene-symbol ->
                     UniProt-accession lookup (see gene_to_uniprot_accession())
  PM2, BA1, BS1, BS2: "show the gnomAD af page"
  PM3:              "use AI to check if in trans previously reported (ClinVar)"
  PP1:              "use AI to check if segregation previously reported (ClinVar)"
  PP2:              "get the gnomAD zscore" (gene-level, not variant-level)
  PP3, BP4:         "show the predictors, follow thresholds" -> NOT URL-able
                     (a static threshold table, not a single reference page)

[Verified reachable, 2026-09-16]
  gnomad_variant_url() and gnomad_gene_url() confirmed HTTP 200 against a
  real variant (MYH7 c.2155C>T) and gene (MYH7). clinvar_search_url()
  confirmed to 302-redirect to a real ClinVar search results page.
  autopvs1_url() confirmed HTTP 200 (a fixed URL, not deep-linkable to a
  specific variant - AutoPVS1 is a form-based tool). gene_to_uniprot_
  accession()/uniprot_page_url() confirmed against MYH7 -> P12883 (the
  real UniProt accession for human beta-myosin heavy chain), via TogoID's
  public REST API (api.togoid.dbcls.jp) - the same conversion the TogoMCP
  connector's togoid_convertId tool uses internally (route "hgnc_symbol,
  hgnc,uniprot"; a direct gene-symbol->uniprot route doesn't exist in
  TogoID, only this 2-hop one, per that tool's own error message).
  TogoMCP itself isn't reachable from acmg_pipeline's own runtime code (it's
  a connector wired into this session, not a server this project's Python
  process can dial into the way it does pubmed.mcp.claude.com) - calling
  the same public REST endpoint TogoMCP itself calls gets the identical
  result without that dependency.

  Genome build: assumes GRCh38 (dataset=gnomad_r4), matching this
  project's demo VCF data's own "##reference=GRCh38" declaration - a
  GRCh37-sourced VariantRecord would need dataset=gnomad_r2_1 instead
  (not handled here; no GRCh37 demo data exists in this project to
  validate against).
"""

from __future__ import annotations

from urllib.parse import quote

import requests

from acmg_pipeline.vcf_record import VariantRecord

AUTOPVS1_URL = "https://autopvs1.bgi.com"
_TOGOID_CONVERT_URL = "https://api.togoid.dbcls.jp/convert"
_TOGOID_TIMEOUT_SEC = 15


def gene_to_uniprot_accession(gene_symbol: str) -> str | None:
    """
    Resolves a gene symbol (e.g. "MYH7") to its UniProt accession (e.g.
    "P12883") via TogoID's public REST API - the 2-hop route "hgnc_symbol,
    hgnc,uniprot" (a direct hgnc_symbol->uniprot route doesn't exist in
    TogoID). Returns None if the gene symbol doesn't resolve (unknown
    symbol, network failure, or a real but symbol-less gene) rather than
    raising - callers building an optional reference link shouldn't crash
    the whole page over one missing cross-reference.
    """
    try:
        resp = requests.get(
            _TOGOID_CONVERT_URL,
            params={"ids": gene_symbol, "route": "hgnc_symbol,hgnc,uniprot",
                    "report": "pair", "format": "json"},
            timeout=_TOGOID_TIMEOUT_SEC,
        )
        resp.raise_for_status()
        data = resp.json()
        pairs = data.get("results", [])
    except (requests.RequestException, ValueError):
        return None
    for source_id, target_id in pairs:
        if source_id == gene_symbol:
            return target_id
    return None


def uniprot_page_url(variant: VariantRecord) -> str | None:
    """For PM1, PM5 - "show variant location/domain" (the domain/hotspot
    judgment itself is still not automatable here, but the UniProt entry
    page a curator or PM1/PM5's real implementer would consult is)."""
    gene = variant.info.get("GENE")
    if not gene:
        return None
    accession = gene_to_uniprot_accession(gene)
    if not accession:
        return None
    return f"https://www.uniprot.org/uniprotkb/{accession}/entry"


def clinvar_search_url(variant: VariantRecord) -> str | None:
    """For PS1, PM3, PP1 - "show/check the ClinVar page" (doc's own wording
    doesn't distinguish a different query per code; all three are asking
    to look at the same variant's ClinVar record)."""
    gene = variant.info.get("GENE")
    hgvsc = variant.info.get("HGVSC")
    if not gene or not hgvsc:
        return None
    term = f"{gene}[gene] AND {hgvsc}"
    return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={quote(term)}"


def gnomad_variant_url(variant: VariantRecord, dataset: str = "gnomad_r4") -> str | None:
    """For PM2, BA1, BS1, BS2 - "show the gnomAD AF page" for this specific
    variant. Returns None if ref/alt aren't concrete alleles (e.g. this
    project's demo data uses ALT="." for some indels/deletions, which
    isn't a real gnomAD-linkable allele representation)."""
    if not variant.ref or not variant.alt or variant.alt in (".", ""):
        return None
    return (f"https://gnomad.broadinstitute.org/variant/"
            f"{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}?dataset={dataset}")


def gnomad_gene_url(variant: VariantRecord, dataset: str = "gnomad_r4") -> str | None:
    """For PP2 - "get the gnomAD zscore" (a gene-level constraint metric,
    shown on the gene page, not a variant page)."""
    gene = variant.info.get("GENE")
    if not gene:
        return None
    return f"https://gnomad.broadinstitute.org/gene/{quote(gene)}?dataset={dataset}"


# code -> the reference-link function to call, for codes this module can
# produce a URL for. PM1/PM5 (needs a gene->UniProt lookup, not string
# assembly) and PP3/BP4 (a threshold table, not a single reference page)
# are deliberately absent - see the module docstring.
_URL_BUILDERS = {
    "PVS1": lambda v: AUTOPVS1_URL,
    "PS1": clinvar_search_url,
    "PM3": clinvar_search_url,
    "PP1": clinvar_search_url,
    "PM1": uniprot_page_url,
    "PM5": uniprot_page_url,
    "PM2": gnomad_variant_url,
    "BA1": gnomad_variant_url,
    "BS1": gnomad_variant_url,
    "BS2": gnomad_variant_url,
    "PP2": gnomad_gene_url,
}


def reference_url_for_criterion(code: str, variant: VariantRecord) -> str | None:
    """
    Single entry point: the reference URL doc/recs for expert board.docx
    asks to show for `code`, given `variant` - or None if this code has no
    URL-able reference page (see _URL_BUILDERS) or the variant lacks the
    fields needed to build one (e.g. ALT="." for gnomad_variant_url).
    """
    builder = _URL_BUILDERS.get(code)
    if builder is None:
        return None
    return builder(variant)


def all_reference_urls(variant: VariantRecord) -> dict[str, str]:
    """Every code -> URL pair this module can build for `variant`, skipping
    codes where the URL couldn't be constructed (rather than including a
    None)."""
    return {
        code: url
        for code in _URL_BUILDERS
        if (url := reference_url_for_criterion(code, variant)) is not None
    }
