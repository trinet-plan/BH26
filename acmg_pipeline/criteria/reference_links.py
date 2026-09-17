"""
acmg_pipeline/criteria/reference_links.py

URL generation helpers for the curator-facing reference pages named in
doc/recs for expert board.docx (per-criterion "what to show the curator"
requirements) - NOT judgment logic. This project does not implement PVS1/
PS1/PM1/PM2/PM3/PM5/PP1/PP2/BA1/BS1/BS2/BP4's own decision logic (see
acmg_pipeline/criteria/stubs.py: PVS1/PS1/PM2/PP2/BA1/BS1/BS2 are Layer-1
codes, another team's responsibility; PP1 was this project's own
literature/segregation judgment for months (acmg_pipeline/criteria/
segregation.py) but was handed off to another team on 2026-09-16 - see
stubs.py's HANDED_OFF_TO_OTHER_TEAM - and is now a stub here too, like the
Layer-1 codes). These functions only build a reference URL a curator (or
an AI assistant acting on the curator's behalf) can open - they never
fetch, parse, or interpret the destination page's content.

Takes VariantRecord (acmg_pipeline.vcf_record) directly as the parameter,
per the same "argument is a pulled-in entity class, not a bespoke wrapper"
convention as acmg_pipeline.criteria.segregation.from_clinical_note_family() -
a later VCF INFO field addition needs no signature change here.

[Where the returned URL ends up in VA-Spec output]
  Every function below just returns a plain `str | None` - none of them
  touch VA-Spec themselves. The one caller today, acmg_pipeline.export.
  build_stub_evidence_line(), puts that string into ONE specific place:
  EvidenceLine.extensions, as Extension(name="referenceLink", value=url) -
  deliberately NOT EvidenceLine.reportedIn (a Document there would be
  spec-valid with just `urls` set, no `pmid` needed, but reportedIn means
  "the source this evidence was reported in", and no evaluation of these
  pages ever happens here - see that function's own docstring, "Why
  extensions, not reportedIn", for the full reasoning; this was a
  deliberate call, don't move it without reading that first). Each
  function's own docstring below repeats this as a one-line "Field:" note
  so it doesn't require jumping to export.py to know where its own output
  lands.

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

# Criteria: PVS1 (doc: "Add autoPVS1 chart data/result"). Site: AutoPVS1
# (autopvs1.bgi.com). Logic: none - AutoPVS1 is a form-based web tool with
# no deep-link/query-string API this project has found, so this is a fixed
# URL, not built per-variant. Referenced directly in _URL_BUILDERS below
# rather than via its own function, since there's no per-variant logic to
# name a function after.
AUTOPVS1_URL = "https://autopvs1.bgi.com"

_TOGOID_CONVERT_URL = "https://api.togoid.dbcls.jp/convert"
_TOGOID_TIMEOUT_SEC = 15


def gene_to_uniprot_accession(gene_symbol: str) -> str | None:
    """
    Helper for uniprot_page_url() (PM1/PM5) below - not itself tied to a
    specific ACMG code.

    Site: TogoID's public REST API (api.togoid.dbcls.jp/convert) - the
      same ID-conversion service the TogoMCP connector's togoid_convertId
      tool calls internally; called directly here since TogoMCP itself
      isn't reachable from this project's own runtime code (it's a
      connector wired into the Claude Code session, not a server this
      project's Python process can dial into the way it does
      pubmed.mcp.claude.com - see acmg_pipeline/pipeline.py's
      connect_pubmed()).
    Logic: resolves a gene symbol (e.g. "MYH7") to its UniProt accession
      (e.g. "P12883") via the 2-hop route "hgnc_symbol,hgnc,uniprot" - a
      direct hgnc_symbol->uniprot route does not exist in TogoID (confirmed
      empirically: requesting it returns HTTP 400 "no route", with the
      2-hop route given in the error's own suggestion). Returns None if the
      gene symbol doesn't resolve (unknown symbol, network failure, or a
      real but symbol-less gene) rather than raising - callers building an
      optional reference link shouldn't crash the whole page over one
      missing cross-reference.
    Field: not applicable here - this returns a bare accession string, not
      a URL. See uniprot_page_url() (its only caller) for where its output
      ultimately lands.
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
    """
    Criteria: PM1, PM5 (doc: "show variant location/domain, use AI to
      check hotspot/nearby benign variants" - the domain/hotspot judgment
      itself is still NOT automated here; this only gets the reference
      page a curator or PM1/PM5's real implementer would open to do it).
    Site: UniProt (uniprot.org), the entry page for the gene's canonical
      protein - shows domain boundaries, active sites, and known variants
      in context, which is what "location/domain" means in practice.
    Logic: variant.info["GENE"] (a symbol, e.g. "MYH7") -> UniProt
      accession (e.g. "P12883") via gene_to_uniprot_accession() above, then
      f"https://www.uniprot.org/uniprotkb/{accession}/entry". Returns None
      if GENE is missing or the accession lookup fails (unknown symbol,
      network error) - no accession means no valid page to link to.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_stub_evidence_line() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    gene = variant.info.get("GENE")
    if not gene:
        return None
    accession = gene_to_uniprot_accession(gene)
    if not accession:
        return None
    return f"https://www.uniprot.org/uniprotkb/{accession}/entry"


def clinvar_search_url(variant: VariantRecord) -> str | None:
    """
    Criteria: PS1 ("show the clinvar page/summary for that codon"), PM3
      ("use AI to check if in trans previously reported (clinvar)"), PP1
      ("use AI to check if segregation previously reported (clinvar)") -
      the doc phrases these three differently but all three point at the
      same action (open this variant's ClinVar record), so one function
      serves all three; the AI-driven "check if X was previously reported"
      part is NOT done here, only the page to check it on.
    Site: ClinVar (ncbi.nlm.nih.gov/clinvar), NCBI's search UI.
    Logic: a term-search URL built as "<GENE>[gene] AND <HGVSc>" (e.g.
      "MYH7[gene] AND c.2155C>T"), URL-encoded and appended to
      ".../clinvar/?term=...". This is a SEARCH, not a direct accession
      link - this project doesn't look up the ClinVar VariationID, so the
      URL lands on ClinVar's own search results page rather than one
      specific record. Returns None if GENE or HGVSC is missing.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_stub_evidence_line() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    gene = variant.info.get("GENE")
    hgvsc = variant.info.get("HGVSC")
    if not gene or not hgvsc:
        return None
    term = f"{gene}[gene] AND {hgvsc}"
    return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={quote(term)}"


def gnomad_variant_url(variant: VariantRecord, dataset: str = "gnomad_r4") -> str | None:
    """
    Criteria: PM2, BA1, BS1, BS2 (doc: "show the gnomAD af page" - all four
      are population-frequency codes read off the same variant-level AF
      figure, just compared against different thresholds elsewhere).
    Site: gnomAD (gnomad.broadinstitute.org), the per-variant page (allele
      frequency broken out by population, plus quality/coverage flags).
    Logic: gnomAD's variant URL is positional, not accession-based:
      f".../variant/{chrom}-{pos}-{ref}-{alt}?dataset={dataset}" straight
      from VariantRecord's own chrom/pos/ref/alt fields (no external
      lookup needed - unlike uniprot_page_url() above). `dataset` defaults
      to "gnomad_r4" (GRCh38); pass "gnomad_r2_1" for a GRCh37-sourced
      VariantRecord (not exercised against real data - this project's own
      demo VCFs all declare "##reference=GRCh38"). Returns None when
      ref/alt aren't concrete alleles (e.g. this project's demo data uses
      ALT="." for some indels/deletions - not a real gnomAD-linkable
      representation, so no URL is better than a broken one).
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_stub_evidence_line() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    if not variant.ref or not variant.alt or variant.alt in (".", ""):
        return None
    return (f"https://gnomad.broadinstitute.org/variant/"
            f"{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}?dataset={dataset}")


def gnomad_gene_url(variant: VariantRecord, dataset: str = "gnomad_r4") -> str | None:
    """
    Criteria: PP2 (doc: "get the gnomAD zscore" - the missense Z-score is a
      GENE-level constraint metric, unlike PM2/BA1/BS1/BS2's variant-level
      AF, so this points at gnomAD's gene page instead of its variant page).
    Site: gnomAD (gnomad.broadinstitute.org), the gene page's "Constraint"
      section (lists the missense Z-score and pLI/LOEUF alongside it).
    Logic: f".../gene/{gene}?dataset={dataset}" from variant.info["GENE"]
      directly - no lookup needed, gnomAD's gene URLs take a bare symbol.
      Same `dataset` genome-build caveat as gnomad_variant_url() above.
      Returns None if GENE is missing.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_stub_evidence_line() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
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

# Public: the codes this module can build a reference URL for at all (a
# per-variant call may still return None for one of these - e.g.
# gnomad_variant_url() on an ALT="." record - but the code itself has a
# defined URL strategy). Callers building a UI/EvidenceLine per code (see
# acmg_pipeline.export.build_stub_evidence_line()) can use this to know
# which of the 25 stub codes are even worth trying.
CODES_WITH_REFERENCE_URL = frozenset(_URL_BUILDERS)


def reference_url_for_criterion(code: str, variant: VariantRecord) -> str | None:
    """
    Single entry point: the reference URL doc/recs for expert board.docx
    asks to show for `code`, given `variant` - or None if this code has no
    URL-able reference page (see _URL_BUILDERS) or the variant lacks the
    fields needed to build one (e.g. ALT="." for gnomad_variant_url).
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_stub_evidence_line() (its actual caller)
      - see this module's own docstring, "Where the returned URL ends up
      in VA-Spec output".
    """
    builder = _URL_BUILDERS.get(code)
    if builder is None:
        return None
    return builder(variant)


def all_reference_urls(variant: VariantRecord) -> dict[str, str]:
    """
    Every code -> URL pair this module can build for `variant`, skipping
    codes where the URL couldn't be constructed (rather than including a
    None).
    Field: same as reference_url_for_criterion() above - each value here
      lands in that code's own EvidenceLine.extensions, not one shared
      bundle (see export.build_stub_evidence_line(), called once per code).
    """
    return {
        code: url
        for code in _URL_BUILDERS
        if (url := reference_url_for_criterion(code, variant)) is not None
    }
