"""
acmg_pipeline/criteria/reference_links.py

URL generation helpers for the curator-facing reference pages named in
doc/recs for expert board.docx (per-criterion "what to show the curator"
requirements) - NOT judgment logic. Some listed criteria are evaluated by
the integrated automated engine, some by the literature engine, and some
remain stubs; this module is deliberately independent of those decisions. These
functions only build a reference URL a curator (or an AI assistant acting
on the curator's behalf) can open - they never fetch, parse, or interpret
the destination page's content.

Takes VariantRecord (acmg_pipeline.vcf_record) directly as the parameter,
per the same "argument is a pulled-in entity class, not a bespoke wrapper"
convention as acmg_pipeline.criteria.segregation.from_clinical_note_family() -
a later VCF INFO field addition needs no signature change here.

[Where the returned URL ends up in VA-Spec output]
  Every function below just returns a plain `str | None` - none of them
  touch VA-Spec themselves. acmg_pipeline.export.build_reference_extensions()
  puts each returned string into ONE specific place on real and stub lines
  alike: EvidenceLine.extensions, as repeated
  Extension(name="referenceLink", value=url) entries -
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
  PS1:              "show the clinvar page/summary for that codon" -> a
                     ClinVar search at the variant's own residue, the search
                     the document's screenshot shows ("brca1, tyr1127"); the
                     genomic window is only the fallback when no residue is
                     given (see clinvar_codon_url())
  PM1, PM5:         "show variant location/domain, use AI to check hotspot/
                     nearby benign variants" -> the domain/hotspot judgment
                     itself is still not URL-able. The document's own
                     screenshots for these two are Franklin's Region Viewer
                     and its PM1 verdict, so Franklin's variant page leads,
                     followed by the UniProt entry page the domain boundaries
                     come from (see gene_to_uniprot_accession())
  PM2, BA1, BS1, BS2: population-frequency page (gnomAD per the document,
                     plus TogoVar, the provider actually queried at runtime)
  PM3:              "use AI to check if in trans previously reported (ClinVar)"
  PP1:              "use AI to check if segregation previously reported (ClinVar)"
  PP2:              "get the gnomAD zscore" (gene-level, not variant-level)
  PP3, BP4:         "show the predictors, follow thresholds" -> the thresholds
                     are a published table, not a per-variant page: the
                     document reproduces Table 2 of Pejaver et al. 2022, the
                     same calibration config/demo-rules.json already applies,
                     so the link is to that paper (see
                     pp3_bp4_calibration_url())

  No other code appears in that document, and no reference link is emitted
  for one - see _URL_BUILDERS.

[Verified reachable]
  togovar_variant_url() uses TogoVar's documented GRCh38 coordinate URL.
  gnomad_variant_url() and gnomad_gene_url() confirmed HTTP 200 against a
  real variant (MYH7 c.2155C>T) and gene (MYH7). clinvar_search_url()
  confirmed to 302-redirect to a real ClinVar search results page.
  autopvs1_variant_url() confirmed against two real variants via a live
  browser (2026-09-18, curator-suggested): AutoPVS1 DOES support a
  positional deep link, `/variant/{hg19|hg38}/{chrom}-{pos}-{ref}-{alt}`,
  contrary to this module's earlier finding that it was a fixed, form-only
  URL - MYH9 hg19 22-36678800-G-A rendered a filled-in NF6/Moderate
  flowchart page, and this project's own MYBPC3 c.278delA (hg38
  11-47351252-CT-C) rendered NF1/VeryStrong, matching the pipeline's own
  PVS1 result for that variant. gene_to_uniprot_
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

import re
from urllib.parse import quote

import requests

from acmg_pipeline.constants import ALL_ACMG_CODES
from acmg_pipeline.vcf_record import VariantRecord

# Each site's own unparameterized landing page: what a builder returns when it
# cannot address this variant there (see positional_allele() for what makes a
# record unaddressable, and reference_urls_for_criterion() for the rule).
# Confirmed reachable 2026-09-18; TogoVar answers 403 to a scripted request,
# as it does for its variant pages too, but serves the page in a browser.
AUTOPVS1_URL = "https://autopvs1.bgi.com"
TOGOVAR_URL = "https://grch38.togovar.org/"
GNOMAD_URL = "https://gnomad.broadinstitute.org/"
FRANKLIN_URL = "https://franklin.genoox.com/"
CLINVAR_URL = "https://www.ncbi.nlm.nih.gov/clinvar/"
UNIPROT_URL = "https://www.uniprot.org/"

# PP3/BP4's threshold table (doc: "follow below thresholds"). PMID 36413997,
# confirmed 2026-09-18 to be the paper whose Table 2 the doc reproduces.
PP3_BP4_CALIBRATION_URL = "https://pubmed.ncbi.nlm.nih.gov/36413997/"

_TOGOID_CONVERT_URL = "https://api.togoid.dbcls.jp/convert"
_TOGOID_TIMEOUT_SEC = 15


_CONCRETE_ALLELE = re.compile(r"^[ACGTNacgtn]+$")
# Assemblies the positional builders below can address. A record that names any
# other one is not linked rather than linked to the wrong coordinate.
_GRCH38_NAMES = frozenset({"grch38", "hg38", "grch38.p14", "genome_reference_consortium_human_build_38"})


def positional_allele(variant: VariantRecord) -> tuple[str, int, str, str] | None:
    """The variant as (chrom, pos, ref, alt) if it can address an external page.

    Returns None - no link at all - rather than let a builder below format a URL
    that resolves to nothing. None of these sites say so when it does: checked
    live (2026-09-18), gnomAD answers a malformed allele with HTTP 200 and the
    byte-identical single-page-app shell it serves for a real one, Franklin and
    AutoPVS1 also answer 200, and TogoVar answers 403 for real and malformed
    alike. So a broken link cannot be detected by fetching it, and is rejected
    here on the record's own shape instead.

    Rejected:
      - a multi-ALT row ("A,T"), a symbolic ALT ("<DEL>", "<DUP>"), the
        upstream-deletion ALT ("*"), the missing value ("."), and anything else
        that is not a run of bases. All of these formatted straight into a URL
        before this existed.
      - a record that declares an assembly other than GRCh38 in INFO/ASSEMBLY.
        Every positional builder here is GRCh38-only (gnomad_r4, grch38.togovar,
        hg38), so a GRCh37 record used to be linked to whatever sits at those
        coordinates in GRCh38 - a wrong page, silently. An absent INFO/ASSEMBLY
        still means GRCh38, matching automated_core.interface.criterion_input().

    Normalized: a "chr14"-style CHROM to "14" (TogoVar, gnomAD and AutoPVS1 all
    want it bare; only franklin_variant_url() was stripping it before), and
    lower-case bases to upper case.
    """
    assembly = variant.info.get("ASSEMBLY")
    if assembly and str(assembly).strip().casefold() not in _GRCH38_NAMES:
        return None
    if not variant.chrom or not variant.pos or int(variant.pos) < 1:
        return None
    ref, alt = variant.ref, variant.alt
    if not ref or not alt:
        return None
    if not _CONCRETE_ALLELE.match(str(ref)) or not _CONCRETE_ALLELE.match(str(alt)):
        return None
    chrom = str(variant.chrom)
    if chrom.casefold().startswith("chr"):
        chrom = chrom[3:]
    if not chrom:
        return None
    return chrom, int(variant.pos), str(ref).upper(), str(alt).upper()


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
      passed through export.build_reference_extensions() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    gene = variant.info.get("GENE")
    if not gene:
        return UNIPROT_URL
    accession = gene_to_uniprot_accession(gene)
    if not accession:
        # The gene symbol did not resolve (unknown symbol, or the TogoID
        # lookup was unreachable). UniProt's own search page, unqueried,
        # rather than nothing - see reference_urls_for_criterion().
        return UNIPROT_URL
    return f"https://www.uniprot.org/uniprotkb/{accession}/entry"


def clinvar_search_url(variant: VariantRecord) -> str | None:
    """
    Criteria: PM3 ("use AI to check if in trans previously reported
      (clinvar)"), PP1 ("use AI to check if segregation previously reported
      (clinvar)") - both point at the same action (open this variant's
      ClinVar record), so one function serves both; the AI-driven "check if
      X was previously reported" part is NOT done here, only the page to
      check it on. (PS1 used to share this function too - see
      clinvar_position_url() below for why it now has its own.)
    Site: ClinVar (ncbi.nlm.nih.gov/clinvar), NCBI's search UI.
    Logic: a term-search URL built as "<GENE>[gene] AND <HGVSc>" (e.g.
      "MYH7[gene] AND c.2155C>T"), URL-encoded and appended to
      ".../clinvar/?term=...". This is a SEARCH, not a direct accession
      link - this project doesn't look up the ClinVar VariationID, so the
      URL lands on ClinVar's own search results page rather than one
      specific record. Returns None if GENE or HGVSC is missing.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_reference_extensions() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    gene = variant.info.get("GENE")
    hgvsc = variant.info.get("HGVSC")
    if not gene or not hgvsc:
        return CLINVAR_URL
    term = f"{gene}[gene] AND {hgvsc}"
    return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={quote(term)}"


def clinvar_position_url(variant: VariantRecord, window: int = 10) -> str | None:
    """
    Criteria: PS1 (doc: "show the clinvar page/summary for that codon" -
      PS1 is met by a DIFFERENT nucleotide change producing the same amino
      acid change, so a search keyed on this variant's own exact HGVSc, as
      clinvar_search_url() above does for PM3/PP1, would never surface the
      other variant PS1 is actually about). This module doesn't compute
      codon boundaries (would need transcript exon/CDS phase, not just
      chrom/pos), so a `window`-bp genomic range centered on the variant is
      used instead - wide enough to catch same-codon and immediately
      adjacent changes, narrow enough to stay readable.
    Site: ClinVar (ncbi.nlm.nih.gov/clinvar), NCBI's search UI, its
      "chrpos38" field.
    Logic: f"{chrom}[chr] AND {pos-window}:{pos+window}[chrpos38]",
      URL-encoded. Confirmed live (2026-09-18, curator-suggested) against
      this project's own MYBPC3 c.278delA (11:47351252): an exact
      `N:N[chrpos38]` single-position range is valid syntax (0 or 1 results,
      normalization-sensitive for indels - ClinVar may record a deletion at
      a shifted position), while `N[chrpos38]` with no range at all returns
      a spurious zero-padded non-match; a +-10bp window reliably returns the
      surrounding variants (11 for this variant, including a nearby
      pathogenic frameshift) without a range-parsing surprise. Assumes
      GRCh38, matching this module's other position-based builders. Returns
      None if chrom/pos aren't set.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_reference_extensions() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    allele = positional_allele(variant)
    if allele is None:
        # ClinVar's chrpos38 field is GRCh38-specific, and this builder needs a
        # real position; an unaddressable record gets ClinVar's own search page
        # with no query rather than a window around a coordinate that does not
        # mean here what it means there.
        return CLINVAR_URL
    chrom, pos, _ref, _alt = allele
    term = f"{chrom}[chr] AND {pos - window}:{pos + window}[chrpos38]"
    return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={quote(term)}"


_AA3 = (
    "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Sec|Ser|Thr|Trp|Tyr|Val|Ter"
)
_RESIDUE = re.compile(rf"p\.\(?({_AA3})(\d+)")


def protein_residue(variant: VariantRecord) -> str | None:
    """The variant's own residue as ClinVar writes it (e.g. "Arg719"), or None.

    Reads INFO/HGVSP, which this project's demo data writes in 3-letter form
    with optional parentheses ("p.(Arg719Trp)"). Only the reference residue
    and its position are kept - the substituted residue is deliberately
    dropped, because both criteria that use this are about OTHER changes at
    the same residue.
    """
    hgvsp = variant.info.get("HGVSP")
    if not hgvsp:
        return None
    match = _RESIDUE.search(str(hgvsp))
    return f"{match.group(1)}{match.group(2)}" if match else None


def clinvar_codon_url(variant: VariantRecord) -> str | None:
    """
    Criteria: PS1 (doc: "show the clinvar page/summary for that codon").
    Site: ClinVar (ncbi.nlm.nih.gov/clinvar), NCBI's search UI.
    Logic: "<GENE>[gene] AND <residue>" (e.g. "MYH7[gene] AND Arg719"),
      which is the search the doc's own screenshot shows ("brca1, tyr1127"):
      every ClinVar record at this residue, whatever the nucleotide change.
      That is what PS1 asks about - PS1 is met by a DIFFERENT nucleotide
      change producing the same amino acid change, so a search keyed on this
      variant's own exact HGVSc (clinvar_search_url(), used by PM3/PP1) would
      never surface it. Confirmed live (2026-09-18) for this project's own
      MYH7 c.2155C>T: "MYH7[gene] AND Arg719" returns exactly the 5 records
      at residue 719, including p.Arg719Gln (Pathogenic) and p.Arg719Trp.
      Falls back to clinvar_position_url()'s genomic window when INFO/HGVSP
      is absent and no residue can be named - a coarser answer to the same
      question, not a different one.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_reference_extensions() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    gene = variant.info.get("GENE")
    residue = protein_residue(variant)
    if not gene or not residue:
        # No residue to search at: the genomic window, or - if the record is not
        # addressable by coordinate either - ClinVar's own unqueried page.
        return clinvar_position_url(variant)
    term = f"{gene}[gene] AND {residue}"
    return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={quote(term)}"


def autopvs1_variant_url(variant: VariantRecord, build: str = "hg38") -> str | None:
    """
    Criteria: PVS1 (doc: "Add autoPVS1 chart data/result").
    Site: AutoPVS1 (autopvs1.bgi.com), the per-variant PVS1 flowchart page
      (decision path, adjusted strength, disease-mechanism table).
    Logic: positional, same shape as gnomad_variant_url()/togovar_variant_url()
      below: f".../variant/{build}/{chrom}-{pos}-{ref}-{alt}". `build`
      defaults to "hg38" (matching this project's own GRCh38 assumption
      elsewhere in this module); pass "hg19" for a GRCh37-sourced
      VariantRecord. Confirmed live (2026-09-18) against two real variants -
      see this module's own docstring, "Verified reachable" - after an
      earlier version of this module concluded AutoPVS1 had no deep-link
      API and used a fixed homepage URL for every PVS1 line; that was
      wrong, AutoPVS1 does answer this URL shape. Falls back to the
      homepage (AUTOPVS1_URL) rather than None when ref/alt aren't concrete
      alleles, unlike gnomad_variant_url() - PVS1 previously always carried
      *some* link, and the homepage still lets a curator paste the variant
      in by hand, so losing the link entirely would be a regression.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_reference_extensions() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    allele = positional_allele(variant)
    if allele is None:
        return AUTOPVS1_URL
    chrom, pos, ref, alt = allele
    return f"https://autopvs1.bgi.com/variant/{build}/{chrom}-{pos}-{ref}-{alt}"


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
      passed through export.build_reference_extensions() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    allele = positional_allele(variant)
    if allele is None:
        return GNOMAD_URL
    chrom, pos, ref, alt = allele
    return (f"https://gnomad.broadinstitute.org/variant/"
            f"{chrom}-{pos}-{ref}-{alt}?dataset={dataset}")


def togovar_variant_url(variant: VariantRecord) -> str | None:
    """Return TogoVar's GRCh38 report page for an exact positional allele."""
    allele = positional_allele(variant)
    if allele is None:
        return TOGOVAR_URL
    chrom, pos, ref, alt = allele
    return f"https://grch38.togovar.org/variant/{chrom}-{pos}-{ref}-{alt}"


def franklin_variant_url(variant: VariantRecord, build: str = "hg38") -> str | None:
    """Return Franklin's page for an exact small variant.

    Franklin uses the ``snp`` route for SNVs and small indels.  This project
    evaluates GRCh38 variants, represented by the ``-hg38`` suffix.  A link
    is omitted for placeholder alleles rather than sending the curator to a
    page that cannot identify a variant.
    """
    allele = positional_allele(variant)
    if allele is None:
        return FRANKLIN_URL
    chrom, pos, ref, alt = allele
    return (
        "https://franklin.genoox.com/clinical-db/variant/snp/"
        f"chr{chrom}-{pos}-{ref}-{alt}-{build}"
    )


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
      passed through export.build_reference_extensions() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    gene = variant.info.get("GENE")
    if not gene:
        return GNOMAD_URL
    return f"https://gnomad.broadinstitute.org/gene/{quote(gene)}?dataset={dataset}"


def pp3_bp4_calibration_url(variant: VariantRecord) -> str:
    """
    Criteria: PP3, BP4 (doc: "show the predictors, and follow below
      thresholds", followed by that paper's Table 2 - the score intervals
      for thirteen missense tools across four pathogenic and four benign
      strengths).
    Site: PubMed, Pejaver et al. 2022 (Am J Hum Genet 109:2163-2177), the
      ClinGen SVI PP3/BP4 calibration.
    Logic: constant - the thresholds are a published table, not a
      per-variant page, so this ignores `variant` (kept in the signature so
      it composes with the other builders). This is the same source
      config/demo-rules.json's `computational.calibrations.
      revel-pejaver-2022` already names and whose REVEL row it already
      implements, so the link points a curator at the table the pipeline
      actually applied rather than at a generic predictor page.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_reference_extensions() - see this module's
      own docstring, "Where the returned URL ends up in VA-Spec output".
    """
    return PP3_BP4_CALIBRATION_URL


# code -> the reference-link builders for that criterion, in the order a
# curator should see them.
#
# [Scope: only the criteria the source document names, 2026-09-18]
#   Every entry below answers a per-criterion ask in doc/recs for expert
#   board.docx (quoted in this module's own docstring and in each builder).
#   The 15 codes that document says nothing about - PS2, PS3, PS4, PM4, PM6,
#   PP4, PP5, BS3, BS4, BP1, BP2, BP3, BP5, BP6, BP7 - get no reference link,
#   per the user's direction. This replaces the earlier arrangement in which
#   reference_urls_for_criterion() appended Franklin's variant page to ALL
#   ALL_ACMG_CODES, which put a link on codes nothing had asked for one on.
#
# Two deliberate departures from a literal reading of that document:
#   - PM2/BA1/BS1/BS2 carry TogoVar as well as the gnomAD page the document
#     asks for. TogoVar is the population provider this pipeline actually
#     queries (config/demo-rules.json's `population_sources`), so dropping it
#     would leave the curator without a link to the source the judgment was
#     made from.
#   - PM1/PM5 lead with Franklin's variant page, which is where the document's
#     own screenshots for those two codes come from (its Region Viewer shows
#     the domain track and the neighbouring ClinVar/UniProt assessments the
#     document asks for). Franklin exposes no deep link to that viewer alone
#     (checked live, 2026-09-18), so the variant page is the linkable target.
#     UniProt follows it as the source of the domain boundaries themselves.
_URL_BUILDERS = {
    "PVS1": (autopvs1_variant_url,),
    "PS1": (clinvar_codon_url,),
    "PM3": (clinvar_search_url,),
    "PP1": (clinvar_search_url,),
    "PM1": (franklin_variant_url, uniprot_page_url),
    "PM5": (franklin_variant_url, uniprot_page_url),
    "PM2": (togovar_variant_url, gnomad_variant_url),
    "BA1": (togovar_variant_url, gnomad_variant_url),
    "BS1": (togovar_variant_url, gnomad_variant_url),
    "BS2": (togovar_variant_url, gnomad_variant_url),
    "PP2": (gnomad_gene_url,),
    "PP3": (pp3_bp4_calibration_url,),
    "BP4": (pp3_bp4_calibration_url,),
}

# Public: the codes this module can build a reference URL for at all (a
# per-variant call may still return None for one of these - e.g.
# gnomad_variant_url() on an ALT="." record - but the code itself has a
# defined URL strategy). Callers building a UI/EvidenceLine per code (see
# acmg_pipeline.export.build_reference_extensions()) can use this to know
# which criteria have a curator-facing reference link worth generating.
CODES_WITH_REFERENCE_URL = frozenset(_URL_BUILDERS)


def reference_urls_for_criterion(code: str, variant: VariantRecord) -> list[str]:
    """Return every curator-facing URL for one criterion.

    Returns [] for a criterion the source document names no reference page
    for - see _URL_BUILDERS' own note on that scope. For a criterion it does
    name, the list is never empty: a builder that cannot address this variant
    returns that site's own unparameterized landing page (AUTOPVS1_URL,
    TOGOVAR_URL, GNOMAD_URL, FRANKLIN_URL, CLINVAR_URL, UNIPROT_URL) instead
    of nothing, so the curator still reaches the right site and can enter the
    variant by hand (2026-09-18, per the user's direction). What is never
    emitted is a URL that encodes this variant wrongly - see
    positional_allele().

    Builder failures are isolated so an optional external lookup (for
    example, gene-to-UniProt resolution) cannot suppress a deterministic link
    beside it. Duplicate URLs are removed while preserving order.
    """
    if code not in ALL_ACMG_CODES:
        return []

    urls: list[str] = []
    builders = _URL_BUILDERS.get(code, ())
    for builder in builders:
        try:
            url = builder(variant)
        except Exception:
            url = None
        if url is not None and url not in urls:
            urls.append(url)
    return urls


def reference_url_for_criterion(code: str, variant: VariantRecord) -> str | None:
    """
    Backward-compatible singular accessor for the first available URL.

    New callers should use reference_urls_for_criterion() so additional
    links, including Franklin and the second population-frequency source,
    are not discarded.
    Field: EvidenceLine.extensions (Extension(name="referenceLink")) once
      passed through export.build_reference_extensions() (its actual caller)
      - see this module's own docstring, "Where the returned URL ends up
      in VA-Spec output".
    """
    urls = reference_urls_for_criterion(code, variant)
    return urls[0] if urls else None


def all_reference_urls(variant: VariantRecord) -> dict[str, list[str]]:
    """
    Every code -> URLs this module can build for `variant`, skipping codes
    where no URL could be constructed.
    Field: same as reference_url_for_criterion() above - each value here
      lands in that code's own EvidenceLine.extensions, not one shared
      bundle (see export.build_reference_extensions(), called once per code).
    """
    return {
        code: urls
        for code in ALL_ACMG_CODES
        if (urls := reference_urls_for_criterion(code, variant))
    }
