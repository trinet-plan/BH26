"""
ps3_bs3_ps4_gate.py
Entry gate for the PS3 / BS3 / PS4 pipeline

[Precondition] The caller has already completed the PVS1 judgment; only cases
where PVS1 does not apply (i.e., where this module has a role to play) reach
this point. The PVS1 judgment logic itself is not included here.

[Responsibilities]
  1. Variant matching gate (shared by PS3/BS3/PS4)
     Checks whether the target variant has already been curated by a VCEP in
     ClinGen ERepo. The gate itself (variant-notation matching + ERepo lookup
     + the production/validation dual mode) is built to accept a generic
     criterion string ("PS3"/"BS3"/"PS4", any of them), so it works for PS4
     unmodified as well (verified with real data: RUNX1 c.601C>T, PS4:Met).
       - production mode: if already curated, return that Met/Not Met as-is
         and do not invoke the downstream stage (PS3/BS3/PS4 direction
         judgment from the literature).
       - validation mode: return only the identification result with the
         ground truth (Met/Not Met) hidden, forcing the downstream stage to
         run regardless. For accuracy measurement (pipeline evaluation) only.
  2. Approved-assay cross-check (PS3/BS3 only, not applicable to PS4)
     Cross-checks the experimental-method description extracted downstream
     (from the literature) against each VCEP's approved-assay list, and
     returns notes on applicability / BS3 inapplicability, etc. The RASopathy
     VCEP table obtained in design doc section 11 (Wilcox et al. 2025, Genet
     Med Open, Supplementary Table 2) is loaded as the initial data.
     [Note] PS4 is based on case-control statistics / segregation analysis,
     not a functional assay, so it is out of scope for this cross-check.
     PS4-specific evidence evaluation (e.g., aggregating case counts across
     multiple papers) is handled separately in design doc section 4 (the
     former PS4-only pipeline).

[Out of scope (the caller's or a higher-level module's responsibility)]
  - The PVS1 judgment itself
  - Generating explanatory text for human curators (handled on the prompt-design side)
  - The literature-reading direction judgment itself (downstream of this gate)
  - PS4-specific evidence aggregation logic (stacking case counts across
    multiple PMIDs, design doc section 4)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

import requests


# ============================================================================
# 1. Variant-notation normalization and matching (organized/generalized from
#    the implementation in design doc section 6-2)
# ============================================================================

_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}


def _protein_equivalents(variant_str: str) -> set[str]:
    """Generate all notational variants of p.Xxx###Yyy (1-letter/3-letter, with/without parentheses)."""
    eqs: set[str] = set()
    for m in re.finditer(r"p\.\(?([A-Za-z]{1,3})(\d+)([A-Za-z\*]{1,3})\)?", variant_str):
        aa1, pos, aa2 = m.groups()
        eqs.add(f"{aa1}{pos}{aa2}")
        eqs.add(f"p.({aa1}{pos}{aa2})")
        a1, a2 = _THREE_TO_ONE.get(aa1), _THREE_TO_ONE.get(aa2)
        if a1 and a2:
            eqs.add(f"{a1}{pos}{a2}")
    return eqs


def _nucleotide_equivalents(variant_str: str) -> set[str]:
    """Absorb whitespace variation in c./m. notation (e.g., c.1012G > A)."""
    eqs: set[str] = set()
    mc = re.search(r"(c\.[0-9_+\-]+(?:[ACGT]>[ACGT]|delins[ACGT]*|del|ins[ACGT]*|dup))", variant_str)
    if mc:
        core = mc.group(1)
        eqs.add(core)
        eqs.add(re.sub(r"([ACGT])>([ACGT])", r"\1 > \2", core))
    mm = re.search(r"(m\.[0-9_]+(?:[ACGT>a-z]+|del|dup))", variant_str)
    if mm:
        eqs.add(mm.group(1))
    return eqs


def _indel_position_range(variant_str: str) -> Optional[tuple[int, int]]:
    """Extract the coordinate range for del/dup/ins-type variants (used for +/- tolerance matching)."""
    m = re.search(r"[cm]\.(-?\d+)(?:[+\-]\d+)?_(\d+)(?:[+\-]\d+)?(?:del|dup|ins)", variant_str)
    if m:
        return int(m.group(1)), int(m.group(2))
    m2 = re.search(r"[cm]\.(-?\d+)(?:[+\-]\d+)?(?:del|dup|ins)", variant_str)
    if m2:
        p = int(m2.group(1))
        return p, p
    return None


class MatchStatus(Enum):
    MATCHED = "matched"          # exact match after absorbing notational variation
    HEURISTIC = "heuristic"      # match within coordinate tolerance (+/-2); needs manual review
    UNSUCCESSFUL = "unsuccessful"  # no match


@dataclass
class VariantMatchResult:
    status: MatchStatus
    matched_text: Optional[str] = None
    tolerance_used: Optional[int] = None  # when heuristic, the actual offset in bases


def match_variant_in_text(variant_str: str, text: str, coordinate_tolerance: int = 2) -> VariantMatchResult:
    """
    Checks whether the target variant (HGVS notation) is present in the paper
    text / ERepo record. Tries, in order: exact match (with notational
    variation absorbed) -> heuristic (coordinate +/- tolerance) -> unsuccessful.

    Finding from design doc section 6-2: exact matching alone misses cases
    due to conventional differences in boundary notation (e.g.,
    m.14512_14513del vs m.14513_14514del). Coordinate-tolerance matching is a
    heuristic meant to catch these; treat it as a flag requiring manual
    review, not a mechanical confirmation.
    """
    equivalents = _protein_equivalents(variant_str) | _nucleotide_equivalents(variant_str)
    for eq in equivalents:
        if eq and eq in text:
            return VariantMatchResult(MatchStatus.MATCHED, matched_text=eq)

    rng = _indel_position_range(variant_str)
    if rng:
        lo, hi = rng
        for cand_lo in range(lo - coordinate_tolerance, lo + coordinate_tolerance + 1):
            for cand_hi in range(hi - coordinate_tolerance, hi + coordinate_tolerance + 1):
                for pat in (f"{cand_lo}_{cand_hi}del", f"{cand_lo}_{cand_hi}dup", f"{cand_lo}_{cand_hi}ins"):
                    if pat in text:
                        drift = max(abs(cand_lo - lo), abs(cand_hi - hi))
                        return VariantMatchResult(MatchStatus.HEURISTIC, matched_text=pat, tolerance_used=drift)

    return VariantMatchResult(MatchStatus.UNSUCCESSFUL)


# ============================================================================
# 2. ClinGen ERepo lookup (the variant-matching gate itself)
# ============================================================================

class CriterionStatus(Enum):
    MET = "met"
    NOT_MET = "not_met"
    UNKNOWN = "unknown"  # the criterion code itself does not exist in ERepo for this variant


@dataclass
class ERepoLookupResult:
    found_in_erepo: bool
    # only meaningful when found_in_erepo=True
    caid: Optional[str] = None
    evidence_code_status: dict[str, CriterionStatus] = field(default_factory=dict)
    outcome: Optional[str] = None  # Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign
    # PMIDs cited as evidence for this variant classification (from ERepo's
    # evidenceLinks field). NOTE: this is the full set of papers cited for
    # the variant's overall classification, not necessarily filtered to
    # papers specific to PS3/BS3 functional evidence alone - a real VCEP
    # classification often cites multiple papers per variant (e.g., RUNX1
    # c.601C>T cited 4 PMIDs in a real ERepo response seen during design
    # doc section 9). Downstream, each cited paper still needs its own
    # variant-identification check (match_status), since not every cited
    # paper necessarily contains functional-assay data relevant to PS3/BS3
    # specifically (some may support other criteria, e.g. PM2 population
    # frequency, PS4 case counts, etc.).
    evidence_pmids: list[str] = field(default_factory=list)
    raw: Optional[dict] = None     # raw API response (kept for audit purposes)


def _extract_pmids_from_evidence_links(evidence_links: list[dict]) -> list[str]:
    """
    ERepo's evidenceLinks field is a list of {"@id": "https://www.ncbi.nlm.nih.gov/pubmed/<PMID>"}
    entries. Extracts the bare PMID string from each URL.
    """
    pmids = []
    for link in evidence_links or []:
        url = link.get("@id", "")
        m = re.search(r"pubmed/(\d+)", url)
        if m:
            pmids.append(m.group(1))
    return pmids


class ERepoClient:
    """
    Handles queries to the ClinGen ERepo API (erepo.clinicalgenome.org).

    Confirmed reachable and working from this environment (2026-09-15),
    superseding the design doc section 9-4 note that this domain was not on
    the sandbox's network allowlist:
      GET https://erepo.clinicalgenome.org/evrepo/api/classifications
          ?gene=RUNX1&hgvs=c.601C>T
      -> HTTP 200, {"variantInterpretations": [...]} with 4 evidenceLinks,
         matching the shape OfflineERepoClient already assumes (gene.label,
         hgvs[], evidenceLinks[], guidelines[].agents[].evidenceCodes[]).
    The actual matching/parsing logic is delegated to OfflineERepoClient
    rather than duplicated here.
    """

    BASE_URL = "https://erepo.clinicalgenome.org"
    TIMEOUT_SEC = 15

    def lookup(self, gene: str, hgvs: str) -> ERepoLookupResult:
        try:
            resp = requests.get(
                f"{self.BASE_URL}/evrepo/api/classifications",
                params={"gene": gene, "hgvs": hgvs},
                timeout=self.TIMEOUT_SEC,
            )
            resp.raise_for_status()
            data = resp.json()
        except requests.RequestException as e:
            raise RuntimeError(
                f"ClinGen ERepo API request failed (gene={gene}, hgvs={hgvs}): {e}"
            ) from e

        return OfflineERepoClient([data]).lookup(gene, hgvs)


class OfflineERepoClient(ERepoClient):
    """
    For testing/backtesting. Uses a list of pre-fetched ERepo responses
    (a list of dicts, in the exact JSON format actually returned by the
    ClinGen API) as the lookup source. Designed around the 36 responses
    obtained during the full verification pass in design doc section 9.
    """

    def __init__(self, cached_responses: list[dict]):
        # build a cache mapping (gene, variant-matching text) -> ERepoLookupResult
        self._entries: list[dict] = []
        for resp in cached_responses:
            for interp in resp.get("variantInterpretations", []):
                self._entries.append(interp)

    def lookup(self, gene: str, hgvs: str) -> ERepoLookupResult:
        hgvs_variants = _protein_equivalents(hgvs) | _nucleotide_equivalents(hgvs) | {hgvs}
        for interp in self._entries:
            interp_gene = (interp.get("gene") or {}).get("label", "")
            if interp_gene.upper() != gene.upper():
                continue
            hgvs_list = interp.get("hgvs", [])
            if not any(any(v in h for v in hgvs_variants if v) for h in hgvs_list):
                continue

            status_map: dict[str, CriterionStatus] = {}
            outcome = None
            for guideline in interp.get("guidelines", []):
                outcome = outcome or (guideline.get("outcome") or {}).get("label")
                for agent in guideline.get("agents", []):
                    for code in agent.get("evidenceCodes", []):
                        label = code.get("label")
                        status = code.get("status", "")
                        status_map[label] = (
                            CriterionStatus.MET if status == "Met" else CriterionStatus.NOT_MET
                        )
            return ERepoLookupResult(
                found_in_erepo=True,
                caid=interp.get("caid"),
                evidence_code_status=status_map,
                outcome=outcome,
                evidence_pmids=_extract_pmids_from_evidence_links(interp.get("evidenceLinks", [])),
                raw=interp,
            )
        return ERepoLookupResult(found_in_erepo=False)


class CsvBackedERepoClient(ERepoClient):
    """
    An implementation that sources lookups directly from either an ERepo bulk
    CSV export (e.g., .../evrepo/api/classifications/all?format=tabbed) or a
    CGBench-style per-variant CSV. Used for testing/regression checks in
    development environments where the live API is unavailable, and also
    stands on its own as a production option: "pre-ingest curated data in
    bulk instead of hitting the API on every call."

    Expected CSV columns: entry_index, variant (HGVS), hgnc_gene, assertion,
    evidence_code, met_status (met/not_met) - exactly the CGBench format used
    in design doc sections 6-9. When multiple rows share the same entry_index
    (i.e., the same variant) with different criterion codes, they are merged
    into a single ERepoLookupResult.
    """

    def __init__(self, csv_path: str):
        import csv as _csv

        groups: dict[str, dict] = {}
        with open(csv_path, newline="", encoding="utf-8") as f:
            for row in _csv.DictReader(f):
                key = row["entry_index"]
                if key not in groups:
                    groups[key] = {
                        "gene": row["hgnc_gene"],
                        "variant": row["variant"],
                        "assertion": row.get("assertion", ""),
                        "codes": {},
                        "pmids": set(),
                    }
                status = CriterionStatus.MET if row["met_status"] == "met" else CriterionStatus.NOT_MET
                groups[key]["codes"][row["evidence_code"]] = status
                # CGBench stores exactly one 'pmid' per (entry_index,
                # evidence_code) row, formatted like "PubMed:12345678". This
                # is a real limitation relative to live ERepo: a variant's
                # true evidenceLinks list often has multiple PMIDs (see the
                # OfflineERepoClient docstring/tests for a real 4-PMID
                # example), but CGBench appears to have already collapsed
                # each evidence code down to one representative citation. As
                # a result, evidence_pmids from this client will typically
                # contain only 1 PMID even when a live ERepo lookup for the
                # same variant would return several.
                pmid_raw = row.get("pmid", "")
                m = re.search(r"(\d+)\s*$", pmid_raw)
                if m:
                    groups[key]["pmids"].add(m.group(1))
        self._groups = list(groups.values())

    def lookup(self, gene: str, hgvs: str) -> ERepoLookupResult:
        hgvs_variants = _protein_equivalents(hgvs) | _nucleotide_equivalents(hgvs) | {hgvs}

        # 1st pass: match only within groups whose gene name matches (the normal path)
        for g in self._groups:
            if g["gene"].upper() != gene.upper():
                continue
            if not any(v and v in g["variant"] for v in hgvs_variants):
                continue
            return ERepoLookupResult(
                found_in_erepo=True,
                evidence_code_status=g["codes"],
                outcome=g["assertion"],
                evidence_pmids=sorted(g["pmids"]),
                raw=g,
            )

        # 2nd pass (fallback): for rows where CGBench's gene column is blank,
        # match on variant notation alone regardless of gene name (a real
        # problem actually encountered in section 11 with the MT-TK/MT-ND6
        # cases). Only accepted when the gene column is empty AND the variant
        # notation matches. If multiple candidates match, do not accept any
        # of them (to avoid false positives) and return unresolved instead.
        candidates = [
            g for g in self._groups
            if not g["gene"].strip() and any(v and v in g["variant"] for v in hgvs_variants)
        ]
        if len(candidates) == 1:
            g = candidates[0]
            return ERepoLookupResult(
                found_in_erepo=True,
                evidence_code_status=g["codes"],
                outcome=g["assertion"],
                evidence_pmids=sorted(g["pmids"]),
                raw=g,
            )

        return ERepoLookupResult(found_in_erepo=False)


# ============================================================================
# 3. Approved-assay lists (structured data, per VCEP)
# ============================================================================

@dataclass
class ApprovedAssay:
    name: str
    description: str
    ps3_criterion: str  # free-text description of "what increase/decrease indicates the PS3 direction"
    applicable_genes: set[str]


@dataclass
class VcepAssaySpec:
    vcep_name: str
    approved_assays: list[ApprovedAssay]
    bs3_applicable: bool = True
    bs3_note: str = ""


# Design doc section 11: Wilcox et al. 2025, Genet Med Open, Supplementary Table 2
RASOPATHY_VCEP_SPEC = VcepAssaySpec(
    vcep_name="RASopathy VCEP",
    bs3_applicable=False,
    bs3_note=(
        "BS3 is not applicable for RASopathy. Known pathogenic variants "
        "exist that show wild-type-like results depending on physiological "
        "conditions, stimulation, or cell line (as stated in the source text)."
    ),
    approved_assays=[
        ApprovedAssay(
            name="RAS Activation Assay",
            description="Measurement of RAS bound and co-immunoprecipitated with RAF1/RBD",
            ps3_criterion="Increase in the RAS/RAF1 or RAS/RBD complex",
            applicable_genes={"MRAS", "HRAS", "KRAS", "SOS1", "SOS2", "NRAS", "LZTR1"},
        ),
        ApprovedAssay(
            name="MEK Activation Assay",
            description="Phosphorylated MEK / unphosphorylated MEK ratio (basal and after RTK stimulation)",
            ps3_criterion="Increase in phosphorylation",
            applicable_genes={
                "MRAS", "BRAF", "HRAS", "KRAS", "MAP2K1", "MAP2K2",
                "PTPN11", "RAF1", "SOS1", "SOS2", "NRAS", "LZTR1",
            },
        ),
        ApprovedAssay(
            name="ERK Activation Assay",
            description="Phosphorylated ERK / unphosphorylated ERK ratio (basal and after stimulation)",
            ps3_criterion="Increase in phosphorylation",
            applicable_genes={
                "MRAS", "BRAF", "HRAS", "KRAS", "MAP2K1", "MAP2K2",
                "PTPN11", "RAF1", "SOS1", "SOS2", "NRAS", "LZTR1",
                # Note (design doc section 11-2): idx 8336 (RIT1, ERK
                # activation assay) was confirmed via a live ClinGen ERepo
                # lookup to have PS3_Supporting: Met, but RIT1 is not listed
                # among this table's applicable genes (the discrepancy is
                # unresolved). Trusting the live ERepo data as authoritative,
                # RIT1 is provisionally included here as well.
                "RIT1",
            },
        ),
        ApprovedAssay(
            name="SHP-2 Phosphatase Activity",
            description="Phosphorylated / dephosphorylated SHP2 ratio",
            ps3_criterion="Increase in dephosphorylation",
            applicable_genes={"PTPN11"},
        ),
        ApprovedAssay(
            name="BRAF Kinase Activity",
            description="Kinase activity that phosphorylates MEK/ERK",
            ps3_criterion="Increase in phosphorylation",
            applicable_genes={"BRAF"},
        ),
        ApprovedAssay(
            name="RAF1 Kinase Activity",
            description="Kinase activity that phosphorylates MEK/ERK",
            ps3_criterion="Increase in phosphorylation",
            applicable_genes={"RAF1"},
        ),
        ApprovedAssay(
            name="LZTR1 Stability/Localization",
            description="Measurement of LZTR1 protein level and localization",
            ps3_criterion="Decreased protein level or abnormal localization",
            applicable_genes={"LZTR1"},
        ),
    ],
)

# Registry mapping VCEP name -> spec. Add more entries here as tables for
# other VCEPs become available.
VCEP_ASSAY_REGISTRY: dict[str, VcepAssaySpec] = {
    "RASopathy VCEP": RASOPATHY_VCEP_SPEC,
}


class AssayApplicability(Enum):
    APPROVED = "approved"              # matches an approved assay and applies to the target gene
    ASSAY_NOT_APPROVED = "assay_not_approved"  # this assay type is not on this VCEP's approved list
    GENE_NOT_LISTED = "gene_not_listed"        # the assay type is approved, but the target gene is not listed
    NO_SPEC_AVAILABLE = "no_spec_available"    # no list is registered for this VCEP at all
    BS3_NOT_APPLICABLE = "bs3_not_applicable"  # BS3 itself cannot be used for this VCEP


@dataclass
class AssayCheckResult:
    applicability: AssayApplicability
    matched_assay: Optional[ApprovedAssay] = None
    note: str = ""


# Common English function words / generic lab-jargon terms that would
# otherwise be extracted as "keywords" from an assay's description and cause
# false matches against unrelated assays that happen to share the same
# common word (e.g., "after", "ratio", "basal" appear in more than one
# assay's description above).
_ASSAY_KEYWORD_STOPWORDS = {
    "and", "the", "for", "with", "via", "after", "before", "ratio",
    "basal", "assay", "measurement", "level", "levels", "activity",
    "increase", "increased", "decrease", "decreased", "measuring",
    "bound", "localization", "stability",
}


def check_approved_assay(
    vcep_name: str,
    gene: str,
    assay_description: str,
    criterion: str = "PS3",
) -> AssayCheckResult:
    """
    Cross-checks a free-text description of the experimental method extracted
    downstream (from the literature) - e.g., "measured ERK1/2 phosphorylation
    by western blot" - against the VCEP's approved-assay list.

    Matching is a simple implementation: it first checks whether the assay's
    own name appears in the description, then falls back to partial keyword
    matches against the assay's description text. This is not strict NLP
    matching, so borderline cases require manual review (this function is
    only a first-pass screen).

    Note: a small stopword list (_ASSAY_KEYWORD_STOPWORDS) filters out common
    English function words extracted from the description text for the
    keyword-fallback pass. Without both the name-first pass and this filter,
    a generic word shared between two different assays' descriptions (e.g.,
    "stimulation", "phosphorylated") could cause a false match against the
    wrong assay. This was discovered in practice: after the approved-assay
    descriptions were translated to English, "ERK Activation Assay
    (... after serum stimulation)" was incorrectly matched against "MEK
    Activation Assay" because both descriptions happen to contain
    "stimulation" (and MEK is checked first in the assay list).
    """
    spec = VCEP_ASSAY_REGISTRY.get(vcep_name)
    if spec is None:
        return AssayCheckResult(
            AssayApplicability.NO_SPEC_AVAILABLE,
            note=f"No approved-assay list is registered for '{vcep_name}'. Needs to be added to VCEP_ASSAY_REGISTRY.",
        )

    if criterion.upper() == "BS3" and not spec.bs3_applicable:
        return AssayCheckResult(AssayApplicability.BS3_NOT_APPLICABLE, note=spec.bs3_note)

    desc_lower = assay_description.lower()

    # Pass 1: match on the assay's own name first. This is the most specific
    # signal and must take priority over description-keyword matching -
    # otherwise a generic word shared between two assays' descriptions (e.g.,
    # "stimulation", "phosphorylated") can cause a match against the wrong
    # assay before the correct one is ever considered (this happened in
    # practice: an "ERK Activation Assay" description matched "MEK
    # Activation Assay" first, since MEK is earlier in the list and both
    # descriptions mention "stimulation").
    for assay in spec.approved_assays:
        if assay.name.lower() in desc_lower:
            if gene.upper() in assay.applicable_genes:
                return AssayCheckResult(AssayApplicability.APPROVED, matched_assay=assay)
            return AssayCheckResult(
                AssayApplicability.GENE_NOT_LISTED,
                matched_assay=assay,
                note=f"{assay.name} is an assay approved by {vcep_name}, but {gene} is not in its applicable-gene list",
            )

    # Pass 2 (fallback): only if no assay name matched, fall back to
    # partial keyword matching against the description text.
    for assay in spec.approved_assays:
        keywords = [
            w for w in re.findall(r"[a-zA-Z0-9\-]{3,}", assay.description.lower())
            if w not in _ASSAY_KEYWORD_STOPWORDS
        ]
        if any(kw in desc_lower for kw in keywords):
            if gene.upper() in assay.applicable_genes:
                return AssayCheckResult(AssayApplicability.APPROVED, matched_assay=assay)
            return AssayCheckResult(
                AssayApplicability.GENE_NOT_LISTED,
                matched_assay=assay,
                note=f"{assay.name} is an assay approved by {vcep_name}, but {gene} is not in its applicable-gene list",
            )

    return AssayCheckResult(
        AssayApplicability.ASSAY_NOT_APPROVED,
        note=f"The described method did not match any assay approved by {vcep_name}",
    )


# ============================================================================
# 4. The gate itself (production / validation dual mode)
# ============================================================================

class GateMode(Enum):
    PRODUCTION = "production"
    VALIDATION = "validation"


@dataclass
class GateDecision:
    mode: GateMode
    variant_match: VariantMatchResult
    erepo_result: ERepoLookupResult
    # only set in production mode when found_in_erepo=True
    resolved_status: Optional[CriterionStatus] = None
    # in validation mode, the downstream pipeline is always run
    should_run_pipeline: bool = True
    # for audit/debug purposes only. Holds the ground truth only in
    # validation mode; never passed downstream.
    _hidden_ground_truth: Optional[CriterionStatus] = field(default=None, repr=False)


def run_gate(
    gene: str,
    hgvs: str,
    criterion: str,  # "PS3" or "BS3" (strength suffixes are handled separately by the caller)
    erepo_client: ERepoClient,
    mode: GateMode = GateMode.PRODUCTION,
) -> GateDecision:
    """
    The variant-matching gate itself.

    production: if ERepo already has Met/Not Met for this criterion, adopt it
                and do not invoke the downstream stage (literature reading)
                (should_run_pipeline=False).
    validation: checks the ERepo result but never passes it downstream (kept
                hidden), always running the pipeline. The caller should only
                use GateDecision._hidden_ground_truth for comparison AFTER
                the downstream judgment result comes back.
    """
    erepo_result = erepo_client.lookup(gene, hgvs)

    match_result = VariantMatchResult(MatchStatus.UNSUCCESSFUL)
    if erepo_result.found_in_erepo:
        match_result = VariantMatchResult(MatchStatus.MATCHED, matched_text=hgvs)

    if not erepo_result.found_in_erepo:
        return GateDecision(
            mode=mode,
            variant_match=match_result,
            erepo_result=erepo_result,
            should_run_pipeline=True,
        )

    status = erepo_result.evidence_code_status.get(criterion, CriterionStatus.UNKNOWN)

    if mode == GateMode.PRODUCTION:
        if status == CriterionStatus.UNKNOWN:
            # listed in ERepo, but this criterion code itself was not
            # evaluated (e.g., classification was already completed using
            # other strong criteria alone - the 3 PI3K-pathway cases in
            # design doc section 8)
            return GateDecision(
                mode=mode,
                variant_match=match_result,
                erepo_result=erepo_result,
                should_run_pipeline=True,
            )
        return GateDecision(
            mode=mode,
            variant_match=match_result,
            erepo_result=erepo_result,
            resolved_status=status,
            should_run_pipeline=False,
        )

    # validation mode: store the ground truth only in _hidden_ground_truth;
    # resolved_status stays None (not shown downstream).
    return GateDecision(
        mode=mode,
        variant_match=match_result,
        erepo_result=erepo_result,
        should_run_pipeline=True,
        _hidden_ground_truth=status,
    )
