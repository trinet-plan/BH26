"""
test_data/full_criteria_ground_truth.py

A ground-truth test dataset spanning ALL 28 ACMG/AMP 2015 codes (not just the
5 this project implements a judgment engine for). This is DATA ONLY - no
judgment logic lives here. It exists for two purposes:

  1. A validation fixture for whoever implements the 16 Layer-1
     (automated/rule-based) codes and PP4 (see acmg_pipeline/criteria/
     stubs.py's AUTOMATED_RULE_BASED_OTHER_TEAM and design doc section
     15-13) - real (gene, hgvsc, criterion) -> met/not_met/strength triples
     to check a new implementation against.
  2. Test data for acmg_pipeline.classification.classify() itself, to see
     how the full 28-code picture compares against each variant's real,
     already-known overall classification (see test_full_criteria_ground_
     truth.py) - useful even though this project only ever supplies 5 of
     the 28 codes' worth of real evidence today.

[Data sources]
  "erepo": Live re-query (2026-09-15) of ClinGen ERepo (acmg_pipeline.gate.
    ERepoClient), unfiltered - i.e. every evidence code ERepo has a status
    for, not just PS3/BS3/PS4/PP1/BS4. Two groups of variants were queried:
      - The 22 unique variants already used in acmg_pipeline/pipeline.py's
        test_cases (re-queried here for their FULL evidenceCodeStatus,
        since the original query there only looked up PMIDs for a single
        criterion per variant).
      - 4 additional variants (HNF1A x3, ATM x1) found by reverse-searching
        ERepo for real PP4 usage across gene-only queries (GET .../
        classifications?gene=<X>, no hgvs - confirmed to return every
        curated variant for that gene) over a candidate gene list. This
        was in direct response to the user's suggestion (2026-09-15) to
        source PP4 examples "backwards" from ClinGen rather than relying
        only on the demo-case document's single non-ERepo PP4 example.
        HNF1A turned out to be the single richest source of real PP4 usage
        (24 hits across many variants - it uses the MODY Probability
        Calculator-driven PP4 specification the ClinGen SVI PP1/PP4
        guidance paper's ACGS follow-up mentions by name).
    The raw snapshot is test_data/erepo_full_requery_2026-09-15.json,
    produced by the query in this docstring's [Regenerating] section below.
    ERepo's own evidence code labels only ever carry an explicit strength
    suffix (e.g. "PP4_Moderate") when a VCEP moved it off that code's
    default tier; a bare label (e.g. "PP4") uses the code's own default
    strength per Richards et al. 2015 Table 3/4 (VeryStrong for PVS1,
    StandAlone for BA1, Strong for PS*/BS*, Moderate for PM*, Supporting
    for PP*/BP*) - see _default_strength() in the regeneration script.
    tier="A" throughout (an ERepo-curated VCEP classification).

  "democase": Hand-transcribed from democase/real_cases_groundtruth_
    integrated_v6_ja.md (already-researched, human-curated ground truth for
    the 4 BH26 demo case variants + a few pedagogically important companion
    variants from the same cases - see inline notes). tier follows that
    document's own A/B/C confidence hierarchy (section 1 there), NOT
    ERepo's - most of these variants are not ClinGen VCEP-reviewed at all,
    which is exactly why the document tracks provenance/confidence
    separately per entry rather than treating every ground truth value as
    equally certain.

[Regenerating the "erepo" entries]
  python3 - <<'EOF'
  import sys, json
  sys.path.insert(0, ".")
  from acmg_pipeline.gate import ERepoClient
  from acmg_pipeline import classification as clsf

  VARIANTS = [...]  # see full_criteria_ground_truth.py's git history / design
                     # doc section 15-14 for the exact list used 2026-09-15
  client = ERepoClient()
  # ... (parse_label / default_strength as in this module) ...
  EOF
  Re-running this will pick up any ERepo updates since 2026-09-15; the
  checked-in JSON is a snapshot, not a live query (matching how pipeline.
  py's own ground_truth strings are dated literals, not live-queried at
  import time - see design doc section 8-9's ground-truth audit history for
  why static, dated snapshots are preferred over re-deriving at test time).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from acmg_pipeline.classification import ALL_ACMG_CODES, Strength
from acmg_pipeline.gate import CriterionStatus

_DATA_DIR = Path(__file__).resolve().parent


@dataclass(frozen=True)
class GroundTruthEntry:
    gene: str
    hgvsc: str
    hgvsp: str
    criterion: str                    # base ACMG code, e.g. "PP4" (never carries a strength suffix)
    status: CriterionStatus
    strength: Optional[Strength]      # required iff status == MET (mirrors classification.CriterionEvidence)
    source: str                       # "erepo" | "democase"
    tier: str                         # "A" | "B" | "C" (democase's confidence hierarchy; erepo entries are always "A")
    variant_outcome: Optional[str] = None  # this variant's overall classification, for classify()-comparison use
    notes: str = ""

    def __post_init__(self):
        if self.criterion not in ALL_ACMG_CODES:
            raise ValueError(f"{self.criterion!r} is not a recognized ACMG/AMP 2015 code")
        if self.status == CriterionStatus.MET and self.strength is None:
            raise ValueError(f"{self.gene} {self.hgvsc} {self.criterion}: status=MET requires a strength")
        if self.source not in ("erepo", "democase"):
            raise ValueError(f"unknown source {self.source!r}")
        if self.tier not in ("A", "B", "C"):
            raise ValueError(f"unknown tier {self.tier!r}")


def _load_erepo_entries() -> list[GroundTruthEntry]:
    raw = json.loads((_DATA_DIR / "erepo_full_requery_2026-09-15.json").read_text())
    return [
        GroundTruthEntry(
            gene=r["gene"], hgvsc=r["hgvsc"], hgvsp=r["hgvsp"], criterion=r["criterion"],
            status=CriterionStatus(r["status"]),
            strength=Strength(r["strength"].lower()) if r["strength"] else None,
            source=r["source"], tier=r["tier"], variant_outcome=r["variant_outcome"],
        )
        for r in raw
    ]


# ============================================================================
# democase entries (hand-transcribed from democase/real_cases_groundtruth_
# integrated_v6_ja.md - see that document's sections 3-6 for full narrative
# context per case). Only variants with an unambiguous c. HGVS in the
# document are included here (a handful of "noise" variants are described
# only by protein change / gene name in prose and are skipped rather than
# guessed at).
# ============================================================================

_DEMOCASE_ENTRIES = [
    # --- Case 1: MYBPC3 + KCNJ5 (phenotype-mismatch teaching pattern) ---
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.278delA", hgvsp="p.Lys93ArgfsTer3", criterion="PVS1",
        status=CriterionStatus.MET, strength=Strength.VERY_STRONG,
        source="democase", tier="B", variant_outcome="Likely Pathogenic",
        notes="Frameshift + premature stop, NMD predicted; MYBPC3 LOF is an established disease mechanism.",
    ),
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.278delA", hgvsp="p.Lys93ArgfsTer3", criterion="PM2",
        status=CriterionStatus.MET, strength=Strength.MODERATE,
        source="democase", tier="B", variant_outcome="Likely Pathogenic",
        notes="Absent from 1000G/ESP6500/ExAC and regional cohorts.",
    ),
    GroundTruthEntry(
        gene="KCNJ5", hgvsc="c.464G>A", hgvsp="p.Arg155Gln", criterion="PM2",
        status=CriterionStatus.MET, strength=Strength.MODERATE,
        source="democase", tier="B", variant_outcome="Uncertain Significance",
        notes="ExAC 8.236e-06 only; absent from 1000G/ESP6500.",
    ),
    GroundTruthEntry(
        gene="KCNJ5", hgvsc="c.464G>A", hgvsp="p.Arg155Gln", criterion="PP3",
        status=CriterionStatus.MET, strength=Strength.SUPPORTING,
        source="democase", tier="B", variant_outcome="Uncertain Significance",
        notes="SIFT/PolyPhen2/MutationTaster all damaging; REVEL=0.934.",
    ),
    GroundTruthEntry(
        gene="KCNJ5", hgvsc="c.464G>A", hgvsp="p.Arg155Gln", criterion="PP4",
        status=CriterionStatus.NOT_MET, strength=None,
        source="democase", tier="B", variant_outcome="Uncertain Significance",
        notes=(
            "Rejected: KCNJ5 is unrelated to the HCM phenotype under investigation "
            "(the son carries only the KCNJ5 variant and has no cardiac phenotype per "
            "the source paper) - the canonical phenotype-mismatch teaching example for "
            "why PP4 needs gene-phenotype specificity, not just a high in-silico score."
        ),
    ),
    # --- Case 2: MYBPC3 compound heterozygote (dual-cause pattern) ---
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.2905+1G>A", hgvsp="N/A", criterion="PVS1",
        status=CriterionStatus.MET, strength=Strength.VERY_STRONG,
        source="democase", tier="B", variant_outcome="Pathogenic",
        notes="Canonical splice-donor-site disruption; LOF established disease mechanism.",
    ),
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.2905+1G>A", hgvsp="N/A", criterion="PS4",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="B", variant_outcome="Pathogenic",
        notes="ClinVar VCV000042666 cites a case-control paper (PMID:27532257); 7 independent labs concordant.",
    ),
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.2905+1G>A", hgvsp="N/A", criterion="PM2",
        status=CriterionStatus.MET, strength=Strength.MODERATE,
        source="democase", tier="B", variant_outcome="Pathogenic",
        notes="gnomAD 0.0014%; 2 UK Biobank carriers; absent from 1000 Genomes.",
    ),
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.836del", hgvsp="p.Gly279ValfsTer21", criterion="PVS1",
        status=CriterionStatus.MET, strength=Strength.VERY_STRONG,
        source="democase", tier="B", variant_outcome="Likely Pathogenic",
        notes="Premature stop, NMD predicted; Western blot confirmed loss of protein (direct functional support).",
    ),
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.836del", hgvsp="p.Gly279ValfsTer21", criterion="PM2",
        status=CriterionStatus.MET, strength=Strength.MODERATE,
        source="democase", tier="B", variant_outcome="Likely Pathogenic",
        notes="Absent from gnomAD/UK Biobank/1000 Genomes/ClinVar/dbSNP.",
    ),
    # --- Case 3: MYH7 (ClinGen-reviewed) + MYBPC3 p.Glu334Lys (Pan-Asian
    # frequency pattern) + a real high-confidence Benign noise variant ---
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.2155C>T", hgvsp="p.Arg719Trp", criterion="PS2",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="A", variant_outcome="Pathogenic",
        notes="De novo confirmed by trio analysis (PMID:10957787). ClinGen Inherited Cardiomyopathy Expert Panel.",
    ),
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.2155C>T", hgvsp="p.Arg719Trp", criterion="PS3",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="A", variant_outcome="Pathogenic",
        notes=(
            "Deleterious cardiac remodeling from reduced hemodynamics (PMID:24829265). "
            "The human curator explicitly noted low self-confidence interpreting the "
            "knock-in/knock-out mouse model literature underlying this call."
        ),
    ),
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.2155C>T", hgvsp="p.Arg719Trp", criterion="PS4",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="A", variant_outcome="Pathogenic",
        notes="Found in 8 affected individuals across 6 unrelated Chinese families; absent in unaffected non-carriers (PMID:19645038).",
    ),
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.2155C>T", hgvsp="p.Arg719Trp", criterion="PM1",
        status=CriterionStatus.MET, strength=Strength.MODERATE,
        source="democase", tier="A", variant_outcome="Pathogenic",
        notes="ClinGen Evidence Repository, formally applied.",
    ),
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.2155C>T", hgvsp="p.Arg719Trp", criterion="PM2",
        status=CriterionStatus.MET, strength=Strength.MODERATE,
        source="democase", tier="A", variant_outcome="Pathogenic",
        notes="ClinGen Evidence Repository, formally applied.",
    ),
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.2155C>T", hgvsp="p.Arg719Trp", criterion="PM5",
        status=CriterionStatus.MET, strength=Strength.MODERATE,
        source="democase", tier="A", variant_outcome="Pathogenic",
        notes="Same residue's other AA change (p.Arg719Gln, Variation ID 14107) is Expert-Panel Pathogenic.",
    ),
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.2155C>T", hgvsp="p.Arg719Trp", criterion="PP1",
        status=CriterionStatus.MET, strength=Strength.SUPPORTING,
        source="democase", tier="A", variant_outcome="Pathogenic",
        notes="ClinGen Evidence Repository, formally applied.",
    ),
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.2155C>T", hgvsp="p.Arg719Trp", criterion="PP3",
        status=CriterionStatus.MET, strength=Strength.SUPPORTING,
        source="democase", tier="A", variant_outcome="Pathogenic",
        notes="ClinGen Evidence Repository, formally applied.",
    ),
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.1000G>A", hgvsp="p.Glu334Lys", criterion="PS3",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="C", variant_outcome="Uncertain Significance",
        notes=(
            "E334K cMyBPC alters UPS function / ion-channel & Ca2+-handling protein levels. "
            "Same low-confidence-in-mouse-model caveat as MYH7's PS3 above; directly "
            "conflicts with the BS1 call below, which is why the overall call is VUS."
        ),
    ),
    GroundTruthEntry(
        gene="MYBPC3", hgvsc="c.1000G>A", hgvsp="p.Glu334Lys", criterion="BS1",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="C", variant_outcome="Uncertain Significance",
        notes=(
            "gnomAD overall AF=0.02368%, East Asian AF=0.3338% (~14x). Only becomes "
            "apparent once regional (not just global) allele frequency is used - global-"
            "only data had this variant called Pathogenic at one point."
        ),
    ),
    GroundTruthEntry(
        gene="MYH7", hgvsc="c.3382G>A", hgvsp="p.Ala1128Thr", criterion="BS1",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="A", variant_outcome="Benign",
        notes=(
            "Noise variant (deliberately included to test correct de-prioritization). "
            "ClinGen Inherited Cardiomyopathy Expert Panel Evidence Repository (Kelly et al. "
            "2018): gnomAD Latino 0.02%, East Asian 0.016% - high in both populations, not "
            "a region-specific-frequency pattern like MYBPC3 c.1000G>A above."
        ),
    ),
    # --- Case 4: DSG2 (Japan/ToMMo, Pan-Asian founder-frequency pattern) ---
    GroundTruthEntry(
        gene="DSG2", hgvsc="c.1592T>G", hgvsp="p.Phe531Cys", criterion="PS3",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="C", variant_outcome=None,
        notes=(
            "Curator's own re-derivation (not a formal VCEP call). Mouse Dsg2 p.Phe536Cys "
            "knock-in model (homozygous): biventricular dilation, reduced LVEF, conduction "
            "abnormalities, fibrosis. Same low-confidence-in-mouse-model caveat as Case 3. "
            "variant_outcome intentionally left None: the curator's re-derivation (PS3+PS4 "
            "-> Pathogenic) does not match the real ClinVar submitter majority (VCV000044283.23: "
            "3 Likely pathogenic / 7 Uncertain significance / 1 Likely benign) - see design doc "
            "section 15-13/15-14 and democase doc section 6 for the full discussion of this "
            "double-uncertainty case."
        ),
    ),
    GroundTruthEntry(
        gene="DSG2", hgvsc="c.1592T>G", hgvsp="p.Phe531Cys", criterion="PS4",
        status=CriterionStatus.MET, strength=Strength.STRONG,
        source="democase", tier="C", variant_outcome=None,
        notes="Family screening: only homozygous carriers show full-penetrance ARVC; heterozygous carriers unaffected.",
    ),
]


GROUND_TRUTH: list[GroundTruthEntry] = _load_erepo_entries() + _DEMOCASE_ENTRIES


def entries_for(gene: str, hgvsc: str) -> list[GroundTruthEntry]:
    return [e for e in GROUND_TRUTH if e.gene == gene and e.hgvsc == hgvsc]


def unique_variants() -> list[tuple[str, str]]:
    seen: list[tuple[str, str]] = []
    for e in GROUND_TRUTH:
        key = (e.gene, e.hgvsc)
        if key not in seen:
            seen.append(key)
    return seen
