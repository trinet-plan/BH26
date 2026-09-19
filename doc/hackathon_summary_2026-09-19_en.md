# ACMG/AMP Classification Pipeline — 6 Days of Development

Summary of what was built during the BH26 hackathon (6 days, with 2026-09-19 as the final
day). Builds on the [previous day's version](hackathon_summary_2026-09-18_en.md), folding
in the full 64-variant run, the 4-demo-case run, and the PVS1 mechanism-gate improvements
that followed from them.

The pipeline combines an LLM (gemma-4 served via vLLM) with the PubMed MCP server to
generate draft classifications across all 28 ACMG/AMP 2015 criteria, accelerating a human
curator's first-pass screening.

## 1. Implementation overview

All 28 ACMG/AMP 2015 criteria are implemented across four layers, grouped by how each
criterion is judged. Three systems are already running: deterministic rule/database-driven
criteria, criteria that require reading the literature (LLM + PubMed MCP), and criteria that
require a phenotype/segregation match. The rest are honestly reported as "not evaluated"
rather than silently omitted.

| Layer | Criteria | Count | Status |
|---|---|---:|---|
| Automated (Layer 1) | `PVS1 PS1 PM1 PM2 PM4 PM5 PP2 PP3 BA1 BS1 BP1 BP3 BP4 BP7` | 14 | Implemented |
| Literature engine | `PS3 BS3 PS4` | 3 | Implemented |
| Phenotype/segregation | `PP1 BS4 PP4` | 3 | Implemented |
| Disabled by policy | `PP5 BP6` | 2 | Always UNKNOWN |
| Not yet implemented | `PS2 PM3 PM6 BS2 BP2 BP5` | 6 | Explicitly NOT_EVALUATED |

The final classification is combined using Tavtigian et al. 2018's Bayesian point-based
framework (Supporting=1 / Moderate=2 / Strong=4 / Very Strong=8; thresholds: ≥10
Pathogenic, 6-9 Likely Pathogenic, 0-5 VUS, -1..-6 Likely Benign, ≤-7 Benign), and returned
from the API (FastAPI + Docker) as a GA4GH VA-Spec EvidenceLine.

## 2. VA-Spec output design work

The area the team invested the most effort in, since it is what a human curator actually
reads.

- **Removing duplication and inconsistency** — Fields bundled into a single
  `bh26AssessmentDetails` object were fully flattened into individually named top-level
  extensions, at the same level as `curatorHints`.
- **Fixing a direction inconsistency** — Benign-direction evidence (BS3) was always being
  reported with a pathogenic-direction outcome.
- **Feedback to GA4GH VA-Spec** — 5 schema-constraint issues were written up in
  [doc/ga4gh_va_spec_feedback_v1_en.md](ga4gh_va_spec_feedback_v1_en.md) and filed upstream
  as [ga4gh/va-spec#448](https://github.com/ga4gh/va-spec/issues/448). GA4GH's maintainer
  (korikuzma) confirmed issue #2 (a generic evidence-item type, `DataItem`) is already
  addressed in the upcoming `1.1.0-ballot.2026-09`; issues #1 (`Method.reportedIn` accepting
  multiple documents) and #3 (`Direction` representing "not evaluated") are not, and the
  `EvidenceLine` class itself is being merged into `Statement` in that same upcoming
  version — a larger structural change worth tracking.

## 3. Extending condition (diagnosis) resolution

PVS1's disease-mechanism gate and PP1/BS4/PP4 all need a patient's diagnosis resolved to a
MONDO disease class. Beyond the existing "read ERepo's own condition field" path, this was
extended with:

- **Retrying MONDO resolution** — `resolve_diagnosis_mondo()`'s final candidate pick is an
  LLM call, and empirically returns NOT_FOUND roughly 1 in 5 times for a borderline
  diagnosis (e.g. "hypertrophic cardiomyopathy with apical ventricular aneurysm") from pure
  sampling variance. Raised from 3 to 5 retries; verified 10/10 successes afterward.
- **Fixing the validation script's own wiring** — `run_integrated_validation_demo.py` had
  fallen out of sync with the real API path (`pipeline_interface.py`) and never called
  `resolve_diagnosis_mondo()` at all. Fixed - MYBPC3 c.278delA and others now correctly
  reach PVS1(very_strong).
- **Literature-based condition search (new)** — For variants ERepo has no record of at all
  (e.g. 3 demo-case-sourced MYBPC3 variants), a new module,
  `acmg_pipeline/condition_from_literature_search.py`, searches PubMed and asks the
  project's own LLM what disease the paper reports the variant causing, then feeds that
  through `resolve_diagnosis_mondo()`. The extraction prompt was tuned to ask for the full
  disease name rather than an abbreviation (e.g. "HCM"), since MONDO/OLS4 indexes full
  names - verified against MYBPC3 c.278delA.

## 4. Extending PVS1's mechanism gate — distinguishing "unknown" from "contradicted"

Investigating why PVS1's MET rate was 0% in the full 64-variant run led to two deliberate
policy changes.

- **Applying through a related (parent/child) disease** — When the case's stated diagnosis
  and the curated disease are only related in the MONDO hierarchy (not an exact match),
  PVS1 previously always stopped at MANUAL_REVIEW. Now, when every record for the related
  disease agrees the mechanism is established, PVS1 is applied and the imprecise disease
  match is disclosed via a curatorHint. Verified: RUNX1 c.601C>T went from UNKNOWN to
  MET(very_strong).
- **Applying when the inheritance mode is merely unstated** — Distinguishes a case whose
  inheritance mode was simply never recorded from one that explicitly contradicts the
  curated mode. In the former case, when the gene/disease has exactly one curated mode and
  every record agrees the mechanism is established, PVS1 is now applied with a disclosure
  that the case was assumed to follow that mode. RPE65 c.495+1dup now passes the mechanism
  gate itself, though it still doesn't reach MET, blocked by an unrelated gap (missing NMD
  prediction data).

Both changes follow the same principle: missing information is not the same as conflicting
information - the former is applied with disclosure, the latter still withheld.

## 5. Accuracy validation

### 4 demo cases (run via the same script path as the real API)

Of the 3 variants with a known answer, **1/3 (33.3%)** matched at the variant level;
**9/15 (60.0%)** matched at the criterion level.

| Code | Match | Code | Match |
|---|---:|---|---:|
| PVS1 | 2/2 (100%) | PS3 | 0/2 (0%) |
| PM1 | 1/1 (100%) | PP1 | 0/1 (0%) |
| PM2 | 3/3 (100%) | PS2 | 0/1 (0%, not implemented) |
| PM5 | 1/1 (100%) | | |
| PP3 | 2/2 (100%) | | |
| PS4 | 1/3 (33.3%) | | |

The automated layer (PVS1, PM1, PM2, PM5, PP3) matched every time. The literature engine
(PS3, PS4) and PP1 remain weak, consistent with the previously-identified pattern: papers
directly discussing the target variant are scarce.

### Full 64-variant ERepo run (first complete run)

Previous days were limited to 12- and 16-variant subsets; this time **all 64 variants ran
overnight** (about 3.5 hours - much faster than the 6-8 hour estimate).

- **59/64 evaluated** (5 skipped: 4 with unresolved coordinates [PTEN, USH2A, TECTA, TPM1],
  1 with no ground-truth classification [DSG2])
- **Variant-level match: 19/59 (32.2%)**
- **Criterion-level (MET→MET recall): 108/224 (48.2%)**

| Code | Match | Code | Match |
|---|---:|---|---:|
| BP4 | 8/8 (100%) | PP4 | 2/11 (18.2%) |
| PM2 | 37/38 (97.4%) | PS4 | 4/23 (17.4%) |
| PP3 | 30/31 (96.8%) | PP2 | 1/6 (16.7%) |
| PM1 | 12/15 (80.0%) | PS3 | 0/21 (0%) |
| BS1 | 4/5 (80.0%) | PVS1 | 0/5 (0%)* |
| PM5 | 4/6 (66.7%) | PS2/PM3/PM6/BP7/BA1/BS2 etc. | 0% (mostly not implemented) |
| PP1 | 5/15 (33.3%) | | |
| BS3 | 1/5 (20.0%) | | |

*Investigating PVS1's 0/5 led to the two policy changes in Section 4 (RUNX1 verified
fixed, RPE65 partially). The full 64-variant set has not yet been re-run with these fixes.

## 6. Other deliverables

- **GA4GH VA-Spec feedback** — filed as [ga4gh/va-spec#448](https://github.com/ga4gh/va-spec/issues/448)
- **`condition_from_literature_search.py` (new)** — literature-based condition resolution for variants ERepo has no record of
- **Live literature search for PP1/BS4 (`pp1_segregation_search.py`)** — searches for family-segregation data when none is supplied; verified to bring some variants to a real PP1 MET
- **64-variant ground-truth dataset** — `test_data/full_criteria_ground_truth.py`
- **Japanese translation of the team's paper**

## 7. Known gaps / handover notes

- **6 criteria not yet implemented** (PS2, PM3, PM6, BS2, BP2, BP5) — future work.
- **RPE65's missing NMD prediction data** — PVS1's mechanism gate now passes, but the
  variant-level decision tree can't reach MET without NMD prediction evidence that isn't
  available for this variant. Root cause found: querying Ensembl VEP directly
  (`numbers=1`) shows a plain splice_donor_variant (e.g. MYBPC3 c.2905+1G>A) gets a real
  intron number back, but RPE65 c.495+1dup - a duplication landing exactly on the
  exon/intron boundary, annotated as frameshift_variant + splice_region_variant - gets
  neither an exon nor an intron number from VEP at all. `acmg_pipeline/providers/nmd.py`
  is designed to report no record rather than guess when numbering doesn't settle the
  question, so this is working as intended, not a bug. A real fix (computing the exon
  position ourselves from CDS coordinates and transcript structure, or adding another
  annotation source) is a real implementation effort - **deferred as a known limitation**
  for now rather than attempted tonight; the current "leave it to a curator"
  (MANUAL_REVIEW/not_met) behavior stays as-is.
- **Re-running the full 64-variant set with today's PVS1 fixes** — not yet done; expected
  to improve the MET rate further.
- **Literature engine accuracy** — the bottleneck is the scarcity of papers that mention
  the target variant, not a defect in the LLM's judgment.
- **AutoPVS1 integration** — currently limited to a reference link; out of scope for this
  project.

---
BH26 ACMG/AMP Classification Pipeline &middot; 6-day hackathon development summary (as of 2026-09-19)
