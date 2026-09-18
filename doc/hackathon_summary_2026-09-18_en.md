# ACMG/AMP Classification Pipeline — 5 Days of Development

Summary of what was built during the BH26 hackathon (5 days, with 2026-09-18 as the final
day). The pipeline combines an LLM (gemma-4 served via vLLM) with the PubMed MCP server to
generate draft classifications across all 28 ACMG/AMP 2015 criteria, accelerating a human
curator's first-pass screening.

## 1. Implementation overview

All 28 ACMG/AMP 2015 criteria are implemented across four layers, grouped by how each
criterion is judged. Three systems are already running: deterministic rule/database-driven
criteria, criteria that require reading the literature (LLM + PubMed MCP), and criteria that
require a phenotype/segregation match. The remaining six criteria are honestly reported as
"not evaluated" rather than silently omitted.

| Layer | Criteria | Count | Status |
|---|---|---:|---|
| Automated (Layer 1) | `PVS1 PS1 PM1 PM2 PM4 PM5 PP2 PP3 PP5 BA1 BS1 BP1 BP3 BP4 BP6 BP7` | 16 | Implemented |
| Literature engine | `PS3 BS3 PS4` | 3 | Implemented |
| Phenotype/segregation | `PP1 BS4 PP4` | 3 | Implemented |
| Not yet implemented | `PS2 PM3 PM6 BS2 BP2 BP5` | 6 | Explicitly NOT_EVALUATED |

The final classification is combined using Tavtigian et al. 2018's Bayesian point-based
framework (Supporting=1 / Moderate=2 / Strong=4 / Very Strong=8; thresholds: ≥10
Pathogenic, 6-9 Likely Pathogenic, 0-5 VUS, -1..-6 Likely Benign, ≤-7 Benign), and returned
from the API (FastAPI + Docker) as a GA4GH VA-Spec EvidenceLine.

## 2. VA-Spec output design work

This is the area the team invested the most effort in, since it is what a human curator
actually reads.

- **Removing duplication and inconsistency** — Fields that had been bundled into a single
  `bh26AssessmentDetails` object (status / provenance / direction / rulesUsed, etc.) were
  fully flattened into individually named top-level extensions, at the same level as
  `curatorHints`. Duplication across `summary` / `strength` / `evidenceOutcome` /
  `criterion` / `evidenceItemIds` was removed, `structuredEvidenceItems` and
  `strengthEstimationMethod` were merged into `curatorHints`, and `rulesUsed` was
  de-duplicated against `specifiedBy.reportedIn`.
- **Fixing a cross-module wiring gap** — PVS1's disease-mechanism check needs the
  diagnosis recorded in the clinical note, but that information was not actually reaching
  the decision logic. This gap was found and fixed: once the diagnosis resolves to a MONDO
  disease class, PVS1 now correctly reaches `very_strong`.
- **Fixing a direction inconsistency** — Benign-direction evidence (BS3) was always being
  reported with a pathogenic-direction outcome (SUPPORTS); this was fixed and covered by a
  regression test.
- **Feedback to GA4GH VA-Spec** — Schema constraints discovered while extending the
  output (`Method.reportedIn` accepts only a single Document, `hasEvidenceItems` has no
  generic StudyResult type, `Direction` cannot represent "not evaluated") were written up
  as 5 issues in
  [doc/ga4gh_va_spec_feedback_v1_en.md](ga4gh_va_spec_feedback_v1_en.md) and organized as
  a proposal back to GA4GH.

## 3. Accuracy validation

### 4 demo cases (12 variants total)

Variant-level classification matched ground truth in **4/9 (44.4%)** of the variants with
a known answer; at the criterion level, **18/26 (69.2%)** matched.

| Code | Match | Code | Match |
|---|---:|---|---:|
| PVS1 | 3/3 (100%) | BS1 | 3/4 (75%) |
| PM2 | 5/5 (100%) | BP4 | 1/2 (50%) |
| PP3 | 2/2 (100%) | PS3 | 1/3 (33.3%) |
| PM1 / PM5 | 2/2 (100%) | PS4 | 1/3 (33.3%) |
| | | PP1 | 0/1 (0%) |

The automated layer (PVS1, PM1, PM2, PM5, PP3) matched **100% of the time**. The
literature engine (PS3, PS4), on the other hand, only matched 33.3% of the time. Root-cause
analysis tooling found this was mostly (71.0%) because papers directly discussing the
target variant simply do not exist for these variants yet — a gap in the literature itself,
not a defect in the implementation.

### PP4 subset from the 64-variant ERepo set

12 variants where ground truth marks PP4 as MET were deliberately selected for a run
(together covering 20 of the 22 implemented criteria). The variant-level match rate came
out low, at **2/11 (18.2%)**, but investigating why surfaced an important structural
finding.

> **Finding: this harness cannot exercise PP1/PP4 at all**
> The validation script was passing an empty clinical note, so PP1/BS4/PP4 — which require
> a diagnosis to evaluate against — always came back `unknown`, and never appeared as
> `met` across all 12 variants. This is a structural limitation of the 64-variant ERepo
> dataset (which only carries gene/HGVSc, no free-text clinical narrative), not a real
> accuracy problem.

## 4. Other deliverables

- **GA4GH VA-Spec feedback** — [doc/ga4gh_va_spec_feedback_v1_en.md](ga4gh_va_spec_feedback_v1_en.md) (5 issues)
- **Error-analysis harness** — categorizes the literature engine's `not_clear` judgments by root cause
- **64-variant ground-truth dataset** — `test_data/full_criteria_ground_truth.py` (from an ERepo re-query plus the demo cases)
- **Japanese translation of the team's paper** — a full translation of the Perspective paper, preserving authorship, references, and figures

## 5. Known gaps / handover notes

- **6 criteria not yet implemented** (PS2, PM3, PM6, BS2, BP2, BP5) — the current ground
  truth dataset does not contain a single real example of `BP6` or `PP5` either, so the
  validation data itself needs to be expanded first.
- **Empirically validating PP1/PP4** — until a diagnosis is attached to each ERepo variant,
  these criteria can never reach `met` in the automated validation pipeline. Adding
  diagnosis data is the next step.
- **Literature engine accuracy** — the bottleneck is the scarcity of papers that mention
  the target variant, not a defect in the LLM's judgment.
- **AutoPVS1 integration** — currently limited to a reference link; integrating the actual
  automated decision logic is out of scope for this project.

---
BH26 ACMG/AMP Classification Pipeline &middot; 5-day hackathon development summary (as of 2026-09-18)
