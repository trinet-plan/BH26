# Accuracy Summary for All 28 Criteria (2026-09-30 edition)

Measurement point: the 2026-09-26 363-variant integrated validation run (the 95 variants affected by an Anthropic credit shortage mid-run were cleanly re-run with `--only-file` and merged back in). 291 of the 363 variants are usable for both the final-classification comparison and the per-criterion accuracy comparison (the remaining 72 are unusable: 71 with unresolved genomic coordinates, plus 1 that hit a VA-Spec output schema error on PS4).

This document supersedes [`doc/criteria_accuracy_2026-09-22_ja.md`](criteria_accuracy_2026-09-22_ja.md) (kept for history), recomputing everything on the same 291-variant real data. The detailed root-cause investigation for each criterion is left to [`doc/current_criteria_logic_ja.md`](current_criteria_logic_ja.md); this document is limited to accuracy numbers, cross-cutting trends, and next steps.

## Met/not_met table for all 28 criteria

| Layer | Code | met | not_met | Notes |
|---|---|---|---|---|
| auto | PVS1 | 44/49 (89.8%) | 62/62 (100%) | MYOC/PTEN's 5 false positives and the Gene2Phenotype mechanism fix are both reflected |
| auto | PS1 | 3/7 (42.9%) | 137/137 (100%) | Known low recall (the ground truth has many splice-affecting variants, which PS1 is designed to leave out of scope). Decided not to address this |
| manual | PS2 | 0/9 (0.0%) | 132/132 (100%) | This validation data has no parental information in the clinical note, so it structurally tends to be unevaluable |
| manual | PS3 | 10/36 (27.8%) | 135/138 (97.8%) | No ERepo-cited papers / the experiment-conflict safety net over-triggering - both fixed |
| manual | PS4 | 5/61 (8.2%) | 114/116 (98.3%) | Added pooled case-series aggregation |
| auto | PM1 | 15/20 (75.0%) | 107/119 (89.9%) | Fixed with VCEP codon-range data |
| auto | PM2 | 111/122 (91.0%) | 90/97 (92.8%) | |
| manual | PM3 | - | - | **Not implemented** (needs trans-phase information against another variant, which ERepo does not carry alongside its final calls - structurally unobtainable) |
| auto | PM4 | 1/1 (100%) | 122/123 (99.2%) | Sample size too small to be more than a reference value |
| auto | PM5 | 11/65 (16.9%) | 95/96 (99.0%) | Known low recall (same reason as PS1). Decided not to address this |
| manual | PM6 | 0/7 (0.0%) | 113/113 (100%) | Same reason as PS2 (missing clinical-note data) |
| manual | PP1 | 2/31 (6.5%) | 115/124 (92.7%) | ERepo-first search plus the query rewrite are both in, but recall remains low (a structural limit of literature search) |
| auto | PP2 | 24/33 (72.7%) | 87/87 (100%) | Fixed with the CSpec applicability gate plus the ClinVar benign-missense fraction |
| auto | PP3 | 52/52 (100%) | 129/141 (91.5%) | |
| semi-auto | PP4 | 3/18 (16.7%) | 81/113 (71.7%) | ERepo-first search and the panel-excluding query are both in, but the cohort-specificity problem remains |
| auto | PP5 | - | - | Not used (always fixed at UNKNOWN, by policy) |
| auto | BA1 | 35/74 (47.3%) | 136/139 (97.8%) | Fixed with per-VCEP thresholds and gnomAD popmax support |
| auto | BS1 | 53/59 (89.8%) | 105/140 (75.0%) | Same as above |
| auto | BS2 | 4/40 (10.0%) | 106/110 (96.4%) | Newly implemented 2026-09-25. The ground truth skews heavily toward autosomal-dominant genes, which fall outside this criterion's automated scope |
| manual | BS3 | 8/56 (14.3%) | 121/123 (98.4%) | Same engine and same fix as PS3 |
| manual | BS4 | 0/24 (0.0%) | 111/111 (100%) | Same engine and same fix as PP1 |
| auto | BP1 | 48/48 (100%) | 103/103 (100%) | Fixed with the CSpec applicability gate |
| manual | BP2 | - | - | **Not implemented** (same as PM3: needs trans/cis phase information, structurally unobtainable) |
| auto | BP3 | 6/14 (42.9%) | 105/105 (100%) | |
| auto | BP4 | 87/90 (96.7%) | 79/109 (72.5%) | |
| manual | BP5 | 0/26 (0.0%) | 109/111 (98.2%) | The structural limit of insufficient ERepo-cited papers remains unresolved |
| auto | BP6 | - | - | Not used (always fixed at UNKNOWN, by policy) |
| auto | BP7 | 49/73 (67.1%) | 90/90 (100%) | |

## Accuracy of the final Classification

Comparing `classify()`'s final category (Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign) against ground truth (291 variants).

- Exact match: 130/291 (44.7%)
- **Within one tier** (a shift to an adjacent category, e.g. Pathogenic → Likely Pathogenic): 262/291 (90.0%)
- **Dangerous flips** (a complete mix-up between the benign side and the pathogenic side): 3/291 (1.0%)

### Confusion matrix (rows = ground truth, columns = engine's call)

| Truth \ Engine | Benign | Likely Benign | VUS | Likely Pathogenic | Pathogenic |
|---|---|---|---|---|---|
| Benign | 47 | 50 | 9 | - | - |
| Likely Benign | 7 | 32 | 10 | 1 | - |
| Uncertain Significance | - | 12 | 26 | 4 | 2 |
| Likely Pathogenic | - | 1 | 8 | 10 | 3 |
| Pathogenic | - | 1 | 15 | 38 | 15 |

Breakdown of the 3 dangerous flips:
- ATM c.7271T>G: Pathogenic → Likely Benign
- TNNI3 c.485G>A: Likely Pathogenic → Likely Benign
- PTEN c.114T>G: Likely Benign → Likely Pathogenic

The exact-match rate is on the low side, but most of the miss is a one-tier shift to an adjacent category; genuinely dangerous mix-ups stay at 3 of 291 (1.0%).

## Cross-cutting trends

Patterns that emerged from this round of development and validation, beyond any single criterion's own root cause.

- Criteria that go through literature search (PubMed × LLM) - PS3/BS3/PS4/PP1·BS4/PP4/BP5 - consistently trail the automated codes in accuracy (the specific cause differs per criterion: insufficient ERepo citations, a search query that isn't specific enough to the right cohort, an over-triggering safety net - but the pattern of "literature-search criteria lag behind" is shared)
- The "UNKNOWN reported as not_met" export design creates a real measurement pitfall: a naive accuracy count reads the not_met side as more accurate than it really is (discovered via BS2 and PP4 - the high not_met figures in the table above may carry some of this same inflation)
- MONDO and other ontology granularity mismatches keep tripping up logic that requires an exact condition match, and this recurs across more than one criterion (BS1, BS2, PVS1)
- Dangerous flips (a complete mix-up between the benign side and the pathogenic side) hold steady at around 1%, regardless of the low exact-match rate - the system behaves conservatively, mostly drifting by one adjacent tier rather than flipping sides outright
- A small, illustrative spot-check tends to overstate real-world performance (BS2's 12-case check found 7/12 matches, but the real figure across all 291 variants is closer to 10%)

## Next steps

- Improve recall for the literature-search criteria (broader search queries and data sources)
- Rework how accuracy itself is measured - today's VA-Spec output files are designed to report "could not be evaluated" as `not_met`, so computing accuracy by reading them makes the real figure look better than it is. This needs a switch to tallying the actual internal judgment (met / not_met / could-not-evaluate as three distinct states) directly, rather than reading it back out of the export
- Structure curatorHints' reason categories - today only a narrow set of categories exists (e.g. "unevaluated", meaning "reported as not_met even though it could not actually be evaluated"). Categorizing the reason behind a not_met result more systematically - "not applicable", "searched but inconclusive", "diluted population", and so on - would let correct accuracy counting be done straight from the exported JSON, which also serves the accuracy-measurement rework above
- Consolidate MONDO normalization into a shared helper, instead of each criterion handling it separately
- Improve BS2's inheritance-mode resolution coverage (broaden Gene2Phenotype / ClinGen Gene-Disease Validity coverage)
- Keep an eye on readiness for the ACMG/AMP 4.0 guidelines
