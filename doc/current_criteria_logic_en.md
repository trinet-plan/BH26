# ACMG/AMP All 28 Criteria — Current Logic Overview

This document summarizes, based on the actual code, how this pipeline evaluates each ACMG/AMP 2015 criterion **as of right now**. It does not include a change history or a comparison with past versions.

Scope: all 28 criteria — 17 automated criteria (PVS1, PS1, PM1, PM2, PM4, PM5, PP2, PP3, PP5, BA1, BS1, BS2, BP1, BP3, BP4, BP6, BP7), 4 literature-LLM-judged criteria (PS3, BS3, PS4, BP5), 3 phenotype/segregation criteria (PP1, BS4, PP4), 2 de-novo criteria (PS2, PM6), 2 unimplemented criteria (PM3, BP2), and the final classification logic (`classify()`).

---

## 1. Automated criteria

### PVS1

Scope is limited to predicted loss-of-function (LoF) variants only: stop-gained, frameshift, canonical splice donor/acceptor, start-loss, or a splice variant with LoF effect confirmed by an RNA assay. When no disease condition is supplied, PVS1 itself is not evaluated; instead, a preliminary assessment of "what PVS1 would conclude if only a condition were known" is computed internally and kept for curators (this never feeds into the final MET/strength).

- **Mandatory gate**: for the gene/disease/inheritance-mode combination, "LoF is an established disease mechanism" must be confirmed from curated records (ClinGen Dosage / Gene2Phenotype / applicability data sourced from ClinGen CSpec). The CSpec applicability snapshot (`config/cspec_applicability.json`) is supplied through the same `gene_disease` category as Gene2Phenotype, so PVS1's own code needs no dedicated branch for it. Even without an exact match, it applies when every related record under a MONDO parent/child disease relationship agrees the mechanism is established (flagged for curator confirmation). When the inheritance mode is not recorded but curation carries only a single mode that agrees the mechanism is established, that mode is assumed and applied (also flagged for confirmation). A gene mismatch, an excluded disease, curation existing only for a different disease, or a mode conflict all stop at UNKNOWN. If Gene2Phenotype's `mechanism_support` is not `"evidence"` (e.g. `"inferred"` — circular reasoning from the variant type alone), that mechanism call is downgraded to "unresolved," and "undetermined" records, which used to be silently dropped, are now always emitted as a record too, so a weak inference-only mechanism claim can never stand unopposed without a newer record able to contest it (fixed 2026-09-25).
- **Decision tree by variant type**:
  - Truncating (stop_gained/frameshift): biological relevance of the transcript → NMD prediction. If NMD is predicted, also check the relevance of the affected exon for **MET (very_strong)**. If NMD is not predicted, disruption of a critical functional region gives **MET (strong)**; failing that, the LoF variant's population frequency (UNKNOWN if too common — that would contradict the premise that LoF causes this disease), the region's biological relevance, and the fraction of protein lost (above the threshold → **MET (strong)**, below → **MET (moderate)**) are evaluated in order. "Critical functional region" refers only to overlap with a UniProt `Domain` / `Active site` / `Coiled coil` / non-disordered `Region`; `Binding site` (a single-residue contact point) and `Motif` (a short trafficking signal, etc.) have been excluded since 2026-09-25 (in the real MYOC/PTEN false-positive cases, these small features — almost always present somewhere near a truncation point — were the cause of the misdetection).
  - Canonical splice / RNA-confirmed splice LoF: check whether alternative splicing rescues it and whether the reading frame is disrupted. If disrupted, it follows the same path as truncating variants; if not, it follows the same region-evaluation path as non-NMD truncating variants.
  - Start-loss: UNKNOWN (rescued) if a sound alternative transcript exists with a different translation start site. Otherwise, **MET (moderate/supporting)** depending on whether a downstream in-frame start site exists and whether there is upstream pathogenic evidence.
- **NOT_MET**: this logic never produces an explicit NOT_MET (not-applicable, missing information, and needs-confirmation all resolve to UNKNOWN).
- **Main UNKNOWN cases**: no condition supplied, an out-of-scope consequence, mechanism not established/not confirmed, missing or conflicting transcript/NMD/region/splice information, or the LoF variant is present at high frequency in the population.
- **Design notes**: internally distinguishes "not applicable (NOT_APPLICABLE)," "not evaluated (NOT_EVALUATED)," and "needs manual review (MANUAL_REVIEW)" (all appear externally as UNKNOWN, but are distinguishable in curator-facing metadata). Any call based on automated inference — UniProt-derived region prediction, the default splice policy, etc. — always carries a flag prompting curator confirmation.

### PS1

Delegated to `comparator.py`'s shared comparison logic. Missense variants only (splice-type consequences are UNKNOWN pending an equivalence check; anything else is not applicable). Searches for a "Pathogenic" candidate in ClinVar with the same residue, the same amino-acid change, an independent variant, and the same transcript.

- **MET**: a candidate satisfies either the "automated-call conditions" (exact match, independence, eligible review status, splice mechanism checked) or the "manually-reviewed conditions" (disease match, everything confirmed, curator notes present). Strength is strong.
- **NOT_MET**: a complete comparison search finished but no eligible pathogenic candidate was found.
- **UNKNOWN**: missing annotation, non-missense, missing protein context, a candidate exists but independence/mechanism confirmation is incomplete, or the search itself is incomplete.
- **Design note**: to avoid circularity (using a ClinVar classification that itself derives from PS1 to justify PS1), independent evidence, mechanism agreement, and a splice-mechanism check are all mandatory.

### PM1

Delegated to `regions.py`'s `evaluate_region`. A choice between two routes: "mutational hotspot" and "critical functional domain." For some genes (RYR1, MECP2, BMPR2, etc.), a dedicated critical-domain call based on codon ranges hand-transcribed from VCEP specifications (`gene_critical_domains`) is evaluated before the ordinary region category.

- **MET**: on the hotspot route, against the configured policy (window width, pathogenic/benign count thresholds), benign variants are at or below the threshold and pathogenic variants are at or above it. On the critical-domain route, both `critical_functional_region` and `benign_depletion` are True. Strength is moderate.
- **NOT_MET**: on the hotspot route, benign variants exceed the threshold; on the critical-domain route, either field is False.
- **UNKNOWN**: missing or inconsistent protein coordinates, region not retrieved, hotspot policy not configured or version mismatch, missing counts, or insufficient density with no benign variants either (insufficient evidence).
- **Design note**: for VCEP-sourced critical-domain range data, `benign_depletion` defaults to True, and this is always surfaced as a review item.

### PM2

Retrieves population-frequency observations via `common.population_context()` and judges rarity.

- **MET**: the highest valid-observation AF is at or below the configured threshold (demo value 5e-05), or every source reports "zero observations" — either gives a weak MET (supporting).
- **NOT_MET**: the highest valid-observation AF exceeds the threshold.
- **UNKNOWN**: threshold/policy not configured, or a real provider error left the search incomplete.
- **Design note**: clearly distinguishes a configured threshold of 0 from "not configured." Set lower than BS1's default threshold.

### PM4

Delegated to `regions.py`'s `evaluate_region`. In-frame insertion, in-frame deletion, and stop-loss only.

- **MET**: functional review complete, the protein-length change is nonzero, and it is not in a non-functional repeat region. Strength is moderate.
- **NOT_MET**: the protein-length change is zero, or it falls inside a non-functional repeat region.
- **UNKNOWN**: out-of-scope consequence, region not retrieved / multiple regions, inconsistent protein coordinates, or functional review incomplete.
- **Design note**: shares the same region category and coordinate-validation logic as PM1.

### PM5

Shares PS1's comparator logic (the same function, the same candidate pool). Differs from PS1 in targeting only candidates with the same residue but a *different* amino-acid change.

- **MET**: a same-residue, different-amino-acid-change, independent pathogenic candidate satisfies either the automated-call or the manually-reviewed condition. Strength is moderate.
- **NOT_MET**: a residue-scoped comparison search finished but no eligible candidate was found.
- **UNKNOWN**: largely the same as PS1 (missing annotation, non-missense, insufficient candidate confirmation, incomplete search).
- **Design note**: to guarantee the validity of "same mechanism even at a different site," this weighs mechanism agreement and the absence of a splice-mechanism mismatch more heavily than PS1 does.

### PP2

Shares `mechanism.py`'s common missense-mechanism logic with BP1. Missense variants only.

- **MET**: a curated gene-disease record has `missense_mechanism_established`, `spectrum_review_complete`, and `low_benign_missense_variation` all True. Strength is supporting. A `CANDIDATE` call from the statistical suggestion (ClinGen Gene-Disease Validity + gnomAD missense-constraint Z-score) may also be adopted as a weak MET (flagged). Since 2026-09-25 (policy v3), this statistical suggestion requires not only misZ exceeding its threshold but also a low (default ≤10%) benign fraction among ClinVar-curated missense variants (real BMPR2 false positive: misZ=3.15 just barely exceeded the threshold, but the real ClinVar pathogenic/benign missense ratio showed 36.8% benign — not "low" at all).
- **NOT_MET**: `missense_mechanism_established` or `low_benign_missense_variation` is False, or the statistical suggestion is `NOT_SUGGESTED`.
- **UNKNOWN**: non-missense, no record retrieved and the suggestion is also insufficient/negative, an expert panel explicitly rules it not applicable, or the ClinGen CSpec applicability snapshot explicitly marks PP2 "Not applicable" for that gene (added 2026-09-25 — for genes where a real VCEP specification has explicitly decided not to use PP2 at all, e.g. the epilepsy sodium-channel panels for SCN2A/SCN1A).
- **Design note**: an expert panel's "not applicable" call is accepted as a record even without a mechanism field. The statistical suggestion produces nothing at all when INSUFFICIENT/CURATED_NEGATIVE, rather than fabricating a result from an uncertain signal. The CSpec "Not applicable" gate is evaluated before the statistical fallback, giving priority to a real VCEP's explicit decision.

### PP3

Delegated to `computational.py`'s calibrated-predictor (REVEL, SpliceAI, etc.) score-band judgment.

- **MET**: the predicted score for either mechanism (protein impact / splicing impact) falls within its calibrated band. When more than one applies, the strongest strength is used.
- **NOT_MET**: a calibration applies, but no score falls within its band.
- **UNKNOWN**: missing/multiple annotation, no calibration matches the target consequence, or the score is unavailable, conflicting, or out of range.
- **Design note**: an explicit calibration approach using only the score bands the source publication defines — not a vote count or a cross-threshold.

### PP5

Always returns a fixed UNKNOWN. A placeholder implementation following ClinGen policy (never score from an external database's assertion alone, without checking the primary evidence); there is no branching logic.

### BA1

`threshold_policy()` determines either a gene-specific VCEP threshold (`ba1_threshold_override`, adopted only when every field is present) or the default threshold (AF > 0.05).

- **MET**: a valid observation exceeds the threshold, and it is not a known exception. Strength is stand_alone.
- **NOT_MET**: every observation is at or below the threshold, or it exceeds the threshold but is a known exception.
- **UNKNOWN**: threshold/statistic/comparison operator not configured, `population_context` itself unresolved, or the exception-list check is missing.
- **Design note**: a gene-specific threshold replaces the default wholesale only when every field is present (no partial override). Whenever `faf95` is used, the record always states that it is this project's own Wilson-score computation, not the value gnomAD itself publishes.

### BS1

Prioritizes a disease-specific threshold (`disease_frequency_threshold`, per-variant or gene-wide); falls back to the default threshold (demo value 0.0001) when condition, provenance, or inheritance mode are not all present.

- **MET**: a valid observation exceeds the threshold. Even when the adopted threshold is still DRAFT (unapproved), MET is still returned, but it is explicitly marked as pending approval.
- **NOT_MET**: every observation is at or below the threshold.
- **UNKNOWN**: neither the disease-specific nor the default threshold is configured or valid, or `population_context` is unresolved.
- **Design note**: treats "is there evidence of exceeding the threshold" and "is the threshold itself approved" as separate questions. Every input the threshold was derived from (prevalence, inheritance, penetrance, etc.) is recorded, and any gap is named explicitly.

### BS2 (implemented 2026-09-25)

Long unimplemented on the grounds that it "needs the patient's own clinical record and is out of scope for the literature pipeline," it is now automated under a different, commonly-used real-world interpretation: a gnomAD genotype count (homozygote count / hemizygote count) stands in as a proxy for "observed in a healthy adult." **Scope is limited to autosomal recessive (AR), semidominant, and X-linked** — autosomal dominant (AD), X-linked dominant, and mitochondrial inheritance are deliberately out of scope (a single heterozygous observation there is essentially the same signal BA1/BS1's allele frequency already captures, and would double-count it).

- **Determining the inheritance mode**: `input_data["inheritance"]` (from the clinical note, unset for most ERepo-sourced records) takes priority. When unset, automated resolution is attempted from Gene2Phenotype (`"gene_disease"`) and ClinGen Gene-Disease Validity (`"gene_disease_validity"`) curated records, **but only when every source agrees** (the same "adopt only on agreement, give up on disagreement" caution as PVS1's own mechanism-candidate presentation). An exact MONDO condition-ID match is not required — in practice, the same disease is often registered under several MONDO IDs at different granularity.
- **Field checked**: `homozygote_count` for AR/semidominant, `hemizygote_count` for X-linked (UNKNOWN if the variant's chromosome is not X). Semidominant (ClinGen Gene-Disease Validity's own MOI notation "SD," e.g. LDLR / familial hypercholesterolemia) is treated the same as AR — the very definition of semidominant (heterozygotes: mild, late-onset; homozygotes: severe, fully penetrant, early-onset) is exactly the homozygous-side signal BS2's own question of "fully penetrant, early onset" is asking about.
- **Threshold**: the same two-tier "gene-specific override → configured default" shape as BA1 (`bs2_threshold_override` takes priority). No gene-specific threshold has been curated yet, so the configured default (provisional value: MET on even a single homozygote/hemizygote) is always used at present.
- **MET**: the genotype count exceeds the threshold. Strength is strong. When the adopted threshold is DRAFT (unapproved) — which is currently always — MET is still returned, but explicitly marked as a prediction pending approval (the same design as BS1).
- **NOT_MET**: every observation is at or below the threshold.
- **UNKNOWN**: the inheritance mode is neither AR/semidominant/X-linked nor resolvable, the mode is X-linked but the chromosome is not X, the threshold is not configured, `population_context` is unresolved, or the genotype-count field itself cannot be retrieved.
- **Design note**: since 2026-09-25 the gnomAD provider fetches `homozygote_count`/`hemizygote_count` (previously AC/AN/AF only). The VA-Spec `methodType` is "Case-Control Enrichment Assessment," not "Population Data Assessment" (following the official GA4GH reference model's own categorization).

### BP1

Shares `mechanism.py`'s common logic with PP2. Missense variants only.

- **MET**: `spectrum_review_complete=True`, `predominantly_truncating=True`, and `missense_mechanism_established=False`, all at once. Strength is supporting.
- **NOT_MET**: those fields are all present but the condition is not met.
- **UNKNOWN**: non-missense, no record retrieved, an expert panel explicitly rules it not applicable, the ClinGen CSpec applicability snapshot explicitly marks BP1 "Not applicable" (the same 2026-09-25 fix as PP2), or review is incomplete.

### BP3

Delegated to `regions.py`'s common region evaluation. In-frame insertion/deletion only.

- **MET**: `repetitive=True` and `functional_importance=False`. Strength is supporting.
- **NOT_MET**: the condition is not met.
- **UNKNOWN**: out-of-scope consequence, region not retrieved, inconsistent protein coordinates, or functional review incomplete.
- **Design note**: the automated provider judges `repetitive` from overlap with a UniProt feature (Repeat / Compositional bias) as well as sequence complexity (the Wootton-Federhen/SEG method, a 24-aa window, a 2.2-bit threshold). Only overlap with a Domain / Binding site / Active site / Motif / non-disordered Region sets `functional_importance` to True.

### BP4

Delegated to `computational.py`'s calibrated-predictor score-band judgment (the counterpart to PP3).

- **MET**: some calibration's BP4-side band is hit, and no calibration hits its PP3-side band.
- **NOT_MET**: some calibration hits its PP3-side band, or no calibration hits its BP4-side band at all.
- **UNKNOWN**: calibration policy not configured, no prediction for the target predictor, an invalid/out-of-range score, or a conflict with `high_confidence_null_or_splice`.
- **Design note**: predictors deliberately have an "abstention band" (a score range that falls in neither band), and abstention is never treated as a counter-argument.

### BP6

Always returns a fixed UNKNOWN. A policy gate for "never score from an external classification alone without checking the primary evidence"; there is no branching logic.

### BP7

Synonymous variants only. Intronic/non-coding/UTR-type consequences are UNKNOWN pending confirmation; any other non-synonymous consequence is not applicable.

- **MET**: `outside_splice_critical_region=True` and `no_predicted_splice_impact=True`. Strength is supporting.
- **NOT_MET**: the condition is not met.
- **UNKNOWN**: out-of-scope consequence, no record retrieved, missing required information, or a conflict with RNA evidence.
- **Design note**: the criterion is judged solely on position outside the splice-critical region and a SpliceAI score threshold, with no requirement of sequence conservation (per ClinGen SVI's 2023 policy). The automated provider reuses the SpliceAI result already bundled with the VEP annotation and makes no new HTTP request of its own.

---

## 2. Literature-LLM-judged criteria (PS3, BS3, PS4, BP5)

**Shared framework**: for the target variant, ERepo's own cited papers (`evidence_pmids`) are checked first, falling back to a live PubMed search for PMID candidates when there are none. Full text (or the abstract, if unavailable) is fetched for every PMID obtained, and each paper is sent to the LLM for structured JSON extraction. Each paper's result feeds into a multi-paper aggregation step: if every paper agrees on the same direction, that direction is adopted; if the direction is split, it falls back safely to `not_clear` (with a hint prompting human confirmation).

### PS3/BS3

First judges whether the target paper actually examined the target variant at all (e.g., in a large-scale saturation mutagenesis screen, if the variant's own specific value never appears in the text, it is treated as not established). Each experiment is then independently judged as functionally abnormal, normal, intermediate, or mixed. As a safety net, a forced override to `not_clear` occurs when the free-text description contradicts the classification, when a confident direction is reached despite zero experiments, or when there is no quantitative numeric support at all. When a confident direction is reached despite abnormal and normal findings coexisting within the same paper, the LLM's own rationale is checked for a reference to the conflicting finding or a contrastive phrase ("although," "despite," etc.); if found, this is treated as "a conclusion actually reached after weighing the conflict," and the original direction is kept with a caution flag; if not, it is overridden to `not_clear` (fixed 2026-09-25 — previously, any conflicting finding triggered an unconditional override, but a real case, PTEN c.112C>T — whose true answer is PS3=MET (moderate) — turned out to have been wrongly overridden by this safety net's own flawed premise).

### PS4

Specialized for extracting case-control data. The study design (case_control / family_cohort / case_series / not_case_control) is judged, and carrier counts for affected/unaffected individuals, odds ratios, and p-values are extracted. A forced override to `not_clear` occurs when the study design is `not_case_control` yet PS4 was concluded, or when a confident direction is reached with absolutely no quantitative support (no counts, no odds ratio, no p-value). Even when a single paper does not meet the case-control design bar, if two or more independent case-series papers (none of which contains a contradicting unaffected-carrier count) can be pooled to reach a total of at least 3 affected carriers, this pooled aggregation is adopted as a fallback — but only when the standard aggregation itself came back `not_clear` (added 2026-09-25 — real ClinGen PS4 determinations are often reached by pooling several case-series papers rather than from a single case-control-design paper).

### BP5

For a patient carrying the target variant, extracts whether a pathogenic variant in a different gene that better explains the phenotype is reported within the same paper (inference from the mere fact that the target variant's own classification is uncertain/benign is explicitly forbidden). A forced override to `not_clear` occurs when BP5 was concluded despite no alternate gene having been extracted at all, or when BP5 was concluded despite the alternate finding never being stated to explain the phenotype.

- **A real limitation of the data source**: ERepo's cited papers are the basis for the variant's overall classification, and do not necessarily point to BP5's own specific basis (the discovery of an alternate diagnosis). The live-search fallback used when ERepo has no cited paper is also a generic query centered on the target gene's name, making it hard to find a case report describing an alternate diagnosis.

---

## 3. Phenotype/family-segregation criteria (PP1, BS4, PP4)

`pp1_bs4_pp4_engine.py` implements this by wiring in ClinGen 2024's Bayesian-points evaluation logic (`pp4_pp1_bs4.evaluate_locus_evidence`).

- **PP4**: the target variant's ERepo-cited papers (`evidence_pmids`) are judged first; if there are none, a live literature search (`pp4_literature_search`) is run using the patient's diagnosis as the query, and evaluation proceeds only if a confirmed diagnostic-yield statistic is found (ERepo-first since 2026-09-25). The search always tries both a query narrowed with `"NOT multigene panel NOT gene panel NOT NGS panel"` and the original generic query; at extraction time the LLM is also asked to judge whether the solved cohort's population is gene-specific or a broad multi-gene panel, and the search keeps going until a gene-specific paper is found (a broad-panel yield is likely diluted by other genes, so it is adopted only as a last-resort fallback, disclosed as a CAUTION). UNKNOWN if neither is found.
- **PP1**: uses the patient's own clinical-note family status (`family.relatives`) directly when present; otherwise supplements it with a live literature search for pedigree data. The search judges the target variant's ERepo-cited papers first, and if still not found, supplements with a query built from the protein notation (both the 1-letter and 3-letter forms) — the raw HGVS c. notation is never included in the search query at all (empirically confirmed to return zero results even for a variant whose real source paper is indexed).
- **PP1's and PP4's points are not floored independently** — when both have evidence, they are converted following Biesecker et al. 2024's Table 4 (a combined +5.0-point cap), by a deterministic rule that apportions the higher strength to whichever side has the larger raw point value.
- **BS4**: judged solely from a separate module, `segregation.py`'s `bs4_met` boolean; fixed at STRONG when MET.

`segregation.py` itself is an independent implementation that applies the same literature-LLM-extraction framework as PS3/BS3 to per-family pedigrees (counts of affected/unaffected × variant-carrying status), judging PP1 when every family shows complete co-segregation and BS4 when some family shows non-co-segregation. A forced override to `not_clear` occurs when the actual counts contradict the concluded direction (the most important safety net), or when a confident direction is reached despite zero families.

---

## 4. De-novo criteria (PS2, PM6)

Rule-based evaluation using only `ClinicalNoteExtraction.de_novo` (parental variant status, paternity/maternity confirmation status) — no literature search or LLM call is made. When both parents are confirmed variant-negative and the patient is affected with a negative family history, a mutually-exclusive call is made: if both paternity and maternity are confirmed, PS2 is MET (STRONG) and PM6 is NOT_MET; if confirmation is incomplete, the reverse — PM6 is MET (MODERATE) and PS2 is NOT_MET.

---

## 5. Unimplemented criteria (PM3, BP2)

`stubs.py` always returns a fixed UNKNOWN (`stub_evidence()`), since these depend on trans/cis phase information against another variant and cannot be judged from literature or from any provider. This is clearly distinguished from NOT_MET (a negative determination). BS2 left this classification on 2026-09-25, implemented as an automated criterion using gnomAD genotype counts as a proxy (see §1).

---

## 6. Final classification logic (the `classify()` function)

`acmg_pipeline/classification.py` takes the 28 codes' worth of `CriterionEvidence` (`code` / `status` / `strength` / `source`) and merges them into a single `ClassificationResult`, based on the Tavtigian et al. 2018 Bayesian point system.

### Bayesian points

| Strength | Pathogenic side | Benign side |
|---|---|---|
| SUPPORTING | +1 | -1 |
| MODERATE | +2 | -2 |
| STRONG | +4 | -4 |
| VERY_STRONG | +8 | -8 |
| STAND_ALONE | +8 (treated the same as VERY_STRONG, for convenience) | -8 (same) |

`score` is the sum of the points of every code that reached MET.

### Category thresholds

| Score range | Category |
|---|---|
| score ≥ 10 | Pathogenic |
| 6 ≤ score ≤ 9 | Likely Pathogenic |
| 0 ≤ score ≤ 5 | Uncertain Significance |
| -6 ≤ score ≤ -1 | Likely Benign |
| score ≤ -7 | Benign |

### The BA1 override

If BA1 is MET, Benign is returned immediately regardless of any other evidence or the score (`ba1_override=True`). The score itself is still computed normally as reference information, and no automatic check exists for conflicting strong pathogenic evidence (since all the evidence stays in the `met` list, a human curator can still review it).

### Handling duplicates and unevaluated codes

Passing more than one `CriterionEvidence` for the same code raises a `ValueError` (the caller must narrow it down to one beforehand). `not_evaluated_codes` collects "codes among the 28 with no evaluation result, or UNKNOWN," and does not distinguish an unimplemented stub code from a code that "is implemented but came back UNKNOWN because this particular input gave it nothing to judge from" (that distinction is `registry.is_implemented()`'s own responsibility).

### `from_aggregated_judgment`

The adapter that converts the literature-LLM pipeline's `AggregatedJudgment` into `CriterionEvidence`. UNKNOWN if the judged direction is `not_clear`; NOT_MET if the direction differs from the code being evaluated (e.g., while evaluating BS4, if evidence pointing toward PS3 is found instead, BS4 is never mistaken for MET); MET if it matches, with strength decided from the number of papers adopted.

### The "UNKNOWN → not_met" conversion in VA-Spec output (a caution for accuracy tallies)

`acmg_pipeline/export.py`'s `build_evidence_line()` / `build_automated_evidence_line()` write the machine-readable `status` extension field as `not_met` whenever an implemented evaluation actually ran but could not reach a verdict (UNKNOWN) — this does not apply to unimplemented stub codes, which `build_stub_evidence_line()` separately and honestly reports as `unknown`. This is a deliberate design (2026-09-18) on the grounds that, for `classify()`'s Bayesian scoring, "not established" and "could not be judged" should both count as zero points; the human-facing `description` text and `curatorHints` (`category: "unevaluated"`) spell out the real reason (that it is being reported as not_met despite being unable to be evaluated).

**This design is correct for `classify()` itself, but if a per-criterion recall/specificity tally is computed by naively reading only the `status` extension field out of the exported JSON, a "could not be evaluated" case gets miscounted as "correctly determined not_met," inflating the apparent accuracy on the not_met side** (discovered 2026-09-26 while re-verifying BS2's and PP4's accuracy — see the "Cross-cutting trends" section of [`doc/criteria_accuracy_2026-09-30_en.md`](criteria_accuracy_2026-09-30_en.md) / [`doc/criteria_accuracy_2026-09-22_ja.md`](criteria_accuracy_2026-09-22_ja.md) for details). Tallying per-criterion accuracy correctly requires either filtering by curatorHints' `unevaluated` flag, or reading `CriterionResult.status` directly.
