# BH26 Demo Case 1-4: Real-Case Overview, ClinGen Registration Status, and Ground Truth (Integrated)

**version: v6.0** (translated from real_cases_groundtruth_integrated_v6_ja.md)

Created by: Trinet, Takahashi / 2026-09-08 (updated 2026-09-10)
Scope: `case1_variants.vcf` through `case4_variants.vcf`, plus the corresponding `case{1-4}_clinical_note.txt` (demo-cases-v3-real)

---

## 1. Definition of the confidence tiers

Rather than simply presenting the ground-truth data the comparison tool references as "the answer" without qualification, we make explicit **how confident that answer itself is**.

| Tier | Definition |
|---|---|
| **A** | Formally verified per criterion by a ClinGen Variant Curation Expert Panel (VCEP) (3-star) |
| **B** | Not formally verified by ClinGen, but the paper's authors explicitly state the rationale (ACMG codes, etc.) |
| **C** | A real, existing disagreement (conflicting interpretations, etc.) is itself the answer. There is no single settled classification |

## 2. ClinGen registration status — overall summary

| Case | Variant | ClinVar registration | ClinGen VCEP registration | Tier |
|---|---|---|---|---|
| 1 | MYBPC3 c.278delA (p.Lys93ArgfsTer3) | Not registered (novel) | ❌ | B |
| 1 | KCNJ5 c.464G>A (p.Arg155Gln) | Not registered (novel) | ❌ (no corresponding VCEP) | B |
| 2 | MYBPC3 c.2905+1G>A | ✅ VCV000042666 | ❌ | B |
| 2 | MYBPC3 c.836del (p.Gly279Valfs*21) | Not registered (novel) | ❌ | B |
| **3** | **MYH7 c.2155C>T (p.Arg719Trp)** | ✅ VCV000014104 | **✅ ClinGen Inherited Cardiomyopathy Expert Panel** | **A** |
| 3 | MYBPC3 c.1000G>A (p.Glu334Lys) | ✅ (conflicting) | ❌ | C |
| 3 | MYH7 c.3382G>A (p.Ala1128Thr, noise) | ✅ (Likely benign) | ✅ ClinGen Inherited Cardiomyopathy Expert Panel (BS1) | **A** |
| 4 | DSG2 c.1592T>G (p.Phe531Cys) | ✅ (conflicting) | ❌ | C |

**Key finding**: The ClinGen Cardiomyopathy VCEP currently has MYH7 in Phase 1 (complete), while MYBPC3, TNNT2, TNNI3, ACTC1, etc. are in Phase 2 (evaluating applicability). **Of the 4 cases, only one variant — the MYH7 variant in Case 3 — has received a formal ClinGen 3-star review.**

**Why we did not restrict ourselves to ClinGen-only variants**: The genes/variants ClinGen covers are only a small fraction, and reproducing the realistic situation of "multiple candidate variants in one case" using only ClinGen-reviewed variants is nearly impossible. In addition, ClinGen's own evidence base tends to be skewed toward data centered on European/American populations, which is a poor fit for the purpose of the Pan-Asian demonstration. We judged that **explicitly stating that the confidence of the ground truth itself has a tier structure is more faithful to the "preservation of provenance and uncertainty" principle emphasized in Dr. Fujiwara's paper** than obscuring it.

---

## 3. Case 1: MYBPC3 + KCNJ5 (phenotype-mismatch pattern)

**Source**: Han S, Zhang Y-Y, Geng J (2026). Case Report: A novel MYBPC3 gene variant in a Chinese patient with hypertrophic cardiomyopathy and apical ventricular aneurysm, with a concurrent novel KCNJ5 gene variant. *Front Cardiovasc Med* 13:1841777. doi:10.3389/fcvm.2026.1841777. Open access, CC BY. **Finalized against human curation (case_report.pptx) on 2026-09-08.**

**Overview**: A Chinese HCM patient with a concurrent apical ventricular aneurysm. In addition to a novel MYBPC3 variant, a novel KCNJ5 variant was incidentally found (a gene causing familial hyperaldosteronism type III / long QT syndrome type 13, unrelated to HCM). The clinicians themselves describe excluding Fabry disease (GLA) and cardiac amyloidosis via actual testing.

**Purpose**: A case designed to test whether the true causal candidate (MYBPC3) can be correctly prioritized without being misled by the phenotype-incongruent incidental finding (KCNJ5). **A real `vep_parser.py` issue (missing single-submitter/novel Pathogenic variants) is actually reproduced with this real data** (see Section 9).

### Ground Truth

| Variant | Criterion | Strength | Accept/Reject | Tier | Rationale |
|---|---|---|---|---|---|
| case1-var1 MYBPC3 c.278delA | PVS1 | Very Strong | Accept | B | Frameshift + premature stop codon, predicted NMD; LOF is an established disease mechanism for MYBPC3 (ClinGen SVI) |
| case1-var1 MYBPC3 c.278delA | PM2 | Moderate | Accept | B | Not registered in 1000G/ESP6500/ExAC/regional cohorts |
| case1-var2 KCNJ5 c.464G>A | PM2 | Moderate | Accept | B | ExAC 8.236e-06 only; undetected in 1000G/ESP6500 |
| case1-var2 KCNJ5 c.464G>A | PP3 | Supporting | Accept | B | SIFT/PolyPhen2/MutationTaster all "D"; REVEL=0.934 |
| case1-var2 KCNJ5 c.464G>A | PP4 | — | **Reject** | B | KCNJ5 is unrelated to the HCM phenotype. The paper reports the proband's son carries only the KCNJ5 variant and has no cardiac phenotype |
| case1-var3 MYH7 p.Ser1491Cys (noise) | BS1, BP4 | Strong, Supporting | Accept | B | ClinVar VCV000043020, Benign per multi-lab consensus |

**Discrepancy between the ACMG guideline table (Table 5) and the Tavtigian point system**: For MYBPC3 (PVS1_VeryStrong + PM2_Moderate), ACMG Table 5 gives **Likely Pathogenic** (1 Very Strong + 1 Moderate), while the Tavtigian point system reaches **+10 points, meeting the Pathogenic threshold**. A real example where the two methods can diverge. KCNJ5 (PM2_Moderate + PP3_Supporting) does not fall into any category in Table 5 and remains **VUS**; the Tavtigian method likewise gives +3, also uncertain.

**Final Classification**: case1-var1 (MYBPC3) = **Likely Pathogenic**. case1-var2 (KCNJ5) = **VUS**; PP4 is rejected due to phenotype mismatch. **Causal variant = case1-var1 (MYBPC3)**

---

## 4. Case 2: MYBPC3 compound heterozygous (dual-answer pattern)

**Source**: Wang J, Hong L, Li Y, Mao Z, Zhu Y, Qi M, Zhou R, Hong X (2025). Case Report: Lethal neonatal hypertrophic cardiomyopathy from compound heterozygous MYBPC3 variants. *Front Cardiovasc Med*. doi:10.3389/fcvm.2025.1726463. CC BY. **Finalized against human curation (case_report.pptx) on 2026-09-08.**

**Overview**: A molecular-autopsy case of lethal neonatal HCM. A paternally inherited splice-site variant (c.2905+1G>A, exon 27, known, ClinVar VCV000042666) and a maternally inherited novel frameshift variant (c.836del, exon 8, novel) — both loss-of-function, compound heterozygous. Confirmed via trio-WES; Western blot failed to detect the protein (null phenotype).

**Purpose**: Rather than the "narrow down to one answer" pattern, this tests a different pattern from Cases 1 and 3: **"both variants are required together as the cause."**

### Ground Truth

| Variant | Criterion | Strength | Accept/Reject | Tier | Rationale |
|---|---|---|---|---|---|
| case2-var1 MYBPC3 c.2905+1G>A (paternal, exon 27) | PVS1 | Very Strong | Accept | B | Disrupts the canonical splice-donor site; LOF is an established disease mechanism |
| case2-var1 | PS4 | Strong | Accept | B | ClinVar (VCV000042666) cites a case-control comparison paper (PMID:27532257). **Seven clinical diagnostic labs independently evaluated it, all in agreement (Pathogenic/Likely Pathogenic)** |
| case2-var1 | PM2 | Moderate | Accept | B | gnomAD 0.0014%; UK Biobank 2 carriers; 1000 Genomes undetected |
| case2-var2 MYBPC3 c.836del (maternal, exon 8) | PVS1 | Very Strong | Accept | B | Premature stop codon, predicted NMD; loss of function directly supported by failed protein detection on Western blot |
| case2-var2 | PM2 | Moderate | Accept | B | Undetected in gnomAD/UK Biobank/1000 Genomes; not registered in ClinVar/dbSNP |
| case2-var3 TNNT2 (noise) | BS1, BP4 | Strong, Supporting | Accept | B | ClinVar VCV000188689 |

**ACMG guideline table vs. Tavtigian point system**: case2-var1 (PVS1_VeryStrong + PS4_Strong + PM2_Moderate) = **Pathogenic** per Table 5, and also Pathogenic at +14 Tavtigian points (the two methods agree). case2-var2 (PVS1_VeryStrong + PM2_Moderate) = **Likely Pathogenic** per Table 5, but reaches the **Pathogenic threshold** at +10 Tavtigian points (the same type of inter-method discrepancy seen for Case 1's MYBPC3 variant).

**Final Classification**: case2-var1 = **Pathogenic**; case2-var2 = **Likely Pathogenic** (not simply "Pathogenic" for both). **Causal variant = both case2-var1 and case2-var2** (this is not a design where a single winner is chosen)

---

## 5. Case 3: MYH7 + MYBPC3 p.Glu334Lys (Pan-Asian frequency pattern) ★ includes a formally ClinGen-reviewed variant

**Source**: MYH7 = ClinVar VCV000014104, ClinGen Inherited Cardiomyopathy Expert Panel (2016-12-15, 3-star). MYBPC3 = "The Need for Inclusive Genomic Research," *Circulation: Genomic and Precision Medicine*, DOI:10.1161/CIRCGEN.122.003736 (citing Tomar et al.'s SG10K analysis); corroborating source (Hong Kong cohort): PMC12123433. **The patient's clinical narrative is constructed** (the variant identification information is real).

**Overview**: The source is a Hong Kong Chinese HCM cohort study (53 cases, 13 P/LP variants, 21 VUS). Combines MYH7 p.Arg719Trp (Pathogenic, ClinGen-approved) as the correct answer with MYBPC3 p.Glu334Lys, a real example previously classified Pathogenic that is now VUS (conflicting) in ClinVar due to East Asian population data. **Finalized against human curation (case_report_1.pptx) on 2026-09-08.**

**Purpose**: The core Pan-Asian variant. Also usable as teaching material for the automation-level distinctions in ① (Evidence Preparation), since a formally ClinGen-approved established answer (MYH7) and a candidate not yet standardized (MYBPC3) coexist in the same case.

### Ground Truth

| Variant | Criterion | Strength | Accept/Reject | Tier | Rationale |
|---|---|---|---|---|---|
| case3-var1 MYH7 c.2155C>T | PS2 | Strong | Accept | **A** | ClinGen Evidence Repository. De novo occurrence confirmed via trio analysis (PMID:10957787) |
| case3-var1 | PS3 | Strong | Accept | **A (confidence caveat)** | Reduced hemodynamics induces adverse myocardial remodeling (PMID:24829265). **The human curator explicitly notes, "I don't have much knowledge of knock-in/knockout mouse models, so I'm not very confident in judging how reliable this is"** |
| case3-var1 | PS4 | Strong | Accept | **A** | Detected in 8 affected individuals across 6 unrelated Chinese families; not detected in unaffected non-carriers (PMID:19645038) |
| case3-var1 | PM1, PM2, PM5, PP1, PP3 | Moderate/Moderate/Moderate/Supporting/Supporting | Accept | **A** | Formally applied by the ClinGen Evidence Repository. PM5: a different amino-acid change at the same residue (p.Arg719Gln, Variation ID 14107) is Pathogenic per Expert Panel |
| case3-var2 MYBPC3 c.1000G>A | PS3 | Strong | Accept | **C (confidence caveat)** | UPS impairment by E334K cMyBPC alters cardiac ion-channel and Ca2+-handling protein levels. **The same lack of confidence regarding knock-in/knockout mouse findings is noted here** |
| case3-var2 | BS1 | — | Accept (when using regional data) | **C** | gnomAD overall AF=0.02368%, East Asian AF=0.3338% (~14x). Given the patient is Chinese, this favors the benign direction |
| case3-var2 | (overall) | — | — | **C** | PS3 (pathogenic-leaning) and BS1 (benign-leaning) conflict, resulting in VUS. **This variant had been classified Pathogenic when only the gnomAD overall AF was known** — a core real-world example of how adding regional data changes the call |
| case3-noise group (1-5) | various | — | Reject | B | Real, high-frequency Benign variants |
| case3-noise6 MYH7 c.3382G>A (p.Ala1128Thr) | BS1 | Strong | Reject | **A** | ClinGen Inherited Cardiomyopathy Expert Panel Evidence Repository (Kelly et al. 2018). BS1 applied based on gnomAD Latino 0.02% and East Asian 0.016%. **This is similarly high-frequency in both populations, and is not an example of the reverse pattern discussed in Section 3-4B** |

**Correction (2026-09-08, resolved via the corrected human curation v2)**: case3-var2's ClinVar breakdown (VCV000177902.53: Uncertain significance 9, Benign 2, Likely benign 3) initially appeared to exactly match Case 4's DSG2 breakdown, raising suspicion of a copy-paste error. **The curator subsequently confirmed, in the corrected version, that Case 4's breakdown is in fact different (Likely pathogenic 3, Uncertain significance 7, Likely benign 1).** The suspicion was resolved.

**Final Classification**: case3-var1 (MYH7) = **Pathogenic** (Tier A; confidence is unshaken even excluding PS3, since PS2 + PS4 alone already give 2 Strong criteria). case3-var2 (MYBPC3) = **VUS** (Tier C; given the confidence caveat on PS3, this may in practice lean somewhat toward BS1). **Causal variant = case3-var1 (MYH7)**

---

## 6. Case 4: DSG2 p.Phe531Cys (Japan/ToMMo, Pan-Asian frequency pattern)

**Source**: Arrhythmogenic right ventricular cardiomyopathy in a Japanese patient with a homozygous founder variant of DSG2 in the East Asian population. *Human Genome Variation*. PMC9360431. **Finalized against human curation (case_report_1.pptx) on 2026-09-08.**

**Overview**: A real Japanese ARVC patient. Homozygous DSG2 p.Phe531Cys variant. Markedly higher frequency than gnomAD overall (up to ~14x) in gnomAD East Asian, ToMMo, and HGVD alike. An East Asian founder variant.

**Purpose**: The Japan-side Pan-Asian case. Connects directly to the live ToMMo/NCBN/GEM-J queries available via TogoMCP.

### Ground Truth

| Variant | Criterion | Strength | Accept/Reject | Tier | Rationale |
|---|---|---|---|---|---|
| case4-var1 DSG2 c.1592T>G (homozygous) | PS3 | Strong | Accept (curator's own re-derivation) | **C (confidence caveat)** | A mouse Dsg2 p.Phe536Cys knock-in model (corresponding to human p.Phe531Cys) confirmed biventricular dilation, reduced LVEF, conduction abnormalities, and fibrosis in the homozygous state. **The human curator explicitly notes, "I don't have much knowledge of knock-in/knockout mouse models, so I'm not very confident in judging how reliable this is"** |
| case4-var1 | PS4 | Strong | Accept | **C (curator's own re-derivation)** | Family screening: only homozygous carriers showed the ARVC phenotype at full penetrance; heterozygous carriers were unaffected |
| case4-var1 | (overall) | — | — | **C** | The curator's own re-derivation (PS3 + PS4, 2 Strong) yields Pathogenic. **The actual ClinVar record (VCV000044283.23) is Likely pathogenic (3) / Uncertain significance (7) / Likely benign (1)** — the discrepancy is less severe than initially suspected before the correction: 3 of 11 submissions already call it Likely pathogenic, though the majority (7) remain Uncertain |
| case4-noise group | various | — | Reject | B | Real Benign variants (some from live TogoMCP queries) |

**Final Classification**: case4-var1 (DSG2) = **no single settled answer, but trending toward Pathogenic** (Tier C). **The curator's own re-derivation (PS3+PS4 → Pathogenic) is directionally consistent with part of the actual ClinVar record (3 of 11 submissions are Likely pathogenic), but the majority (7) remain Uncertain. Given the low confidence in PS3 (caveat regarding the mouse-model interpretation), the most accurate description of this variant is "trending toward Pathogenic but not yet settled." Causal variant = case4-var1 (DSG2), reported with an explicit note that confidence is limited**

---

## 7. Summary of Tier distribution

| Case | Tier A | Tier B | Tier C |
|---|---|---|---|
| 1 | — | var1 (MYBPC3), var3 (noise) | — |
| 2 | — | var1 & var2 (MYBPC3 x2), var3 (noise) | — |
| 3 | var1 (MYH7), noise6 (MYH7 BS1) | noise group | var2 (MYBPC3) |
| 4 | — | noise group | var1 (DSG2) |

**Only the two MYH7-related variants in Case 3 fall under Tier A.** When the comparison tool reports results, the design is to record not just a plain concordance rate but also **which Tier of ground truth each match/mismatch was against** (see Section 3-1 of the main plan document).

**A further breakdown of confidence (important)**: Even within Tier A and Tier C, **there are multiple places where the human curator explicitly states low confidence regarding PS3 (functional studies, especially those using knock-in/knockout mouse models)** (both MYH7 and MYBPC3 in Case 3; DSG2 in Case 4). Since the Tier system alone cannot capture this, we have individually annotated the relevant Evidence with a "confidence caveat" note. Case 4 in particular carries a double uncertainty: low confidence in PS3, plus a mismatch between the curator's own re-derivation (PS3+PS4 → Pathogenic) and the actual consensus recorded in ClinVar (leaning Uncertain).

---

## 8. Overview of the demo data (file composition, triage verification results)

| File | Content |
|---|---|
| `case{1-4}_clinical_note.txt` | Clinical narrative for each case (English). #1, #2, and #4 are recorded/summarized from real case reports; for #3 the variants are real but the patient narrative is constructed |
| `case{1-4}_variants.vcf` | Variant list for each case. In addition to the target variants, real benign noise variants gathered via live TogoMCP queries and ClinVar have been added, expanding each case to a scale (6-7 variants per case) also usable for testing the ① triage layer (`vep_parser.py`) |

**Triage measurement results** (`vep_parser.py`, unmodified main branch):

| Case | Number of variants (target + noise) | Measured result |
|---|---|---|
| 1 | 2 targets + 5 noise = 7 | Final verification using the correct data confirmed via human curation (case_report.pptx) (MYBPC3 = PVS1+PM2/Likely Pathogenic, KCNJ5 = PM2+PP3/VUS). **The `vep_parser.py` issue actually manifests**: KCNJ5 (VUS, Rank 5) ranks above the true causal candidate MYBPC3 (Likely Pathogenic, but novel so review status is unestablished, Rank 100) |
| 2 | 2 targets + 5 noise = 7 | Final verification reflecting the CLNREVSTAT values confirmed via human curation. **A difference emerged even within the same compound-heterozygous pair**: case2-var1 (known, 7-lab consensus) correctly ranks first at Rank 2, while case2-var2 (novel, review status unestablished) is stuck at Rank 100 due to the issue. Both, however, correctly rank above the noise (200) |
| 3 | 2 targets + 6 noise = 8 | MYH7 (the correct answer) correctly ranks first (Rank 2); MYBPC3 (the Pan-Asian candidate) also correctly ranks above the noise. The added MYH7 BS1 variant (Tier A) also correctly ranks lowest among the noise (Rank 200) |
| 4 | 1 target + 5 noise = 6 | DSG2 correctly ranks above the noise (Rank 9 vs. 200) |

**On verifying the response to the issue (finalized, 2026-09-08)**: This item went through multiple corrections in a single day. ① Initially, before human curation was available, MYBPC3 was recorded as PVS1+PM2/Likely Pathogenic based on inference, and the issue's reproduction was confirmed → ② A misreading of a search snippet (mistaking a sentence about KCNJ5 for one about MYBPC3) led to an incorrect correction to PM2+PP3/VUS, and an incorrect report that "the issue does not manifest" → ③ Human curation (case_report.pptx, based on reading the full text of the paper) definitively confirmed that ① had been correct. **Final conclusion: when verified with the correct data — MYBPC3 = Likely Pathogenic (PVS1+PM2), KCNJ5 = VUS (PM2+PP3) — the `vep_parser.py` issue does in fact manifest.** In Case 2 as well, we confirmed that only the novel variant within the same compound-heterozygous pair is affected by this issue. The response version (`vep_parser_bugfixed.py`) is kept in the reference materials as a candidate to propose to Francis on the day.

---

## 9. AlphaMissense / AlphaGenome annotation (added 2026-09-10)

The user independently annotated all variants with AlphaMissense (missense pathogenicity prediction) and AlphaGenome (splicing/regulatory-region impact prediction, raw score and PHRED score) (`annotation_alphamissense_alphagenome_v1.xlsx`). In this process, **the genomic coordinates (POS/REF/ALT) for all 28 variants across the 4 cases were actually verified, and at least 5 variants (case2-var1, case3-var1, case3-var2, case4-var1, case3-noise6) were found to differ from the ClinVar records and were corrected.** **This is exactly the concern I had previously and honestly flagged — that these were "approximate values for demo purposes, without rigorous re-mapping." This coordinate verification ended up covering all variants, so the corresponding open item in Section 11 is considered resolved.**

After the correction, annotation was successfully obtained for every variant. The variants for which annotation could not be obtained (the MYBPC3 deletion in case1-var1, and the MYBPC3 deletion in case2-var2) were, as expected, simply out of scope for AlphaMissense (which targets missense variants), since they are frameshifts/deletions.

**A notable finding**: for case1-var2 (KCNJ5 c.464G>A, designed as the phenotype-mismatch trap), AlphaMissense independently computed a high pathogenicity score of **am_pathogenicity=0.9485 (likely_pathogenic)**. Under the ACMG/ClinVar-based evaluation this remains VUS (PM2+PP3), but taken in isolation, AlphaMissense makes it look like "a variant that is likely pathogenic." **This actually strengthens the educational value of this case**: "how high a variant's own pathogenicity score is" and "whether it is the correct gene to explain this patient's phenotype (HCM)" are separate questions — no matter how high KCNJ5's AlphaMissense score is, it should not be adopted as the causal candidate, because it is a gene unrelated to HCM. This contrast is a concrete example reinforcing both the deterministic-vs-LLM comparison experiments (Section 3-4 of the plan) and the design philosophy of the Evidence Copilot (the need to independently evaluate consistency with the phenotype).

Re-verification (`vep_parser.py`) confirmed that the coordinate corrections do not affect the Rank output (since `vep_parser.py` uses only CLNSIG/CLNREVSTAT/frequency/IMPACT for its determination, not the genomic coordinates themselves).

---

## 10. Notes on demo-data creation

In the course of creating this document, we learned that **not just the ACMG evaluation values (Criterion/Strength/Classification), but the case's clinical narrative (clinical note) also requires the same level of verification.** The specific errors found are as follows.

| Location | Nature of the error | Cause | How it was discovered |
|---|---|---|---|
| Case 1 MYBPC3 ACMG evaluation | Incorrectly "corrected" PVS1+PM2/Likely Pathogenic to PM2+PP3/VUS | Misread a search snippet (mistook a sentence about KCNJ5 for one about MYBPC3) | Cross-check against human curation (case_report.pptx) |
| Case 2 MYBPC3 (maternal) ACMG evaluation | Incorrectly applied PM3 and PM6 | Independently supplemented criteria not stated in the paper as a "standard extension" | Cross-check against human curation |
| Case 2 MYBPC3 (paternal) ACMG evaluation | Overlooked PS4 | Insufficient reading of search snippets | Cross-check against human curation |
| **Case 1 clinical note** | **Fabricated the patient's age and clinical course (34 years old, 6-month history of exertional dyspnea → actually 63 years old, 7-day history of palpitations)** | Filled in a "plausible-sounding" clinical picture without verification | Direct search of the original paper's full text |
| **Case 2 clinical note** | **Content from an entirely different, unrelated paper (Alsters et al. 2019, a Dutch case) was mixed in** (maternal diabetes history, negative family history) | Confusion arising while searching multiple similar papers | Direct search of the original paper's full text |
| Case 3 MYH7/MYBPC3 and Case 4 DSG2 ACMG evaluation | The detailed breakdown of the PS-series (PS2/PS3/PS4/PM5, etc.) was insufficient (my own search alone only surfaced ClinVar's Pathogenic/Conflicting label) | Search snippets did not capture the detailed breakdown of the ACMG evaluation (which criteria applied and why) | Cross-check against human curation (case_report_1.pptx) |
| Case 3 MYBPC3 / Case 4 DSG2 gnomAD figures | Had used illustrative approximate values; the actual precise figures (gnomAD global 0.02368%, East Asian 0.3338%, etc.) came to light | Had not confirmed the precise figures via search | Cross-check against human curation |
| **Case 4 clinical note** | **The patient's age/sex (58-year-old man) and the parental consanguineous marriage (the very reason for the homozygous presentation) were completely missing** | Only summary-level information had been gathered from search snippets; no verification against the full text had been done | Direct search of the original paper's full text |

**Additional lessons (Cases 3-4, 2026-09-08)**:
5. **Evaluating functional studies (PS3), especially the validity of studies using knock-in/knockout mouse models, cannot be done with confidence without domain expertise.** The human curator explicitly noted, for Case 3 (both MYH7 and MYBPC3) and Case 4 (DSG2) alike, "I don't have much knowledge of knock-in/knockout mouse models, so I'm not very confident in judging how reliable this is" — this reveals a "variability in confidence per criterion" that the Tier system alone cannot express. It may be worth leaving room in the Evidence Record to record, separately from the Tier (provenance rank), **the inherent difficulty of evaluating a specific criterion type (such as functional studies)**.
6. **When the same breakdown numbers appear for two different variants, suspect a copy-paste error.** For Case 3's MYBPC3 and Case 4's DSG2, the ClinVar submitter breakdowns initially matched exactly, raising suspicion of a copy-paste error. **The curator confirmed in the corrected version (case_report_v2.pptx) that the breakdowns actually differ (DSG2: Likely pathogenic 3, Uncertain significance 7, Likely benign 1), and the suspicion was resolved.** Noticing the suspicious match and requesting confirmation led directly to the correct data.
7. **A human curator's own re-derivation (a calculation from ACMG criteria) and the consensus actually recorded in a database by the expert community do not necessarily agree** (this surfaced in Case 4). Rather than adopting one or the other as "the answer," we record both and treat the discrepancy itself as useful information.

**Lessons learned**:
1. **ACMG evaluation figures appear "conspicuously important" and so are handled carefully, but a narrative clinical description can seem fine as long as it "sounds plausible," which invites carelessness. However, since both are used as ground truth, both require the same level of verification.**
2. **Before adopting a search snippet, confirm — from the surrounding context — exactly which variant or paper it is describing.** Judging solely from a fragmentary match risks confusing it with a description of an adjacent, different subject.
3. **Because multiple papers address the same disease and the same gene (MYBPC3), when several similar cases appear mixed together in search results, cross-check the author names, publication year, and DOI each time to prevent conflation.**
4. Cross-checking against human review (in this case, case_report.pptx) was decisively effective in uncovering these errors. **On Day 1 as well, it is worth building in an opportunity for practitioners like Ruth and Francis to visually verify content rather than taking auto-generated content at face value.**

---

## 11. Open items / future work

- ClinGen registration status has not been individually checked for the noise variants (benign variants) — only the target variants have been confirmed
- How to handle the disease axis across the demo cases as a whole (HCM and ARVC are mixed) remains undecided
- The comparison tool's output specification (e.g., displaying concordance rates by Tier) is design-only and not yet implemented
