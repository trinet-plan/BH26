# Human Curation Results: ACMG Evaluation for Case 1-4

**version: v2** (translated from human_curation_v2_ja.md, originally human_curation_v2.pptx)

---

## Case 1: A novel MYBPC3 gene variant in a Chinese patient with hypertrophic cardiomyopathy and apical ventricular aneurysm, with a concurrent novel KCNJ5 gene variant

### Overview

Chinese HCM patient with concurrent apical ventricular aneurysm. In addition to a novel MYBPC3 variant, a novel KCNJ5 variant was incidentally found (a gene causing familial hyperaldosteronism type III / long QT syndrome type 13, unrelated to HCM).

Question: Can the true causal candidate (MYBPC3) be correctly identified without being misled by the phenotype-incongruent incidental finding (KCNJ5)?

### Clinical features

- Hypertrophic cardiomyopathy (HCM)
- MYBPC3 and MYH7 are the most common causal genes in HCM patients
- Left ventricular apical aneurysm in the absence of coronary artery disease is reported in <5% of HCM cases
- KCNJ5 variants often cause familial hyperaldosteronism type III or long QT syndrome type 13
- HCM is primarily characterized by left ventricular hypertrophy (LVH) in the absence of other cardiac, systemic, or metabolic disease
- MYH7 and MYBPC3 are most commonly implicated, accounting for the majority of variant-positive cases
- Other genes (TNNI3, TNNT2, TPM1, MYL2, MYL3, ACTC1, etc.) each contribute to 1-5% of patients

### Two variants detected

- Novel heterozygous variant in MYBPC3 exon 2: NM_000256.3(MYBPC3):c.278delA (p.Lys93ArgfsTer3)
- Novel heterozygous variant in KCNJ5 exon 2: NM_000256.3(KCNJ5):c.464G>A (p.Arg155Gln)

### About MYBPC3 c.278delA (p.Lys93ArgfsTer3)

- Population frequency database review shows this variant is rare and not reported in 1000 Genomes, ESP6500, or ExAC
- Also absent from a regional population database including both cardiomyopathy cases and controls
- This single-base deletion replaces the positively charged polar residue lysine at position 93, causing a frameshift during translation. A premature stop codon is introduced at the next amino acid position, predicted to trigger NMD (Nonsense-Mediated mRNA Decay)
- Downstream frameshift or nonsense variants at the same locus are consistently reported in ClinVar as pathogenic for hypertrophic cardiomyopathy and related cardiac conditions
- Per ACMG guidelines, this variant is classified as "Likely Pathogenic" for HCM, meeting the PVS1 and PM2 criteria

**ACMG evaluation:**

| Criterion | Rationale | ACMG strength | Tavtigian points |
|---|---|---|---|
| PVS1 (null variant in a gene where LOF is a known disease mechanism) | Single-base deletion causes frameshift + premature stop codon → NMD. ClinGen MYBPC3-specific recommendations state PVS1 is applicable since LOF is an established disease mechanism | Very Strong | +8 |
| PM2 (absent from population/control databases) | Not reported in 1000 Genomes, ESP6500, or ExAC | Moderate | +2 |

Per ACMG Table 5, 1 Very Strong + 1 Moderate = Likely Pathogenic. Per the Tavtigian point system, +10 points = Pathogenic. (ClinGen recommends downgrading PM2 from Moderate to Supporting.)

Reference: ACMG guidelines https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/ / ClinGen MYBPC3 specifications https://zenodo.org/records/21434332

### About KCNJ5 c.464G>A (p.Arg155Gln)

- Population frequency search confirms this is a rare variant (1000 Genomes: absent; ESP6500: absent; ExAC: 8.236e-06)
- Cross-prediction using multiple bioinformatics tools (including SIFT and Polyphen-2) consistently indicated a deleterious effect (SIFT: "D"; Polyphen2: "D"; MutationTaster_pred: "D"; VEST4 score: 0.953; REVEL score: 0.934; other tools: 7 D / 1 H)
- Suggests the resulting amino acid substitution may affect protein function
- Not found in ClinVar or HGMD
- Neighboring missense variants (e.g., c.451G>A p.Gly151Arg, c.452G>A p.Gly151Glu, c.470T>G p.Ile157Ser, c.473C>G p.Thr158Arg) are repeatedly recorded in ClinVar as pathogenic/likely pathogenic variants associated with long QT syndrome or aldosteronism

**ACMG evaluation:**

| Criterion | Rationale | ACMG strength | Tavtigian points |
|---|---|---|---|
| PM2 (absent from population/control databases) | Reported only in ExAC at 8.236e-06 | Moderate | +2 |
| PP3 (multiple computational lines of evidence support a deleterious effect) | SIFT/Polyphen2/MutationTaster all "D"; VEST4=0.953; REVEL=0.934; other tools 7D/1H → judged to adversely affect protein function | Supporting | +1 |

Per ACMG Table 5, Moderate 1 + Supporting 1 does not fall into a defined category, so this remains VUS. Per the Tavtigian point system, +3 points = uncertain (VUS). (Although reported in ExAC, the original paper's authors adopted PM2 based on rarity.)

### Conclusion

Adopt **Likely Pathogenic** MYBPC3 exon 2: c.278delA (p.Lys93ArgfsTer3).
Reject **VUS** KCNJ5 exon 2: c.464G>A (p.Arg155Gln).

---

## Case 2: Lethal neonatal hypertrophic cardiomyopathy from compound heterozygous MYBPC3 variants

### Overview

A molecular autopsy case of lethal neonatal HCM. A paternally inherited splice-site variant (c.2905+1G>A) and a maternally inherited novel frameshift variant (c.836del; p.Gly279Valfs*21) — both loss-of-function, compound heterozygous. Confirmed via trio-WES, interpreted per ACMG criteria. Both parents are asymptomatic carriers.

Question: Rather than "narrowing down to a single answer," can the system recognize that "both variants are required together as the cause"?

### Clinical features

A case report highlighting the clinical utility of postmortem molecular diagnosis. A two-month-old infant died of sudden-onset acute heart failure. Forensic autopsy was performed, and WES was conducted on the proband and parents. Autopsy revealed severe HCM, atrial septal defect (ASD), and extensive myocardial necrosis and fibrosis. WES identified variants in MYBPC3: a known paternal splice-site variant (c.2905+1G>A) and a novel maternal truncating frameshift variant (c.836del; p.Gly279Valfs*21). Both variants are predicted to result in complete loss of protein function.

### Two variants detected

- MYBPC3, exon 27, known, paternal, splice-site variant, c.2905+1G>A
- MYBPC3, exon 8, novel, maternal, frameshift deletion, c.836del (p.Gly279Valfs*21)

### About MYBPC3 c.2905+1G>A

Predicted to result in complete loss of function. Disrupts the canonical splice donor site, causing exon skipping or intron retention. LOF is a known disease mechanism for MYBPC3.

**ACMG evaluation:**

| Criterion | Rationale | ACMG strength | Tavtigian points |
|---|---|---|---|
| PVS1 (null variant) | Disrupts the canonical splice donor site. LOF is a known disease mechanism for MYBPC3 | Very Strong | +8 |
| PS4 (significantly increased prevalence in cases vs. controls) | ClinVar (VCV000042666, https://www.ncbi.nlm.nih.gov/clinvar/variation/42666/) cites a case-control comparison study for cardiomyopathy (PMID:27532257, https://pubmed.ncbi.nlm.nih.gov/27532257/) | Strong | +4 |
| PM2 (absent from population/control databases) | gnomAD: 0.0014%; UK Biobank: 2 carriers; 1000 Genomes: undetected. Has an ID in ClinVar/dbSNP (VCV000042666, rs397515991) | Moderate | +2 |

Per ACMG Table 5, 1 Very Strong + 1 Strong + 1 Moderate = Pathogenic. Per the Tavtigian point system, +14 points = Pathogenic. (ClinGen recommends downgrading PM2 to Supporting; although reported, the original authors adopted PM2 based on rarity.)

### About MYBPC3 c.836del (p.Gly279Valfs*21)

Predicted to result in complete loss of function. Introduces a premature termination codon (PTC) early in the coding sequence. Transcripts containing such a PTC are typically targeted for degradation by nonsense-mediated decay (NMD). This compound heterozygous state likely results in a near-total absence of functional cMyBP-C protein (null phenotype). Western blot confirmation of the protein failed.

**ACMG evaluation:**

| Criterion | Rationale | ACMG strength | Tavtigian points |
|---|---|---|---|
| PVS1 (null variant) | Frameshift/PTC predicted to trigger NMD. LOF is a known disease mechanism for MYBPC3 | Very Strong | +8 |
| PM2 (absent from population/control databases) | gnomAD: undetected; UK Biobank: undetected; 1000 Genomes: undetected. No ClinVar/dbSNP ID | Moderate | +2 |

Per ACMG Table 5, 1 Very Strong + 1 Moderate = Likely Pathogenic. Per the Tavtigian point system, +10 points = Pathogenic. (ClinGen recommends downgrading PM2 to Supporting.)

### Conclusion

Adopt both variants.
**Pathogenic**: MYBPC3 exon 27: c.2905+1G>A
**Likely Pathogenic**: MYBPC3 exon 8: c.836del (p.Gly279Valfs*21)

---

## Case 3: Genetic landscape of hypertrophic cardiomyopathy in Hong Kong Chinese population

### Overview

MYBPC3 p.Glu334Lys, previously classified Pathogenic, was found in 13/4810 (0.03%) of the SG10K cohort; ClinVar currently lists it as having conflicting interpretations. An example of ancestry-matched reference data revealing a variant's true frequency. The same variant has also been reported as VUS in 3 patients in a Hong Kong Chinese cohort study, corroborating the finding across multiple East Asian cohorts.

Finding: "Appears rare in gnomAD overall, but is actually non-negligible in frequency among Asian populations."

(The clinical-features slide is essentially the same as Case 1 and is skipped.)

### About the target variants

This paper sequenced 53 HCM patients and classified detected variants per ACMG guidelines: 13 P/LP variants and 21 VUS were found. Of these, we evaluate the P/LP variant MYH7 p.Arg719Trp (ClinGen-approved Pathogenic) and the VUS MYBPC3 p.Glu334Lys per ACMG guidelines.

### MYH7 p.Arg719Trp

Evaluated per guidelines with reference to ClinGen (https://erepo.clinicalgenome.org/evrepo/ui/classification/7b17a8bd-b169-46a1-8efa-b7388c4ecdef?version=1.0)

**ACMG evaluation:**

| Criterion | Rationale | ACMG strength |
|---|---|---|
| PS2 (de novo, no family history) | Trio analysis in the cited paper detected a de novo MYH7 p.Arg719Trp variant (PMID:10957787, https://pubmed.ncbi.nlm.nih.gov/10957787/) | Strong |
| PS3 (well-established in vitro/in vivo functional study showing a deleterious effect) | The cited paper reports that reduced hemodynamics induces adverse myocardial remodeling via both sympathetic nervous system and renin-angiotensin-aldosterone system stimulation (PMID:24829265, https://pubmed.ncbi.nlm.nih.gov/24829265/). **⚠️ We do not have much knowledge of knock-in/knockout mouse models, so we are not very confident in judging how reliable this is...** | Strong |
| PS4 (significantly increased prevalence vs. controls) | The cited paper studied 6 unrelated Chinese families with autosomal dominant HCM and detected the variant in 8 affected individuals; it was not detected in unaffected non-carriers (PMID:19645038, https://pubmed.ncbi.nlm.nih.gov/19645038/) | Strong |

Also applicable: PP3 (tools predict damaging), PM2 (absent from ExAC), PM5 (c.2156G>A / p.Arg719Gln — Variation ID 14107 — Pathogenic by Expert Panel), PP1, PM1.

Per ACMG guidelines, 2 or more Strong criteria among PS1-4 → Pathogenic.

### MYBPC3 p.Glu334Lys

ClinVar: "Conflicting classifications of pathogenicity" — Uncertain significance (9); Benign (2); Likely benign (3) (VCV000177902.53, https://www.ncbi.nlm.nih.gov/clinvar/variation/177902/)

**ACMG evaluation:**

| Criterion | Rationale | ACMG strength |
|---|---|---|
| PS3 (deleterious effect on gene product) | UPS impairment by E334K cMyBPC alters protein levels of cardiac ion channels and Ca2+-handling proteins, increasing calcium transient amplitude and causing electrophysiological dysfunction (https://www.sciencedirect.com/science/article/abs/pii/S002228361100996X ). **⚠️ We do not have much knowledge of knock-in/knockout mouse models, so we are not very confident in judging how reliable this is...** | Strong |
| BS1 (allele frequency greater than expected for the disorder) | gnomAD global AF=0.0002368; gnomAD East Asian AF=0.003338. Given the relatively high frequency in East Asians and that the affected patient is Chinese, this is taken to favor the benign direction | — |

Per ACMG guidelines, the pathogenic-leaning PS3 and benign-leaning BS1 conflict, resulting in a VUS classification. **Note: when only the gnomAD global AF was known, this variant had been classified Pathogenic.**

---

## Case 4: Arrhythmogenic right ventricular cardiomyopathy in a Japanese patient with a homozygous founder variant of DSG2 in the East Asian population

### Overview

A case report of a Japanese ARVC (arrhythmogenic right ventricular cardiomyopathy) patient carrying a homozygous DSG2 p.Phe531Cys variant. HGMD classifies it as "disease-causing," while the actual population frequency is clearly higher in East Asian populations, as detailed below. This can be queried via ToMMo/NCBN/GEM-J through TogoMCP.

### Clinical features

Arrhythmogenic right ventricular cardiomyopathy (ARVC) is a hereditary cardiomyopathy causing fatal arrhythmias and heart failure. It is one of the leading causes of arrhythmic cardiac arrest in young people and athletes, with an estimated prevalence of 1 in 1000-5000. The most common clinical presentation is palpitations or effort-induced syncope in adolescents or young adults, with T-wave inversion in the right precordial leads on ECG, ventricular arrhythmias, and right ventricular abnormalities on imaging. A homozygous DSG2 variant presumed to be an East Asian founder variant, c.1592T>G (p.Phe531Cys), was identified.

### Variant detected

Homozygous DSG2 variant c.1592T>G (p.Phe531Cys)

### DSG2 p.Phe531Cys

ClinVar: "Conflicting classifications of pathogenicity" — Likely pathogenic (3); Uncertain significance (7); Likely benign (1) (VCV000044283.23, https://www.ncbi.nlm.nih.gov/clinvar/variation/44283/)

**ACMG evaluation:**

| Criterion | Rationale | ACMG strength |
|---|---|---|
| PS3 (deleterious effect on gene product) | A mouse Dsg2 p.Phe536Cys knock-in model (corresponding to human p.Phe531Cys) was generated, comparing wild-type, heterozygous, and homozygous mice. Homozygous mice showed biventricular dilation, reduced left ventricular systolic function, ECG conduction abnormalities, myocyte dropout, and replacement/interstitial/perivascular fibrosis in both ventricles (https://link.springer.com/article/10.1186/s12916-024-03593-8 ). **⚠️ We do not have much knowledge of knock-in/knockout mouse models, so we are not very confident in judging how reliable this is...** | Strong |
| PS4 (significantly increased prevalence vs. controls) | Family screening showed that only homozygous carriers exhibited a definite ARVC phenotype with 100% penetrance, while heterozygous carriers were unaffected (Chen L et al. 2019, https://doi.org/10.1016/j.ijcard.2018.06.105 ) | Strong |

Per ACMG guidelines, 2 or more Strong criteria among PS1-4 → Pathogenic.

---

## Appendix

### ACMG guideline text (Richards et al. 2015)

**PVS1** (Very strong evidence of pathogenicity): Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multi-exon deletion) in a gene where loss of function (LOF) is a known mechanism of disease

Caveats:
- Beware of genes where LOF is not a known disease mechanism (**e.g. GFAP, MYH7**) — **notably, the guideline itself explicitly warns against using PVS1 for MYH7**
- Use caution interpreting LOF variants at the extreme 3' end of a gene
- Use caution with splice variants that are predicted to lead to exon skipping but leave the remainder of the protein intact
- Use caution in the presence of multiple transcripts

**PM2**: Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes, or ExAC

Caveat: Population data for indels may be poorly called by next generation sequencing

**PP3**: Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc)

Caveat: As many in silico algorithms use the same or very similar input for their predictions, each algorithm should not be counted as an independent criterion. PP3 can be used only once in any evaluation of a variant.

### Table 5: Rules for Combining Criteria to Classify Sequence Variants

**Pathogenic**
1. 1 Very Strong (PVS1) AND
   a. ≥1 Strong (PS1-PS4), OR
   b. ≥2 Moderate (PM1-PM6), OR
   c. 1 Moderate (PM1-PM6) and 1 Supporting (PP1-PP5), OR
   d. ≥2 Supporting (PP1-PP5)
2. ≥2 Strong (PS1-PS4), OR
3. 1 Strong (PS1-PS4) AND
   a. ≥3 Moderate (PM1-PM6), OR
   b. 2 Moderate (PM1-PM6) AND ≥2 Supporting (PP1-PP5), OR
   c. 1 Moderate (PM1-PM6) AND ≥4 Supporting (PP1-PP5)

**Likely Pathogenic**
1. 1 Very Strong (PVS1) AND 1 Moderate (PM1-PM6), OR
2. 1 Strong (PS1-PS4) AND 1-2 Moderate (PM1-PM6), OR
3. 1 Strong (PS1-PS4) AND ≥2 Supporting (PP1-PP5), OR
4. ≥3 Moderate (PM1-PM6), OR
5. 2 Moderate (PM1-PM6) AND ≥2 Supporting (PP1-PP5), OR
6. 1 Moderate (PM1-PM6) AND ≥4 Supporting (PP1-PP5)

**Benign**
1. 1 Stand-Alone (BA1), OR
2. ≥2 Strong (BS1-BS4)

**Likely Benign**
1. 1 Strong (BS1-BS4) and 1 Supporting (BP1-BP7), OR
2. ≥2 Supporting (BP1-BP7)

### Tavtigian point-based (Bayesian) system

**Point values by ACMG/AMP evidence strength category**

| Evidence Strength | Pathogenic | Benign |
|---|---|---|
| Indeterminate | 0 | 0 |
| Supporting | 1 | -1 |
| Moderate | 2 | -2 |
| Strong | 4 | -4 |
| Very Strong | 8 | -8 |

**Point-based classification categories**

| Category | Point ranges |
|---|---|
| Pathogenic | ≥10 |
| Likely Pathogenic | 6-9 |
| Uncertain | 0-5 |
| Likely Benign | -1 to -6 |
| Benign | ≤-7 |
