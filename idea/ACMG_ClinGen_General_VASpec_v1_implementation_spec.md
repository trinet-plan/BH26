# ACMG / ClinGen General Guidance → GA4GH VA-Spec v1 実装仕様

## 1. 適用範囲

本仕様は、疾患・遺伝子固有の ClinGen VCEP Criteria Specification（CSpec）を適用せず、
以下の優先順位で ACMG criterion を評価する。

```text
ACMG/AMP 2015
    ↓
ClinGen General Guidance
    ↓
Implementation Rule（閾値・DB・自動処理）
```

- 分類体系: ACMG/AMP 2015 の 28 criteria
- ClinGen: 疾患非依存の General Guidance
- CSpec: **適用しない**
- 出力: **GA4GH VA-Spec v1.0 系 ACMG 2015 aligned profile**
- 最終出力: `Variant Pathogenicity Statement (ACMG 2015)`
- Criterion単位: `Variant Pathogenicity Evidence Line (ACMG 2015)`

## 2. VA-Spec v1 出力規則

VA-Spec v1 の ACMG EvidenceLine では `methodType` に criterion code
（`PVS1` ～ `BP7`）を直接設定する。

- Pathogenic criterion が成立: `directionOfEvidenceProvided = "supports"`
- Benign criterion が成立: `directionOfEvidenceProvided = "disputes"`
- 評価済みだが不成立: `directionOfEvidenceProvided = "none"`
- `none` の場合は `strengthOfEvidenceProvided` を出力しない

`evidenceOutcome`:

```text
成立 + default strength        → <criterion>
成立 + strength変更             → <criterion>_<strength>
評価済み + 不成立                → <criterion>_not_met
```

例:

```text
PVS1
PVS1_strong
PM2_supporting
PS3_moderate
BP4_strong
PP3_not_met
```

`NOT_EVALUATED / NOT_APPLICABLE / MANUAL_REVIEW` は内部 workflow status とし、
VA-Spec の `evidenceOutcome` に独自コードを追加しない。必要に応じて Extension / local metadata に保持する。

## 3. 28 criteria 実装マトリクス

| Criterion | ACMG default / General | ClinGen General | 具体的処理条件 | 主なDB/API | 自動化 | VA-Spec v1 evidenceOutcome | 推奨 Evidence Item |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PVS1 | very strong | Use the ClinGen PVS1 decision tree; consider transcript relevance, NMD, exon relevance, alternative initiation, and the importance/extent of the truncated region. Strength may be downgraded. | Eligible LoF consequence → confirm gene–disease LoF mechanism → select relevant transcript → predict NMD → if NMD escape, assess critical region/protein loss → assign VS/S/M/P per PVS1 decision tree. If mechanism is unknown, do not auto-fire. | VEP; MANE/RefSeq; ClinGen Gene-Disease Validity; UniProt/InterPro | Semi-automatic | PVS1 if default VS; otherwise PVS1_<strength>; PVS1_not_met when evaluated and not met | Generic StudyResult for annotation/NMD/mechanism evidence; Statement may be used for curated gene-disease mechanism. No dedicated PVS1 StudyResult profile in VA-Spec v1.0. |
| PS1 | strong | For splice-related variants, ClinGen splicing guidance allows comparison of equivalent predicted/observed splice effects; avoid double counting with PVS1/PP3. | Missense: same residue + same alternate amino acid + different nucleotide change + comparator is sufficiently established Pathogenic; exclude cases where comparator pathogenicity is mainly due to a different splice effect. Splicing: evaluate equivalent splice consequence under ClinGen guidance. | ClinVar; NCBI Variation; ClinGen Evidence Repository; VEP; PubMed | Semi-automatic | PS1 if default strong; otherwise PS1_<strength>; PS1_not_met | Statement describing comparator pathogenicity plus generic StudyResult for amino-acid/splice equivalence. |
| PS2 | strong | Use ClinGen PS2/PM6 de novo point framework; points depend on parental relationship confirmation and phenotype consistency; aggregate independent probands to determine strength. | Proband carries ALT; both parents do not; parental relationships confirmed; genotype QC passes → assign proband points based on phenotype consistency → sum independent probands → map total to evidence strength. | Internal family/case data; HPO; MONDO; ClinGen Gene-Disease Validity | Semi-automatic | PS2 if default strong; otherwise PS2_<strength>; PS2_not_met | Generic StudyResult for trio/genotype data plus Statement/StudyResult for phenotype assessment. |
| PS3 | strong | ClinGen functional-evidence guidance evaluates assay relevance, controls, validation variants, replication/statistics and allows strength adjustment. | Require curated assay record → confirm assay measures a disease-relevant mechanism → assess positive/negative controls and validation → damaging result → assign strength from assay validity/quality. Do not infer assay validity from a PMID alone. | MaveDB; PubMed; ClinGen Evidence Repository | Curated / semi-automatic | PS3 if default strong; otherwise PS3_<strength>; PS3_not_met | ExperimentalVariantFunctionalImpactStudyResult and/or a functional-impact Statement supported by that StudyResult. |
| PS4 | strong | No disease-agnostic single threshold is adopted here; case-control statistics and/or independent-case evidence must be evaluated in disease context. | If case-control data exist, calculate/ingest effect estimate and confidence interval; require enrichment supporting association. For rare variants without adequate case-control data, use a documented case-count rule only when justified by the applicable general/local method. | PubMed; case-control datasets; internal case database | Semi-automatic / manual | PS4 if default strong; otherwise PS4_<strength>; PS4_not_met | Generic StudyResult containing case-control/case-count statistics; source Document for study. |
| PM1 | moderate | No universal disease-agnostic hotspot window is specified; domain/hotspot evidence must be justified from curated biology/variant distribution. | Variant lies in a curated critical domain/hotspot + pathogenic variation is enriched + benign variation is absent/limited → PM1. If only generic domain overlap is known, route to review rather than auto-fire. | UniProt; InterPro; ClinVar; VEP | Semi-automatic | PM1 if default moderate; otherwise PM1_<strength>; PM1_not_met | Generic StudyResult for domain/hotspot annotation and local variant spectrum; supporting Statements/Documents. |
| PM2 | moderate (ACMG); generally supporting under ClinGen recommendation | ClinGen recommends using PM2 at Supporting strength in general and emphasizes reliable population data/rarity rather than naïve absence. | Use reliable ancestry-aware population frequency/FAF. If absent or sufficiently rare under the configured general implementation rule, emit PM2_supporting. If coverage/AN/region quality is insufficient, mark internal status NOT_EVALUATED rather than not-met. | gnomAD; NCBI ALFA | Automatic with quality gates | PM2_supporting (because strength is adjusted from ACMG default moderate); PM2_not_met | CohortAlleleFrequencyStudyResult; use ancillaryResults/qualityMeasures for FAF and coverage/quality metadata as needed. |
| PM3 | moderate | ClinGen PM3 point framework weights phase confirmation, classification of the other allele, homozygous observations and independent probands; total points map to strength. | For each proband: determine phase and classification of second allele → assign points → cap/reconcile repeated evidence as required → sum independent probands → map to Supporting/Moderate/Strong/Very Strong. | Internal family data; ClinVar; ClinGen Evidence Repository | Semi-automatic | PM3 if default moderate; otherwise PM3_<strength>; PM3_not_met | Generic StudyResult for phase/genotypes plus Statement for classification of the second allele. |
| PM4 | moderate | Avoid double counting with PVS1; interpret protein-region relevance and repetitive sequence context. | Consequence is in-frame insertion/deletion or stop-loss + changes protein length + not solely in non-functional repetitive region + same molecular consequence is not already counted under PVS1 → PM4. | VEP; UniProt; InterPro; RepeatMasker | Automatic / semi-automatic | PM4 if default moderate; otherwise PM4_<strength>; PM4_not_met | Generic StudyResult for consequence, length change, repeat/domain annotation. |
| PM5 | moderate | Comparator pathogenicity and mechanism must be sufficiently established; same amino-acid substitution belongs under PS1, not PM5. | Evaluated variant is missense + same residue has a different Pathogenic missense + comparator has adequate evidence + not same AA substitution → PM5. | ClinVar; ClinGen Evidence Repository; PubMed | Semi-automatic | PM5 if default moderate; otherwise PM5_<strength>; PM5_not_met | Statement describing comparator pathogenic variant plus generic StudyResult for residue equivalence. |
| PM6 | moderate | Evaluated jointly with PS2 using the ClinGen de novo point framework; parental relationship status and phenotype specificity control points/strength. | Proband carries ALT + available parent(s) do not + parental relationship not fully confirmed → PM6 candidate → phenotype/relationship points → aggregate independent probands → strength. | Internal family data; HPO; MONDO | Semi-automatic | PM6 if default moderate; otherwise PM6_<strength>; PM6_not_met | Generic StudyResult for trio/family observations and phenotype evidence. |
| PP1 | supporting | ClinGen PP1/BS4/PP4 guidance recommends quantitative/structured segregation assessment and integration with phenotype/locus specificity; strength may increase. | Build informative pedigree observations → account for inheritance, penetrance, phenocopy and locus heterogeneity → calculate segregation LR/points or equivalent structured measure → map to evidence strength. | Internal pedigree data; HPO; MONDO | Semi-automatic / curator | PP1 if default supporting; otherwise PP1_<strength>; PP1_not_met | Generic StudyResult for segregation observations; Statement/StudyResult for phenotype specificity. |
| PP2 | supporting | Gene-level constraint alone is insufficient; disease mechanism and pathogenic variant spectrum must support missense causality. | Variant is missense + established gene-disease mechanism includes pathogenic missense + benign missense variation is low relative to expectation → PP2. If disease mechanism is unknown, do not auto-fire. | ClinGen Gene-Disease Validity; gnomAD constraint; ClinVar | Semi-automatic | PP2; PP2_not_met | Statement for disease mechanism plus generic StudyResult for gene-level missense constraint/variant spectrum. |
| PP3 | supporting | Use calibrated predictors rather than simple majority voting; ClinGen provides calibrated evidence strengths for missense prediction and separate splicing guidance. | Select predictor appropriate to variant mechanism → compare score with a documented calibrated interval → emit PP3 at the calibrated strength. Keep protein-impact and splice-impact evidence distinct and reconcile double counting. | REVEL; SpliceAI; VEP; other calibrated predictor source | Automatic | PP3 if supporting; PP3_moderate / PP3_strong when adjusted; PP3_not_met | Generic StudyResult for computational predictor score(s); no dedicated computational StudyResult profile in VA-Spec v1.0. |
| PP4 | supporting | ClinGen guidance explicitly links PP4 phenotype specificity with locus heterogeneity and segregation evidence (PP1/BS4); disease-specificity must be demonstrated. | Normalize phenotype → compare with candidate disease/gene phenotype → assess specificity and genetic heterogeneity → combine with segregation context as required → assign evidence strength or route to curator. | HPO; MONDO; ClinGen Gene-Disease Validity | Semi-automatic / curator | PP4 if default supporting; otherwise PP4_<strength>; PP4_not_met | Statement or generic StudyResult capturing phenotype-to-disease/gene match and specificity assessment. |
| PP5 | supporting | ClinGen recommends that PP5 not be used; retrieve and evaluate primary evidence instead. | Never produce MET under the General ClinGen policy. Use external classifications only to discover primary evidence. | ClinVar for discovery only | Disabled | Do not emit a MET EvidenceLine. If the implementation explicitly evaluates it as an ACMG criterion and records not-met, PP5_not_met is schema-compatible; deprecation itself should be kept in Extension/local workflow metadata. | Normally none. |
| BA1 | stand alone | ClinGen BA1 recommendations emphasize ancestry/population context, adequate allele number, data quality, and an exception list. | Use reliable population AF/FAF → require adequate observed alleles → if frequency >0.05 and not a BA1 exception → BA1. If frequency quality is unreliable, NOT_EVALUATED internally. | gnomAD; ClinGen BA1 exception information | Automatic with quality gates | BA1; BA1_not_met | CohortAlleleFrequencyStudyResult with sourceDataSet/cohort and relevant quality metadata. |
| BS1 | strong | Requires comparison with the maximum credible allele frequency expected for the disease; depends on prevalence, inheritance, penetrance and allelic/genetic contribution. | Compute/obtain disease-specific maximum credible AF → compare ancestry-aware observed AF/FAF → if observed exceeds threshold, BS1. Without disease-level assumptions/threshold, do not claim a general automatic BS1. | gnomAD; epidemiologic/curated disease data | Semi-automatic | BS1 if default strong; otherwise BS1_<strength>; BS1_not_met | CohortAlleleFrequencyStudyResult plus Statement/StudyResult carrying disease-frequency assumptions. |
| BS2 | strong | Interpretation depends on expected penetrance, age of onset, genotype/zygosity and reliable phenotype ascertainment. | Confirm disease expected to be highly penetrant by observed age → identify genotype incompatible with disease in well-phenotyped healthy individual(s) → assess zygosity/inheritance → BS2. Population counts alone require caution. | Well-phenotyped cohorts; PubMed; population datasets as supportive context | Semi-automatic | BS2 if default strong; otherwise BS2_<strength>; BS2_not_met | Generic StudyResult for healthy-individual genotype/phenotype observations; population StudyResult may be supporting context. |
| BS3 | strong | Same ClinGen functional-evidence framework as PS3; assay validity and controls determine whether/at what strength benign evidence is appropriate. | Require curated disease-relevant validated assay → normal/non-damaging result with appropriate controls/validation → assign strength from assay validity. | MaveDB; PubMed; ClinGen Evidence Repository | Curated / semi-automatic | BS3 if default strong; otherwise BS3_<strength>; BS3_not_met | ExperimentalVariantFunctionalImpactStudyResult and/or functional-impact Statement. |
| BS4 | strong | ClinGen PP1/BS4/PP4 guidance requires structured segregation interpretation and consideration of penetrance, phenocopy and locus heterogeneity. | Identify genotype-phenotype inconsistency → test whether reduced penetrance, phenocopy, age-dependent onset, heterogeneity or genotype error can explain it → if not, quantify non-segregation evidence → BS4. | Internal pedigree data; HPO/MONDO as context | Semi-automatic / curator | BS4 if default strong; otherwise BS4_<strength>; BS4_not_met | Generic StudyResult for segregation/non-segregation observations. |
| BP1 | supporting | Disease mechanism and known pathogenic variant spectrum are required; gene constraint alone should not define the rule. | Variant is missense + established disease mechanism is predominantly LoF/truncating + pathogenic missense mechanism is not established → BP1. | ClinGen Gene-Disease Validity; ClinVar | Semi-automatic | BP1; BP1_not_met | Statement for disease mechanism plus generic StudyResult for pathogenic variant spectrum. |
| BP2 | supporting | Interpret phase in the context of inheritance and penetrance; require reliable phasing and comparator classification. | Obtain two variants + phase + inheritance → if configuration satisfies BP2 logic and second variant is sufficiently established pathogenic → BP2. Unknown phase → do not auto-fire. | Internal family data; ClinVar; ClinGen Evidence Repository | Semi-automatic | BP2; BP2_not_met | Generic StudyResult for phase/genotypes plus Statement for second-variant pathogenicity. |
| BP3 | supporting | Confirm that the repeat region lacks established functional importance; avoid conflict with PM4. | In-frame insertion/deletion + overlaps repetitive region + no established functional domain/critical region → BP3. If domain/repeat interpretation is uncertain, review. | VEP; RepeatMasker; UniProt; InterPro | Automatic / semi-automatic | BP3; BP3_not_met | Generic StudyResult for in-frame consequence, repeat and domain annotation. |
| BP4 | supporting | Use calibrated computational predictors instead of majority vote; benign calibrated intervals can support adjusted strengths. | Select calibrated predictor → compare score to benign calibration interval → assign BP4 strength → reconcile with splice/functional evidence and avoid applying benign prediction where incompatible with high-confidence null/splice evidence. | REVEL; SpliceAI; VEP | Automatic | BP4 if supporting; BP4_moderate / BP4_strong / BP4_very_strong when adjusted; BP4_not_met | Generic StudyResult for computational score(s); no dedicated computational StudyResult profile in VA-Spec v1.0. |
| BP5 | supporting | Alternative diagnosis must adequately explain the phenotype; dual diagnosis must be considered. | Identify alternate causal variant/diagnosis → assess whether it fully or near-fully explains phenotype → if evaluated variant lacks independent phenotype contribution, BP5; if dual diagnosis remains plausible, MANUAL_REVIEW internally. | Internal patient data; ClinVar; HPO; MONDO | Manual / semi-automatic | BP5; BP5_not_met | Statement describing alternative molecular diagnosis plus phenotype-match evidence. |
| BP6 | supporting | ClinGen recommends that BP6 not be used; retrieve and evaluate primary evidence instead. | Never produce MET under the General ClinGen policy. External benign classifications are evidence-discovery aids only. | ClinVar for discovery only | Disabled | Do not emit a MET EvidenceLine. If explicitly evaluated and recorded not-met, BP6_not_met is schema-compatible; deprecation belongs in Extension/local workflow metadata. | Normally none. |
| BP7 | supporting | ClinGen splicing guidance refines location/prediction/RNA considerations and warns against double counting across PVS1/PS1/PP3/BP4/BP7. | Synonymous/noncoding candidate appropriate for BP7 → outside splice-critical positions as defined by adopted rule → calibrated splice prediction supports no impact → conservation/RNA evidence does not suggest function → BP7. Reconcile with other splice criteria. | VEP; SpliceAI; PhyloP/PhyloCons; PubMed/RNA evidence | Automatic / semi-automatic | BP7; BP7_not_met | Generic StudyResult for splice prediction/conservation; Statement or functional StudyResult for RNA evidence where available. |

## 4. Criterion別 詳細

### PVS1

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Null variant (nonsense, frameshift, canonical ±1/2 splice, initiation codon, exon deletion, etc.) in a gene where loss of function is an established disease mechanism.
- **ClinGen General**: Use the ClinGen PVS1 decision tree; consider transcript relevance, NMD, exon relevance, alternative initiation, and the importance/extent of the truncated region. Strength may be downgraded.
- **具体的処理**: Eligible LoF consequence → confirm gene–disease LoF mechanism → select relevant transcript → predict NMD → if NMD escape, assess critical region/protein loss → assign VS/S/M/P per PVS1 decision tree. If mechanism is unknown, do not auto-fire.
- **Implementation default / 注意**: NMD rule may use PTC >50 nt upstream of final exon–exon junction as a practical rule. Gene constraint (e.g., pLI/LOEUF) is supporting context only, not a substitute for disease mechanism.
- **必要入力**: normalized variant; consequence; transcript/exon/CDS position; NMD context; gene–disease mechanism; protein region/domain; alternative transcript/initiation
- **DB/API**: VEP; MANE/RefSeq; ClinGen Gene-Disease Validity; UniProt/InterPro
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `PVS1`
- **MET時 direction**: `supports`
- **strength**: very strong / strong / moderate / supporting
- **`evidenceOutcome`**: PVS1 if default VS; otherwise PVS1_<strength>; PVS1_not_met when evaluated and not met
- **Evidence Item**: Generic StudyResult for annotation/NMD/mechanism evidence; Statement may be used for curated gene-disease mechanism. No dedicated PVS1 StudyResult profile in VA-Spec v1.0.
- **Curator確認点**: LoF disease mechanism, biologically relevant transcript, exon/region relevance.

### PS1

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Same amino-acid change as an established pathogenic variant, caused by a different nucleotide change.
- **ClinGen General**: For splice-related variants, ClinGen splicing guidance allows comparison of equivalent predicted/observed splice effects; avoid double counting with PVS1/PP3.
- **具体的処理**: Missense: same residue + same alternate amino acid + different nucleotide change + comparator is sufficiently established Pathogenic; exclude cases where comparator pathogenicity is mainly due to a different splice effect. Splicing: evaluate equivalent splice consequence under ClinGen guidance.
- **Implementation default / 注意**: Use primary evidence behind the comparator where possible; do not rely solely on an unreviewed database label.
- **必要入力**: protein/codon consequence; comparator variant; comparator classification/evidence; splice prediction/RNA evidence
- **DB/API**: ClinVar; NCBI Variation; ClinGen Evidence Repository; VEP; PubMed
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `PS1`
- **MET時 direction**: `supports`
- **strength**: strong or adjusted strength
- **`evidenceOutcome`**: PS1 if default strong; otherwise PS1_<strength>; PS1_not_met
- **Evidence Item**: Statement describing comparator pathogenicity plus generic StudyResult for amino-acid/splice equivalence.
- **Curator確認点**: Quality/independence of comparator pathogenic evidence and splice-effect equivalence.

### PS2

- **方向**: Pathogenic
- **ACMG/AMP 2015**: De novo variant in a patient with the disease and no family history, with maternity and paternity confirmed.
- **ClinGen General**: Use ClinGen PS2/PM6 de novo point framework; points depend on parental relationship confirmation and phenotype consistency; aggregate independent probands to determine strength.
- **具体的処理**: Proband carries ALT; both parents do not; parental relationships confirmed; genotype QC passes → assign proband points based on phenotype consistency → sum independent probands → map total to evidence strength.
- **Implementation default / 注意**: Typical genotype QC can require DP/GQ thresholds, but these are implementation policy, not ACMG/ClinGen universal thresholds.
- **必要入力**: trio genotypes; parental relationship confirmation; genotype QC; phenotype; family history; independent proband count
- **DB/API**: Internal family/case data; HPO; MONDO; ClinGen Gene-Disease Validity
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `PS2`
- **MET時 direction**: `supports`
- **strength**: supporting / moderate / strong / very strong per point framework
- **`evidenceOutcome`**: PS2 if default strong; otherwise PS2_<strength>; PS2_not_met
- **Evidence Item**: Generic StudyResult for trio/genotype data plus Statement/StudyResult for phenotype assessment.
- **Curator確認点**: Phenotype consistency and independence of observations.

### PS3

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Well-established functional studies supportive of a damaging effect on the gene or gene product.
- **ClinGen General**: ClinGen functional-evidence guidance evaluates assay relevance, controls, validation variants, replication/statistics and allows strength adjustment.
- **具体的処理**: Require curated assay record → confirm assay measures a disease-relevant mechanism → assess positive/negative controls and validation → damaging result → assign strength from assay validity/quality. Do not infer assay validity from a PMID alone.
- **Implementation default / 注意**: Prefer curated functional evidence or MAVE data with explicit assay metadata. Keep raw assay result separate from ACMG interpretation.
- **必要入力**: assay type/system/readout; controls; validation variants; result; statistics; replication; PMID/MaveDB ID
- **DB/API**: MaveDB; PubMed; ClinGen Evidence Repository
- **自動化**: Curated / semi-automatic
- **VA-Spec `methodType`**: `PS3`
- **MET時 direction**: `supports`
- **strength**: supporting / moderate / strong / very strong when justified
- **`evidenceOutcome`**: PS3 if default strong; otherwise PS3_<strength>; PS3_not_met
- **Evidence Item**: ExperimentalVariantFunctionalImpactStudyResult and/or a functional-impact Statement supported by that StudyResult.
- **Curator確認点**: Assay validity and relevance to disease mechanism.

### PS4

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Prevalence of the variant in affected individuals is significantly increased compared with controls.
- **ClinGen General**: No disease-agnostic single threshold is adopted here; case-control statistics and/or independent-case evidence must be evaluated in disease context.
- **具体的処理**: If case-control data exist, calculate/ingest effect estimate and confidence interval; require enrichment supporting association. For rare variants without adequate case-control data, use a documented case-count rule only when justified by the applicable general/local method.
- **Implementation default / 注意**: Do not treat ClinVar review status or number of submissions as PS4 evidence.
- **必要入力**: case/control AC/AN; OR/RR; CI; p-value; independent probands; phenotype; study design
- **DB/API**: PubMed; case-control datasets; internal case database
- **自動化**: Semi-automatic / manual
- **VA-Spec `methodType`**: `PS4`
- **MET時 direction**: `supports`
- **strength**: strong or adjusted strength
- **`evidenceOutcome`**: PS4 if default strong; otherwise PS4_<strength>; PS4_not_met
- **Evidence Item**: Generic StudyResult containing case-control/case-count statistics; source Document for study.
- **Curator確認点**: Study design, population matching, ascertainment, independence of cases.

### PM1

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Located in a mutational hotspot and/or critical and well-established functional domain without benign variation.
- **ClinGen General**: No universal disease-agnostic hotspot window is specified; domain/hotspot evidence must be justified from curated biology/variant distribution.
- **具体的処理**: Variant lies in a curated critical domain/hotspot + pathogenic variation is enriched + benign variation is absent/limited → PM1. If only generic domain overlap is known, route to review rather than auto-fire.
- **Implementation default / 注意**: fastVEP-style ±5 aa / ≥3 pathogenic / 0 benign is an implementation heuristic, not a ClinGen universal rule.
- **必要入力**: protein position; domain/hotspot coordinates; pathogenic/benign variants in region; domain function
- **DB/API**: UniProt; InterPro; ClinVar; VEP
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `PM1`
- **MET時 direction**: `supports`
- **strength**: moderate or adjusted strength
- **`evidenceOutcome`**: PM1 if default moderate; otherwise PM1_<strength>; PM1_not_met
- **Evidence Item**: Generic StudyResult for domain/hotspot annotation and local variant spectrum; supporting Statements/Documents.
- **Curator確認点**: Whether the region is truly critical/hotspot and sufficiently depleted of benign variation.

### PM2

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Absent from controls, or at extremely low frequency for a recessive disorder, in population databases.
- **ClinGen General**: ClinGen recommends using PM2 at Supporting strength in general and emphasizes reliable population data/rarity rather than naïve absence.
- **具体的処理**: Use reliable ancestry-aware population frequency/FAF. If absent or sufficiently rare under the configured general implementation rule, emit PM2_supporting. If coverage/AN/region quality is insufficient, mark internal status NOT_EVALUATED rather than not-met.
- **Implementation default / 注意**: No universal disease-agnostic nonzero AF cutoff is asserted. If adopting fastVEP defaults (e.g., AD/unknown ≤4e-5; AR ≤7e-5), label them IMPLEMENTATION thresholds.
- **必要入力**: AC; AN; AF; ancestry-specific AF; FAF; coverage/quality; inheritance context
- **DB/API**: gnomAD; NCBI ALFA
- **自動化**: Automatic with quality gates
- **VA-Spec `methodType`**: `PM2`
- **MET時 direction**: `supports`
- **strength**: supporting
- **`evidenceOutcome`**: PM2_supporting (because strength is adjusted from ACMG default moderate); PM2_not_met
- **Evidence Item**: CohortAlleleFrequencyStudyResult; use ancillaryResults/qualityMeasures for FAF and coverage/quality metadata as needed.
- **Curator確認点**: Rarity threshold policy and reliability of population observation.

### PM3

- **方向**: Pathogenic
- **ACMG/AMP 2015**: For recessive disorders, detected in trans with a pathogenic variant.
- **ClinGen General**: ClinGen PM3 point framework weights phase confirmation, classification of the other allele, homozygous observations and independent probands; total points map to strength.
- **具体的処理**: For each proband: determine phase and classification of second allele → assign points → cap/reconcile repeated evidence as required → sum independent probands → map to Supporting/Moderate/Strong/Very Strong.
- **Implementation default / 注意**: A practical mapping often uses 0.5/1/2/4 points for Supporting/Moderate/Strong/Very Strong; preserve proband-level points.
- **必要入力**: proband genotype; second variant; second-variant classification; phase; parental genotypes; phenotype; homozygous status
- **DB/API**: Internal family data; ClinVar; ClinGen Evidence Repository
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `PM3`
- **MET時 direction**: `supports`
- **strength**: supporting / moderate / strong / very strong
- **`evidenceOutcome`**: PM3 if default moderate; otherwise PM3_<strength>; PM3_not_met
- **Evidence Item**: Generic StudyResult for phase/genotypes plus Statement for classification of the second allele.
- **Curator確認点**: Second-allele classification, phase certainty and duplicate probands.

### PM4

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Protein length changes due to in-frame deletions/insertions in a non-repeat region or stop-loss variants.
- **ClinGen General**: Avoid double counting with PVS1; interpret protein-region relevance and repetitive sequence context.
- **具体的処理**: Consequence is in-frame insertion/deletion or stop-loss + changes protein length + not solely in non-functional repetitive region + same molecular consequence is not already counted under PVS1 → PM4.
- **Implementation default / 注意**: If repeat/domain annotation is unavailable, do not auto-fire.
- **必要入力**: consequence; amino-acid length change; repeat annotation; functional region/domain
- **DB/API**: VEP; UniProt; InterPro; RepeatMasker
- **自動化**: Automatic / semi-automatic
- **VA-Spec `methodType`**: `PM4`
- **MET時 direction**: `supports`
- **strength**: moderate or adjusted strength
- **`evidenceOutcome`**: PM4 if default moderate; otherwise PM4_<strength>; PM4_not_met
- **Evidence Item**: Generic StudyResult for consequence, length change, repeat/domain annotation.
- **Curator確認点**: Functional relevance of affected region and overlap with PVS1.

### PM5

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Novel missense change at an amino-acid residue where a different missense change is established as pathogenic.
- **ClinGen General**: Comparator pathogenicity and mechanism must be sufficiently established; same amino-acid substitution belongs under PS1, not PM5.
- **具体的処理**: Evaluated variant is missense + same residue has a different Pathogenic missense + comparator has adequate evidence + not same AA substitution → PM5.
- **Implementation default / 注意**: Prefer comparator with reviewed primary evidence; exclude comparator whose pathogenicity is mainly splice-mediated if not mechanistically comparable.
- **必要入力**: residue; evaluated AA change; comparator AA change; comparator classification/evidence
- **DB/API**: ClinVar; ClinGen Evidence Repository; PubMed
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `PM5`
- **MET時 direction**: `supports`
- **strength**: moderate or adjusted strength
- **`evidenceOutcome`**: PM5 if default moderate; otherwise PM5_<strength>; PM5_not_met
- **Evidence Item**: Statement describing comparator pathogenic variant plus generic StudyResult for residue equivalence.
- **Curator確認点**: Quality and mechanistic relevance of comparator evidence.

### PM6

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Assumed de novo, but without confirmation of maternity and paternity.
- **ClinGen General**: Evaluated jointly with PS2 using the ClinGen de novo point framework; parental relationship status and phenotype specificity control points/strength.
- **具体的処理**: Proband carries ALT + available parent(s) do not + parental relationship not fully confirmed → PM6 candidate → phenotype/relationship points → aggregate independent probands → strength.
- **Implementation default / 注意**: Genotype QC thresholds are implementation policy and must be provenance-tracked.
- **必要入力**: proband/parent genotypes; relationship confirmation status; genotype QC; phenotype; independent probands
- **DB/API**: Internal family data; HPO; MONDO
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `PM6`
- **MET時 direction**: `supports`
- **strength**: supporting / moderate / strong / very strong when justified by point framework
- **`evidenceOutcome`**: PM6 if default moderate; otherwise PM6_<strength>; PM6_not_met
- **Evidence Item**: Generic StudyResult for trio/family observations and phenotype evidence.
- **Curator確認点**: Whether the observation qualifies as confirmed vs assumed de novo.

### PP1

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Co-segregation with disease in multiple affected family members in a gene definitively known to cause the disease.
- **ClinGen General**: ClinGen PP1/BS4/PP4 guidance recommends quantitative/structured segregation assessment and integration with phenotype/locus specificity; strength may increase.
- **具体的処理**: Build informative pedigree observations → account for inheritance, penetrance, phenocopy and locus heterogeneity → calculate segregation LR/points or equivalent structured measure → map to evidence strength.
- **Implementation default / 注意**: Do not use simple affected-carrier count alone when pedigree structure materially changes evidence.
- **必要入力**: pedigree; genotypes; affected status; age; penetrance assumptions; inheritance; phenotype/locus specificity
- **DB/API**: Internal pedigree data; HPO; MONDO
- **自動化**: Semi-automatic / curator
- **VA-Spec `methodType`**: `PP1`
- **MET時 direction**: `supports`
- **strength**: supporting / moderate / strong as justified
- **`evidenceOutcome`**: PP1 if default supporting; otherwise PP1_<strength>; PP1_not_met
- **Evidence Item**: Generic StudyResult for segregation observations; Statement/StudyResult for phenotype specificity.
- **Curator確認点**: Informative meioses, penetrance, phenocopy, locus heterogeneity.

### PP2

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Missense variant in a gene with low rate of benign missense variation and where missense variants are a common disease mechanism.
- **ClinGen General**: Gene-level constraint alone is insufficient; disease mechanism and pathogenic variant spectrum must support missense causality.
- **具体的処理**: Variant is missense + established gene-disease mechanism includes pathogenic missense + benign missense variation is low relative to expectation → PP2. If disease mechanism is unknown, do not auto-fire.
- **Implementation default / 注意**: Missense constraint Z-score can be a heuristic, but must not replace gene-disease mechanism evidence.
- **必要入力**: variant type; gene-disease mechanism; pathogenic/benign missense spectrum; gene constraint
- **DB/API**: ClinGen Gene-Disease Validity; gnomAD constraint; ClinVar
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `PP2`
- **MET時 direction**: `supports`
- **strength**: supporting
- **`evidenceOutcome`**: PP2; PP2_not_met
- **Evidence Item**: Statement for disease mechanism plus generic StudyResult for gene-level missense constraint/variant spectrum.
- **Curator確認点**: Whether missense is genuinely a common mechanism for the evaluated disease.

### PP3

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Multiple lines of computational evidence support a deleterious effect on the gene/gene product.
- **ClinGen General**: Use calibrated predictors rather than simple majority voting; ClinGen provides calibrated evidence strengths for missense prediction and separate splicing guidance.
- **具体的処理**: Select predictor appropriate to variant mechanism → compare score with a documented calibrated interval → emit PP3 at the calibrated strength. Keep protein-impact and splice-impact evidence distinct and reconcile double counting.
- **Implementation default / 注意**: fastVEP uses REVEL ≥0.644 Supporting, ≥0.773 Moderate, ≥0.932 Strong; SpliceAI ≥0.2 as Supporting. Treat these as adopted implementation thresholds with source/version.
- **必要入力**: predictor name/version; raw score; variant type; calibration source/threshold; splice context if applicable
- **DB/API**: REVEL; SpliceAI; VEP; other calibrated predictor source
- **自動化**: Automatic
- **VA-Spec `methodType`**: `PP3`
- **MET時 direction**: `supports`
- **strength**: supporting / moderate / strong
- **`evidenceOutcome`**: PP3 if supporting; PP3_moderate / PP3_strong when adjusted; PP3_not_met
- **Evidence Item**: Generic StudyResult for computational predictor score(s); no dedicated computational StudyResult profile in VA-Spec v1.0.
- **Curator確認点**: Calibration source/version and double-counting with splice/functional evidence.

### PP4

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Patient phenotype/family history is highly specific for a disease with a single genetic etiology.
- **ClinGen General**: ClinGen guidance explicitly links PP4 phenotype specificity with locus heterogeneity and segregation evidence (PP1/BS4); disease-specificity must be demonstrated.
- **具体的処理**: Normalize phenotype → compare with candidate disease/gene phenotype → assess specificity and genetic heterogeneity → combine with segregation context as required → assign evidence strength or route to curator.
- **Implementation default / 注意**: Do not impose a universal HPO similarity cutoff as a ClinGen rule.
- **必要入力**: HPO terms; diagnosis; age/sex; labs; family history; candidate disease/gene; locus heterogeneity
- **DB/API**: HPO; MONDO; ClinGen Gene-Disease Validity
- **自動化**: Semi-automatic / curator
- **VA-Spec `methodType`**: `PP4`
- **MET時 direction**: `supports`
- **strength**: supporting or adjusted strength when justified
- **`evidenceOutcome`**: PP4 if default supporting; otherwise PP4_<strength>; PP4_not_met
- **Evidence Item**: Statement or generic StudyResult capturing phenotype-to-disease/gene match and specificity assessment.
- **Curator確認点**: Phenotype specificity and alternative genetic etiologies.

### PP5

- **方向**: Pathogenic
- **ACMG/AMP 2015**: Reputable source recently reports variant as pathogenic, but evidence is not available to the laboratory.
- **ClinGen General**: ClinGen recommends that PP5 not be used; retrieve and evaluate primary evidence instead.
- **具体的処理**: Never produce MET under the General ClinGen policy. Use external classifications only to discover primary evidence.
- **Implementation default / 注意**: Internal workflow status = NOT_APPLICABLE/DEPRECATED.
- **必要入力**: None for scoring; external classifications may be used for evidence discovery
- **DB/API**: ClinVar for discovery only
- **自動化**: Disabled
- **VA-Spec `methodType`**: `PP5`
- **MET時 direction**: `supports`
- **strength**: supporting (legacy only)
- **`evidenceOutcome`**: Do not emit a MET EvidenceLine. If the implementation explicitly evaluates it as an ACMG criterion and records not-met, PP5_not_met is schema-compatible; deprecation itself should be kept in Extension/local workflow metadata.
- **Evidence Item**: Normally none.
- **Curator確認点**: N/A; deprecated.

### BA1

- **方向**: Benign
- **ACMG/AMP 2015**: Allele frequency is >5% in population databases.
- **ClinGen General**: ClinGen BA1 recommendations emphasize ancestry/population context, adequate allele number, data quality, and an exception list.
- **具体的処理**: Use reliable population AF/FAF → require adequate observed alleles → if frequency >0.05 and not a BA1 exception → BA1. If frequency quality is unreliable, NOT_EVALUATED internally.
- **Implementation default / 注意**: A minimum observed allele count/AN gate may be adopted; document it. The user source proposes at least 2,000 observed alleles based on ClinGen guidance.
- **必要入力**: population; AC; AN; AF/FAF; homozygotes; quality; BA1 exception status
- **DB/API**: gnomAD; ClinGen BA1 exception information
- **自動化**: Automatic with quality gates
- **VA-Spec `methodType`**: `BA1`
- **MET時 direction**: `disputes`
- **strength**: stand alone
- **`evidenceOutcome`**: BA1; BA1_not_met
- **Evidence Item**: CohortAlleleFrequencyStudyResult with sourceDataSet/cohort and relevant quality metadata.
- **Curator確認点**: Exception-list status and unreliable frequency regions.

### BS1

- **方向**: Benign
- **ACMG/AMP 2015**: Allele frequency is greater than expected for the disorder.
- **ClinGen General**: Requires comparison with the maximum credible allele frequency expected for the disease; depends on prevalence, inheritance, penetrance and allelic/genetic contribution.
- **具体的処理**: Compute/obtain disease-specific maximum credible AF → compare ancestry-aware observed AF/FAF → if observed exceeds threshold, BS1. Without disease-level assumptions/threshold, do not claim a general automatic BS1.
- **Implementation default / 注意**: A fixed fallback such as AF>1% is an implementation heuristic only, not a universal ClinGen rule.
- **必要入力**: AF/FAF; prevalence; inheritance; penetrance; genetic contribution; allelic contribution; threshold assumptions
- **DB/API**: gnomAD; epidemiologic/curated disease data
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `BS1`
- **MET時 direction**: `disputes`
- **strength**: strong or adjusted
- **`evidenceOutcome`**: BS1 if default strong; otherwise BS1_<strength>; BS1_not_met
- **Evidence Item**: CohortAlleleFrequencyStudyResult plus Statement/StudyResult carrying disease-frequency assumptions.
- **Curator確認点**: Maximum credible AF assumptions.

### BS2

- **方向**: Benign
- **ACMG/AMP 2015**: Observed in a healthy adult individual for a recessive, dominant, or X-linked disorder where full penetrance is expected at an early age.
- **ClinGen General**: Interpretation depends on expected penetrance, age of onset, genotype/zygosity and reliable phenotype ascertainment.
- **具体的処理**: Confirm disease expected to be highly penetrant by observed age → identify genotype incompatible with disease in well-phenotyped healthy individual(s) → assess zygosity/inheritance → BS2. Population counts alone require caution.
- **Implementation default / 注意**: fastVEP-style AC or homozygote-count heuristics are implementation defaults, not universal ClinGen thresholds.
- **必要入力**: genotype/zygosity; age; phenotype; penetrance; expected onset; inheritance; healthy-individual count
- **DB/API**: Well-phenotyped cohorts; PubMed; population datasets as supportive context
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `BS2`
- **MET時 direction**: `disputes`
- **strength**: strong or adjusted
- **`evidenceOutcome`**: BS2 if default strong; otherwise BS2_<strength>; BS2_not_met
- **Evidence Item**: Generic StudyResult for healthy-individual genotype/phenotype observations; population StudyResult may be supporting context.
- **Curator確認点**: Penetrance, age of onset and quality of phenotype ascertainment.

### BS3

- **方向**: Benign
- **ACMG/AMP 2015**: Well-established functional studies show no damaging effect on protein function or splicing.
- **ClinGen General**: Same ClinGen functional-evidence framework as PS3; assay validity and controls determine whether/at what strength benign evidence is appropriate.
- **具体的処理**: Require curated disease-relevant validated assay → normal/non-damaging result with appropriate controls/validation → assign strength from assay validity.
- **Implementation default / 注意**: Keep raw assay result separate from ACMG interpretation and avoid double counting the same functional evidence.
- **必要入力**: assay; system; controls; validation variants; result; statistics; replication; PMID/MaveDB ID
- **DB/API**: MaveDB; PubMed; ClinGen Evidence Repository
- **自動化**: Curated / semi-automatic
- **VA-Spec `methodType`**: `BS3`
- **MET時 direction**: `disputes`
- **strength**: supporting / moderate / strong / very strong when justified
- **`evidenceOutcome`**: BS3 if default strong; otherwise BS3_<strength>; BS3_not_met
- **Evidence Item**: ExperimentalVariantFunctionalImpactStudyResult and/or functional-impact Statement.
- **Curator確認点**: Assay validity and disease-mechanism relevance.

### BS4

- **方向**: Benign
- **ACMG/AMP 2015**: Lack of segregation in affected members of a family.
- **ClinGen General**: ClinGen PP1/BS4/PP4 guidance requires structured segregation interpretation and consideration of penetrance, phenocopy and locus heterogeneity.
- **具体的処理**: Identify genotype-phenotype inconsistency → test whether reduced penetrance, phenocopy, age-dependent onset, heterogeneity or genotype error can explain it → if not, quantify non-segregation evidence → BS4.
- **Implementation default / 注意**: Do not auto-fire from a single discordant relative without contextual review.
- **必要入力**: pedigree; genotypes; phenotypes; ages; penetrance; inheritance; locus heterogeneity
- **DB/API**: Internal pedigree data; HPO/MONDO as context
- **自動化**: Semi-automatic / curator
- **VA-Spec `methodType`**: `BS4`
- **MET時 direction**: `disputes`
- **strength**: strong or adjusted
- **`evidenceOutcome`**: BS4 if default strong; otherwise BS4_<strength>; BS4_not_met
- **Evidence Item**: Generic StudyResult for segregation/non-segregation observations.
- **Curator確認点**: Phenocopy, reduced penetrance, age-dependent onset, locus heterogeneity.

### BP1

- **方向**: Benign
- **ACMG/AMP 2015**: Missense variant in a gene for which primarily truncating variants are known to cause disease.
- **ClinGen General**: Disease mechanism and known pathogenic variant spectrum are required; gene constraint alone should not define the rule.
- **具体的処理**: Variant is missense + established disease mechanism is predominantly LoF/truncating + pathogenic missense mechanism is not established → BP1.
- **Implementation default / 注意**: Counts of known pathogenic missense variants may be used as a heuristic but must be provenance-tracked.
- **必要入力**: variant type; disease mechanism; known pathogenic variant spectrum
- **DB/API**: ClinGen Gene-Disease Validity; ClinVar
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `BP1`
- **MET時 direction**: `disputes`
- **strength**: supporting
- **`evidenceOutcome`**: BP1; BP1_not_met
- **Evidence Item**: Statement for disease mechanism plus generic StudyResult for pathogenic variant spectrum.
- **Curator確認点**: Whether missense disease mechanism is truly absent/rare.

### BP2

- **方向**: Benign
- **ACMG/AMP 2015**: Observed in trans with a pathogenic variant for a fully penetrant dominant disorder, or observed in cis with a pathogenic variant in any inheritance pattern.
- **ClinGen General**: Interpret phase in the context of inheritance and penetrance; require reliable phasing and comparator classification.
- **具体的処理**: Obtain two variants + phase + inheritance → if configuration satisfies BP2 logic and second variant is sufficiently established pathogenic → BP2. Unknown phase → do not auto-fire.
- **Implementation default / 注意**: Phased family/read-based evidence should be preferred to inferred phase.
- **必要入力**: two variants; phase; inheritance; penetrance; second-variant classification/evidence
- **DB/API**: Internal family data; ClinVar; ClinGen Evidence Repository
- **自動化**: Semi-automatic
- **VA-Spec `methodType`**: `BP2`
- **MET時 direction**: `disputes`
- **strength**: supporting
- **`evidenceOutcome`**: BP2; BP2_not_met
- **Evidence Item**: Generic StudyResult for phase/genotypes plus Statement for second-variant pathogenicity.
- **Curator確認点**: Phase certainty and inheritance/penetrance assumptions.

### BP3

- **方向**: Benign
- **ACMG/AMP 2015**: In-frame deletion/insertion in a repetitive region without a known function.
- **ClinGen General**: Confirm that the repeat region lacks established functional importance; avoid conflict with PM4.
- **具体的処理**: In-frame insertion/deletion + overlaps repetitive region + no established functional domain/critical region → BP3. If domain/repeat interpretation is uncertain, review.
- **Implementation default / 注意**: RepeatMasker overlap alone is insufficient if the region has known function.
- **必要入力**: indel length/frame; repeat annotation; functional domain/region annotation
- **DB/API**: VEP; RepeatMasker; UniProt; InterPro
- **自動化**: Automatic / semi-automatic
- **VA-Spec `methodType`**: `BP3`
- **MET時 direction**: `disputes`
- **strength**: supporting
- **`evidenceOutcome`**: BP3; BP3_not_met
- **Evidence Item**: Generic StudyResult for in-frame consequence, repeat and domain annotation.
- **Curator確認点**: Functional relevance of repeat region.

### BP4

- **方向**: Benign
- **ACMG/AMP 2015**: Multiple lines of computational evidence suggest no impact on gene/gene product.
- **ClinGen General**: Use calibrated computational predictors instead of majority vote; benign calibrated intervals can support adjusted strengths.
- **具体的処理**: Select calibrated predictor → compare score to benign calibration interval → assign BP4 strength → reconcile with splice/functional evidence and avoid applying benign prediction where incompatible with high-confidence null/splice evidence.
- **Implementation default / 注意**: fastVEP uses REVEL ≤0.290 Supporting, ≤0.183 Moderate, ≤0.016 Strong, ≤0.003 Very Strong; SpliceAI ≤0.1 Supporting. Record predictor/version and calibration source.
- **必要入力**: predictor/version; score; variant type; calibration thresholds; splice context
- **DB/API**: REVEL; SpliceAI; VEP
- **自動化**: Automatic
- **VA-Spec `methodType`**: `BP4`
- **MET時 direction**: `disputes`
- **strength**: supporting / moderate / strong / very strong
- **`evidenceOutcome`**: BP4 if supporting; BP4_moderate / BP4_strong / BP4_very_strong when adjusted; BP4_not_met
- **Evidence Item**: Generic StudyResult for computational score(s); no dedicated computational StudyResult profile in VA-Spec v1.0.
- **Curator確認点**: Calibration source/version and conflicts/double counting.

### BP5

- **方向**: Benign
- **ACMG/AMP 2015**: Variant found in a case with an alternate molecular basis for disease.
- **ClinGen General**: Alternative diagnosis must adequately explain the phenotype; dual diagnosis must be considered.
- **具体的処理**: Identify alternate causal variant/diagnosis → assess whether it fully or near-fully explains phenotype → if evaluated variant lacks independent phenotype contribution, BP5; if dual diagnosis remains plausible, MANUAL_REVIEW internally.
- **Implementation default / 注意**: Do not auto-fire solely because another pathogenic variant exists.
- **必要入力**: patient phenotype; alternate variant/gene; alternate classification; phenotype match; dual-diagnosis assessment
- **DB/API**: Internal patient data; ClinVar; HPO; MONDO
- **自動化**: Manual / semi-automatic
- **VA-Spec `methodType`**: `BP5`
- **MET時 direction**: `disputes`
- **strength**: supporting
- **`evidenceOutcome`**: BP5; BP5_not_met
- **Evidence Item**: Statement describing alternative molecular diagnosis plus phenotype-match evidence.
- **Curator確認点**: Completeness of alternative explanation and dual diagnosis.

### BP6

- **方向**: Benign
- **ACMG/AMP 2015**: Reputable source recently reports variant as benign, but evidence is not available to the laboratory.
- **ClinGen General**: ClinGen recommends that BP6 not be used; retrieve and evaluate primary evidence instead.
- **具体的処理**: Never produce MET under the General ClinGen policy. External benign classifications are evidence-discovery aids only.
- **Implementation default / 注意**: Internal workflow status = NOT_APPLICABLE/DEPRECATED.
- **必要入力**: None for scoring; external classifications may be used for evidence discovery
- **DB/API**: ClinVar for discovery only
- **自動化**: Disabled
- **VA-Spec `methodType`**: `BP6`
- **MET時 direction**: `disputes`
- **strength**: supporting (legacy only)
- **`evidenceOutcome`**: Do not emit a MET EvidenceLine. If explicitly evaluated and recorded not-met, BP6_not_met is schema-compatible; deprecation belongs in Extension/local workflow metadata.
- **Evidence Item**: Normally none.
- **Curator確認点**: N/A; deprecated.

### BP7

- **方向**: Benign
- **ACMG/AMP 2015**: Synonymous variant with no predicted splice impact and low nucleotide conservation.
- **ClinGen General**: ClinGen splicing guidance refines location/prediction/RNA considerations and warns against double counting across PVS1/PS1/PP3/BP4/BP7.
- **具体的処理**: Synonymous/noncoding candidate appropriate for BP7 → outside splice-critical positions as defined by adopted rule → calibrated splice prediction supports no impact → conservation/RNA evidence does not suggest function → BP7. Reconcile with other splice criteria.
- **Implementation default / 注意**: A practical implementation may use SpliceAI ≤0.1 and a documented conservation threshold; any intronic distance limits are implementation policy, not universal ClinGen thresholds.
- **必要入力**: consequence; exon/intron offset; splice predictor/version/score; conservation; RNA evidence if available
- **DB/API**: VEP; SpliceAI; PhyloP/PhyloCons; PubMed/RNA evidence
- **自動化**: Automatic / semi-automatic
- **VA-Spec `methodType`**: `BP7`
- **MET時 direction**: `disputes`
- **strength**: supporting
- **`evidenceOutcome`**: BP7; BP7_not_met
- **Evidence Item**: Generic StudyResult for splice prediction/conservation; Statement or functional StudyResult for RNA evidence where available.
- **Curator確認点**: Splicing edge cases, RNA evidence, and double counting.

## 5. Evidence Item の使い分け

VA-Spec v1.0 で専用 profile があるものは優先して使う。

### Population evidence
`CohortAlleleFrequencyStudyResult`

主対象:
- BA1
- BS1
- PM2
- BS2 の population observation を補助的に表現する場合

### Functional evidence
`ExperimentalVariantFunctionalImpactStudyResult`

主対象:
- PS3
- BS3

### その他
VA-Spec v1.0 には de novo、segregation、computational prediction、phenotype specificity 等について、
ACMG criterion ごとの専用 StudyResult profile は標準では用意されていない。
このため以下を Evidence Item として使用する。

```text
StudyResult（generic）
Statement
EvidenceLine
iriReference
```

Criterion 固有の追加データは、ACMG EvidenceLine の標準項目を壊さず、
Evidence Item / ancillaryResults / qualityMeasures / Extension 等へ保持する。

## 6. 内部 Rule Engine と VA-Spec Export の境界

```text
Raw / Curated Evidence
  ↓
Observation / Evidence Store
  ↓
ACMG + ClinGen General Rule Engine
  ↓
MET / NOT_MET / NOT_EVALUATED /
NOT_APPLICABLE / MANUAL_REVIEW
  ↓
VA-Spec Export Mapper
  ↓
EvidenceLine[]
  ↓
ACMG Combination Rule
  ↓
Variant Pathogenicity Statement
(P / LP / VUS / LB / B)
```

## 7. 最終 Statement 例

```json
{
  "id": "ex:statement001",
  "type": "Statement",
  "proposition": {
    "id": "ex:proposition001"
  },
  "classification": {
    "primaryCoding": {
      "code": "likely pathogenic",
      "system": "ACMG Guidelines, 2015"
    }
  },
  "specifiedBy": {
    "type": "Method",
    "name": "ACMG/AMP 2015 with ClinGen General Guidance"
  },
  "hasEvidenceLines": [
    "ex:PVS1_001",
    "ex:PM2_001",
    "ex:PP3_001"
  ]
}
```

最終分類:
- `pathogenic`
- `likely pathogenic`
- `uncertain significance`
- `likely benign`
- `benign`

## 8. Provenance 最低要件

- VA-Spec version
- criterion / methodType
- evidenceOutcome
- direction / strength
- ACMG/ClinGen guidance 名・version
- implementation rule version
- variant representation / reference genome
- transcript
- gene
- disease/condition context
- inheritance context
- input DB / dataset / version
- predictor / tool / version
- retrieval date / evaluation date
- evaluator
- curator review status

`IMPLEMENTATION` の固定閾値は ClinGen General Guidance と混同しない。

## 9. 主な根拠資料

- Richards S, et al. ACMG/AMP 2015.
- ClinGen Variant Classification Guidance  
  https://www.clinicalgenome.org/tools/clingen-variant-classification-guidance/
- GA4GH VA-Spec v1.0  
  https://va-spec.ga4gh.org/en/1.0/
- ACMG Variant Pathogenicity Statement with Evidence example  
  https://va-spec.ga4gh.org/en/1.0/examples/acmg-variant-pathogenicity-statement-with-evidence.html
- fastVEP ACMG implementation  
  https://github.com/Huang-lab/fastVEP/blob/master/docs/ACMG.md

## 10. 適用上の注意

- CSpec/VCEP-specific rule は本仕様では適用しない。
- fastVEP の固定閾値は、ClinGen universal rule でないものを `IMPLEMENTATION` として区別する。
- VA-Spec v1 の ACMG EvidenceLine profile は draft maturity と記載されるため、実装時は利用する v1.x JSON Schema を固定する。
