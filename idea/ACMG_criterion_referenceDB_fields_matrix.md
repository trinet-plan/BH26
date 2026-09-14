# ACMG / ClinGen General Guidance
## ACMG criterion × 参照DB × 取得フィールド マトリクス

## 1. 目的

本資料は、ACMG/AMP 2015 の各 criterion について、ClinGen General Guidance を前提として、
判定に必要な外部・内部データソースと取得フィールドを整理したものである。

本資料では、API や配布ファイル固有の物理カラム名ではなく、実装上必要となる
**論理フィールド**を定義する。

想定フロー:

```text
External / Internal Source
        ↓
Normalized Evidence Store
        ↓
Criterion Rule Engine
        ↓
VA-Spec EvidenceLine
```

---

## 2. ACMG criterion × 参照DB × 取得フィールド

| Criterion | 参照DB / Source | 取得フィールド | 主な用途 | 取得形態 |
|---|---|---|---|---|
| **PVS1** | VEP | consequence, gene, transcript, exon, intron, CDS position, protein position, HGVS.c, HGVS.p | LoF variant判定、位置判定 | Direct |
| PVS1 | MANE / RefSeq | transcript ID, MANE status, exon structure, CDS start/end | biologically relevant transcript選択、NMD判定 | Direct |
| PVS1 | ClinGen Gene-Disease Validity | gene, disease, validity classification, inheritance / disease mechanism context | LoFが疾患機序として妥当か | Direct + Curated |
| PVS1 | UniProt | protein length, functional region, domain, active/binding site | NMD escape時の残存protein評価 | Direct |
| PVS1 | InterPro | domain ID, domain coordinates, domain description | 欠失領域の重要性評価 | Direct |
| **PS1** | VEP | codon, amino-acid position, reference AA, alternate AA, HGVS.p | 同一AA置換の比較 | Direct |
| PS1 | ClinVar | variant ID, HGVS, clinical significance, review status, condition, assertion source | comparator pathogenic variant探索 | Direct |
| PS1 | ClinGen Evidence Repository | variant, classification, criterion evidence, expert-panel assertion | comparatorの根拠確認 | Direct |
| PS1 | PubMed | PMID, study metadata, functional/splicing evidence | comparator primary evidence | Curated |
| **PS2** | Internal family data | proband genotype, mother genotype, father genotype, DP/GQ等QC, pedigree | de novo判定 | Internal |
| PS2 | Internal family data | maternity confirmed, paternity confirmed, family history | PS2 / PM6区別 | Internal |
| PS2 | HPO | patient HPO terms | phenotype consistency判定 | Internal + Ontology |
| PS2 | MONDO | disease ID, disease term | disease正規化 | Direct |
| PS2 | ClinGen Gene-Disease Validity | gene, disease, validity | phenotype/gene-disease妥当性 | Direct |
| **PS3** | MaveDB | variant, assay ID, assay type, score, score direction, experimental system, controls, replication/metadata | functional damaging evidence | Direct |
| PS3 | PubMed | PMID, assay method, system, readout, controls, statistics, result | assay validity評価 | Curated |
| PS3 | ClinGen Evidence Repository | functional evidence, criterion assignment, strength | expert-curated PS3 evidence | Direct |
| **PS4** | Publication / cohort DB | case AC, case AN/N, control AC, control AN/N | case-control enrichment | Direct / Curated |
| PS4 | Publication / cohort DB | OR/RR, CI, p-value | association strength評価 | Direct / Derived |
| PS4 | Internal case DB | independent proband count, phenotype, variant | rare disease case accumulation | Internal |
| PS4 | PubMed | study design, ancestry, ascertainment, cohort definition | study validity確認 | Curated |
| **PM1** | UniProt | domain/region name, start/end, functional annotation | critical domain判定 | Direct |
| PM1 | InterPro | domain ID, start/end, domain function | domain overlap判定 | Direct |
| PM1 | ClinVar | pathogenic variants in region, benign variants in region, AA positions | hotspot / variant clustering | Direct + Derived |
| PM1 | VEP | protein position, consequence | regionとのintersection | Direct |
| **PM2** | gnomAD | AC, AN, AF, ancestry-specific AC/AN/AF, FAF, homozygote count | rarity / absence判定 | Direct |
| PM2 | gnomAD | filter / quality, available allele number / coverage context | absenceが評価可能か | Direct |
| PM2 | NCBI ALFA | population, AC, AN, AF | population frequency補助 | Direct |
| **PM3** | Internal family data | proband genotype, second variant, parental genotypes, phase, homozygous status | trans / homozygous判定 | Internal |
| PM3 | ClinVar | second variant clinical significance, review status, disease | second allele classification | Direct |
| PM3 | ClinGen Evidence Repository | second variant classification, supporting evidence | pathogenic性確認 | Direct |
| **PM4** | VEP | consequence, in-frame insertion/deletion, stop-loss, AA length change | protein length変化判定 | Direct |
| PM4 | RepeatMasker | repeat type, repeat coordinates, overlap | repeat領域除外 | Direct |
| PM4 | UniProt | functional region/domain coordinates | affected region評価 | Direct |
| PM4 | InterPro | domain coordinates/function | affected region評価 | Direct |
| **PM5** | VEP | amino-acid position, ref AA, alt AA | same residue判定 | Direct |
| PM5 | ClinVar | variants at same residue, clinical significance, review status, condition | different pathogenic missense探索 | Direct |
| PM5 | ClinGen Evidence Repository | comparator classification/evidence | comparator evidence確認 | Direct |
| PM5 | PubMed | comparator primary evidence | pathogenic性確認 | Curated |
| **PM6** | Internal family data | proband genotype, parental genotypes, parent availability | assumed de novo判定 | Internal |
| PM6 | Internal family data | maternity/paternity confirmation status, genotype QC | PS2との差別化 | Internal |
| PM6 | HPO / MONDO | phenotype terms, disease ID | phenotype specificity | Direct + Internal |
| **PP1** | Internal pedigree | pedigree structure, genotype, affected status, age, meioses | segregation評価 | Internal |
| PP1 | Internal / curated metadata | penetrance assumption, inheritance pattern | LR/point計算 | Curated |
| PP1 | HPO / MONDO | phenotype, disease | phenocopy/locus specificity評価 | Internal + Direct |
| **PP2** | ClinGen Gene-Disease Validity | gene, disease, gene-disease validity, disease mechanism context | missense mechanism確認 | Direct + Curated |
| PP2 | gnomAD | missense constraint / observed vs expected context | benign missense variation評価 | Direct |
| PP2 | ClinVar | pathogenic missense count/spectrum, benign missense spectrum | variant spectrum評価 | Direct + Derived |
| **PP3** | REVEL | REVEL score | missense deleterious prediction | Direct |
| PP3 | SpliceAI | splice score(s), predicted splice effect | splicing impact prediction | Direct |
| PP3 | VEP | consequence, variant type, transcript | predictor適用対象判定 | Direct |
| PP3 | Rule configuration | predictor version, calibration source, PP3 thresholds | strength決定 | Configuration |
| **PP4** | Internal patient data | HPO terms, age, sex, labs, family history | phenotype定義 | Internal |
| PP4 | HPO | HPO ID, label, hierarchy | phenotype正規化 | Direct |
| PP4 | MONDO | disease ID, disease hierarchy | disease正規化 | Direct |
| PP4 | ClinGen Gene-Disease Validity | gene, disease, validity | candidate gene-disease関係 | Direct |
| **PP5** | ClinVar | external classification, source | primary evidence探索のみ | Direct |
| PP5 | — | scoring fieldなし | **ClinGen Generalでは使用しない** | Disabled |
| **BA1** | gnomAD | AC, AN, AF, ancestry-specific AF, FAF, homozygote count | AF > 5%判定 | Direct |
| BA1 | gnomAD | quality/filter, allele number | frequency reliability確認 | Direct |
| BA1 | ClinGen BA1 exception source | variant, exception status, rationale | BA1 exception除外 | Curated |
| **BS1** | gnomAD | AC, AN, AF, ancestry-specific AF, FAF | observed AF | Direct |
| BS1 | disease evidence source | prevalence | maximum credible AF計算 | Curated |
| BS1 | disease evidence source | penetrance | maximum credible AF計算 | Curated |
| BS1 | disease evidence source | inheritance | maximum credible AF計算 | Curated |
| BS1 | disease evidence source | genetic contribution, allelic contribution | maximum credible AF計算 | Curated |
| **BS2** | Well-phenotyped cohort | genotype, zygosity, age, healthy status | healthy carrier observation | Direct / Internal |
| BS2 | Population DB | homozygote count, allele count | population-level supporting context | Direct |
| BS2 | PubMed / disease source | onset age, penetrance, inheritance | healthy observationの解釈 | Curated |
| **BS3** | MaveDB | assay ID, functional score, normal/damaging direction, experimental system, controls | non-damaging functional evidence | Direct |
| BS3 | PubMed | assay method, controls, result, statistics, replication | assay validity | Curated |
| BS3 | ClinGen Evidence Repository | BS3 evidence, strength | expert-curated evidence | Direct |
| **BS4** | Internal pedigree | genotype, phenotype, age, pedigree | non-segregation | Internal |
| BS4 | Internal / curated | penetrance, inheritance, phenocopy, locus heterogeneity | discordance解釈 | Curated |
| BS4 | HPO / MONDO | phenotype/disease | phenotype consistency | Direct |
| **BP1** | ClinGen Gene-Disease Validity | gene, disease, disease mechanism context | truncating mechanism確認 | Direct + Curated |
| BP1 | ClinVar | pathogenic LoF spectrum, pathogenic missense spectrum | known pathogenic variant spectrum | Direct + Derived |
| **BP2** | Internal family data | variant1, variant2, phase, parental genotype | cis/trans判定 | Internal |
| BP2 | ClinVar | second variant classification, review status, disease | comparator pathogenic性 | Direct |
| BP2 | ClinGen Evidence Repository | second variant evidence | comparator quality | Direct |
| BP2 | Disease/inheritance metadata | inheritance, penetrance | BP2 logic適用 | Curated |
| **BP3** | VEP | consequence, frame status, insertion/deletion length | in-frame indel判定 | Direct |
| BP3 | RepeatMasker | repeat coordinates/type | repeat overlap | Direct |
| BP3 | UniProt | functional region/domain | repeat regionに機能があるか | Direct |
| BP3 | InterPro | domain coordinates/function | functional region確認 | Direct |
| **BP4** | REVEL | REVEL score | benign missense prediction | Direct |
| BP4 | SpliceAI | splice score(s) | benign splice prediction | Direct |
| BP4 | VEP | consequence / variant type | predictor選択 | Direct |
| BP4 | Rule configuration | predictor version, calibration source, BP4 thresholds | evidence strength | Configuration |
| **BP5** | Internal patient data | phenotype/HPO | phenotype definition | Internal |
| BP5 | Internal / ClinVar | alternative variant, gene, classification | alternative molecular diagnosis | Internal + Direct |
| BP5 | HPO / MONDO | alternate disease phenotype | phenotype match | Direct |
| BP5 | Curator assessment | complete explanation / partial explanation / dual diagnosis plausible | BP5判定 | Curated |
| **BP6** | ClinVar | external benign classification | primary evidence探索のみ | Direct |
| BP6 | — | scoring fieldなし | **ClinGen Generalでは使用しない** | Disabled |
| **BP7** | VEP | consequence, exon/intron offset, transcript | synonymous/noncoding候補判定 | Direct |
| BP7 | SpliceAI | splice score(s) | no predicted splice impact | Direct |
| BP7 | PhyloP / PhyloCons | conservation score | nucleotide conservation | Direct |
| BP7 | PubMed / RNA evidence | RNA assay result, splicing observation | predictionとの矛盾確認 | Curated |

---

## 3. 参照リソース別の共通取得フィールド

| Source | 共通取得フィールド群 | 利用Criteria |
|---|---|---|
| **gnomAD** | AC, AN, AF, ancestry AF, FAF, homozygote, quality, constraint | PM2, PP2, BA1, BS1, BS2 |
| **ClinVar** | variant ID, HGVS, classification, review status, condition, evidence/source | PS1, PM1, PM3, PM5, PP2, PP5, BP1, BP2, BP5, BP6 |
| **ClinGen Evidence Repository** | variant assertion, criterion evidence, strength, expert classification | PS1, PS3, PM3, PM5, BS3, BP2 |
| **ClinGen Gene-Disease Validity** | gene, disease, validity, disease-mechanism context | PVS1, PS2, PP2, PP4, BP1 |
| **VEP** | consequence, transcript, exon/intron, HGVS, codon, AA position | PVS1, PS1, PM1, PM4, PM5, PP3, BP3, BP4, BP7 |
| **MANE / RefSeq** | transcript, transcript status, exon/CDS structure | PVS1中心 |
| **UniProt / InterPro** | domain/region, coordinates, protein function | PVS1, PM1, PM4, BP3 |
| **RepeatMasker** | repeat type, coordinates | PM4, BP3 |
| **REVEL** | score | PP3, BP4 |
| **SpliceAI** | splice score / predicted consequence | PS1補助, PP3, BP4, BP7、PVS1補助 |
| **MaveDB** | assay, functional score, experimental metadata | PS3, BS3 |
| **HPO / MONDO** | phenotype ID, disease ID、階層 | PS2, PM6, PP1, PP4, BS4, BP5 |
| **PubMed** | PMID + primary evidence | PS1, PS3, PS4, PM5, BS2, BS3, BP7等 |
| **Internal case/family DB** | genotype, pedigree, phenotype, phase, QC | PS2, PS4, PM3, PM6, PP1, BS4, BP2, BP5 |

---

## 4. Evidence Store 推奨データモデル

### 4.1 Population evidence

```text
variant
source
source_version
population
AC
AN
AF
FAF
homozygote_count
quality_status
reference_assembly
retrieved_at
```

主な利用先:

- PM2
- BA1
- BS1
- BS2
- PP2（constraintを別途保持する場合）

---

### 4.2 Computational evidence

```text
variant
transcript
predictor
predictor_version
score
calibration_source
calibration_version
threshold
interpretation
reference_assembly
retrieved_at
```

主な利用先:

- PP3
- BP4
- BP7
- PVS1 / PS1 の splice 補助判定

---

### 4.3 Functional evidence

```text
variant
source
source_version
assay_id
assay_type
experimental_system
readout
result
score
score_direction
controls
validation_variants
statistics
replication
source_document
retrieved_at
```

主な利用先:

- PS3
- BS3

---

### 4.4 Variant classification / comparator evidence

```text
variant
source
source_version
classification
review_status
condition
gene
assertion_source
evidence_summary
criterion_evidence
source_document
retrieved_at
```

主な利用先:

- PS1
- PM1
- PM3
- PM5
- PP2
- BP1
- BP2
- BP5

---

### 4.5 Gene-Disease evidence

```text
gene
disease
source
source_version
validity_classification
inheritance
disease_mechanism_context
evidence_summary
retrieved_at
```

主な利用先:

- PVS1
- PS2
- PP2
- PP4
- BP1

---

### 4.6 Family / Case evidence

```text
case_id
family_id
variant
genotype
zygosity
phase
father_genotype
mother_genotype
relationship_confirmation
phenotype
affected_status
age
family_history
genotype_qc
independent_proband_flag
```

主な利用先:

- PS2
- PS4
- PM3
- PM6
- PP1
- BS4
- BP2
- BP5

---

## 5. 実装上の基本方針

DBから取得した値を直接 ACMG criterion の判定結果として保存するのではなく、
**Raw / Normalized Evidence と criterion 判定を分離する**。

```text
Raw / Curated Evidence
        ↓
Normalized Evidence Store
        ↓
ACMG + ClinGen General Rule Engine
        ↓
MET / NOT_MET / NOT_EVALUATED /
NOT_APPLICABLE / MANUAL_REVIEW
        ↓
VA-Spec Export Mapper
        ↓
Variant Pathogenicity EvidenceLine
```

この構造にすることで、例えば:

- gnomAD の1つの population record を PM2 / BA1 / BS1 で共有できる
- REVEL の1つの score を PP3 / BP4 の双方で利用できる
- ClinVar の comparator variant を PS1 / PM5 / BP2 等で共通利用できる
- Evidence source の version 更新時に criterion 判定を再実行できる

---

## 6. VA-Specとの対応

### Population evidence

VA-Spec:

```text
CohortAlleleFrequencyStudyResult
```

主対象:

- PM2
- BA1
- BS1
- BS2 の population observation

### Functional evidence

VA-Spec:

```text
ExperimentalVariantFunctionalImpactStudyResult
```

主対象:

- PS3
- BS3

### その他の evidence

VA-Spec v1.0 では、de novo、segregation、computational prediction、
phenotype specificity などについて criterion 専用 StudyResult profile は標準では用意されていないため、
以下を利用する。

```text
StudyResult
Statement
EvidenceLine
iriReference
```

必要な criterion 固有データは、Evidence Item / ancillaryResults /
qualityMeasures / Extension 等に保持する。

---

## 7. Provenanceとして最低限保持する項目

各DB・予測ツール・内部データについて、少なくとも以下を保持する。

```text
source_name
source_type
source_version
dataset_version
reference_assembly
retrieval_date
evaluation_date
tool_name
tool_version
predictor_version
transcript
gene
disease
condition
inheritance
criterion
implementation_rule_version
curator_review_status
```

固定閾値を使用する場合は、ClinGen General Guidance と区別して
**Implementation Rule** としてversion管理する。

---

## 8. 次の物理実装レイヤ

本資料は論理フィールドの定義である。

次段階では以下を定義する。

```text
Source
× 実フィールド名
× API / Download endpoint
× dataset version
× 更新頻度
× reference assembly
× local cache可否
× license / redistribution条件
```

優先して物理マッピングする候補:

1. gnomAD
2. ClinVar
3. VEP
4. ClinGen Evidence Repository
5. ClinGen Gene-Disease Validity
6. REVEL
7. SpliceAI
8. MANE / RefSeq
9. UniProt / InterPro
10. MaveDB
