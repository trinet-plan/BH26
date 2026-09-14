# ACMG / ClinGen Variant Interpretation Tool
## 期待される入力データ設計

## 1. 基本方針

本ツールでは、**VCFのvariant 1件を1評価単位とする**。

ただし、ACMG/AMP criterionの多くはvariant情報だけでは評価できないため、
実際の入力モデルは以下のように考える。

```text
必須:
  Variant

任意:
  Disease / Condition
  Phenotype
  Case information
  Family / Pedigree
  Genotype
  Functional evidence
```

つまり、本ツールは単なる「VCF判定ツール」ではなく、

> **1 Variantを起点にEvidenceを集約し、ACMG criterionを評価するツール**

として設計する。

---

## 2. 最小入力: Variant

最小入力はvariantそのものとする。

必要項目:

```text
CHROM
POS
REF
ALT
reference assembly
```

VCF例:

```vcf
chr13   32316461   .   C   T   .   PASS   .
```

ツール内部では、例えば以下のような正規化済みvariant表現に変換する。

```json
{
  "assembly": "GRCh38",
  "chrom": "13",
  "pos": 32316461,
  "ref": "C",
  "alt": "T"
}
```

---

## 3. Variantを起点とした処理

```text
Variant
  ↓
Normalization
  ↓
VEP
  ├ consequence
  ├ transcript
  ├ exon / intron
  ├ HGVS
  ├ codon
  └ amino-acid position
  ↓
External Evidence
  ├ gnomAD
  ├ ClinVar
  ├ ClinGen
  ├ REVEL
  ├ SpliceAI
  ├ UniProt / InterPro
  └ MaveDB / PubMed
  ↓
Evidence Store
  ↓
ACMG / ClinGen General Rule Engine
  ↓
VA-Spec EvidenceLine
```

---

## 4. Variantのみで評価可能なcriterion

Variantと外部DB・annotationから評価可能な主なcriterion:

```text
PVS1
PS1
PM1
PM2
PM4
PM5
PP2
PP3
BA1
BP1
BP3
BP4
BP7
```

ただし、PVS1やPP2などはgene-disease mechanismの情報が不足する場合、
完全自動ではなく `MANUAL_REVIEW` または `NOT_EVALUATED` とする。

### 例

```text
Variant
  ↓
VEP
  ├ consequence
  ├ transcript
  ├ HGVS
  └ protein position

  ↓
gnomAD
  ├ AF
  ├ FAF
  ├ AC
  └ AN

  ↓
ClinVar
  └ comparator variants

  ↓
REVEL / SpliceAI
  └ prediction score
```

---

## 5. Disease / Conditionを追加した場合

入力例:

```json
{
  "variant": {
    "assembly": "GRCh38",
    "chrom": "13",
    "pos": 32316461,
    "ref": "C",
    "alt": "T"
  },
  "condition": {
    "mondo_id": "MONDO:xxxxxxx"
  },
  "inheritance": "autosomal_dominant"
}
```

Disease contextがあることで、以下のようなcriterionの評価精度が上がる。

```text
PVS1
PP2
PP4
BS1
BP1
```

特にPVS1では、

```text
LoF consequence
  ↓
そのgene-disease relationshipで
LoFが疾患メカニズムとして成立しているか
```

を確認する必要がある。

そのため、

```text
stop_gained
```

だけを理由にPVS1を自動成立させない。

---

## 6. Phenotypeを追加した場合

入力例:

```json
{
  "case": {
    "patient_id": "CASE001",
    "hpo_terms": [
      "HP:0001250",
      "HP:0001263"
    ]
  }
}
```

Phenotype contextがあることで、主に以下を評価可能になる。

```text
PS2 / PM6 の phenotype consistency
PP1
PP4
BS4
BP5
```

HPOを用いてphenotypeを正規化し、
MONDO等のdisease ontologyと組み合わせて評価する。

---

## 7. Family / Trioデータを追加した場合

入力例:

```json
{
  "family": {
    "proband_genotype": "0/1",
    "father_genotype": "0/0",
    "mother_genotype": "0/0",
    "maternity_confirmed": true,
    "paternity_confirmed": true
  }
}
```

family / trio情報があることで、以下を評価可能になる。

```text
PS2
PM3
PM6
PP1
BS4
BP2
```

---

## 8. 推奨入力レベル

ツールとしては、入力を3段階に分けると扱いやすい。

| Level | 入力 | 主に評価できる内容 |
|---|---|---|
| **Level 1** | Variant | DB・annotation・prediction中心の自動判定 |
| **Level 2** | Variant + Disease | gene-disease mechanism、疾患頻度、inheritanceを含む判定 |
| **Level 3** | Variant + Disease + Case/Family | de novo、segregation、phenotype、phaseまで評価 |

### 入力モデル

```text
Required
  └ Variant

Optional
  ├ Disease / Condition
  ├ Inheritance
  ├ Phenotype
  ├ Case
  ├ Family / Pedigree
  ├ Genotype
  └ Functional evidence
```

---

## 9. VCFを入力にする場合

ユーザー操作としては、VCFアップロードを基本としてよい。

VCF:

```vcf
#CHROM POS      REF ALT
1      123456   A   G
2      345678   C   T
3      456789   G   A,C
```

内部ではALT allele単位に分解する。

```text
Variant 1: 1:123456 A>G
Variant 2: 2:345678 C>T
Variant 3: 3:456789 G>A
Variant 4: 3:456789 G>C
```

つまり、

> **1 ALT allele = 1 Variant Evaluation**

とする。

multi-allelic siteは評価前にsplit / normalizeする。

---

## 10. VCF FORMAT情報の利用

VCFにsample genotypeが含まれている場合は、
GT / GQ / DPなどをfamily evidenceとして利用できる。

例:

```vcf
FORMAT     Proband      Father       Mother
GT:GQ:DP   0/1:99:50    0/0:99:47    0/0:99:55
```

内部では、

```text
Proband ALTあり
Father ALTなし
Mother ALTなし
```

まで取得可能。

ただし、VCFだけでは、

```text
どのsampleが
proband / father / motherなのか
```

という関係は通常判別できない。

そのためfamily evaluationを行う場合は、

```text
VCF
+
PED
```

または同等のfamily metadataが必要。

---

## 11. 推奨入力セット

### 最小構成

```text
VCF
```

### Disease-aware構成

```text
VCF
+
Disease / MONDO
+
Inheritance
```

### Full case-aware構成

```text
VCF
+
PED
+
Phenotype / HPO
+
Disease / MONDO
+
Inheritance
```

---

## 12. 初期実装として推奨する範囲

最初の実装では、以下の構成が最も実装しやすい。

```text
Input
  ↓
VCF 1 variant
  ↓
Normalize
  ↓
VEP
gnomAD
ClinVar
ClinGen
REVEL
SpliceAI
  ↓
Automatic / Candidate criteria
  ↓
VA-Spec
```

### 初期入力

```json
{
  "assembly": "GRCh38",
  "chrom": "13",
  "pos": 32316461,
  "ref": "C",
  "alt": "T"
}
```

### 初期出力イメージ

```text
PVS1    MANUAL_REVIEW
PM2     MET / Supporting
PP3     MET / Moderate
BA1     NOT_MET
BP4     NOT_MET
```

併せて、criterion判定に用いたEvidenceも返す。

```text
gnomAD AF
gnomAD FAF
ClinVar classification
ClinVar review status
REVEL score
SpliceAI score
VEP consequence
Transcript
Protein position
```

---

## 13. 推奨API入力モデル

将来的には、以下のようなJSONを基本入力モデルとする。

```json
{
  "variant": {
    "assembly": "GRCh38",
    "chrom": "13",
    "pos": 32316461,
    "ref": "C",
    "alt": "T"
  },
  "condition": {
    "mondo_id": "MONDO:xxxxxxx"
  },
  "inheritance": "autosomal_dominant",
  "case": {
    "patient_id": "CASE001",
    "hpo_terms": [
      "HP:0001250",
      "HP:0001263"
    ]
  },
  "family": {
    "proband_genotype": "0/1",
    "father_genotype": "0/0",
    "mother_genotype": "0/0",
    "maternity_confirmed": true,
    "paternity_confirmed": true
  }
}
```

`variant`のみ必須とし、
その他の要素はoptionalとする。

---

## 14. ツールの定義

本ツールは、

> **VCFのvariantを起点として、外部DB・予測器・case/family情報からEvidenceを収集し、ACMG/ClinGen General Guidanceに基づいてcriterionを評価し、VA-Spec形式で出力するツール**

として定義する。

入力情報が不足するcriterionについては、
無理に `MET / NOT_MET` を決定せず、

```text
NOT_EVALUATED
NOT_APPLICABLE
MANUAL_REVIEW
```

として内部workflow statusを保持する。
