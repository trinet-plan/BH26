# Ensembl蛋白注釈・予測Evidence

更新: 2026-09-15。対象はEnsembl REST/VEP release 116、GRCh38、入力で指定された
version付きRefSeq transcriptである。

## 取得とtranscript選択

VEP HGVS endpointを `refseq=1`, `protein=1`, `transcript_version=1` で呼び、
`transcript_consequences[].transcript_id` が入力のRefSeq accession/versionと完全一致する
1要素だけを採用する。別isoformのconsequenceや蛋白座標は混合しない。

annotation Evidenceには次を保存する。

- `protein_id`: VEP `protein_id`（RefSeq入力ではNP accession/version）
- `protein_start`, `protein_end`: VEPの1-based蛋白座標
- `ref_aa`, `alt_aa`: VEP `amino_acids`。synonymousの単一値はref/alt双方に保存
- `protein_length_change`: missense/synonymousは0、in-frame indelはVEPの変更peptide長差
- `hgvsp`: VEPが返したtranscript固有HGVS protein表現

stop-loss、frameshift、start-loss、stop-gainについて、変更tokenだけから完全な蛋白長差を
推測しない。この場合の `protein_length_change` はnullとする。

## 予測値

同じtranscript consequenceから以下を独立した `category=computational` Evidenceとして保存する。

| predictor | 保存値 | VEP指定 | 現在の判定利用 |
|---|---|---|---|
| AlphaMissense | `am_pathogenicity`, `am_class` | `AlphaMissense=1` | 不可 |
| REVEL | score（返る場合のみ） | `REVEL=1` | 不可 |
| SpliceAI | 4 delta scoreとその最大値 | `SpliceAI=2` | 不可 |
| Ensembl Compara conservation | RESTの数値 | `Conservation=1` | 不可 |

`source`, Ensembl release、取得時刻、transcript、利用データセット説明を保持する。
ただしEnsembl REST応答はREVEL backing file、SpliceAI model、conservation track/methodの版を
十分に識別しない。AlphaMissenseの公開class閾値もACMG evidence strengthの校正ではない。
このため全レコードを `calibration_eligible=false` とし、PP3/BP4 evaluatorは明示的に除外する。
値を取得できたことと、臨床判定に利用できることを同一視しない。

## 校正policyの要件

PP3/BP4へ渡すには、少なくともpredictor/model/data版、対象機序（proteinまたはsplicing）、
対象consequence、score domain、criterion別interval・strength、校正文献またはVCEP仕様を
固定したEvidence adapterが必要である。相関する複数predictorの票決は行わない。

一般的な候補資料は、missenseにはPejaver et al. 2022（PMID: 36413997）、splicingには
Walker et al. 2023（PMID: 37352859）である。ただし疾患・遺伝子別VCEP仕様がある場合は
その適用範囲と版を優先して選ぶ。

BP7はさらに、synonymous/noncoding variantが標準splice region外であることをtranscript exon
境界に対して判定するversion付きposition ruleと、採用する保存性trackの校正が必要である。
VEPの今回の応答だけではRefSeq exon境界とconservation methodを十分に特定できないため、
自動的な `synonymous_assessment` 生成は行わない。現在取得するSpliceAI・保存性値は、後続の
校正済みBP7 adapterへ渡すraw Evidenceである。

## PS1 ClinVar comparator

PS1は疾患入力を必須とせず、次のprotein-level Evidenceを自動評価する。

1. version付きRefSeq protein HGVSをClinVar ESearchのexact phraseで検索する。
2. VCVのcoding/protein HGVS、GRCh38 variant、集約GermlineClassificationとreview statusを取得する。
3. comparatorのcoding HGVSをEnsembl VEPで再注釈する。
4. 元変異とは異なるnucleotide variantで、同一protein上の変更残基が1個かつ同じ置換か確認する。
5. ClinVar集約分類がPathogenicで、criteria provided/expert panel/practice guideline、かつ
   conflictingでなく、両変異のSpliceAIが0.1以下ならPS1をMETとする。

疾患未指定時のMETは `assessment_scope=protein_level`,
`condition_assessment=NOT_EVALUATED` とし、疾患整合性とClinVar分類の循環参照確認を
`review_points` に残す。疾患が入力され不一致ならMANUAL_REVIEWとする。検索が完全で適格な
別variantがなければNOT_MET、検索・mapping・splice確認が不完全ならNOT_EVALUATEDまたは
MANUAL_REVIEWとし、陰性へ変換しない。
