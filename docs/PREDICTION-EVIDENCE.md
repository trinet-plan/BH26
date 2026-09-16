# Ensembl蛋白注釈・予測Evidence

更新: 2026-09-16。対象はEnsembl REST/VEP release 116とMyVariant.info経由のdbNSFP 4.8a、
GRCh38、入力で指定されたversion付きRefSeq transcriptである。

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

| predictor | 保存値 | VEP指定 | 判定利用 |
|---|---|---|---|
| AlphaMissense | `am_pathogenicity`, `am_class` | `AlphaMissense=1` | 不可 |
| REVEL | score（返る場合のみ） | `REVEL=1` | 不可 |
| SpliceAI | 4 delta scoreとその最大値 | `SpliceAI=2` | 校正側の版宣言がある場合のみ |
| Ensembl Compara conservation | RESTの数値 | `Conservation=1` | 不可 |

`source`, Ensembl release、取得時刻、transcript、利用データセット説明を保持する。
ただしEnsembl REST応答はREVEL backing file、SpliceAI model、conservation track/methodの版を
十分に識別しない。AlphaMissenseの公開class閾値もACMG evidence strengthの校正ではない。
このためVEP由来の全レコードを `calibration_eligible=false` とする。
値を取得できたことと、臨床判定に利用できることを同一視しない。

## dbNSFP（MyVariant.info経由）

PP3/BP4に使う予測値はこちらから取得する。`/v1/metadata` の `src.dbnsfp.version` で
dbNSFP release（現在4.8a）を確定し、各スコアの `predictor_version` に固定したうえで
`calibration_eligible=true` とする。

- 採用: REVEL、AlphaMissense（いずれもmechanism=protein）
- 不採用: SIFT、PolyPhen-2。ClinGenの校正推奨対象ではないため基準入力に渡さない
- dbNSFPはnonsynonymous SNV対象のため、indelには問い合わせない
- 404は「値なし」として扱い、通信失敗とは区別する
- transcript単位のリスト値は先頭値を採用し、欠損を0へ変換しない

## 校正policy

校正区間は `config/demo-rules.json` の `computational` に置く。1つのcalibrationは
predictor/version、出典と版、対象機序、対象consequence、score domain、criterion別の
interval・strengthを固定する。相関する複数predictorの票決は行わない。

`selected_calibrations` に複数のcalibrationを並べられる。現在の設定は次の2つである。

| calibration | 機序 | 出典 | 区間 |
|---|---|---|---|
| `revel-pejaver-2022` | protein | Pejaver et al. 2022（PMID: 36413997） | PP3 >=0.644 / 0.773 / 0.932、BP4 <=0.290 / 0.183 / 0.016 / 0.003 |
| `spliceai-svi-2023` | splicing | Walker et al. 2023（PMID: 37352859） | PP3 >=0.2 supporting、BP4 <=0.1 supporting |

機序は加算しない。PP3はいずれかの機序が区間を満たせば成立し、最も強い区間を採用する。
BP4は適用対象の全機序が良性側を示す必要があり、splice影響が予測される場合は成立させず、
阻止した機序を `benign_blocked_by` に記録する。低い蛋白スコアはsplice部位の破壊について
何も言わないためである。

Ensembl VEPはSpliceAIのモデル版を公開しないため、calibrationが `version_assertion`
（source、unreported_version、asserted_version、asserted_by、justification）を明示した
場合に限り使用する。宣言は結果のprovenanceへ `version_assertions` として残り、
宣言がなければそのスコアは判定に使わない。疾患・遺伝子別VCEP仕様がある場合は、
その適用範囲と版を優先して選ぶ。

BP7はさらに、synonymous/noncoding variantが標準splice region外であることをtranscript exon
境界に対して判定するversion付きposition ruleと、採用する保存性trackの校正が必要である。
VEPの今回の応答だけではRefSeq exon境界とconservation methodを十分に特定できないため、
自動的な `synonymous_assessment` 生成は行わない。現在取得するSpliceAI・保存性値は、後続の
校正済みBP7 adapterへ渡すraw Evidenceである。

## PM5 ClinVar residue comparator

PM5は同一残基の別アミノ酸変化を対象とするため、PS1のexact protein change検索では
「他の置換が報告されていない」ことを示せない。遺伝子単位のmissense検索とesummaryで
同一残基の候補を絞り、候補ごとにefetchとEnsembl照合で蛋白参照・残基・参照アミノ酸の一致と
置換の相違を確認する。検索記録には `search_scope`（`exact_protein_change` / `residue`）を
持たせ、PM5は `residue` の完了記録だけを不成立の根拠として受け付ける。

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
