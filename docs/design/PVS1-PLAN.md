# PVS1 General decision tree

更新: 2026-09-16。ACMG/AMP 2015、ClinGen SVI PVS1 2018（PMID: 30192042）、ClinGen SVI
Splicing 2023（PMID: 37352859）のGeneral Guidanceを対象とする。VCEP/CSpec、tool固有score、
gene固有例外は適用しない。

## 対象範囲

現行のGRCh38 SNV・小規模indel入力で表現できるstop gained、frameshift、canonical splice
donor/acceptor、RNA assayで実証されたsplice LoF、start lostを評価する。single/multi-exon
deletion、duplication、その他CNV/SVは対象外である。NMD、MANE、RNA、region情報を外部から
自動取得せず、取得層と判定層を分離した正規化Evidenceとして受け取る。

## Evidence契約

全Evidenceはvariant key、evidence ID、source/version、retrieved time、PASS qualityに加え、
人手reviewのcurator/reviewed_at、またはversion付きautomated policyを持つ。

| category | 主なフィールド |
| --- | --- |
| `gene_disease` | `gene`, optional `condition`, `lof_mechanism_established` |
| `transcript_assessment` | `transcript`, `relevance`, `exon_relevance` |
| `nmd_prediction` | `predicted`, `exon`, `distance_to_final_junction`, `rule_source` |
| `splice_assessment` / `rna_assay` | `splice_outcome`, `reading_frame_disrupted`, `alternative_rescue` |
| `protein_region` | `critical_region_disrupted`, `region_biologically_relevant`, `lost_residues`, `total_protein_length` |
| `population_lof` | `lof_variants_frequent` |
| `initiation_assessment` | `intact_alternative_transcript`, `downstream_in_frame_start`, `upstream_pathogenic_evidence` |

conditionありでは一致するmechanismを優先し、それが存在しない場合だけgene-levelへfallbackする。
明示的なfalse/unknownをfallbackで上書きしない。conditionなしではcondition-specific Evidenceを
流用しない。複数の矛盾やgene/transcript不一致はMANUAL_REVIEWとする。

## 判定

- LoF mechanismなし/unknownはNOT_EVALUATED、明示的な非LoF機序はNOT_APPLICABLE。
- NMDを受けるtruncating variantでtranscript/exonがrelevantならPVS1 very strong。
- NMD escapeでcritical regionを破壊すればstrong。population LoF頻発または非関連regionなら
  NOT_APPLICABLE。それ以外はprotein lossが10%超ならstrong、10%以下ならmoderate。
- spliceはrescue、frame、NMDを明示評価する。RNA実証Evidenceには`used_by: [PVS1]`を記録する。
- start-lossはintact alternative transcriptがあればNOT_APPLICABLE。downstream start評価後、
  上流の病原性Evidenceがあればmoderate、なければsupporting。
- missingはNOT_EVALUATED、矛盾はMANUAL_REVIEWとし、absenceやfalseへ変換しない。

閾値、ruleset名・版・一次資料は`config/demo-rules.json`の`PVS1`に置く。各結果は
`evaluation_context`、`decision_trace`、`rules_used`、`warnings`、`unresolved_requirements`を持つ。
condition未指定でもgene-level mechanismがあればMETを許容するが、その評価範囲をwarningと
mechanism scopeへ必ず記録する。

## 出力と検証

METはPVS1、PVS1_strong、PVS1_moderate、PVS1_supportingとしてVA-Spec 1.0.1 Evidence Lineへ
変換する。workflow状態は標準Evidence Lineへ出力しない。合成テストでcondition解決、全強度、
全停止分岐、10%境界、splice/RNA、start-loss、trace/provenanceを検証する。固定キャッシュの
demo-dataはPVS1用curated Evidenceを含まないため、25件NOT_APPLICABLE・3件NOT_EVALUATEDである。
