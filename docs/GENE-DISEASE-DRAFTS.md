# Gene–disease機序のDRAFT生成

PP2・BP1・PVS1に使う最終`gene_disease` Evidenceは、疾患別にreview済みでなければならない。
一方、レビュー対象を探す段階では、InterVarの固定リストやfastVEP型のconstraint・変異分布集計が
有用である。この実装では両者を混同せず、`GeneDiseaseDraftProvider`がレビュー専用の
`gene_disease_draft`を生成する。

## 入力

- ClinGen Gene-Disease Validity相当の、遺伝子・疾患・classification・版・取得日時
- gnomAD gene constraint相当のmisZ・pLI・LOEUF・版・取得日時
- ClinVar gene-wide spectrum相当の病的missense、病的truncating、良性missense件数
- 版付きの候補抽出policy。閾値はProviderコードへ固定しない

ClinVar検索が途中で打ち切られた場合は`complete=false`とする。0件を「missense機序なし」と
解釈せず、BP1候補を`INSUFFICIENT`にする。

## 出力と安全境界

出力は常に次の属性を持つ。

```json
{
  "category": "gene_disease_draft",
  "quality_status": "DRAFT",
  "assessment_status": "DRAFT",
  "criterion_eligible": false
}
```

PP2・BP1・PVS1は`category=gene_disease`かつprovenanceが完全なEvidenceだけを読むため、DRAFTは
criterionを成立させない。DRAFTの`CANDIDATE`はレビュー順序を示すだけである。

## レビュー後の昇格

レビュー担当者は原資料と疾患文脈を確認し、少なくとも以下を明示する。

- `missense_mechanism_established`
- `low_benign_missense_variation`
- `predominantly_truncating`
- `lof_mechanism_established`
- `spectrum_review_complete`
- `curator`と`reviewed_at`

これらを別の`gene_disease` Evidenceとして保存する。DRAFTレコードのcategoryを書き換えて流用せず、
レビュー済み成果物を新しく作成する。

## デモ39判定への投入（2026-09-16）

実デモでは次の版付き入力を追加した。

- `config/gene-disease-source-data.json`: ClinGen GDV、gnomAD v4.1.1 constraint、
  2026-09-16取得のClinVar gene-wide件数
- `config/gene-disease-review-decisions.json`: 疾患別の機序・criterion applicability判断
- `config/curated-context.json`: 28 ALTレコードをrecord単位でHCMまたはARVCへ対応付け

同じvariant keyがHCM症例とARVC症例に現れるため、variant単位のcontextよりrecord単位の
`record_contexts`を優先する。これにより、TNNI3–HCMのmissense機序をARVCへ流用しない。

```powershell
python -c "import sys; sys.path.insert(0, 'src'); from acmg.cli import main; main()" `
  build-gene-disease-drafts `
  --input config/gene-disease-source-data.json `
  --output work/gene-disease/drafts.json

python -c "import sys; sys.path.insert(0, 'src'); from acmg.cli import main; main()" `
  build-gene-disease-evidence `
  --input config/gene-disease-review-decisions.json `
  --base-evidence work/check-20260915/prepared/evidence.json `
  --output work/gene-disease/evidence.json
```

生成されるcriterion用Evidenceは`assessment_method=automated`、`human_signoff=false`を明示する。
これは出典固定されたデモ評価であり、人間による臨床承認を装わない。ClinGen VCEP仕様がある場合は
数値triageよりVCEPを優先する。特にMYH7 VCEP v2.0.0はPP2・BP1・PVS1をNot Applicableとするため、
gnomAD missense Zが高くてもPP2を成立させない。

再評価した元の39件は`gene_disease`不足が0件となり、内訳はMET 1、NOT_MET 19、
NOT_APPLICABLE 16、NOT_EVALUATED 3である。残る3件はすべてMYBPC3 PVS1で、次の入力待ちへ進んだ。

- case1 frameshift: `transcript_assessment`
- case2 canonical splice: `splice_assessment`
- case2 frameshift: `transcript_assessment`
