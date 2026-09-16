# ACMG criterion判定ロジック・入力Evidenceガイド

更新: 2026-09-16。本書は、現在の実装が各criterionについて「何を入力として見て、どの順序で、
どのstatus・strengthを返すか」をコード準拠で説明する。ACMG/ClinGen文書の一般的な解説ではなく、
このリポジトリの実際の評価動作を監査するための資料である。

## 1. 評価範囲

対象はGRCh38上のSNV・小規模indelで、次の16 criterionを評価する。

| 系統 | criterion |
| --- | --- |
| 病原性・非常に強い | PVS1 |
| 病原性・強い | PS1 |
| 病原性・中程度 | PM1、PM2、PM4、PM5 |
| 病原性・支持 | PP2、PP3、PP5 |
| 良性・単独 | BA1 |
| 良性・強い | BS1 |
| 良性・支持 | BP1、BP3、BP4、BP6、BP7 |

PP5とBP6はClinGen General Guidanceに従い、常に`DEPRECATED`である。外部機関の分類そのものを
証拠として加点せず、その分類を支える一次Evidenceを各criterionへ入力する。

本処理はcriterion単位のEvidence評価であり、Pathogenic/Likely pathogenic/VUS/Likely benign/
Benignの最終5段階分類は行わない。CNV/SV、exon deletion/duplicationも対象外である。

## 2. statusの意味

| status | 意味 |
| --- | --- |
| `MET` | 必要なEvidenceと規則が揃い、criterionが成立した |
| `NOT_MET` | 必要な検索・評価が完了し、成立条件を満たさなかった |
| `NOT_EVALUATED` | Evidence、設定、検索完了性、または必要な文脈が不足している |
| `NOT_APPLICABLE` | consequenceや明示的な機序評価から、このcriterionの対象外と判断できる |
| `MANUAL_REVIEW` | 複数Evidenceの競合、対応関係の不一致、RNA矛盾など、人による解決が必要 |
| `DEPRECATED` | 現行方針で使用しないcriterion |

重要なのは、欠測を`NOT_MET`へ変換しないことである。`NOT_MET`は原則として、検索が完了している、
または成立を否定する信頼できる観測がある場合にだけ返す。

## 3. 全criterion共通の入力・品質ゲート

### 3.1 variantとannotation

全criterionの公開`evaluate()`は、共通契約の
`acmg_pipeline.vcf_record.VariantRecord`と
`acmg_pipeline.clinical_note.ClinicalNoteExtraction`を最初の2引数として受け取る。
CLIのprepared JSONも、criterionを呼ぶ前に必ずこの2クラスへ変換される。

入力variantの`chrom/pos/ref/alt`を内部の正規化variantへ変換する。assemblyはGRCh38のみを
サポートし、`VariantRecord`自体にはassembly欄がないため、省略時はGRCh38、明示時は
`info["ASSEMBLY"]`を使う。REF/ALTは空でないACGT配列で、symbolic ALTやSVは受け付けない。
`ClinicalNoteExtraction`はproband、family、de novo情報を保持し、現在の16 criterionでは
familyの`inheritance_pattern`を遺伝形式の補助入力として利用できる。記載のない臨床情報を
推測して補わない。

多くのcriterionは最初に`annotation` Evidenceを取得する。annotationには少なくとも次が必要である。

- 評価variantと一致する`variant_key`
- version付き`transcript`
- `consequences`
- `gene`、`protein_id`、蛋白座標、アミノ酸など、criterion固有の注釈
- `evidence_id`、`source`、`source_version`、`retrieved_at`
- `quality_status: PASS`

該当annotationが0件なら`NOT_EVALUATED`、同じ評価文脈に複数あればtranscriptを勝手に選ばず
`MANUAL_REVIEW`とする。annotationは主にEnsembl VEPから取得するが、VEPのconsequenceだけで
遺伝子疾患機序、NMD、疾患関連transcript、critical regionなどを推測しない。

### 3.2 Evidenceの選択

Evidence serviceは次を満たすrecordだけを候補にする。

1. categoryが要求されたものと一致する。
2. `variant_key`が完全一致する。
3. `quality_status`が`PASS`である。
4. `evidence_id/source/source_version/retrieved_at`が揃っている。
5. 通常の取得ではtranscriptとconditionも評価文脈に一致する。condition付きの
   Evidenceが入力conditionと異なる場合は除外し、conditionを持たないgene/protein-level
   Evidenceは候補に残す。

人手評価を要求するcategoryでは、次のいずれかが必要である。

- 人手: `curator`と`reviewed_at`
- 自動: `assessment_method: automated`、`method`、`policy_version`

### 3.3 conditionの扱い

conditionは基準ごとに扱いが異なる。

- BS1は疾患別頻度閾値を使うためcondition必須。
- PVS1はcondition-specific LoF機序を優先し、なければgene-levelへfallbackする。
- PS1/PM5、PM1、PP2/BP1は蛋白または遺伝子レベルで評価可能。別疾患に明示的に
  紐づくEvidenceを流用はしないが、condition-agnostic Evidenceによる成立は許し、
  `condition_assessment`とreview pointに疾患文脈の未確定を残す。
- computational evidenceやBA1は原則としてvariant/populationレベルで評価する。

### 3.4 Evidence category辞書

| category | 利用criterion | 判定に使う主な内容 |
| --- | --- | --- |
| `annotation` | deprecated以外の多く | consequence、gene、version付きtranscript、protein ID/座標、AA変化、蛋白長変化 |
| `population` | PM2、BA1、BS1 | provider/cohort別AC・AN・AF、callability、リリース、取得状態 |
| `comparator` | PS1、PM5 | 既知variant、protein ID/残基/置換、classification、review状態、splice機序照合 |
| `comparator_search` | PS1、PM5 | 蛋白全体または同一残基検索が完了したか |
| `region` | PM1、PM4、BP3 | protein範囲、`hotspot`/`critical_domain`、反復領域、functional importance、病的/良性集積 |
| `gene_disease` | PVS1、PP2、BP1 | gene、condition、LoF/missense/truncating機序、variant spectrum |
| `computational` | PP3、BP4 | predictor、version、mechanism、score、calibration eligibility |
| `synonymous_assessment` | BP7 | splice-critical位置、splice影響、保存性、RNA矛盾とそれぞれの出典 |
| `transcript_assessment` | PVS1 | transcript relevance、exon relevanceと判定根拠 |
| `nmd_prediction` | PVS1 | NMD予測boolean、exon/最終junction情報、適用rule source |
| `splice_assessment` | PVS1 | splice outcome、frame disruption、alternative rescue |
| `rna_assay` | PVS1 | RNAで確認されたLoF、assay provenance、`used_by` |
| `protein_region` | PVS1 | criticality、biological relevance、lost residues、total protein length |
| `population_lof` | PVS1 | 対象exon/regionでLoF variantが頻発するか |
| `initiation_assessment` | PVS1 | intact alternative transcript、downstream in-frame start、upstream pathogenic Evidence |

どのcategoryでも、上表の値だけでなく共通provenanceが必要である。例えば`score: 0.9`
だけのrecordはPP3のEvidenceにはならず、predictor、version、出典、校正policyとの一致が要る。

## 4. 一覧表

| criterion | 主に見る情報 | 対象 | METの中心条件 | strength |
| --- | --- | --- | --- | --- |
| PVS1 | LoF機序、transcript/exon、NMD、splice/RNA、region | truncating、splice、start-loss | decision treeの成立終端 | very strong～supporting |
| PS1 | 同一蛋白置換の病原性comparator | missense | 別塩基変化が同じAA置換を生じる | strong |
| PM1 | hotspot密度またはcritical domain | 主に蛋白変化 | 病原性集積＋良性枯渇、またはreview済みcritical domain | moderate |
| PM2 | population AF、AN、callability | 全variant | 全有効観測が設定AF以下、検索失敗なし | supporting |
| PM4 | 蛋白長変化、repeat機能 | in-frame indel、stop-loss | 長さが変化し、nonfunctional repeatではない | moderate |
| PM5 | 同一残基の異なる病原性置換 | missense | 同一残基に異なるAAへの病原性comparator | moderate |
| PP2 | missense疾患機序とvariant spectrum | missense | missense機序あり＋良性missenseが少ない | supporting |
| PP3 | 校正済み計算予測 | 校正scope内 | protein/spliceのいずれかが閾値を満たす | supporting～strong等 |
| PP5 | 使用しない | ― | 常にDEPRECATED | ― |
| BA1 | population AF、BA1例外 | 全variant | AF>5%かつ例外でない | stand-alone |
| BS1 | 疾患別最大credible AF | condition必須 | いずれかの有効AFが疾患閾値を超える | strong |
| BP1 | 疾患variant spectrum | missense | truncating優位＋missense機序でない | supporting |
| BP3 | repeatと機能的重要性 | in-frame indel | repetitiveかつfunctional importanceなし | supporting |
| BP4 | 校正済み計算予測 | 校正scope内 | 適用済み全機序が良性側閾値を満たす | supporting～very strong |
| BP6 | 使用しない | ― | 常にDEPRECATED | ― |
| BP7 | splice位置、splice予測、保存性、RNA | synonymous | splice-critical外＋splice影響なし＋非保存 | supporting |

### 4.1 適用strengthと出力outcome

| criterion | MET時の適用strength | VA-Spec outcome |
| --- | --- | --- |
| PVS1 | `very_strong` / `strong` / `moderate` / `supporting` | `PVS1` / `PVS1_strong` / `PVS1_moderate` / `PVS1_supporting` |
| PS1 | `strong` | `PS1` |
| PM1 | `moderate` | `PM1` |
| PM2 | `supporting` | `PM2_supporting` |
| PM4 | `moderate` | `PM4` |
| PM5 | `moderate` | `PM5` |
| PP2 | `supporting` | `PP2` |
| PP3 | 校正bandが指定するstrength | 既定strengthなら`PP3`、調整時はstrength suffix付き |
| BA1 | `stand_alone` | `BA1` |
| BS1 | `strong` | `BS1` |
| BP1 / BP3 / BP7 | `supporting` | 各criterion code |
| BP4 | 校正bandが指定するstrength | 既定strengthなら`BP4`、調整時はstrength suffix付き |

PM2はACMG原文名の`PM`に関わらず、現在の実装policyではsupportingである。本書は基準名から
strengthを推測せず、実際の出力を記載する。

### 4.2 criterionごとの到達可能status

| criterion | MET | NOT_MET | NOT_EVALUATED | NOT_APPLICABLE | MANUAL_REVIEW | DEPRECATED |
| --- | :---: | :---: | :---: | :---: | :---: | :---: |
| PVS1 | ○ | — | ○ | ○ | ○ | — |
| PS1 | ○ | ○ | ○ | ○ | ○ | — |
| PM1 | ○ | ○ | ○ | — | ○ | — |
| PM2 | ○ | ○ | ○ | — | — | — |
| PM4 | ○ | ○ | ○ | ○ | ○ | — |
| PM5 | ○ | ○ | ○ | ○ | ○ | — |
| PP2 | ○ | ○ | ○ | ○ | ○ | — |
| PP3 | ○ | ○ | ○ | ○ | ○ | — |
| PP5 | — | — | — | — | — | ○ |
| BA1 | ○ | ○ | ○ | — | ○ | — |
| BS1 | ○ | ○ | ○ | — | ○ | — |
| BP1 | ○ | ○ | ○ | ○ | ○ | — |
| BP3 | ○ | ○ | ○ | ○ | ○ | — |
| BP4 | ○ | ○ | ○ | ○ | ○ | — |
| BP6 | — | — | — | — | — | ○ |
| BP7 | ○ | ○ | ○ | ○ | ○ | — |

PVS1に`NOT_MET`がないのは意図的である。明示的な非LoF機序、非関連transcript、rescueなどは
`NOT_APPLICABLE`、必要情報不足は`NOT_EVALUATED`とし、不完全な情報から「PVS1不成立」
という陰性結論を作らない。

### 4.3 判定順序の読み方

各モジュールはおおむね「適用対象→必要Evidence→整合性→成立条件」の順に評価する。
そのため、例えば非missense変異のPP2は、gene-disease Evidenceがなくても先に
`NOT_APPLICABLE`となる。一方、対象consequenceでEvidenceがなければ`NOT_EVALUATED`である。
詳細節の分岐は記載順がその優先順序を表す。

## 5. criterion別の詳細

### 5.1 PVS1 — predicted loss of function

#### 使用情報

| Evidence category | 主なフィールド |
| --- | --- |
| `annotation` | gene、transcript、consequences |
| `gene_disease` | gene、condition、`lof_mechanism_established` |
| `transcript_assessment` | transcript、`relevance`、`exon_relevance` |
| `nmd_prediction` | `predicted`、exon、final junction距離、`rule_source` |
| `splice_assessment` / `rna_assay` | rescue、splice outcome、reading frame、RNA LoF実証 |
| `protein_region` | criticality、biological relevance、lost/total residues |
| `population_lof` | affected exon/regionでLoFが頻発するか |
| `initiation_assessment` | intact alternative transcript、downstream start、上流病原性Evidence |

設定`PVS1.ruleset`には名称、版、出典一覧が必要で、蛋白喪失割合閾値は
`PVS1.rules.protein_loss_threshold`から読む。現行demo設定は0.10である。

#### 共通入口

1. annotationからvariant typeを決める。
   - `stop_gained`、`frameshift_variant`
   - canonical `splice_donor_variant` / `splice_acceptor_variant`
   - `start_lost`
   - 非canonicalでもRNAでLoF実証済みなら`SPLICE_LOF_CONFIRMED`
2. 上記以外は`NOT_APPLICABLE`。
3. `transcript_ablation`はLoF候補だが実装済み経路がないため`NOT_EVALUATED`。
4. LoF機序Evidenceを解決する。
   - conditionあり: 完全一致を優先し、存在しない場合だけgene-levelを使う。
   - conditionなし: gene-levelだけを使い、他疾患のcondition-specific Evidenceは借用しない。
   - gene不一致または複数の競合は`MANUAL_REVIEW`。
   - Evidenceなし・unknownは`NOT_EVALUATED`。
   - `lof_mechanism_established=false`は`NOT_APPLICABLE`。

#### nonsense / frameshift

```text
transcript relevance
├─ NOT_RELEVANT → NOT_APPLICABLE
├─ unknown/missing → NOT_EVALUATED
└─ RELEVANT
   └─ NMD prediction
      ├─ unknown/missing → NOT_EVALUATED
      ├─ true
      │  ├─ exon RELEVANT → MET / very_strong
      │  ├─ exon NOT_RELEVANT → NOT_APPLICABLE
      │  └─ exon unknown → NOT_EVALUATED
      └─ false → NMD escape region path
```

NMD escape region pathは次の順に評価する。

1. `critical_region_disrupted=true`なら`MET / strong`。
2. それ以外で`lof_variants_frequent=true`なら`NOT_APPLICABLE`。
3. regionが明示的に非関連なら`NOT_APPLICABLE`、不明なら`NOT_EVALUATED`。
4. `lost_residues / total_protein_length`を計算する。
   - 不正値は`MANUAL_REVIEW`。
   - 割合が閾値より大きい: `MET / strong`。
   - 割合が閾値以下: `MET / moderate`。

`critical_region_disrupted`と`lof_variants_frequent`は、明示的にtrueの場合にそれぞれの
分岐を確定する。現在のコードでは、これらのrecordがないことだけで分岐を止めず、
biological relevanceとprotein lossを評価できれば先へ進む。

境界は`> 0.10`がstrongであり、ちょうど0.10はmoderateである。

#### splice

1. RNAでLoFが実証されていればRNA assessmentを優先し、`used_by`にPVS1を記録する。
2. `alternative_rescue=true`なら`NOT_APPLICABLE`。
3. rescue不明、splice outcome不明は`NOT_EVALUATED`。
4. outcomeと`reading_frame_disrupted`が矛盾すれば`MANUAL_REVIEW`。
5. out-of-frame/PTCは上記NMD経路へ進む。
6. in-frameはtranscript relevance確認後、NMD escape region pathへ進む。

#### start-loss

1. biologically relevantなintact alternative transcriptあり → `NOT_APPLICABLE`。
2. alternative transcript status不明 → `NOT_EVALUATED`。
3. downstream in-frame startがbooleanでなければ`NOT_EVALUATED`。
4. upstream pathogenic Evidenceあり → `MET / moderate`、なし → `MET / supporting`。

現在の実装では、`downstream_in_frame_start`は「評価済みbooleanであること」を要求するが、true/falseの
値自体はstrength分岐に使用していない。これは現在のコード動作を示す重要な注意点であり、今後の
PVS1 start-loss仕様精緻化候補である。

PVS1だけは`evaluation_context`、C01/G01/G02/V01/NFxx/SPxx/ICxxの`decision_trace`、
`rules_used`、`warnings`、`unresolved_requirements`を詳細出力する。

### 5.2 PS1 — same amino acid change

#### 使用情報

- annotation: missense、transcript、protein ID、AA position、ref/alt AA
- `comparator`: 別genomic variant、同一protein change、Pathogenic分類、review status、
  independence、splice影響、condition、一次Evidence
- `comparator_search`: 検索scopeと完了性

#### 判定

1. missense以外は原則`NOT_APPLICABLE`。splice consequenceは等価性を自動判定せず
   `MANUAL_REVIEW`。
2. protein ID/position/ref/alt AA不足は`NOT_EVALUATED`。
3. 同じprotein、残基、ref AAで、別塩基variantが同じalt AAを生じる候補だけを残す。
4. comparator自身と同一variant、transcript不一致、Pathogenic以外を除外する。
5. 次のどちらかを満たすcomparatorをqualifiedとする。
   - 自動: exact protein match、different nucleotide variant、eligible review status、
     splice checked、splice conflictなし。
   - review済み: condition一致、pathogenic Evidence確認、独立Evidence、機序一致、splice確認、
     curator/review日/一次Evidenceあり。
6. qualified comparatorあり → `MET / strong`。
7. conditionが指定され、qualified comparatorの疾患が一致しない → `MANUAL_REVIEW`。
8. 病原性候補はあるが資格確認が不十分 → `MANUAL_REVIEW`。
9. 完了したexact protein-change検索でeligible候補なし → `NOT_MET`。検索不完了 →
   `NOT_EVALUATED`。

condition未指定でもprotein-levelではMETになり得るが、疾患関連性と循環的PS1利用の確認を
review pointへ残す。

### 5.3 PM1 — hotspot / critical functional domain

#### 使用情報

- annotation: transcript、protein ID、altered protein interval
- `region`: interval、`region_type`、hotspot countsまたはcriticality、良性変異枯渇
- `PM1.hotspot`設定: window、病原性/良性count閾値、method、policy版・出典

regionがaltered protein interval全体を包含する必要がある。部分重複は`MANUAL_REVIEW`。

#### mutational hotspot route

現行demo policyは±5 aa window、Pathogenic/Likely pathogenic 3件以上、Benign/Likely benign 0件以下、
policy version `PM1-hotspot-v1`である。

- policy/method不一致、count不正 → `NOT_EVALUATED`。
- benign countが閾値内だがpathogenic count不足 → 情報不足として`NOT_EVALUATED`。
- pathogenic enrichmentあり、benign depletionあり → `MET / moderate`。
- それ以外で評価が完了している → `NOT_MET`。
- record内の派生booleanとcount計算が矛盾 → `MANUAL_REVIEW`。

#### critical functional domain route

自動生成Evidenceだけではcriticalityを確立できず`NOT_EVALUATED`。人手review済みEvidenceで
`critical_functional_region=true`かつ`benign_depletion=true`なら`MET / moderate`、それ以外は
`NOT_MET`。

condition一致は記録するが必須ではない。protein-levelでMETかつcondition未確認ならreview pointを
付ける。

### 5.4 PM2 — population rarity

#### 使用情報

- population Evidence: AF、AC、AN、population、callability、quality、provider/release
- `PM2.minimum_an`
- `PM2.max_af`
- policy source/version

population観測は次を満たす場合だけ使用する。

- variant一致、provenance完備、quality PASS
- ANが整数でminimum AN以上
- `0 <= AC <= AN`、`0 <= AF <= 1`
- AC/ANとAFの差が1e-6以内
- AC=0の場合は`callable=true`

判定は重複cohortを合算せず、有効観測中の最大AFを使う。

- 最大AF `> max_af` → `NOT_MET`。他providerが失敗していても反例が優先される。
- 全AFが閾値以下だがprovider失敗あり → `NOT_EVALUATED`。
- 全AFが閾値以下で検索完了 → `MET / supporting`。
- 有効観測なし → `NOT_EVALUATED`。

比較は`<=`を希少側とする。現行demoの`max_af=0`はデモ設定であり、疾患別臨床閾値ではない。

### 5.5 PM4 — protein length change

#### 使用情報

- annotation: in-frame insertion/deletionまたはstop-loss、protein interval、
  `protein_length_change`
- review済み`region`: altered interval、`nonfunctional_repeat`、`functional_review_complete`

#### 判定

- in-frame indel/stop-loss以外 → `NOT_APPLICABLE`。
- regionがaltered intervalを包含しない → `MANUAL_REVIEW`。
- functional review未完了 → `MANUAL_REVIEW`。
- protein length change不明 → `NOT_EVALUATED`。
- length changeが0以外、かつnonfunctional repeatではない → `MET / moderate`。
- それ以外 → `NOT_MET`。

### 5.6 PM5 — different missense change at same residue

PS1と同じcomparator evaluatorを使うが、同じalt AAではなく、同じprotein残基における異なるalt AAを
要求する。自動Evidenceでは`residue_match=true`が必要である。

- qualifiedな異なる病原性AA置換あり → `MET / moderate`。
- 病原性候補はあるが資格確認不十分 → `MANUAL_REVIEW`。
- 完了したresidue-scoped検索でeligible候補なし → `NOT_MET`。
- exact protein-change検索だけでは他のAA変化がないことを証明できないため、PM5では検索不完了扱い。
- conditionの扱い、splice確認、independence要件はPS1と同じ。

### 5.7 PP2 — missense mechanism

#### 使用情報

- annotation: missense、gene、transcript
- review済み`gene_disease`:
  - `missense_mechanism_established`
  - `spectrum_review_complete`
  - `low_benign_missense_variation`

#### 判定

- missense以外 → `NOT_APPLICABLE`。
- gene不一致 → `MANUAL_REVIEW`。
- 必須boolean不足 → `NOT_EVALUATED`。
- spectrum review未完了 → `MANUAL_REVIEW`。
- missense機序が確立し、良性missense variationが少ない → `MET / supporting`。
- それ以外 → `NOT_MET`。

pLI、LOEUF、missense Z scoreなどのconstraint値だけではPP2を成立させない。condition未確認でも
gene-levelでMETになり得るが、疾患文脈確認をreview pointに残す。

### 5.8 PP3 — calibrated pathogenic computational evidence

#### 使用情報

- annotation consequence
- `computational` Evidence: predictor、predictor version、mechanism、score、source
- version固定されたcalibration: consequence scope、score domain、strength bands

現行設定の主なinclusive閾値は次のとおり。

| predictor | mechanism | supporting | moderate | strong |
| --- | --- | ---: | ---: | ---: |
| REVEL / dbNSFP 4.8a | protein | >=0.644 | >=0.773 | >=0.932 |
| SpliceAI / asserted Ensembl source | splicing | >=0.20 | ― | ― |

1. consequenceが選択calibrationのscope外なら`NOT_APPLICABLE`。
2. predictor/version/mechanismが一致するEvidenceだけを使用する。VEPがモデル版を公開しない
   SpliceAIは、設定内のversion assertionが完全な場合だけ使用する。
3. 同一calibrationに複数の異なるscoreがあれば`MANUAL_REVIEW`。
4. score domain外、band設定なし、適合scoreなしは`NOT_EVALUATED`。
5. 適用できたprotein/splicing機序のいずれかがbandを満たせば、最強strengthで`MET`。
6. 適用できたものがすべてband外なら`NOT_MET`。

異なる機序のEvidence strengthは加算しない。未取得predictorはprovenanceに残るが、別の適用済み
predictorがPP3を満たせばMETを妨げない。

### 5.9 PP5 — deprecated

常に`DEPRECATED`。外部sourceの「Pathogenic」という主張を直接加点しない。

### 5.10 BA1 — high population frequency

#### 使用情報

- PM2と同じ品質条件を満たすpopulation Evidence
- `BA1.minimum_an`とpolicy provenance
- version付き・review済みBA1 exception assessment

#### 判定

- 有効なAFが1件も5%を超えない:
  - provider失敗あり → `NOT_EVALUATED`
  - 検索完了 → `NOT_MET`
- AF `> 0.05`だが例外確認なし → `NOT_EVALUATED`
- BA1例外である → `NOT_MET`
- 例外status不明 → `MANUAL_REVIEW`
- AF `> 0.05`かつ例外でない → `MET / stand-alone`

5%比較はstrictな`>`で、ちょうど5%はBA1成立ではない。現在この0.05はコード固定値である。

### 5.11 BS1 — AF above disease threshold

#### 使用情報

- conditionとinheritance
- review済み`disease_frequency_threshold`: condition、inheritance、`max_credible_af`、source/version
- 品質条件を満たすpopulation Evidence

#### 判定

- conditionまたは適合する疾患閾値なし → `NOT_EVALUATED`。
- threshold provenance不足 → `NOT_EVALUATED`。
- inheritance不一致 → `MANUAL_REVIEW`。
- いずれかの有効AFが閾値を超える → `MET / strong`。
- 超えずprovider失敗あり → `NOT_EVALUATED`。
- 超えず検索完了 → `NOT_MET`。

比較はstrictな`>`で、ちょうど閾値と同値ならBS1成立ではない。信頼できる高頻度観測があれば、
別providerの失敗があってもMETになる。

### 5.12 BP1 — missense in truncating-mechanism gene

PP2と同じgene-disease spectrum evaluatorを使う。

- missense以外 → `NOT_APPLICABLE`。
- spectrum review完了後、`predominantly_truncating=true`かつ
  `missense_mechanism_established=false` → `MET / supporting`。
- それ以外 → `NOT_MET`。
- gene不一致、review未完了、欠測、condition review pointの扱いはPP2と同じ。

### 5.13 BP3 — in-frame indel in nonfunctional repeat

#### 使用情報

- annotation: in-frame insertion/deletion、protein interval
- review済み`region`: `repetitive`、`functional_importance`、`functional_review_complete`

#### 判定

- in-frame indel以外 → `NOT_APPLICABLE`。
- regionがaltered interval全体を包含しない → `MANUAL_REVIEW`。
- functional review未完了 → `MANUAL_REVIEW`。
- `repetitive=true`かつ`functional_importance=false` → `MET / supporting`。
- それ以外 → `NOT_MET`。

### 5.14 BP4 — calibrated benign computational evidence

PP3と同じcalibration engineを使う。現行設定のinclusive上限は次のとおり。

| predictor | mechanism | supporting | moderate | strong | very strong |
| --- | --- | ---: | ---: | ---: | ---: |
| REVEL / dbNSFP 4.8a | protein | <=0.290 | <=0.183 | <=0.016 | <=0.003 |
| SpliceAI / asserted Ensembl source | splicing | <=0.10 | ― | ― | ― |

PP3との違いは、適用済みの全機序が良性側を支持する必要があることである。例えばREVELが良性側でも、
適用済みSpliceAIが良性band外なら`NOT_MET`。ただし、取得できなかったpredictorはprovenanceに残るが、
現在はそれだけでMETを止めない。

- 適用済み全機序が良性bandを満たす → 最強strengthで`MET`。
- いずれかの適用済み機序がband外 → `NOT_MET`。
- annotationがhigh-confidence null/spliceなのに良性予測が成立 → `MANUAL_REVIEW`。
- version不一致、score競合、scope外等の扱いはPP3と同じ。

### 5.15 BP6 — deprecated

常に`DEPRECATED`。外部sourceの「Benign」という主張を直接加点しない。

### 5.16 BP7 — synonymous with no splice impact

#### 使用情報

- annotation: synonymous consequence、transcript
- review済み`synonymous_assessment`:
  - `outside_splice_critical_region`
  - `no_predicted_splice_impact`
  - `not_conserved`
  - `contradictory_rna_evidence`
  - splice prediction、calibration、conservation、position ruleのprovenance

#### 判定

- synonymous以外:
  - intronic/noncoding/UTR → General BP7 extensionを自動適用せず`MANUAL_REVIEW`
  - その他 → `NOT_APPLICABLE`
- 必須boolean不足またはprediction/position/conservation provenance不足 → `NOT_EVALUATED`
- contradictory RNA Evidenceあり → `MANUAL_REVIEW`
- splice-critical領域外、予測splice影響なし、非保存の3条件すべてtrue → `MET / supporting`
- いずれかfalse → `NOT_MET`

## 6. population Evidenceの共通注意

PM2、BA1、BS1は、population間・provider間のAC/ANを合算しない。重複個体を含み得るcohortを
合算すると見かけのANを水増しするためである。各観測を独立に品質確認し、criterionに応じて最大AF
または閾値超過の有無を使う。

APIやDBにvariantが見つからないことはAF=0ではない。ANとcallabilityがない未登録応答は
`NO_OBSERVATION`として扱い、原則`NOT_EVALUATED`へつながる。

## 7. criterion間の重複・競合表示

各criterionは独立に評価する。同時METを自動的に取り消したりstrengthを変更したりせず、次の組合せを
`REVIEW_OVERLAP:<left>:<right>`として両resultへ付与する。

- PVS1–PM4、PVS1–PP3
- PS1–PP3、PS1–PM5
- PM4–BP3
- PP3–BP4、PP3–BP7、BP4–BP7
- PP2–BP1
- PM2–BA1、PM2–BS1、BA1–BS1

このflagは「どちらかが誤り」という意味ではなく、Evidenceの二重計上、異なる機序、設定矛盾を人が
確認するための通知である。

## 8. 出力で確認する場所

### `results.json`

全criterionのstatus、strength、summary、使用Evidence、missing input、review point、conflict flag、
provenanceを保持する。PVS1は追加でevaluation contextとdecision traceを持つ。

### `evidence-lines.json`

監査用envelopeである。

- `criterion_assessments`: MET以外を含む全criterionの判断
- `evidence_lines`: VA-Specへ変換できるMET/NOT_MET
- `referenced_evidence`: Evidence IDから実体を解決するcatalog

`MET`の病原性criterionは`supports`、良性criterionは`disputes`に変換する。`NOT_MET`は
`neutral`かつ`<criterion>_not_met`であり、反対方向のEvidenceに置き換えない。

### `va-spec-1.0.1/*.json`

1ファイル1 Evidence Lineの公式出力である。`bh26AssessmentDetails` extensionからstatus、規則版、
Evidence ID、警告、PVS1 traceを確認できる。`NOT_EVALUATED`、`NOT_APPLICABLE`、
`MANUAL_REVIEW`、`DEPRECATED`はEvidenceの方向・強度を確定できないworkflow状態なので、公式
Evidence Lineには変換せず監査envelopeへ保持する。

## 9. 情報源と評価結果を混同しないための原則

- ClinVar/ClinGenの最終classificationやACMG codeを、そのまま自己証明として加点しない。
- comparatorは独立した一次Evidence、review status、splice機序等を確認する。
- computational scoreはpredictor名だけでなく、versionと校正intervalまで一致させる。
- curated region/gene-disease assessmentは、取得annotationと別Evidenceとして管理する。
- 欠測を良性・病原性の陰性Evidenceへ変換しない。
- 自動判定でもmethodとpolicy versionがなければreview済みEvidenceとして扱わない。
- 入力資料の`CLNSIG`、`ACMG_CODES`、`NOTE`は照合・監査用のsource labelであり、
  criterionのstatusを直接決めない。

## 10. 設定値、キュレーション、取得Evidenceの分離

| 種類 | 代表ファイル | 役割 |
| --- | --- | --- |
| rule/config | `config/demo-rules.json` | AF、AN、計算予測band、PM1 hotspot、PVS1 protein-loss閾値と出典 |
| curated context | `config/curated-context.json` | BA1例外、疾患別AF閾値などのreview済み文脈 |
| prepared Evidence | `work/.../prepared/evidence.json` | VEP、gnomAD、ClinVar、dbNSFP等から取得・正規化した観測 |
| evaluation output | `work/.../evaluated/` | criterion判定、監査envelope、VA-Spec Evidence Line |

`demo-rules.json`の数値は固定・版管理された現行demo policyであり、全疾患に適用できる
万能な臨床閾値ではない。コード、policy、Evidenceのどれが結果に影響したかを
分離して再現できるよう、出力にversionとEvidence IDを残す。

## 11. 実装・テスト対応表

| 領域 | 判定実装 | 主なテスト |
| --- | --- | --- |
| 共通model・Evidence gate | `src/acmg/core/models.py`、`src/acmg/criteria/common.py`、`src/acmg/services/evidence.py` | `tests/test_criteria.py`、`tests/test_curated.py` |
| PVS1 | `src/acmg/criteria/pvs1.py` | `tests/test_pvs1.py` |
| PS1・PM5 | `src/acmg/criteria/comparator.py` | `tests/test_criteria.py`、`tests/test_curated.py` |
| PM1・PM4・BP3 | `src/acmg/criteria/regions.py` | `tests/test_pm1_hotspot.py`、`tests/test_curated.py` |
| PM2・BA1・BS1 | `src/acmg/criteria/pm2.py`、`ba1.py`、`bs1.py`、`src/acmg/services/population.py` | `tests/test_population.py`、`tests/test_curated.py` |
| PP2・BP1 | `src/acmg/criteria/mechanism.py` | `tests/test_curated.py` |
| PP3・BP4 | `src/acmg/criteria/computational.py` | `tests/test_computational.py` |
| BP7 | `src/acmg/criteria/bp7.py` | `tests/test_curated.py` |
| PP5・BP6 | `src/acmg/criteria/deprecated.py` | `tests/test_criteria.py` |
| overlap flag | `src/acmg/engine.py` | `tests/test_criteria.py`、`tests/test_pvs1.py` |
| VA-Spec変換 | `src/acmg/va_spec/mapper.py` | `tests/test_va_spec.py`、`tests/test_pvs1.py` |

## 12. 実装上の既知の制約

- PVS1のNMD、transcript/exon relevance、protein region、population LoF等の判定核は実装済みだが、
  外部DBからの自動Evidence生成は未完了である。
- PVS1 start-lossのdownstream startは現在、既知booleanであることのみを確認している。
- BA1の5%は現状コード固定である。
- PP3/BP4で取得不能predictorがある場合、別の適用済みpredictorによるMETを止めない。
- BP7のnoncoding extensionは自動判定せずmanual reviewにする。
- PM1のhotspot routeは現行demo proxy policyで、gene/disease-specific CSpecではない。
- condition未指定でもprotein/gene-levelで成立可能なcriterionがあるため、METは最終疾患分類の完成を
  意味しない。review pointと`condition_assessment`を必ず併読する。

## 13. 関連資料

- `config/demo-rules.json`: 実際に読み込む閾値とpolicy provenance
- `docs/design/PVS1-PLAN.md`: PVS1 decision treeとEvidence契約の設計詳細
- `docs/PREDICTION-EVIDENCE.md`: 計算予測の取得・version・校正
- `docs/VA-SPEC-OUTPUT.md`: 内部結果からVA-Spec 1.0.1へのマッピング
- `docs/CLINGEN-POSITIVE-REFERENCES.md`: PP5/BP6を除く14 criterionの独立正例
- `docs/MET-RESULTS.md`: 現行demoでのMET実績と根拠数値
