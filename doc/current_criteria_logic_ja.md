# ACMG/AMP 28基準・現状ロジック一覧

本ドキュメントは、このパイプラインが**現時点で**各ACMG/AMP 2015基準をどう評価しているかを、コードの実装に基づいてまとめたものである。変更履歴や過去のバージョンとの比較は含まない。

対象は28基準すべて:自動化17基準(PVS1, PS1, PM1, PM2, PM4, PM5, PP2, PP3, PP5, BA1, BS1, BS2, BP1, BP3, BP4, BP6, BP7)、文献LLM判定4基準(PS3, BS3, PS4, BP5)、表現型・家系分離3基準(PP1, BS4, PP4)、デノボ2基準(PS2, PM6)、未実装2基準(PM3, BP2)、および最終分類ロジック(`classify()`)。

---

## 1. 自動化コード

### PVS1

対象は予測LoF(loss-of-function)バリアントのみ:ストップ獲得・フレームシフト・正準スプライスドナー/アクセプター・開始コドン損失、またはRNAアッセイでLoF効果が確認されたスプライス変異。疾患条件(condition)が未指定の場合はPVS1自体を評価せず、代わりに「疾患条件さえ分かればPVS1がどう判定するか」という予備評価(preliminary_assessment)を内部的に計算してキュレーター向けに保持する(これは最終的なMET/strengthには反映されない)。

- **必須ゲート**: 遺伝子・疾患・遺伝形式の組み合わせについて「LoFが確立した疾患機序である」ことをキュレーション済みレコード(ClinGen Dosage/Gene2Phenotype/ClinGen CSpec由来の適用性データ)から確認する。CSpecの適用性スナップショット(`config/cspec_applicability.json`)はGene2Phenotypeと同じ`gene_disease`カテゴリで供給されるため、PVS1自身のコード側に専用分岐は不要。完全一致がなくても、MONDO上の親子疾患関係で関連レコードが全て機序確立に同意していれば適用(キュレーターへの確認フラグ付き)。遺伝形式が未記載でキュレーション側が単一モードのみを持ち機序確立に同意していれば、そのモードを仮定して適用(こちらも確認フラグ付き)。遺伝子不一致・疾患除外・別疾患のみキュレーション済み・モード矛盾などは全てUNKNOWNで停止。Gene2Phenotypeの`mechanism_support`が`"evidence"`でない(=`"inferred"`など、バリアント種別からの循環推論)場合はその機序を「未解決」に格下げし、以前は黙って捨てられていた"undetermined"レコードも必ずレコードとして出力するため、弱い推論のみの機序主張が対抗する新しいレコードなしに一人歩きすることを防ぐ(2026-09-25修正)。
- **バリアント種別ごとの判定木**:
  - 切断型(stop_gained/frameshift): 転写産物の生物学的関連性→NMD予測。NMDが予測されるなら、影響エクソンの関連性も確認して**MET(very_strong)**。NMDが予測されないなら、重要機能領域の破壊で**MET(strong)**、それもなければLoFバリアントの集団頻度(高頻度ならUNKNOWN=LoFが疾患原因という前提と矛盾)・領域の生物学的関連性・蛋白質喪失割合(閾値超えで**MET(strong)**、以下で**MET(moderate)**)を順に評価。「重要機能領域」はUniProtの`Domain`/`Active site`/`Coiled coil`/非disordered`Region`の重なりのみを指し、`Binding site`(単一残基の接触点)と`Motif`(短い輸送シグナル等)は2026-09-25以降対象外(MYOC/PTENの実偽陽性5件で、切断点近傍にほぼ必ず存在するこれらの小さな特徴量が誤検出の原因だったため)。
  - 正準スプライス/RNA確認済みスプライスLoF: 選択的スプライシングによるレスキューの有無、リーディングフレーム破壊の有無を確認。破壊されていれば切断型と同じ経路へ、されていなければ非NMD切断型と同じ領域評価経路へ。
  - 開始コドン損失: 別の翻訳開始点を持つ健全な代替転写産物があればUNKNOWN(レスキューされる)。なければ下流インフレーム開始点の有無と上流の病的証拠有無を見て**MET(moderate/supporting)**。
- **NOT_MET**: このロジック上、明示的なNOT_METは発生しない(適用外・情報欠損・要確認は全てUNKNOWN)。
- **UNKNOWN主要ケース**: 条件未指定、対象外のconsequence、機序未確立/未確認、転写産物・NMD・領域・スプライス関連情報の欠損またはコンフリクト、集団中にLoFバリアントが高頻度に存在。
- **設計上の特徴**: 「適用外(NOT_APPLICABLE)」「未評価(NOT_EVALUATED)」「要人手確認(MANUAL_REVIEW)」を内部的に区別して記録する(外部にはすべてUNKNOWNとして見えるが、curator向けメタデータで区別可能)。UniProt由来の領域予測やスプライスのデフォルトポリシー適用など、自動推定に基づく判定には常にキュレーター確認を促すフラグが付く。

### PS1

`comparator.py`の共通比較ロジックに委譲。ミスセンスバリアントのみ対象(スプライス系consequenceは「等価性要確認」としてUNKNOWN、それ以外は適用外)。同一残基・同一アミノ酸変化・独立バリアント・同一トランスクリプトで、ClinVar上「Pathogenic」の候補を検索する。

- **MET**: 候補が「自動判定条件」(完全一致・独立性・レビューステータス適格・スプライス機序チェック済み)または「人手レビュー済み条件」(疾患一致・全項目確認済み・キュレーター記載あり)のいずれかを満たす。strengthはstrong。
- **NOT_MET**: 完全な比較検索が完了しているのに適格な病的候補が見つからない場合。
- **UNKNOWN**: 注釈欠損、非ミスセンス、蛋白質コンテキスト欠損、候補はあるが独立性・機序確認が未完了、検索自体が未完了。
- **特徴**: 循環参照(PS1判定にPS1由来のClinVar分類を使ってしまう)を避けるため、独立エビデンス・機序一致・スプライス機序チェックを必須化している。

### PM1

`regions.py`の`evaluate_region`に委譲。「変異ホットスポット」と「重要機能ドメイン」の2ルートの択一。RYR1/MECP2/BMPR2など一部遺伝子については、VCEP仕様書から手動転記したコドン範囲(`gene_critical_domains`)に基づく専用の重要ドメイン判定が、通常のregionカテゴリより先に評価される。

- **MET**: ホットスポットルートは設定済みポリシー(窓幅・病的/良性カウント閾値)に対し良性変異が閾値以下かつ病的変異が閾値以上。重要ドメインルートは`critical_functional_region`と`benign_depletion`が両方True。strengthはmoderate。
- **NOT_MET**: ホットスポットルートで良性変異が閾値超過、または重要ドメインルートでいずれかの項目がFalse。
- **UNKNOWN**: 蛋白座標の欠損・不整合、region未取得、ホットスポットポリシー未設定/バージョン不一致、カウント欠損、密度不足かつ良性変異もない(証拠不足)。
- **特徴**: VCEP由来の重要ドメイン範囲データでは`benign_depletion`をデフォルトTrueとし、その旨を常にレビュー項目として明示する。

### PM2

`common.population_context()`経由で人口頻度観測を取得し、稀少性を判定。

- **MET**: 有効観測の最大AFが設定閾値(デモ値5e-05)以下、または全ソースが「観測ゼロ」応答なら弱いMET(supporting)。
- **NOT_MET**: 有効観測の最大AFが閾値超過。
- **UNKNOWN**: 閾値・ポリシー未設定、一部プロバイダに実エラーがあり探索不完全。
- **特徴**: 閾値のデフォルト値0と「未設定」を明確に区別する。BS1のデフォルト閾値より小さく設定されている。

### PM4

`regions.py`の`evaluate_region`に委譲。インフレーム挿入・欠失・stop_lostのみ対象。

- **MET**: 機能レビュー完了かつ蛋白長変化が0でなく、非機能的リピート領域でもない場合。strengthはmoderate。
- **NOT_MET**: 蛋白長変化が0、または非機能的リピート領域内。
- **UNKNOWN**: 対象外のconsequence、region未取得/複数、蛋白座標不整合、機能レビュー未完了。
- **特徴**: PM1と同じregionカテゴリ・座標検証ロジックを共有。

### PM5

`comparator.py`の共通比較ロジックをPS1と共有(同一関数・同一候補プール)。同一残基だが異なるアミノ酸変化の候補のみを対象とする点がPS1と異なる。

- **MET**: 同一残基・異なるアミノ酸変化・独立バリアントの病的候補が自動判定または人手レビュー済み条件を満たす。strengthはmoderate。
- **NOT_MET**: 残基スコープの比較検索が完了しているのに適格な候補が見つからない。
- **UNKNOWN**: PS1と概ね同様(注釈欠損、非ミスセンス、候補確認不足、検索未完了)。
- **特徴**: 「別部位変異でも同一機序」の妥当性を保証するため、機序一致・スプライス機序不一致でないことをPS1より重視する。

### PP2

`mechanism.py`の共通ミスセンス機序ロジックをBP1と共有。ミスセンスバリアントのみ対象。

- **MET**: 遺伝子-疾患のキュレーション済みレコードで、`missense_mechanism_established`・`spectrum_review_complete`・`low_benign_missense_variation`が全てTrue。strengthはsupporting。統計的示唆(ClinGen Gene-Disease Validity + gnomAD missense制約Z-score)による`CANDIDATE`判定でも弱いMETとして採用可(フラグ付き)。この統計示唆は2026-09-25(policy v3)以降、misZが閾値を超えるだけでなく、ClinVarキュレーション済みのミスセンス系統(benign missense fraction)も低い(デフォルト10%以下)ことを両方満たす必要がある(BMPR2の実偽陽性: misZ=3.15は閾値をわずかに超えるが、実際のClinVar病的/良性ミスセンス比では36.8%が良性で「低頻度」とは言えなかった)。
- **NOT_MET**: `missense_mechanism_established`または`low_benign_missense_variation`がFalse、または統計示唆が`NOT_SUGGESTED`。
- **UNKNOWN**: 非ミスセンス、レコード未取得かつ示唆も不十分/否定的、専門家パネルが明示的に適用外と判定、ClinGen CSpecの適用性スナップショットが当該遺伝子でPP2を「Not applicable」と明記(2026-09-25追加。実VCEP仕様が明示的にPP2を使わないと決めている遺伝子、例: SCN2A/SCN1Aのてんかんナトリウムチャネルパネル)。
- **特徴**: 専門家パネルの「適用外」判定は機序フィールドなしでもレコードとして採用される。統計的示唆はINSUFFICIENT/CURATED_NEGATIVEのときは何も生成せず、不確実なシグナルから結果を捏造しない。CSpecの「Not applicable」ゲートは統計的フォールバックより先に評価され、実VCEPの明示的判断を優先する。

### PP3

`computational.py`の較正済み予測ツール(REVEL、SpliceAI等)スコア区間判定に委譲。

- **MET**: いずれかの機構(蛋白質影響/スプライシング影響)の予測スコアが較正区間に入る。複数該当時は最強のstrengthを採用。
- **NOT_MET**: 該当する較正はあるがどのスコアも区間に入らない。
- **UNKNOWN**: 注釈欠損/複数、対象consequenceに合う較正なし、スコア取得不可・競合・範囲外。
- **特徴**: 投票数やクロス閾値ではなく、論文が定めたスコア区間のみで判定する明示的キャリブレーション方式。

### PP5

常に固定でUNKNOWNを返す。ClinGen方針(外部データベースの主張だけでは判定しない)に基づくプレースホルダー実装で、分岐ロジックは存在しない。

### BA1

`threshold_policy()`が遺伝子特異的VCEP閾値(`ba1_threshold_override`、全フィールド揃っている場合のみ採用)またはデフォルト閾値(AF>0.05)を決定する。

- **MET**: 有効観測が閾値を超え、かつ既知の除外例外リストに該当しない。strengthはstand_alone。
- **NOT_MET**: 全観測が閾値以下、または閾値超えだが既知の除外対象。
- **UNKNOWN**: 閾値/統計量/比較演算子未設定、`population_context`自体が未解決、除外リスト評価が欠落。
- **特徴**: 遺伝子特異的閾値は全フィールドが揃っている場合のみ全面的にデフォルトを置き換える(部分上書き不可)。`faf95`使用時はgnomAD公表値ではなくWilson score法による独自計算である旨を必ず記録する。

### BS1

疾患特異的閾値(`disease_frequency_threshold`、per-variantまたは遺伝子ワイド)を最優先し、条件・provenance・遺伝形式が揃わなければデフォルト閾値(デモ値0.0001)にフォールバックする。

- **MET**: 有効観測が閾値を超える。採用閾値がDRAFT(未承認)扱いの場合でもMETは返すが、承認待ちである旨を明記する。
- **NOT_MET**: 全観測が閾値以下。
- **UNKNOWN**: 疾患特異的・デフォルトいずれの閾値も未設定/不正、`population_context`未解決。
- **特徴**: 「閾値超えの証拠があるか」と「その閾値自体が承認済みか」を分離して扱う。閾値の由来根拠(prevalence, inheritance, penetrance等)をすべて記録し、欠落があれば名指しする。

### BS2(2026-09-25実装)

長らく「患者本人の臨床記録が必要で文献パイプラインの対象外」として未実装だったが、gnomADの遺伝子型カウント(ホモ接合体数/ヘミ接合体数)を「健常成人での観察」の代理指標として使う、実務でよく使われる別解釈で自動化。**常染色体劣性(AR)・半優性(semidominant)・X連鎖のみが対象**で、常染色体優性(AD)・X連鎖優性・ミトコンドリアは意図的にスコープ外(ヘテロ接合体1件の観察はBA1/BS1の頻度シグナルと本質的に同じもので、二重計上になるため)。

- **遺伝形式の決定**: `input_data["inheritance"]`(臨床ノート由来、ほとんどのERepo由来レコードでは未設定)を最優先。未設定ならGene2Phenotype("gene_disease")とClinGen Gene-Disease Validity("gene_disease_validity")のキュレーション済みレコードから、**全ソースが一致した場合のみ**自動解決を試みる(PVS1の機序候補提示と同じ「一致した場合のみ採用、不一致なら諦める」慎重さ。MONDO条件IDの完全一致は要求しない — 同一疾患が粒度違いの複数MONDO IDで登録されているケースが実際に多いため)。
- **チェックするフィールド**: AR/semidominantは`homozygote_count`、X連鎖は`hemizygote_count`(対象バリアントの染色体がXでない場合はUNKNOWN)。semidominant(ClinGen Gene-Disease Validityの実MOI表記"SD"、例: LDLR/家族性高コレステロール血症)はARと同じ扱い — ヘテロ接合体は軽症・遅発、ホモ接合体が重症・完全浸透・早発という定義そのものが、BS2の問う「完全浸透・早期発症」に対応するホモ接合体側のシグナルであるため。
- **閾値**: BA1と同じ「遺伝子特異的override→configured default」の2段構成(`bs2_threshold_override`優先)。現時点で遺伝子特異的な閾値は一件もキュレーションされておらず、常にconfigured default(暫定値: ホモ/ヘミ接合体1件でもMET)を使う。
- **MET**: 遺伝子型カウントが閾値を超える。strengthはstrong。採用閾値がDRAFT(未承認)扱いの場合(=現状は常に)、METは返すが承認待ちの予測である旨を明記する(BS1と同じ設計)。
- **NOT_MET**: 全観測が閾値以下。
- **UNKNOWN**: 遺伝形式がAR/semidominant/X連鎖のいずれでもない・未解決、X連鎖なのに染色体がXでない、閾値未設定、`population_context`未解決、遺伝子型カウントフィールド自体が取得できない。
- **特徴**: gnomADプロバイダは2026-09-25以降`homozygote_count`/`hemizygote_count`を取得する(以前はAC/AN/AFのみ)。VA-Specの`methodType`は"Population Data Assessment"ではなく"Case-Control Enrichment Assessment"(GA4GHの公式リファレンスモデルの分類に準拠)。

### BP1

`mechanism.py`の共通ロジックをPP2と共有。ミスセンスバリアントのみ対象。

- **MET**: `spectrum_review_complete=True`かつ`predominantly_truncating=True`かつ`missense_mechanism_established=False`。strengthはsupporting。
- **NOT_MET**: 上記フィールドが揃うが条件を満たさない。
- **UNKNOWN**: 非ミスセンス、レコード未取得、専門家パネルが明示的に適用外と判定、ClinGen CSpecの適用性スナップショットがBP1を「Not applicable」と明記(PP2と同じ2026-09-25の修正)、レビュー未完了。

### BP3

`regions.py`の共通region評価に委譲。インフレーム挿入・欠失のみ対象。

- **MET**: `repetitive=True`かつ`functional_importance=False`。strengthはsupporting。
- **NOT_MET**: 条件を満たさない。
- **UNKNOWN**: 対象外のconsequence、region未取得、蛋白座標不整合、機能レビュー未完了。
- **特徴**: 自動プロバイダがUniProt feature(Repeat/Compositional bias)の重なりに加え、配列複雑度(Wootton-Federhen/SEG法、窓幅24aa・閾値2.2bit)でも`repetitive`を判定する。Domain/Binding site/Active site/Motif/非disordered Regionの重なりのみが`functional_importance`をTrueにする。

### BP4

`computational.py`の較正済み予測ツールスコア区間判定に委譲(PP3と対をなす)。

- **MET**: いずれかの較正でBP4側バンドがヒットし、かつどの較正もPP3側バンドをヒットしていない。
- **NOT_MET**: いずれかの較正がPP3側バンドにヒット、またはBP4側バンドが一つもヒットしない。
- **UNKNOWN**: 較正ポリシー未設定、対象predictorの予測なし、スコア無効/範囲外、`high_confidence_null_or_splice`との矛盾。
- **特徴**: 予測ツールには意図的な「棄権帯」(どちらのバンドにも入らないスコア区間)があり、棄権は反証として扱わない。

### BP6

常に固定でUNKNOWNを返す。「一次エビデンスを確認せずに外部分類だけでスコアしない」というポリシーゲートで、分岐ロジックは存在しない。

### BP7

シノニマス(同義)バリアントのみ対象。イントロン・ノンコーディング・UTR系consequenceは「要確認」としてUNKNOWン、それ以外の非シノニマスは適用外。

- **MET**: `outside_splice_critical_region=True`かつ`no_predicted_splice_impact=True`。strengthはsupporting。
- **NOT_MET**: 条件を満たさない。
- **UNKNOWN**: 対象外のconsequence、レコード未取得、必須情報欠損、RNAエビデンスとの矛盾。
- **特徴**: 判定条件はスプライス臨界領域外の位置とSpliceAIスコア閾値のみで、配列保存性は要求しない(ClinGen SVI 2023方針)。自動プロバイダはVEP注釈に同梱のSpliceAI結果を再利用し、新規HTTP取得は行わない。

---

## 2. 文献LLM判定コード(PS3・BS3・PS4・BP5)

**共通枠組み**: 対象バリアントについて、まずERepoの引用論文(`evidence_pmids`)を確認し、無ければライブPubMed検索にフォールバックしてPMID候補を得る。得られた全PMIDについて全文(取得不可ならアブストラクト)を取得し、論文ごとにLLMへ構造化JSON抽出を依頼する。各論文の結果は複数論文集約ロジックに渡され、全論文が同じ方向で一致すればその方向を採用、方向が割れれば安全側に倒して`not_clear`とする(人間確認を促すヒント付き)。

### PS3/BS3

対象論文が実際に対象バリアントを検証しているかをまず判定し(大規模飽和変異スキャン等で対象バリアント固有の値が本文にない場合は不成立扱い)、実験ごとに機能的異常/正常/中間/混合を独立に判定させる。安全弁として、自由文の記述と分類が矛盾する場合、実験0件なのに確定方向が出た場合、定量的な数値裏付けが一切ない場合は、いずれも強制的に`not_clear`へ上書きする。同一論文内で異常/正常所見が混在するのに確信的な方向が出た場合は、LLM自身のrationaleに矛盾する所見への言及・対比的な言い回し("although"/"despite"等)があるかを確認し、あれば「実際に矛盾を検討した上での結論」とみなして元の方向をcautionフラグ付きで維持、無ければ`not_clear`へ上書きする(2026-09-25修正。以前は矛盾所見があれば無条件に上書きしていたが、実データのPTEN c.112C>Tケース — 正解はPS3=MET(moderate) — がこの安全弁自体の誤った前提で不当に上書きされていたことが判明したため)。

### PS4

症例対照(ケース・コントロール)データの抽出に特化。研究デザイン(case_control/family_cohort/case_series/not_case_control)を判定させ、罹患者・非罹患者の保因者数やオッズ比・p値を抽出する。研究デザインが`not_case_control`なのにPS4と結論した場合、定量的裏付け(件数もオッズ比もp値も)が皆無なのに確定方向が出た場合は`not_clear`へ強制上書きする。単一論文がcase-controlデザインの基準を満たさなくても、2件以上の独立した症例シリーズ論文(いずれも反証となる非罹患保因者数を含まない)を積み上げて罹患保因者数の合計が3件以上に達する場合、標準集約が`not_clear`のときに限りこの積み上げ集計をフォールバックとして採用する(2026-09-25追加。実際のClinGen PS4判定は単一のcase-controlデザイン論文ではなく複数の症例シリーズの積み上げで下されることが多いため)。

### BP5

対象バリアントを持つ患者について、その表現型をより良く説明する別遺伝子の病的バリアントが同一論文内で報告されているかを抽出する(対象バリアント自身の分類が不確実/良性であること自体からの推論は明示的に禁止)。代替遺伝子が全く抽出されていないのにBP5と結論した場合、代替所見が表現型を説明すると明記されていないのにBP5と結論した場合は`not_clear`へ強制上書きする。

- **データソースの実態的な限界**: ERepoの引用論文は変異の全体分類に対する根拠であり、BP5固有の根拠(代替診断の発見)を必ずしも指していない。ERepoに引用論文がない場合のライブ検索フォールバックも、対象遺伝子名を軸にした汎用クエリのため、代替診断を報告する症例報告を見つけにくい。

---

## 3. 表現型・家系分離コード(PP1・BS4・PP4)

`pp1_bs4_pp4_engine.py`が、ClinGen 2024のBayesianポイント評価ロジック(`pp4_pp1_bs4.evaluate_locus_evidence`)を接続する形で実装されている。

- **PP4**: 対象バリアントのERepo引用論文(`evidence_pmids`)を先に判定し、なければ患者の診断名をクエリにライブ文献検索(`pp4_literature_search`)を実行、確認された診断的浸透率(diagnostic yield)が見つかった場合のみ評価する(2026-09-25、ERepo優先化)。検索クエリは`"NOT multigene panel NOT gene panel NOT NGS panel"`で絞り込んだものと従来の汎用クエリの両方を常に試し、抽出時にLLMへ「解決済みコホートの母集団が遺伝子特異的(gene_specific)か広範なマルチ遺伝子パネル(broad_panel)か」も判定させ、gene_specificな論文が見つかるまで探索を続ける(broad_panelはyield値が他の遺伝子で希釈されている可能性が高いため、最終手段のフォールバックとしてのみ採用し、その旨をCAUTIONとして開示する)。いずれも見つからなければUNKNOWN。
- **PP1**: 患者自身の臨床ノートの家族ステータス(`family.relatives`)があればそれを直接使い、なければライブ文献検索でペディグリーデータを補う。検索は対象バリアントのERepo引用論文を先に判定し、それでも見つからなければ蛋白質表記(1文字/3文字両方)を使ったクエリで補う(生のHGVS c.表記は検索クエリに一切含めない — 実在論文でも0件になることが実測で確認されているため)。
- **PP1とPP4のポイントは独立にフロアリングせず**、両方に証拠がある場合はBiesecker et al. 2024のTable 4(合算+5.0ポイント上限)に従い、大きい生ポイントを持つ側に高いstrengthを按分する決定的ルールで変換する。
- **BS4**: 別モジュール`segregation.py`の`bs4_met`ブール値のみで判定され、MET時はSTRONG固定。

`segregation.py`自体は、PS3/BS3と同じ文献LLM抽出の枠組みを家系ごとのペディグリー(罹患/非罹患×変異有無のカウント)に適用する独立実装で、全家系が完全共分離ならPP1、非共分離を示す家系があればBS4と判定させる。実カウントと結論方向が矛盾する場合(最重要の安全弁)や、家系0件なのに確定方向が出た場合は`not_clear`へ強制上書きする。

---

## 4. デノボコード(PS2・PM6)

文献検索・LLM呼び出しを行わない、`ClinicalNoteExtraction.de_novo`(両親の変異ステータス、親子鑑定確認状況)のみを使ったルールベース評価。両親とも変異陰性が確認され、患者が罹患・家族歴陰性の場合、父性・母性の両方が確認済みならPS2をMET(STRONG)としPM6はNOT_MET、確認が揃わなければ逆にPM6をMET(MODERATE)としPS2はNOT_METとする、相互排他的な判定を行う。

---

## 5. 未実装コード(PM3・BP2)

`stubs.py`が、別バリアントとのtrans/cis位相情報に依存し、文献からもいかなるプロバイダからも判定できないため、常に固定でUNKNOWN(`stub_evidence()`)を返す。NOT_MET(否定的判定)とは明確に区別される。BS2は2026-09-25にこの分類から外れ、gnomADの遺伝子型カウントを代理指標とする自動化コードとして実装された(§1参照)。

---

## 6. 最終分類ロジック(`classify()`関数)

`acmg_pipeline/classification.py`が、28コード分の`CriterionEvidence`(`code`/`status`/`strength`/`source`)を受け取り、単一の`ClassificationResult`へ統合する。Tavtigian et al. 2018のBayesian点数システムがベース。

### Bayesian点数

| Strength | 病的側 | 良性側 |
|---|---|---|
| SUPPORTING | +1 | -1 |
| MODERATE | +2 | -2 |
| STRONG | +4 | -4 |
| VERY_STRONG | +8 | -8 |
| STAND_ALONE | +8(便宜上VERY_STRONGと同じ扱い) | -8(同上) |

METとなった全コードの点数を合計したものが`score`。

### カテゴリ閾値

| スコア範囲 | カテゴリ |
|---|---|
| score ≥ 10 | Pathogenic |
| 6 ≤ score ≤ 9 | Likely Pathogenic |
| 0 ≤ score ≤ 5 | Uncertain Significance |
| -6 ≤ score ≤ -1 | Likely Benign |
| score ≤ -7 | Benign |

### BA1オーバーライド

BA1がMETであれば、他の証拠やスコアに関わらず即座にBenignを返す(`ba1_override=True`)。スコア自体は参考情報として通常通り計算され、矛盾する強い病的証拠の有無は自動チェックされない(`met`リストに全証拠が残るため人間キュレーターがレビュー可能)。

### 重複・未評価の扱い

同一コードに複数の`CriterionEvidence`が渡されると`ValueError`(呼び出し側で事前に1つへ絞り込む必要がある)。`not_evaluated_codes`は「28コードのうち評価結果が無い、またはUNKNOWNのコード」をまとめたもので、未実装スタブコードと「実装はされているが今回の入力では判断材料が無くUNKNOWNになった」コードを区別しない(区別は`registry.is_implemented()`側の責務)。

### from_aggregated_judgment

文献LLMパイプラインの`AggregatedJudgment`を`CriterionEvidence`に変換するアダプタ。判定方向が`not_clear`ならUNKNOWN、評価対象のコードと異なる方向ならNOT_MET(例:BS4を評価中にPS3方向の証拠が見つかった場合、BS4がMETと誤認されない)、一致すればMETとし、採用論文数からstrengthを決定する。
