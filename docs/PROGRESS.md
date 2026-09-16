# 実装進捗

更新: 2026-09-16。初期目標の完了条件を満たした。以下の拡張残件は別工程。
この文書は工程順の記録である。現時点の全体像・16基準の状況・残件は
[俯瞰](OVERVIEW.md) を参照すること。実デモ28件×16基準の内訳は工程ごとに変わるため、
本文中の件数はその工程時点の値である。

## 実装・検証済み

- 計画をdocs/PLAN.mdへ保存。元のidea/demo-dataは変更していない。
- 同梱Python 3.12.14から `.venv/Scripts/python.exe` を利用可能。
- audit-demo: 4症例、28 ALTレコードを保持するaudit.jsonを生成。
- 注釈Excel 28行を標準ライブラリだけで読み、case/variant_idでVCF全行と1対1照合。
  VCFとExcelの識別情報不一致は0件。分類ラベルは評価境界へ渡さない。
- ALT欠損2件、LMNAの同一HGVSに異なる座標がある不整合を検出。
- ドメイン契約、REF検証、indel左寄せ、同定候補照合、ラベルを除去した入力境界。
- prepare-demo: indexed FASTAと出典付きローカル同定候補から検証済み入力を生成。
- 16コードの独立評価モジュールと一括engine。機序/領域/比較根拠は出典付きreviewed Evidenceを要求。
- PVS1は当初、分岐履歴・推奨強度だけを返す暫定実装として開始した。
- 基準間の重複候補フラグ、ローカルEvidence/Population Service。
- HTTPキャッシュ（内容ハッシュ検証、offline、タイムアウト、限定リトライ）とVEP adapterを実装。
- Ensembl release 116で全17固有HGVSを同定し、GRCh38参照配列検証と応答キャッシュを固定。
  全28件の内訳はVERIFIED 24件、CORRECTED 4件、PENDING 0件。
- 補正はcase1-var1、case1-var2、case1-noise5、case2-var2。元表現と補正後表現をaudit.jsonに併記。
- 同じ固定キャッシュから完全オフラインで全28件を再処理し、audit.json/variants.jsonが
  オンライン実行と同一SHA-256になることを確認。
- Ensembl由来の17注釈Evidenceをground truthラベルから分離してevidence.jsonへ出力。
- VEPをversion付きRefSeq transcriptへ限定し、17固有変異すべてのprotein_id、15変異の蛋白座標、
  missense/synonymousのアミノ酸、導出可能な蛋白長差を保存。別isoformのconsequence混入を禁止。
- VEPからAlphaMissense 10件、SpliceAI 15件、保存性17件をraw computational Evidenceとして保存。
  RESTでmodel/data版を確定できないものはcalibration_eligible=falseとし、PP3/BP4へ不使用。
- ClinVar exact HGVSp検索によるPS1 comparator adapterを実装。疾患未指定でもprotein-levelを
  評価し、condition assessmentを分離。10固有missense検索は全件完了し、適格Pathogenic comparator
  は0件。MYBPC3 p.Ser236Glyには別MNVのBenign/Likely benign候補1件があり要確認として保存。
- evaluate: 準備済みJSONから16コードのresults.json、summary.tsv、run-manifest.jsonと
  VA-Spec Evidence Lineを出力。VA-Spec 1.0.1/GKS-Core 1.0.0に合わせ、not-met=`neutral`、
  methodType=criterion、MappableConceptのtypeなしで書込み前に検証。
- 実デモ28件×16基準=448結果を生成。NOT_APPLICABLE 113、MANUAL_REVIEW 11、
  NOT_EVALUATED 268、DEPRECATED 56。独立Evidence不足のためMET/NOT_METは推測していない。
- Ensembl注釈のみの実デモ評価ではMET/NOT_METがないためEvidence Lineは0件。
  合成fixtureではMET/NOT_METのEvidence Line生成とオフライン同一性を検証。
- gnomAD 4.1.1 adapterを実装。17固有座位を1バッチで取得し、15座位・293観測を固定。
  未登録2座位はNO_OBSERVATIONとして扱い、AF=0へ変換しない。
- ClinVar E-utilities adapterを実装。13 VCVを取得し、全件で正規化済みGRCh38座標と一致。
  現行VCV版と入力版を併記し、集約分類はPP5/BP6に使用しない。
- gnomAD Evidenceを使うデモ評価はPM2がNOT_MET 26件/NOT_EVALUATED 2件、BA1が
  NOT_MET 15件/NOT_EVALUATED 13件。VA-Spec Evidence Lineは41件を生成・検証。
- 監査用の独自envelopeとは別に、`va-spec-1.0.1`へ1 Evidence Line/1 JSONを出力。
- 公式1.0.1例に合わせ、gnomAD人口頻度をCohortAlleleFrequencyStudyResultとして埋め込み、
  DataSet・StudyGroup・Method・DocumentとQC provenanceを保持。詳細仕様を文書化。
- PM2（AC=0、AN=10000、AF=0）とBA1（AC=501、AN=10000、AF=0.0501、例外なし）の
  MET合成fixtureを追加し、supports/supportingとdisputes/standaloneのVA-Spec出力をE2E検証。
- case2-var1と固定gnomAD 4.1.1キャッシュを使う回帰テストを追加。テスト限定の
  max AF=0.000014ではPM2 METとなり、VA-Specのsupports/supporting出力まで検証。
  gnomAD未登録のcase2-var2はcallability/AN不明のためAF=0にせずNOT_EVALUATEDを維持。
- pytest 79テスト成功、Ruff成功。3 Provider固定キャッシュによる全28件E2Eを含む。
- work/synthetic-run-1/2で内部CLI実行成功、入力エラー0。結果は合成variantでありdemo同定結果ではない。

- PM1をACMG定義どおり2ルート（mutational hotspot / critical functional domain）で評価。
  共通条件はbenign_depletionのみ。critical domainにpathogenic enrichmentを要求しない。
  conditionなしでもprotein-levelで判定し、condition_assessmentを分離して記録。
  自動Evidenceはhotspotルート限定で、critical domainは人手のreviewed Evidenceのみ受け付ける。
- ClinVar missense density Provider（esearch+esummary、遺伝子単位検索、残基±N aa）を実装。
  閾値はconfigのversion付きpolicy（PM1-hotspot-v1: ±5 aa, P/LP>=3, B/LB=0）から取得し、
  コードに数値を持たない。数えたVCVとconditionを保存し、conflictingはどちらにも数えない。
  検索が打ち切られた場合はdensity Evidenceを出さない。P/LP報告不足はNOT_METではなく未評価。
- 実デモ28件でPM1はMET 2件、NOT_MET 6件、NOT_EVALUATED 20件。VA-Spec Evidence Lineは65件。
  MET 2件はconditionが未指定のためcondition_assessment=NOT_EVALUATEDと確認事項を併記。
- pytest 98テスト成功、Ruff成功。固定キャッシュによるPM1 hotspotのオフライン再現を含む。

- conditionを全体で任意入力とし、PM5もPS1・PM1と同じくprotein-levelで評価するよう変更。
  疾患関連性はcondition_assessmentとして分離記録する。
- PM5用のClinVar residue検索を実装。遺伝子単位のmissense検索とesummaryで同一残基の候補を
  絞り、候補ごとにefetchとEnsembl照合で蛋白参照・残基・参照アミノ酸の一致と置換の相違を
  確認する。exact検索は「他の置換がない」ことを示せないため、PM5はsearch_scope=residueの
  完了記録のみを不成立の根拠にする。
- 比較候補がVUS・良性・conflictingの場合はMANUAL_REVIEWにせず、検索済みの記録として
  NOT_METのEvidenceに残す。人手レビューで成立しうるP/LP候補だけを確認待ちにする。
- 実デモではPM5がMET 1件（MYH7 p.Arg719、比較候補p.Arg719Gln Pathogenic）、NOT_MET 17件、
  NOT_APPLICABLE 10件。PS1は上記変更によりMANUAL_REVIEW 3件がNOT_MET 2件＋splice 1件へ移動。
  VA-Spec Evidence Lineは79件。

- PP2・BP1もconditionを任意化。gene_disease Evidenceを疾患非依存でも受け付け、
  assessment_scope=gene_levelとcondition_assessmentを分離記録する。MET時に疾患未確認なら
  確認事項を併記。conditionを与えた場合に別疾患のEvidenceを流用しない点は従来どおり。
  PVS1も同様に任意化。LoF機序は遺伝子単位のキュレーションとして扱う（autoPVS1もvariantと
  genome版だけで動き、疾患入力を取らずPVS1.level等の遺伝子リストで機序を解決する）。
  これでBS1以外の全基準がcondition任意。BS1のみ疾患別閾値が疾患なしに定義できないため必須。
- 実デモではPP2・BP1の18件がcondition待ちからgene_disease Evidence待ちへ移動。

- PM1 hotspot密度のisoform numbering取り違えを解消。窓に入った候補を1件ずつefetchし、
  対象protein_id上のClinVar表記から残基を確定してから数える。対象蛋白に表記がない候補と
  確定位置が窓外の候補は数えずにunplaced_excluded / outside_window_excludedへ記録する。
  RUNX1 K110Eでは16件中、自身1件とA134の別isoform表記3件が除かれ12件になる。
  デモ10領域の件数は変化しない（対象遺伝子ではisoform表記の混入がなかった）。

- dbNSFP Provider（MyVariant.info経由）を実装。metadataからdbNSFP版4.8aを取得して各スコアの
  predictor_versionに固定し、版が確定したスコアだけcalibration_eligible=trueとする。
  REVELとAlphaMissenseのみ採用し、SIFT/PolyPhen-2は基準入力に渡さない。
  dbNSFPはnonsynonymous SNV対象のためindelには問い合わせず、404は「値なし」として扱う。
- PP3/BP4の校正区間をconfigへ追加（REVEL: Pejaver et al. 2022のClinGen SVI校正。
  PP3 0.644/0.773/0.932、BP4 0.290/0.183/0.016/0.003）。強度別の複数区間に対応。
- 実デモでPP3がMET 6件（strong 1・moderate 1・supporting 4）、NOT_MET 12件、N/A 10件。
  BP4がMET 8件（moderate 2・supporting 6）、NOT_MET 10件、N/A 10件。
  KCNJ5 p.Arg155Gln（REVEL 0.934）でPP3 strong、MYH7 p.Arg719Trp（0.814）でPP3 moderate。
- 実デモ全体はMET 17、NOT_MET 98、NOT_EVALUATED 112、NOT_APPLICABLE 161、
  MANUAL_REVIEW 4、DEPRECATED 56。VA-Spec Evidence Lineは115件。
- pytest 119テスト成功、Ruff成功。5 Provider固定キャッシュによる全28件E2Eを含む。

- PP3/BP4を複数校正の併用に対応（`selected_calibrations`）。蛋白機序（REVEL）と
  splicing機序（SpliceAI: Walker et al. 2023 ClinGen SVI、PP3>=0.2、BP4<=0.1）を併用する。
  機序は加算せず、PP3はいずれか成立で最強区間を採用。BP4は適用対象の全機序が良性側を示す
  ことを要求し、splice影響が予測される場合は成立させない（benign_blocked_byに記録）。
- Ensembl VEPが版を公開しないSpliceAIは、校正側の`version_assertion`（出典・宣言者・理由）を
  明示した場合に限り使用し、宣言をprovenanceの`version_assertions`に残す。宣言がなければ不使用。
- 実デモでPP3 MET 6件、BP4 MET 12件（うち3件はsynonymousでsplicing校正のみ成立）。
  全体はMET 21、NOT_MET 108、NOT_EVALUATED 112、NOT_APPLICABLE 147、MANUAL_REVIEW 4、
  DEPRECATED 56。VA-Spec Evidence Lineは129件。
- pytest 127テスト成功、Ruff成功。

- `evaluate --context` を追加。キュレーション済み臨床文脈（condition、inheritance、
  疾患別頻度閾値、BA1例外リスト）を取得Evidenceと分離して読み込む。Evidenceに使う
  フィールドを文脈として混入させることは拒否し、キーはvariant keyに限定する。
  読み込んだ文脈の版・件数はrun-manifestの`curated_context`に記録する。
- BA1例外リストは`complete`が真のときだけ判定に使う。未転記の空リスト（complete=false）は
  「リストに無い」と主張せず未評価を維持する。config/curated-context.jsonを口として同梱。
- BA1例外評価とBS1閾値は、IRIを持つ場合のみEvidence itemとして出力し、そうでなければ
  provenanceに記録する。BA1が例外評価を受け取るとVA-Spec出力がIRI要件で失敗する不具合を修正。
- pytest 138テスト成功、Ruff成功。

- ClinGen SVI BA1例外リスト（Ghosh et al. 2018 Table 1の9変異）を**手入力で転記**し、
  complete=trueで同梱。機械可読形式での公開もAPIも存在しないため自動取得はできない。
  転記者・転記元・機械処理した範囲・未実施の再確認をentry_method/transcriptionとして記録し、
  run-manifestにも毎回出力する。第三者による再確認は未実施であり、臨床利用前に原典と要照合。
  各変異はCAID経由でClinGen Allele RegistryのGRCh38 alleleへ解決し、遺伝子一致と転写産物HGVS
  一致を確認。ACAD9 c.-44_-41dupのみregistryのgenomic interval
  （NC_000003.12:g.128879648_128879651dup）からアンカー付き左寄せVCF表現を導出した。
  entryごとにcaid・hgvs_c・resolved_byを必須とし、誤った座標で別変異を免除しないようにする。
- 実デモではAF>5%の11件がBA1 MET（stand_alone）に確定。全体はMET 32、NOT_MET 108、
  NOT_EVALUATED 101、NOT_APPLICABLE 147、MANUAL_REVIEW 4、DEPRECATED 56。
  VA-Spec Evidence Lineは140件。pytest 139テスト成功、Ruff成功。

- ClinGen Evidence RepositoryのAPI（erepo）から専門家パネル判定556件を走査し、consequence別に
  適用ゲートを確認するfixtureを追加。PVS1（splice acceptor / frameshift）、PM4とBP3
  （in-frame deletion）、BP7（synonymous）の5変異をCAID・適用コード・最終分類つきで取り込み、
  各基準が対象consequenceで開き、それ以外では閉じることをテストで固定した。
  キュレーション済みEvidenceは与えていないため、開いたゲートはNOT_EVALUATEDで止まる。
- BP3はerepoの53遺伝子556件で適用例が0件だった。VCEPが実運用でほとんど使っていないため、
  実データでの成立例は得られていない。ゲートはin-frame deletionのfixtureで確認している。
- pytest 146テスト成功、Ruff成功。

- PVS1の暫定MANUAL_REVIEW制約を廃止し、ACMG/AMP 2015、ClinGen SVI PVS1 2018、
  ClinGen SVI Splicing 2023に基づくGeneral decision treeを実装。nonsense/frameshift、
  canonical splice、RNA実証splice LoF、start-lossからvery strong/strong/moderate/supportingを
  正式適用し、VA-Spec Evidence Lineへ変換する。
- condition-specific LoF mechanismを優先し、未提供時または該当Evidenceがない場合はgene-levelへ
  fallbackする。condition status、mechanism scope、decision trace、rule source、warning、
  unresolved requirementをresults.json schema 1.1に記録する。
- transcript、NMD、splice/RNA、protein region、population LoF、initiationの正規化Evidence契約を追加。
  外部自動取得とCNV/SVは対象外で、不足を陰性へ変換しない。
- 実デモの集計は従来どおりMET 32、NOT_MET 108、NOT_EVALUATED 101、NOT_APPLICABLE 147、
  MANUAL_REVIEW 4、DEPRECATED 56。PVS1はNOT_APPLICABLE 25、NOT_EVALUATED 3。
  pytest 159テスト成功、Ruff成功、全28件VA-Spec検証成功。

- ClinGen Evidence Repository 2.5.6のsummary APIを13,265件全走査し、PP5/BP6を除く
  実装対象14 criterionすべてにMET例があることを確認。BP3は現行23件で、以前の限定556件走査に
  基づく「0件」を更新した。
- 13変異のClinGen正例fixtureを追加。ERepoの適用コードは期待値manifestだけに隔離し、評価入力は
  GRCh38 variantと出典付き独立Evidenceだけで構成した。PVS1/PS1/PM1/PM2/PM4/PM5/PP2/PP3/
  BA1/BS1/BP1/BP3/BP4/BP7が各1件以上METになり、VA-Spec 1.0.1出力も再読込検証する。
- 参照値は公開解釈に照合した。例としてPAX6 REVEL=0.967、SLC6A8 REVEL=0.079、ITGB3
  1342/24024、GUCY2D 1/1613704・PM2閾値0.0004、MYH7 10/34232・BS1閾値0.0002を使用。
  GCK/FOXG1の旧assemblyまたは誤った
  VCF anchorを修正し、indel anchorはEnsembl GRCh38参照配列で確認した。
- 正例6テストを含む全165テスト成功、Ruff成功。既存28件の固定キャッシュE2Eは入力エラー0、
  VA-Spec VALIDATEDを維持した。

## 未完了（次工程）

1. BA1例外リストは手入力のため、原典Table 1との第三者照合が未実施。
2. SpliceAIは版がVEP経由の宣言に依存する。版を確定できる配布元からの取得が望ましい。
3. PM1のgene/disease-specific対応（強度可変、PM1不使用geneのNOT_APPLICABLE表明）と、
   critical functional domainのreviewed Evidence入力経路。BS1疾患閾値のキュレーション。
4. BP7のRefSeq exon境界position adapterと校正済み保存性policy、PVS1 NMD等のEvidence取得自動化。
   現在はSpliceAI/保存性raw Evidenceとreviewed assessment入力経路まで。
5. evaluateへのVCF直接入力・オンライン取得、症例文脈補足、準備済み入力の厳密なschema検証。
6. run-manifestに全Provider・実装ルールのハッシュを統合、ground truth比較レポート。
7. Windows/Linux CI。

## 確認済み環境制約

- 通常サンドボックスでは `.git/index.lock` が拒否されるが、明示依頼に基づく限定昇格で
  ローカルコミット `f25e38e` を作成済み。
- Python 3.12環境と依存は導入済み。editable installはTemporaryDirectoryへの書込み拒否があるため、
  現環境ではPYTHONPATH=srcでCLIを実行する。
- 外部pushは行っていない。
- 上記は目標の達成扱いにはしない。環境に依存しない実装は継続可能。

## 注意する設計残件

- prepare-demoに参照FASTAのassembly/出典/ハッシュのmanifest検証を追加すること。
- identity候補はローカル契約とEnsemblライブProviderに対応し、全件取得・固定済み。
- demo用の品質・校正設定は出典付きで導入済みだが、疾患・遺伝子別の臨床承認済みpolicyではない。
  VEP予測Evidenceの取得・版不足による採点除外条件は実装済み。
- 評価モデル/Provider契約にschema検証を追加し、未信頼な型・範囲値を早期に拒否すること。
- ground truthのACMGコードや最終分類は、評価結果と混ぜない。
