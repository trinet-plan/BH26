# 実装進捗

更新: 2026-09-15。初期目標の完了条件を満たした。以下の拡張残件は別工程。
実デモ28件×16基準の内訳は工程ごとに変わる。最新値は末尾のPM1節を参照。

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
- PVS1は分岐履歴・推奨強度のみ。複雑な未自動化分岐はcurated branchで補足する。
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
  これでPS1・PM5・PM1・PP2・BP1・BP7がcondition任意。conditionを要求するのは
  BS1（疾患別閾値が疾患なしに定義できないため必須）とPVS1（未変更）のみ。
- 実デモではPP2・BP1の18件がcondition待ちからgene_disease Evidence待ちへ移動。

- PM1 hotspot密度のisoform numbering取り違えを解消。窓に入った候補を1件ずつefetchし、
  対象protein_id上のClinVar表記から残基を確定してから数える。対象蛋白に表記がない候補と
  確定位置が窓外の候補は数えずにunplaced_excluded / outside_window_excludedへ記録する。
  RUNX1 K110Eでは16件中、自身1件とA134の別isoform表記3件が除かれ12件になる。
  デモ10領域の件数は変化しない（対象遺伝子ではisoform表記の混入がなかった）。

## 未完了（次工程）

1. PM1のgene/disease-specific対応（強度可変、PM1不使用geneのNOT_APPLICABLE表明）と、
   critical functional domainのreviewed Evidence入力経路。conditionとBS1疾患閾値の入力経路。
2. BP7のRefSeq exon境界position adapterと校正済み保存性policy、PVS1 NMD等の計算の自動化拡張。
   現在はSpliceAI/保存性raw Evidenceとreviewed assessment入力経路まで。
3. evaluateへのVCF直接入力・オンライン取得、症例文脈補足、準備済み入力の厳密なschema検証。
3. run-manifestに全Provider・実装ルールのハッシュを統合、ground truth比較レポート。
4. Windows/Linux CI。

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
- 品質・校正設定は現在テスト内の合成プロファイルのみ。本番の出典付き設定は未導入。
  VEP予測Evidenceの取得・版不足による採点除外条件は実装済み。
- 評価モデル/Provider契約にschema検証を追加し、未信頼な型・範囲値を早期に拒否すること。
- ground truthのACMGコードや最終分類は、評価結果と混ぜない。
