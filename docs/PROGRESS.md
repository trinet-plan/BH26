# 実装進捗

更新: 2026-09-15。初期目標の完了条件を満たした。以下の拡張残件は別工程。

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
- pytest 71テスト成功、Ruff成功。3 Provider固定キャッシュによる全28件E2Eを含む。
- work/synthetic-run-1/2で内部CLI実行成功、入力エラー0。結果は合成variantでありdemo同定結果ではない。

## 未完了（次工程）

1. BP7の予測/保存性・PVS1 NMD等の計算の自動化拡張。現在はreviewed assessmentを入力する経路。
2. evaluateへのVCF直接入力・オンライン取得、症例文脈補足、準備済み入力の厳密なschema検証。
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
- 評価モデル/Provider契約にschema検証を追加し、未信頼な型・範囲値を早期に拒否すること。
- ground truthのACMGコードや最終分類は、評価結果と混ぜない。
