# 実装進捗

更新: 2026-09-15。目標全体は未完了。

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
- Ensembl HGVS同定・GRCh38参照配列検証・キャッシュmanifestを実装。単一HGVSの実応答で形式を確認。
  全件照会はdemo-data識別子の外部送信に対する明示許可がないため安全審査で停止。
- evaluate: 準備済みJSONから16コードのresults.json、summary.tsv、run-manifest.jsonと
  VA-Spec Evidence Lineを出力。ga4gh.va-spec 0.5.0a4のモデルで書込み前に検証。
- pytest 57テスト成功、Ruff成功。入力統合、ルール分岐、実デモ監査、VA-Spec、再現性を検証。
- work/synthetic-run-1/2で内部CLI実行成功、入力エラー0。結果は合成variantでありdemo同定結果ではない。

## 未完了（次工程）

1. 外部送信の明示許可後、実装済みEnsembl経路で28件を同定し、キャッシュを固定する。
2. gnomAD/ClinVar adapter、同定からEvidence取得へのオンライン接続。
3. BP7の予測/保存性・PVS1 NMD等の計算の自動化拡張。現在はreviewed assessmentを入力する経路。
4. evaluateへのVCF直接入力・オンライン取得、症例文脈補足、準備済み入力の厳密なschema検証。
5. run-manifestに全Provider・実装ルールのハッシュを統合、ground truth比較レポート。
6. Windows/Linux CI、実デモ全経路の再現テスト。
7. 工程単位のローカルコミット。

## 確認済み環境制約

- git add/commitは `.git/index.lock: Permission denied` で失敗。変更は作業ツリーにある。
- Python 3.12環境と依存は導入済み。editable installはTemporaryDirectoryへの書込み拒否があるため、
  現環境ではPYTHONPATH=srcでCLIを実行する。
- 外部pushは行っていない。
- 上記は目標の達成扱いにはしない。環境に依存しない実装は継続可能。

## 注意する設計残件

- prepare-demoに参照FASTAのassembly/出典/ハッシュのmanifest検証を追加すること。
- identity候補はローカル契約とEnsemblライブProviderに対応。全件取得・固定は外部送信許可待ち。
- 品質・校正設定は現在テスト内の合成プロファイルのみ。本番の出典付き設定は未導入。
- 評価モデル/Provider契約にschema検証を追加し、未信頼な型・範囲値を早期に拒否すること。
- ground truthのACMGコードや最終分類は、評価結果と混ぜない。
