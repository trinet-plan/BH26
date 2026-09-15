# 2026 Biohackathon — ACMG Evidence CLI

demo-dataの4症例を対象とする、Evidence単位のACMG基準評価ツール。実装途中です。
正式な対象範囲と完了条件は [計画](docs/PLAN.md)、現状と残作業は
[進捗](docs/PROGRESS.md)、JSONのフィールド対応と制約は
[VA-Spec 1.0.1 JSON出力仕様](docs/VA-SPEC-OUTPUT.md)、蛋白注釈と予測値は
[予測Evidence仕様](docs/PREDICTION-EVIDENCE.md)、PM1の2ルート設計と
hotspot policyは [PM1設計](docs/PM1-PLAN.md) を参照してください。

## 現在動く機能

- 4症例・28レコードのVCFと注釈Excelを行単位で照合する入力監査（元ファイルは変更しません）。
- REF検証、indel左寄せ、同定Evidenceの照合、補正/保留の記録。
- 全16コードの独立評価モジュール、一括評価、競合候補の表示。
- Ensembl HGVS/GRCh38参照検証の内容ハッシュ付きキャッシュとオフライン再生。
- version付きRefSeqに一致するVEP蛋白注釈と、AlphaMissense・SpliceAI・保存性raw Evidence。
- ClinVar exact protein comparator検索と、疾患評価を分離したPS1 protein-level判定。
- 同一残基を対象とするPM5 residue検索と、疾患評価を分離したPM5 protein-level判定。
- PM1のhotspot/critical domain 2ルート評価と、ClinVar missense densityによるhotspot proxy。
- gnomAD 4.1.1集団頻度とClinVar VCVの内容ハッシュ付きキャッシュ、オフライン再生。
- 準備済みJSONとローカルEvidenceによる内部評価JSON・TSV・manifest出力。
- GA4GH VA-Spec 1.0.1 ACMG Evidence Line互換JSONと、監査用envelope出力。

実デモ全件の同定とEnsembl転写産物注釈の取得は完了しています。
疾患別閾値・キュレーション・校正済み予測Evidenceの取得は未完了です。
監査結果のPENDINGを、変異の同定完了や基準の不成立と解釈しないでください。

## 開発環境

Python 3.12以上を使用します。通常の環境では以下でセットアップできます。

```powershell
py -3.12 -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements.lock
.venv/Scripts/python.exe -m pip install --no-deps -e .
```

現セッションの `.venv` には固定済み依存を導入済みです。一時ディレクトリの権限制約により
editable installだけ失敗する環境では、次のように `PYTHONPATH` を指定して実行できます。

```powershell
$env:PYTHONPATH = 'src'
.venv/Scripts/python.exe -m unittest discover -s tests -v
.venv/Scripts/python.exe -m acmg audit-demo --input-dir demo-data --output-dir work/demo
```

`work/demo/audit.json` は全行の原表現・識別情報・問題点を保持します。
Excel由来のAlphaMissense/AlphaGenome値も出典行とともに保持します。
`CLNSIG`、`ACMG_CODES`、NOTE等を含む原行・ラベルは監査用であり評価入力には渡しません。

## 同定Evidenceが用意されている場合

```powershell
.venv/Scripts/python.exe -m acmg prepare-demo --input-dir demo-data --reference GRCh38.fa --identity-evidence identity.json --output-dir work/demo
```

参照配列は未圧縮GRCh38 FASTAと対応する `.fai` が必要です。
同定Evidenceは監査結果のrecord_idをキーにした候補配列です。候補には `variant`、
`matched_identifiers`（TRANSCRIPT/HGVSCまたはCLNVARIATIONID）、`source`、
`source_version`、`retrieved_at` が必要です。HGVSを持つ入力にはHGVSの確認が必須です。
手入力した未検証座標で同定完了を代用しません。

Ensemblへの外部送信が許可された環境では、HGVSと参照配列を取得し、応答を固定できます。

```powershell
$env:PYTHONPATH = 'src'
.venv/Scripts/python.exe -m acmg prepare-demo-online --input-dir demo-data --cache-dir tests/fixtures/ensembl-cache --output-dir work/demo-live
```

このコマンドはTRANSCRIPT/HGVSCとindel正規化に必要なゲノム座標をEnsemblへ送信します。
臨床データの取り扱い方針を確認してから
使用してください。初回成功後は同じ引数に `--offline --ensembl-release <release>` を加えると、
保存済み応答だけで再生できます。

gnomADとClinVarも含める場合は次のように実行します。送信するのは正規化済みGRCh38座位、
VCFに記載されたVCV accession、およびPM1 hotspot検索時の遺伝子記号だけです。

```powershell
$env:PYTHONPATH = 'src'
.venv/Scripts/python.exe -m acmg prepare-demo-online --input-dir demo-data --cache-dir tests/fixtures/ensembl-cache --evidence-cache-dir tests/fixtures/external-cache --output-dir work/demo-external --ensembl-release 116 --with-gnomad --gnomad-release 4.1.1 --with-clinvar --clinvar-release 2026-09-15 --with-pm1-hotspot --rules config/demo-rules.json --offline
.venv/Scripts/python.exe -m acmg evaluate --input work/demo-external/variants.json --evidence work/demo-external/evidence.json --config config/demo-rules.json --criteria all --offline --output-dir work/demo-evaluated
```

`--with-pm1-hotspot` は各missense変異の残基±N aaにあるClinVar missense変異を数え、PM1の
hotspot Evidenceを生成します。閾値は `config/demo-rules.json` の `PM1.hotspot` にversion付きで
置き、コードには持ちません。これはhotspotの近似指標であり、critical functional domainの
根拠には使いません（domainルートは人手のreviewed Evidenceのみ）。固定キャッシュでは10件の
hotspot Evidenceが生成され、PM1はMET 2件・NOT_MET 0件・NOT_EVALUATED 26件です。
評価対象の変異自身は密度に数えません。P/LPの報告が
閾値未満の場合はNOT_METにせず、報告不足として未評価に区別します。

固定キャッシュではgnomAD 17座位中15座位に観測があり、2座位は未登録です。未登録をAF=0とは
扱いません。ClinVarは13 VCVを取得し、GRCh38座標一致を確認します。ClinVarの集約分類は同定・
比較候補の監査情報に限定し、PP5/BP6の判定には使用しません。
PS1は同一アミノ酸変化のexact検索、PM5は同一残基のresidue検索を使います。exact検索では
「同一残基の別のアミノ酸変化が報告されていない」ことを示せないため、PM5はresidue検索の
完了記録だけを不成立の根拠として受け付けます。residue検索の候補は1件ずつ蛋白参照と残基を
照合してから採用します。

準備後の `variants.json` には同定済み入力のみが入り、保留レコードは `audit.json` に残ります。
`evidence.json` にはラベルから分離したEnsembl注釈が入ります。現在の固定キャッシュでは全28件が
解決し、24件は元座標を確認、4件はHGVS/GRCh38に基づいて補正されます。

## 評価の動作確認（合成データ）

```powershell
.venv/Scripts/python.exe -m acmg evaluate --input tests/fixtures/synthetic-prepared.json --evidence tests/fixtures/synthetic-evidence.json --config tests/fixtures/synthetic-rules.json --criteria all --offline --output-dir work/example-run
```

出力先は新しいディレクトリを指定してください。既存実行結果は上書きしません。
16コードを評価し、results.json・summary.tsv・evidence-lines.json・run-manifest.jsonを生成します。
`evidence-lines.json` は監査情報を含む独自envelopeです。各VA-Specオブジェクトは
`va-spec-1.0.1/*.json` に1 Evidence Line/1ファイルで出力します。
gnomAD集団頻度は公式例に合わせ、各Evidence Lineの `hasEvidenceItems` に
`CohortAlleleFrequencyStudyResult` としてDataSet・StudyGroup・Methodを含めて埋め込みます。
VA-Specを意図的に省略する場合だけ `--internal-only` を指定します。
manifestにはVA-Spec 1.0.1のcanonical schema IDを記録し、固定した自己完結型の制限schemaと
ACMG cross-field検査で書込み前に検証します。構造正規化に使った`ga4gh.va-spec`版も記録します。
1.0.1のmachine-readable schemaに従い、not-met方向は`neutral`、`methodType`はACMG criterion、
`MappableConcept`にはGKS-Core 1.0.0に存在しない`type`を出力しません。
上のfixtureの閾値は合成テスト用であり、実変異評価向けの推奨設定ではありません。
`config/demo-rules.json` のminimum ANはプロジェクト用QC設定であり、疾患別に校正された臨床閾値
ではありません。各基準のNOT_EVALUATED、NOT_APPLICABLE、MANUAL_REVIEW、DEPRECATEDも
`results.json` と `summary.tsv` に記録されます。

## Gitとデータ

工程単位でローカルコミットします。remoteが設定されている場合は区切りのよい時点でpushします。
依存・実行出力はGit対象外とし、再現テスト用キャッシュはGit管理、demo-data原本は保全します。
