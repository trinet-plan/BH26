# 2026 Biohackathon — ACMG Evidence CLI

demo-dataの4症例を対象とする、Evidence単位のACMG基準評価ツール。実装途中です。
正式な対象範囲と完了条件は [計画](docs/PLAN.md)、現状と残作業は
[進捗](docs/PROGRESS.md) を参照してください。

## 現在動く機能

- 4症例・28レコードのVCFと注釈Excelを行単位で照合する入力監査（元ファイルは変更しません）。
- REF検証、indel左寄せ、同定Evidenceの照合、補正/保留の記録。
- 全16コードの独立評価モジュール、一括評価、競合候補の表示。
- Ensembl HGVS/GRCh38参照検証の内容ハッシュ付きキャッシュとオフライン再生。
- 準備済みJSONとローカルEvidenceによる内部評価JSON・TSV・manifest出力。
- GA4GH VA-Spec Pythonモデルで検証したEvidence Line JSON出力。

実デモ全件のオンライン同定と外部Evidence取得はまだ未完了です。
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

このコマンドはTRANSCRIPT/HGVSCをEnsemblへ送信します。臨床データの取り扱い方針を確認してから
使用してください。初回成功後は同じ引数に `--offline --ensembl-release <release>` を加えると、
保存済み応答だけで再生できます。

準備後の `variants.json` には同定済み入力のみが入り、保留レコードは `audit.json` に残ります。

## 評価の動作確認（合成データ）

```powershell
.venv/Scripts/python.exe -m acmg evaluate --input tests/fixtures/synthetic-prepared.json --evidence tests/fixtures/synthetic-evidence.json --config tests/fixtures/synthetic-rules.json --criteria all --offline --output-dir work/example-run
```

出力先は新しいディレクトリを指定してください。既存実行結果は上書きしません。
16コードを評価し、results.json・summary.tsv・evidence-lines.json・run-manifest.jsonを生成します。
VA-Specを意図的に省略する場合だけ `--internal-only` を指定します。
出力には検証に使った `ga4gh.va-spec` バージョンとモデルschema IDを記録します。
上のfixtureの閾値は合成テスト用であり、実変異評価向けの推奨設定ではありません。

## Gitとデータ

外部pushはしません。工程単位のローカルコミットを試行していますが、現権限では
`.git/index.lock` の作成が拒否されるため、まだ初回コミット以降の変更は未コミットです。
依存・出力・キャッシュはGit対象外とし、demo-data原本を保全します。
