# 2026 Biohackathon — ACMG Evidence CLI

demo-dataの4症例を対象とする、Evidence単位のACMG基準評価ツール。実装途中です。
全体像（データ源・16基準の状況・判定方針・残件）は [俯瞰](docs/OVERVIEW.md)、
正式な対象範囲と完了条件は [計画](docs/PLAN.md)、工程単位の記録は
[進捗](docs/PROGRESS.md)、JSONのフィールド対応と制約は
[VA-Spec 1.0.1 JSON出力仕様](docs/VA-SPEC-OUTPUT.md)、蛋白注釈と予測値は
[予測Evidence仕様](docs/PREDICTION-EVIDENCE.md)、PM1の2ルート設計と
hotspot policyは [PM1設計](docs/design/PM1-PLAN.md)、PVS1のGeneral decision treeと
Evidence契約は [PVS1設計](docs/design/PVS1-PLAN.md) を参照してください。
全16 criterionの入力、判定順序、全status終端、strengthは
[判定ロジック・入力Evidenceガイド](docs/CRITERION-DECISION-GUIDE.md)にまとめています。
全criterionのClinGen正例対応は [ClinGen正例リファレンス](docs/CLINGEN-POSITIVE-REFERENCES.md)
にまとめています。
MET 52件の変異・強度・根拠数値は [MET一覧](docs/MET-RESULTS.md) にあります。

## 現在動く機能

- 4症例・28レコードのVCFと注釈Excelを行単位で照合する入力監査（元ファイルは変更しません）。
- REF検証、indel左寄せ、同定Evidenceの照合、補正/保留の記録。
- 全16コードの独立評価モジュール、一括評価、競合候補の表示。
- Ensembl HGVS/GRCh38参照検証の内容ハッシュ付きキャッシュとオフライン再生。
- version付きRefSeqに一致するVEP蛋白注釈と、AlphaMissense・SpliceAI・保存性raw Evidence。
- ClinGen SVI BA1例外リスト9変異による、AF>5%変異のBA1確定。
- ClinVar exact protein comparator検索と、疾患評価を分離したPS1 protein-level判定。
- 同一残基を対象とするPM5 residue検索と、疾患評価を分離したPM5 protein-level判定。
- PM1のhotspot/critical domain 2ルート評価と、ClinVar missense densityによるhotspot proxy。
- ClinGen General Guidanceに基づくPVS1 decision tree、4段階強度、condition/gene-level文脈、
  decision trace、VA-Spec出力。
- dbNSFP releaseを固定したREVEL/AlphaMissenseと、校正区間によるPP3/BP4判定。
- 蛋白機序と splicing機序の2校正を併用するPP3/BP4（加算せず、BP4は全機序の一致を要求）。
- gnomAD 4.1.1集団頻度とClinVar VCVの内容ハッシュ付きキャッシュ、オフライン再生。
- 準備済みJSONとローカルEvidenceによる内部評価JSON・TSV・manifest出力。
- GA4GH VA-Spec 1.0.1 ACMG Evidence Line互換JSONと、監査用envelope出力。
- PP5/BP6を除く14 criterionすべてについて、ClinGen ERepo由来の独立正規化Evidenceから
  最低1件のMETを再現する正例E2Eスイート。

実デモ全件の同定・注釈取得と、校正済み予測Evidence（dbNSFP REVEL）による
PP3/BP4評価は完了しています。遺伝子-疾患機序（PP2/BP1/PVS1）と疾患別頻度閾値（BS1）の
キュレーションは未完了です。現時点の16基準の状況は [俯瞰](docs/OVERVIEW.md) を参照してください。
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
.venv/Scripts/python.exe -m unittest tests.test_clingen_positive -v
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
.venv/Scripts/python.exe -m acmg prepare-demo-online --input-dir demo-data --cache-dir tests/fixtures/ensembl-cache --evidence-cache-dir tests/fixtures/external-cache --output-dir work/demo-external --ensembl-release 116 --with-gnomad --gnomad-release 4.1.1 --with-clinvar --clinvar-release 2026-09-15 --with-pm1-hotspot --rules config/demo-rules.json --with-dbnsfp --offline
.venv/Scripts/python.exe -m acmg evaluate --input work/demo-external/variants.json --evidence work/demo-external/evidence.json --config config/demo-rules.json --criteria all --offline --output-dir work/demo-evaluated
```

`--with-dbnsfp` はMyVariant.info経由でdbNSFPのREVEL/AlphaMissenseを取得します。送信するのは
正規化済みGRCh38座位だけです。metadataエンドポイントからdbNSFP版（現在4.8a）を取得して
各スコアに記録し、版が確定した場合に限り `calibration_eligible` とします。Ensembl VEP RESTは
REVEL/SpliceAIの由来版を公開しないため、VEP由来の予測値はPP3/BP4に使いません。
校正区間は `config/demo-rules.json` の `computational` に出典付きで置きます。
`selected_calibrations` に複数の校正を並べられ、蛋白機序（REVEL: Pejaver et al. 2022、
PP3 supporting/moderate/strong、BP4 supporting/moderate/strong/very strong）と
splicing機序（SpliceAI: Walker et al. 2023 ClinGen SVI Splicing Subgroup、
PP3 >=0.2 supporting、BP4 <=0.1 supporting）を併用します。
機序は加算しません。PP3はいずれかの機序が区間を満たせば成立し、最も強い区間を採ります。
BP4は適用対象の全機序が良性側を示す必要があり、splice影響が予測される場合は成立しません
（低い蛋白スコアはsplice部位の破壊について何も言わないため）。
dbNSFPはnonsynonymous SNV用のため、indelには問い合わせません。

Ensembl VEPはSpliceAIのモデル版を公開しないため、校正側で `version_assertion` として
「どのreleaseが配信したスコアか」を出典・宣言者・理由つきで明示した場合に限り使用します。
宣言は結果のprovenanceに `version_assertions` として残ります。宣言がなければそのスコアは
判定に使いません。

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
監査envelope schema 1.1は全criterionのworkflow statusを `criterion_assessments` に保持し、
MET/NOT_MET以外の未評価・対象外・要確認・非推奨も欠落させません。個別Evidence Lineには
`bh26AssessmentDetails` extensionを付け、判定status、Evidence ID、規則版、警告、PVS1 traceを
単独ファイルでも確認できます。IRI参照の実体はenvelopeの `referenced_evidence` から解決でき、
provider/version、取得日時、品質、method、curator、観測値を監査できます。
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

キュレーション済みの臨床文脈は `--context config/curated-context.json` で渡します。取得した
Evidenceには混ぜず、疾患・遺伝形式・疾患別頻度閾値・BA1例外リストだけを受け付けます。
BA1例外リストは `complete` が真のときだけ判定に使い、未転記（偽）の間は「リストに無い」と
主張せず未評価のままにします。同梱の `config/curated-context.json` には
ClinGen SVI（Ghosh et al. 2018）の9変異を収録済みです。**このリストは論文Table 1からの
手入力**であり（機械可読形式での公開もAPIも存在しないため）、第三者による再確認は
受けていません。`entry_method: manual_transcription` と転記の経緯をリストに記録し、
run-manifestの `curated_context` にも毎回出力します。臨床利用の前には原典の表と
照合してください。各変異はClinGen Allele RegistryのCAIDを
経由してGRCh38 alleleへ解決し、解決方法・遺伝子・転写産物HGVS・登録IDを各entryに記録しています
（ACAD9の重複のみ、registryのgenomic intervalからアンカー付き左寄せ表現を導出）。読み込んだ文脈の版と件数はrun-manifestに記録します。

`condition`（疾患）は任意入力です。PS1・PM5・PM1・PP2・BP1・PVS1は蛋白/遺伝子レベルで判定し、
疾患関連性は判定に混ぜず `condition_assessment`（MATCHED / NOT_EVALUATED）として分離記録し、
未確認のままMETになった場合は確認事項を併記します。`condition` を与えた場合、別疾患の
キュレーション済みEvidenceは流用されません。BS1は疾患別閾値そのものが疾患なしに定義できない
ため、引き続き`condition`が必須です。PVS1はcondition-specific LoF機序を優先し、存在しない場合に
gene-level機序へfallbackします。condition未指定時もgene-level Evidenceがあれば評価を継続し、
`evaluation_context`とwarningに評価範囲を残します。NMD・transcript relevance・RNA assay・
protein region等は出典付きの正規化Evidenceを要求し、未取得値を推測しません。

## Gitとデータ

工程単位でローカルコミットします。remoteが設定されている場合は区切りのよい時点でpushします。
依存・実行出力はGit対象外とし、再現テスト用キャッシュはGit管理、demo-data原本は保全します。
