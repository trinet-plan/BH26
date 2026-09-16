# demo-dataを入力とするACMG優先16コード評価CLI

## 目的と不変条件

ideaの概要設計v2、入力設計、VA-Spec実装仕様を基に、GRCh38のSNV・小規模indelを評価する。
入力はdemo-dataの4症例VCF全行（noiseを含む）。原本は変更しない。
全レコードの同定・補正・保留状況を記録し、検証済みvariantに16コードの結果を返す。
PVS1はClinGen General decision treeから適用可否と4段階強度を返し、PP5/BP6はDEPRECATED。
通常対象はPS1, PM1, PM2, PM4, PM5, PP2, PP3, BA1, BS1, BP1, BP3, BP4, BP7。
最終5段階分類、ランキング、Web UI、CSpec、残り12基準、GRCh37/CNV/SVは対象外。

## 入力同定

- case1_variants_v2.vcfからcase4_variants_v2.vcfを明示的に列挙する。
- 原本ファイル・SHA256・行・症例ID・variant ID・原表現を保持する。
- 概算座標、ALT欠損、HGVSC/HGVSP不整合を監査する。
- 転写産物付きHGVS、ClinVar、GRCh38参照配列を照合して一意に同定する。
- 原表現を無条件に優先せず、不一致を解消できない場合は保留する。
- 補正は別成果物に出典と理由を記録。REF検証・最小化・indel左寄せを行う。
- 症例間でEvidenceは共有可能だが、症例との対応・文脈は維持する。

## Evidenceと評価

- VEP REST、gnomAD GraphQL、ClinVar E-utilitiesのProviderとオフラインキャッシュ。
- その他の予測・機序・機能領域・保存性・BA1例外・疾患閾値は版付きローカルEvidence。
- API障害、未検出、品質不足、真のabsenceを区別する。AF/FAF・集団別観測を保持する。
- CLNSIG、ACMG_CODES、点数、NOTEの判定記述、ground truthは評価入力に使わない。
- AlphaGenome/AlphaMissenseにSpliceAI/REVELの閾値を流用しない。
- Excelと臨床ノートの自動解析は必須経路外。疾患等は任意の構造化入力。
- criterion/service/provider/resolver/VA-Spec mapperを分離する。
- statusはMET, NOT_MET, NOT_EVALUATED, NOT_APPLICABLE, MANUAL_REVIEW, DEPRECATED。
- PM2はSupporting。非ゼロAFの普遍閾値は設定しない。BS1には疾患別閾値が必須。
- 領域・機序・比較変異の根拠不足は確認待ち/未評価。予測器の校正は版管理する。
- PVS1はversion付きreviewed Evidenceが揃った分岐だけMETとし、不足はNOT_EVALUATED、
  明示的な非適用はNOT_APPLICABLE、矛盾はMANUAL_REVIEWとする。
- 基準間の競合と共有Evidenceを出す。最終分類/加点/自動重複解消はしない。

## インターフェースと成果物

Python 3.12、srcレイアウト。各基準はevaluate(input_data, services, config)を実装する。
CLI: prepare-demo（監査・同定）、evaluate（JSON/VCF、選択基準、offline）。
成果物: audit.json、variants.json、同定履歴、results.json、evidence-lines.json、summary.tsv、run-manifest.json。
VA-Spec 1.0.1 schemaに合わせてACMG EvidenceLineを検証する。MET/NOT_METのみ標準出力へ変換する。
独自envelopeと、1 Evidence Line/1 JSONのschema準拠成果物を分離する。
公式例に合わせ、構造化できる集団頻度EvidenceはCohortAlleleFrequencyStudyResultとして
EvidenceLineへ埋め込み、DataSet・StudyGroup・Method・Documentの由来を保持する。
未同定variantからEvidenceLineを作らない。暫定・非推奨・未評価は内部出力へ保持する。

## 工程と検証

1. 入力監査・同定、元データ保全と全件追跡。
2. 共通モデル、ローカルEvidence、PM2からVA-Spec出力への縦断経路。
3. 実データProvider、キャッシュ、BA1/BS1、オフライン再現（gnomAD/ClinVarまで完了）。
4. 予測/領域、比較/機序、PVS1 General decision tree、非推奨コード、競合表示。
5. 全症例デモ、README、依存ロック、pytest/Ruff、Windows/Linux CI。

各基準の成立/不成立/欠損/対象外/確認待ち、境界値、API障害、同定矛盾、
多ALT、indel正規化、REF不一致、症例間共有、schema検証、オフライン再現を検証する。
デモで不足する分岐は合成fixtureを使う。正解ラベル変更で評価が変わらないことも検証する。
ground truthとの相違はルール版・強度・対象外基準・不足Evidenceに分けて説明する。

## 権限とGit

現在のworkspace権限で進める。追加承認を要する操作は基本的に避ける。
区切りごとに検証しローカルコミットする。remote設定後は区切りごとにpushする。
権限/依存/通信で未検証の部分は完了扱いにしない。
