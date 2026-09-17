# 依頼: BS1の疾患別頻度閾値をClinGen公開仕様から取り込み、DRAFTキューを廃止する

**このドキュメントの使い方**: そのまま別の開発者(またはその開発者が使うAIコーディングエージェント)への作業指示として渡せるように書いています。「現状」「依頼する作業」「具体的な変更箇所」「期待するスキーマ」「テストへの影響」「未決定事項」の順に必要な情報を揃えています。

**2026-09-17 追記(このドキュメントの後に入った方針変更)**: `bs1.py`は元々「未承認/既定の閾値を超過していても、レビュー未承認ならMETにせずUNKNOWN+draft候補として返す」という設計だったが、この後の方針転換で「未承認でもMET/strengthは出し、未承認である旨をcuratorHint(`review_points`)と`provenance.draft_threshold_candidate`に残す」という設計に変わった(`policy_status()`自体は残っているが、その結果はもうMET/UNKNOWNの分岐には使われず、curatorHintの文言にのみ反映される)。以下の本文はこの変更前に書かれているため、「MET(良性方向の確定判定)には絶対にならない」等の記述は現状と合わなくなっている箇所がある。該当箇所には注記を追加した。

---

## 1. 現状

BS1(集団頻度が高すぎることを示す良性方向の根拠)は [`acmg_pipeline/criteria/bs1.py`](../acmg_pipeline/criteria/bs1.py) が評価する。判定は疾患別に承認された「最大許容頻度(maximum credible allele frequency)」を必要とするが、これは今のところ**どこからも自動供給されていない**。

- `bs1.py`の`disease_specific_threshold()`(21-48行目)は、`input_data.get("disease_frequency_threshold", {})`という辞書を読むだけで、**自分でファイルを読みには行かない**。
- この辞書は`acmg_pipeline/automated_core/interface.py`の`criterion_input()`が、VCFのINFO列(または`inputs_from_prepared_record()`のprepared JSON)にある`disease_frequency_threshold`キーをそのまま転記しているだけ(96-98行目)。つまり**呼び出し側が事前に埋めておかない限り、常に空**。
- 埋まっていない場合、BS1は`config`側の汎用デフォルト閾値([`config/demo-rules.json`](../config/demo-rules.json)の`BS1.default_*`)にフォールバックする。このデフォルトは`policy_status: DRAFT`扱いになる(`bs1.py`の`policy_status()`, 120-130行目)。**(2026-09-17追記)** 元々はDRAFT扱いだとMETには絶対にならない設計だったが、方針転換によりMET+strengthは出るようになり、未承認である旨は`review_points`(curatorHint)と`provenance.draft_threshold_candidate`に残す形になった。
- [`config/bs1_thresholds_draft.json`](../config/bs1_thresholds_draft.json)には、ClinGen Criteria Specification Registry(<https://cspec.clinicalgenome.org/>)掲載の7件のVCEP仕様(LDLR/家族性高コレステロール血症、PAH/フェニルケトン尿症、GAA/ポンペ病、難聴AR/AD、PTEN、MYOC/緑内障)がDRAFTとして転記済みだが、**これを読み込んで`disease_frequency_threshold`に変換するコードが存在しない**。つまり現状、この7件はキュレーター向けの参考資料止まりで、実行時には一切効いていない。

過去の実インシデント(重要な背景): このDRAFT/APPROVEDという区別自体は思いつきではない。以前、実APIでMYH7 c.4472C>Gがレビュー前のPLACEHOLDER閾値だけで"likely benign"に倒れたことがあり、それを受けて2026-09-17に「未承認の閾値は最終判定に使わない(MET/NOT_METを出さずUNKNOWN+draft候補とする)」という安全弁が一度追加された(詳細: [`tests/test_bs1_threshold_policy.py`](../tests/test_bs1_threshold_policy.py)冒頭のコメント、[`doc/external_data_coverage_ja.md`](external_data_coverage_ja.md) 27節)。**その同じ日のうちに方針が再転換され**、「未承認でも判定自体は出す(LLM等による予測として扱う)。未承認である/レビュー待ちであることはcuratorHint・descriptionに残す」という現在の形になった。過去の実インシデント自体の重みは変わらないので、新しいファイルを取り込む際も「未承認/暫定の値であることを必ずcuratorHintに残す」という担保だけは省略しないこと。

---

## 2. 依頼する作業

### (1) ClinGenの疾患別頻度閾値をダウンロードしてJSON化する

- 出典: ClinGen Criteria Specification Registry (<https://cspec.clinicalgenome.org/>) に掲載されている、各VCEPが正式承認したACMG/AMP仕様書(PDF)。
- `config/bs1_thresholds_draft.json`の既存7件は、このレジストリから手動転記されたものなので、フォーマットの参考にできる(ただし新しいファイルでは`status: "DRAFT"`は使わない。下記(3)参照)。
- 保存先: `config/`配下に新規ファイル(例: `config/bs1_disease_frequency_thresholds.json`)。
- 1エントリに必要なフィールドは「3. 期待するスキーマ」参照。

### (2) `bs1.py`がそのファイルを参照するようにする

現状`bs1.py`はファイルI/Oを一切しない設計なので、どこかに「(gene, condition) → 閾値」のルックアップを追加する必要がある。**参考にできる既存パターンが同じリポジトリ内にある**: [`acmg_pipeline/criteria/pp1_bs4_pp4_engine.py`](../acmg_pipeline/criteria/pp1_bs4_pp4_engine.py)の`load_reference_records()`(83-112行目)が、`config/pp4_reference_records.json`を遺伝子名でルックアップして返す、ほぼ同じ形の実装になっている。

2つの選択肢がある(どちらを取るかは受け取った開発者の判断に委ねる。「4. 未決定事項」参照):

- **案A**: `bs1.py`の`evaluate()`内に、pp1_bs4_pp4_engineと同様の自己完結したローダーを追加し、`disease_specific_threshold()`を呼ぶ前に`input_data["disease_frequency_threshold"]`を埋める。
- **案B**: `bs1.py`は今の「外から渡された辞書を読むだけ」という設計を保ち、`acmg_pipeline/services/resolve.py`の`ProviderEvidenceResolver`(他の外部データもここで解決している)か、`acmg_pipeline/automated_core/interface.py`の`criterion_input()`より手前の層に、gene/condition単位のルックアップを追加する。

### (3) DRAFTキューを廃止する

- `config/bs1_thresholds_draft.json`を削除し、(1)で作った新しいファイルに一本化する。
- 新しいファイルのエントリは、ロード時点で`status: DRAFT/APPROVED`のようなゲートを設けず、**そのまま最終判定の入力として使ってよいものとして扱う**(ClinGen自身のVCEP承認プロセスが、このプロジェクトが今までやっていた「人間キュレーターによる二重チェック」の代わりになる、という判断)。
- **(2026-09-17追記)** この項目は元々「`bs1.py`の`policy_status()`が未承認をMETにしない」という前提の上で、「じゃあ新しいファイルのエントリにstatusを持たせる意味は薄い」という文脈で書いた。しかしその後の方針転換で`policy_status()`自体は残るものの、**その結果はMET/UNKNOWNの分岐にはもう使われず、curatorHintの文言(「未承認/レビュー待ち」)にのみ反映される**ようになった。したがって新しいファイルに`status: DRAFT/APPROVED`を持たせるかどうかは、「MET/NOT_METを出すかどうかのゲート」としてではなく、**「curatorHintに何と書くか(承認済みVCEP仕様か、まだ内部で確認中か)」を決めるための情報**として位置づけ直す必要がある。実装する開発者は、`status`フィールドを完全に削るのではなく、この用途で残すかどうかを判断すること。
- **代替の安全策(必須、変更なし)**: 実行時にMET/NOT_METを止めるゲートが無くなった以上、このJSONを変更するPRは必ず人がレビューする、という運用ルールをファイル自体のコメント欄(`purpose`等)か`CLAUDE.md`に明記すること。過去のインシデントの再発を防ぐための最低限の担保として、これは省略しない。
- `bs1.py`の`policy_status()`(120-130行目)は、`config`側デフォルト閾値のフォールバック用ゲート(常にDRAFT)、および新しいファイルの`status`をcuratorHintの文言に反映するための判定として引き続き使われる。

---

## 3. 期待するスキーマ

`bs1.py`の`disease_specific_threshold()`(21-48行目)が`disease_frequency_threshold`辞書に要求するフィールド:

| フィールド | 必須 | 内容 |
|---|---|---|
| `condition` | ○ | このvariantの`input_data["condition"]`と**完全一致**する疾患名文字列(一致しないと閾値が適用されない) |
| `inheritance` | ○ | `input_data["inheritance"]`と一致する遺伝形式(例: `"autosomal_dominant"`) |
| `source` | ○ | 出典(VCEP名など) |
| `source_version` | ○ | 出典の版 |
| `reviewed_at` | ○ | レビュー日(ClinGenの公開日、または取り込み日) |
| `max_credible_af` | ○ | 0〜1の実数。最大許容頻度 |
| `frequency_statistic` | 任意 | `"af"`または`"faf95"`。省略時は`"af"`とみなされる(その旨が`review`に記録される) |
| `comparison` | 任意 | `">"`または`">="`。省略時は`">"`とみなされる(同上) |
| `policy_status` | 任意 | `"DRAFT"`にすると、閾値超え自体はMETとして出るが、curatorHintに「未承認/レビュー待ち」の注記が付く(120-130行目)。**(2026-09-17追記)** 以前は「DRAFTならMETを止める」ゲートだったが、方針転換によりcuratorHintへの注記のみに変わった。新しい取り込みフローでは基本的に設定しない想定 |

`gene`と`condition`の組でルックアップする設計にする場合、新JSONの各エントリにも`gene`キーを持たせること(`bs1_thresholds_draft.json`の既存7件も`gene`キーを持っている)。

---

## 4. 未決定事項(受け取った開発者が判断すること)

- 案A・案Bどちらの実装位置にするか。
- ルックアップキーを`gene`単独にするか、`(gene, condition)`の組にするか。1遺伝子が複数疾患に関わる場合(例: 表現型多様性のある遺伝子)を想定するなら後者が必須。
- `condition`文字列の一致は完全一致(現状の`bs1.py`の実装のまま)でよいか、正規化(大文字小文字・表記ゆれ吸収)を入れるか。
- ClinGenレジストリの更新(VCEPが仕様を改訂した場合)をどう追従するか — 定期的な手動更新か、何らかの自動チェックを入れるか。

---

## 5. テストへの影響

- [`tests/test_bs1_threshold_policy.py`](../tests/test_bs1_threshold_policy.py): 現在は`disease_frequency_threshold`を手動で組み立てたdictをテスト入力に直接渡している。ファイルベースのルックアップに変える場合、このテスト自体は変更不要(`bs1.evaluate()`に渡す`input_data`の作り方は変わらない)だが、新しいローダー関数自体の単体テストを追加すること。
- [`doc/external_data_coverage_ja.md`](external_data_coverage_ja.md) 331-333行目: `config/bs1_thresholds_draft.json`への言及があるので、ファイル名変更・DRAFT廃止に合わせて更新すること。
