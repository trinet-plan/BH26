# ACMG判定API 実装方針 v1

このドキュメントは、`acmg_pipeline` のLLM判定パイプラインをREST API化して
Docker上でホストするための実装方針を定義する。すでに `acmg_pipeline/api_input.py`
や `acmg_pipeline/vcf_record.py` のdocstringから section 0-1 が参照されているため、
本ドキュメントの節番号はそれらの参照と一致させてある。

## 0. 前提・スコープ

- **このドキュメントのスコープ**: REST APIの外形(エンドポイント、リクエスト/レスポンス、
  ジョブ管理・ポーリング、パイプラインとのインターフェース契約、va_spec出力の組み立て、
  Docker配置)を定義する。
- **スコープ外(別担当が実装)**: `acmg_pipeline/criteria/*` や `gate.py` の中身、
  vLLM/PubMed MCP呼び出しを含む判定ロジックそのもの。本ドキュメントでは
  「1関数呼び出し」として抽象化し、そのインターフェース(引数・戻り値の型)だけを
  契約として定義する(4章参照)。パイプライン担当者はこの契約さえ満たせば、
  内部実装(LLMのプロンプト、ERepo参照、判定アルゴリズム等)を自由に変更できる。
- **既存資産の再利用**:
  - `acmg_pipeline/api_input.py`(`ApiCaseInput`) … リクエストの生入力パース
  - `acmg_pipeline/classification.py`(`classify`, `ClassificationResult`) … 最終判定
  - `acmg_pipeline/export.py`(`build_evidence_line`) … 基準1件分のVA-Spec EvidenceLine化
  - これらは変更せず、その上にAPI層を薄く被せる。

## 0-1. リクエスト契約: 1リクエストにつき1件POST

`ApiCaseInput.parse_vcf()` がすでに「VCF文字列のデータ行は必ず1件」を強制している
(`acmg_pipeline/vcf_record.py` の `parse_vcf` が0件/2件以上で `ValueError`)。
これに合わせ、API側も **1つのPOSTリクエスト = 1バリアントの判定依頼** を契約として固定する。
複数バリアントの一括投入(バッチAPI)は現時点でスコープ外とし、必要になった時点で
別エンドポイント(例: `POST /v1/variants/batch`)として追加を検討する。

## 0-2. 入力スキーマ(TODO)

`democase/case*_api_input_*.json` の `{"vcf": "...", "clinical_note": "..."}` は
あくまで検証用のデモ入力(VEPパーサ・LLM抽出のテストフィクスチャ)であり、
「正規化済み構造化データ」としての最終的なリクエストボディスキーマは未確定のTODOである。

暫定方針:
- **v1(暫定)**: `ApiCaseInput` と同じ形 `{"vcf": str, "clinical_note": str}` を
  そのままリクエストボディとして受け付ける。パース・型変換は既存の
  `ApiCaseInput.parse_vcf()` / `ClinicalNoteExtraction.from_json()` に委譲する。
  これによりデモデータ(`democase/case1_api_input_case1-noise2.json` 等)を
  そのままAPIに投げて動作確認できる。
- **v2(将来)**: gene/HGVSc/HGVSp/zygosity等を個別フィールドに正規化した
  構造化スキーマへ移行する。移行時は `ApiCaseInput` に
  `from_normalized_dict()` のような別コンストラクタを追加する形にし、
  `parse_vcf()` 依存のVCF文字列パースを廃止しない(v1入力との後方互換を保つ)。
  この移行の要否・スキーマ詳細は別途要検討。

## 1. 全体アーキテクチャ

パイプラインは複数論文のPubMed全文取得 + LLM呼び出し(1論文1呼び出し、基準ごと)を
直列に行うため、1バリアントあたり数十秒〜数分かかりうる(`acmg_pipeline/pipeline.py`
の `judge_variant` が典型: PMID数 × 基準数 だけLLM呼び出しが発生する)。
同期HTTPリクエストでこれを待たせるのは非現実的なため、**非同期ジョブ + ポーリング**方式を採る。

パイプライン(`run_pipeline`, 4章)は別サービスではなく、このリポジトリ内の
`acmg_pipeline` モジュールとしてAPIサーバーと**同一プロセス**にimportされて動く想定
(別コンテナ・別APIとして切り出す予定はない)。そのため登場人物は
Client / API server / Job store の3者のみで、パイプラインはAPI server内部の
一関数呼び出し(バックグラウンドワーカースレッド上でのin-process呼び出し)として表す。
Job storeを直接読み書きするのは常にAPI serverであり、パイプライン自身がJob storeに
触れることはない(戻り値をAPI serverに返すだけ)。

```
Client                         API server                              Job store
                          (uvicornプロセス。バックグラウンド
                           ワーカーが acmg_pipeline.run_pipeline()
                           を直接呼び出す)

  | POST /v1/variant          |                                            |
  |--------------------------->| create job (queued)                       |
  |                            |------------------------------------------>|
  |<---- 202 {job_id} ---------|                                            |
  |                            | job -> running へ更新                      |
  |                            |------------------------------------------>|
  |                            | run_pipeline(case)  # 同一プロセス内の関数呼び出し
  |                            | (worker thread / executorで実行、戻り値を待つ)
  |                            | job -> succeeded, va_spec を書き込み        |
  |                            |------------------------------------------>|
  |          ...(ポーリング継続)...                                          |
  | GET /v1/variant/{id}       |                                            |
  |--------------------------->| job読み取り                                 |
  |                            |<-------------------------------------------|
  |<-- 200 {status: succeeded, va_spec: {...}} --|                          |
```

## 2. エンドポイント設計

| メソッド | パス | 説明 |
|---|---|---|
| `POST` | `/v1/variant` | 1バリアント分の判定ジョブを投入。202 Acceptedで `job_id` を返す |
| `GET` | `/v1/variant/{job_id}` | ジョブの状態・(完了していれば)結果を取得(ポーリング用) |
| `DELETE` | `/v1/variant/{job_id}` | (任意)ジョブ・結果の破棄。ストレージ肥大対策 |
| `GET` | `/health` | ヘルスチェック(vLLM/PubMed MCP疎通確認は含めない。あくまでAPIプロセス生存確認) |

`GET /v1/variant/{job_id}/result` のように結果取得を別エンドポイントに分ける案もあるが、
ポーリング対象が1つで済む方がクライアント実装が単純なため、v1は `GET /v1/variant/{job_id}`
のレスポンス内に状態と結果を同居させる方針とする。

### POST /v1/variant

リクエスト: 0-2節の入力スキーマ(v1暫定形)。

レスポンス(202):
```json
{
  "job_id": "3e1b7c2a-...-uuid",
  "status": "queued",
  "poll_url": "/v1/variant/3e1b7c2a-...-uuid"
}
```

### GET /v1/variant/{job_id}

```json
{
  "job_id": "3e1b7c2a-...-uuid",
  "status": "succeeded",
  "created_at": "2026-09-16T10:00:00Z",
  "started_at": "2026-09-16T10:00:01Z",
  "completed_at": "2026-09-16T10:04:32Z",
  "variant": {"gene": "MYBPC3", "hgvsc": "c.278delA", "hgvsp": "p.(Lys93ArgfsTer3)"},
  "va_spec": { "...": "5章参照" },
  "error": null
}
```

`status` は3章の状態遷移に従う。`succeeded` 以外では `va_spec` は `null`、
`failed` では `error` に理由を格納する。

## 3. ジョブモデル・状態遷移

```
queued --> running --> succeeded
                    \-> failed
```

| 状態 | 意味 |
|---|---|
| `queued` | 受理済み、実行待ち |
| `running` | パイプライン実行中 |
| `succeeded` | 完了、`va_spec` が利用可能 |
| `failed` | 例外・タイムアウト等で失敗。`error` に詳細 |

途中経過(「今どの基準/論文を処理中か」)を返す `progress` フィールドは
v1では**含めない**(パイプライン内部の進捗をAPI層に伝播する仕組みが
別実装側にまだ無いため)。必要になった場合、パイプライン側から
コールバックで進捗を書き込めるインターフェースを別途追加する
(7章「未決事項」参照)。

## 4. パイプラインとのインターフェース契約

パイプライン担当者(別の人)が実装する関数は、以下の1関数に集約する。
API層はこの関数をバックグラウンドで呼び出すだけで、中身(MCP/vLLM呼び出し、
ERepoゲート、基準ごとの判定等)には関与しない。

```python
# 別担当が実装する関数のシグネチャ契約
def run_pipeline(case: ApiCaseInput) -> PipelineOutput:
    ...

@dataclass
class PipelineOutput:
    classification: ClassificationResult          # acmg_pipeline.classification
    evidence_lines: dict[str, dict]                # {"PS3": {...EvidenceLine dict...}, ...}
```

- `classification` は `acmg_pipeline.classification.classify()` の戻り値そのもの
  (28コード全件分の met/not_met/not_evaluated を含む)。
- `evidence_lines` は、実際にLLMパイプラインが走った基準(`IMPLEMENTED_CODES` の
  部分集合)についてのみ、`acmg_pipeline.export.build_evidence_line()` が返す
  dict(paper単位のネストEvidenceLineを含む)をキー=ACMGコードで格納したもの。
  `pipeline.py` の `main()` が現在ローカルファイルに書き出しているのと同じ形を、
  ファイルに書く代わりにこのdictとして返してもらう想定。
- `classification.met` の各 `CriterionEvidence.source` は
  `"llm_pipeline:PMID..."` のような要約文字列に過ぎず、論文単位の詳細
  (方向・rationale・curator hint)を失っている。API側のva_spec出力で
  その詳細を保持するために `evidence_lines` を別途要求している
  (`ClassificationResult` 単体からは論文レベルの情報を復元できないため)。
- 例外・タイムアウトはこの関数呼び出しの外側(API層)でハンドリングする
  (7章)。`run_pipeline` 自身はPython例外を投げて失敗を伝えてよい。
- 同期関数として定義する(`pipeline.py` の現状が `asyncio.run(main())` で
  完結する同期プロセスであるため)。非同期化(`async def run_pipeline`)は
  パイプライン担当者の裁量に委ねるが、いずれの場合もAPI層の
  ワーカー実行(6章)は対応できるようにする。

## 5. ClassificationResult → va_spec JSON への変換

`GET /v1/variant/{job_id}` が返す `va_spec` フィールドの中身。
新規モジュール `acmg_pipeline/va_spec_statement.py` を追加し、
`PipelineOutput` から1つのJSON文書を組み立てる関数を実装する
(この変換ロジック自体は本APIの実装範囲内であり、パイプライン担当者の
スコープ外)。

```python
def build_variant_statement(output: PipelineOutput, gene: str, hgvsc: str) -> dict:
    ...
```

出力イメージ(v1・暫定形。`export.py` が最初 hand-built dict → 後に
`ga4gh.va_spec` の実Pydanticモデルへ移行した前例に倣い、まずは手組みの
dictで実装し、`VariantPathogenicityStatement` 等の実モデルが要件に
合うか確認でき次第移行する):

```json
{
  "id": "stmt:MYBPC3_c_278delA",
  "subjectVariant": {"gene": "MYBPC3", "hgvsc": "c.278delA", "hgvsp": "p.(Lys93ArgfsTer3)"},
  "classification": {
    "category": "Likely Pathogenic",
    "score": 8,
    "ba1Override": false
  },
  "hasEvidenceLines": [
    { "...": "evidence_lines['PVS1'] の中身(build_evidence_lineの出力そのもの)" }
  ],
  "criteriaNotEvaluated": ["PS1", "PM1", "..."],
  "criteriaMet": [{"code": "PVS1", "strength": "very_strong"}],
  "criteriaNotMet": [{"code": "BS3"}],
  "generatedAt": "2026-09-16T10:04:32Z"
}
```

- `hasEvidenceLines` には `evidence_lines` dictの値(基準ごとのEvidenceLine)を
  そのまま列挙する。`classification.met`/`not_met` に含まれるが
  `evidence_lines` に無いコード(=stub、またはLLM以外の手段で判定されたコード)は
  `criteriaMet`/`criteriaNotMet` の簡易表現(code + strengthのみ)として残し、
  詳細なEvidenceLineは無い旨が分かるようにする。
- `criteriaNotEvaluated` は `ClassificationResult.not_evaluated_codes` をそのまま使う。

## 6. 実行基盤

- フレームワークは **FastAPI**(既存コードがPython/dataclass中心のため親和性が高い)。
- ジョブ実行:
  - v1(最小構成): FastAPIの `BackgroundTasks` + プロセス内 `dict` によるジョブストア。
    単一ワーカープロセスでのデモ運用を想定。デプロイをシングルコンテナ・
    シングルワーカーに限定する前提で成立する。
  - 将来(スケール要件が出た場合): Redis等をブローカーにした
    タスクキュー(RQ/Celery/arq)へ移行。マルチワーカー・複数コンテナ構成でも
    ジョブ状態を共有できるようにする。`run_pipeline` のインターフェース
    (4章)自体はどちらの構成でも変わらない。
  - `run_pipeline` の実処理は同期・ブロッキングI/O(vLLM/PubMed MCP待ち)を
    含むため、FastAPIのイベントループを塞がないよう
    `run_in_executor` またはワーカースレッド/プロセスで実行する。
- ジョブストアはv1ではプロセス内メモリで良いが、コンテナ再起動で
  ジョブ状態が消える点は許容する(デモ用途のため)。永続化が必要になれば
  SQLite/Redis等に差し替える。

## 7. エラーハンドリング・タイムアウト

- 入力バリデーション(`ApiCaseInput.parse_vcf()` の `ValueError` 等)は
  ジョブ作成前に検査し、`400 Bad Request` を即時返す(ジョブを作らない)。
- `run_pipeline` 実行中の例外はジョブを `failed` にし、`error` に
  例外の型・メッセージを格納する(スタックトレースはログにのみ残し、
  レスポンスには含めない)。
- タイムアウト: パイプラインの想定所要時間に対して十分マージンを取った
  上限(例: 30分)を設け、超過したジョブは `failed`(`error.reason=timeout`)
  とする。具体的な上限値はパイプライン担当者の実測値を踏まえて決定する
  (未決事項)。

## 8. Docker化

- 単一の `Dockerfile` でAPIプロセスを起動(`uvicorn` 経由)。
- vLLM(`VLLM_BASE_URL`/`VLLM_API_KEY`)・PubMed MCPのエンドポイントは
  `pipeline.py` 同様、環境変数 + `.env` で注入する。
- `docker-compose.yml`(将来、タスクキュー導入時): `api` / `redis` /
  (必要なら)`worker` サービスに分割する。v1(BackgroundTasks構成)では
  `api` サービス単体で完結させる。
- vLLM本体・PubMed MCPサーバーはこのAPIコンテナの外(既存の
  `113.43.212.98:8000` 等)にあり、コンテナ化の対象外。

## 9. テスト方針

pytest等のテストフレームワークによる単体テストは組まない。
サーバーを実際に起動し、**curlで叩いて期待通りのJSONが返ってくるかを確認する**
だけのシンプルな確認方法とする。

```bash
# 1. サーバー起動(別ターミナル、またはバックグラウンド)
uvicorn api.main:app --reload

# 2. ジョブ投入(democaseのJSONファイルをそのままボディにできる = 0-2節のv1入力スキーマ)
JOB_ID=$(curl -s -X POST http://localhost:8000/v1/variant \
  -H "Content-Type: application/json" \
  -d @democase/case1_api_input_case1-noise2.json | jq -r .job_id)
echo "job_id=$JOB_ID"

# 3. ポーリング(succeeded/failedになるまでループ)
until [ "$(curl -s http://localhost:8000/v1/variant/$JOB_ID | jq -r .status)" != "running" ] 2>/dev/null; do
  sleep 5
done

# 4. 結果のva_specを確認
curl -s http://localhost:8000/v1/variant/$JOB_ID | jq .
```

確認する観点:
- `status` が最終的に `succeeded`(または想定通り `failed`)になるか
- `va_spec.classification.category` / `.score` が想定した分類(例:
  `case1-var1` = MYBPC3 c.278delA なら Likely Pathogenic寄り)になっているか
- `va_spec.hasEvidenceLines` に実行した基準分のEvidenceLineが入っているか
- `va_spec.criteriaNotEvaluated` に stub コードが正しく列挙されているか

実パイプライン未実装の間は、`run_pipeline` を4章の契約(`PipelineOutput` を
返す)だけ満たすダミー関数に差し替えて同じcurl手順でAPI層(ジョブ管理・
ポーリング・va_spec変換)の動作を先に確認する。実パイプライン結合後は
同じcurlコマンドをそのまま使い回して、返ってくる値が本物の判定結果に
置き換わっていることを確認する。

## 10. 未決事項・TODO

- 0-2節: 正規化済み入力スキーマの詳細(v2)
- 3章: 実行中の進捗(`progress`)をどう伝播するか
- 5章: `va_spec` の最終形を `ga4gh.va_spec` の実モデル
  (`VariantPathogenicityStatement` 等)で構築できるか要検証
  (`export.py` のEvidenceLineと同様、実在するフィールド構成を
  ライブラリ側で確認してから移行する)
- 6章: マルチワーカー化が必要になるタイミング(タスクキュー移行の判断基準)
- 7章: タイムアウト上限値(パイプライン担当者の実測を待つ)
