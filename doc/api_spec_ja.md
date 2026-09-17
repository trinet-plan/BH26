# ACMG判定API エンドポイント仕様

[`README.md`](../README.md) の「APIサーバーの起動」でサーバーを立ち上げた後に参照する、
`api/main.py` が公開するエンドポイントのリファレンス。実装そのものは
[`api/main.py`](../api/main.py)・ジョブの状態遷移は [`api/job_store.py`](../api/job_store.py) を参照。

いずれのエンドポイントも入力は「1バリアント分のVCF(`vcf`)+ 任意の臨床所見(`clinical_note`)」
共通で、リクエストボディはJSON。

ACMG/AMP 2015の28基準は、判定ロジックの実装状況で3グループに分かれる
(詳細は [`acmg_pipeline/constants.py`](../acmg_pipeline/constants.py) 参照)。

| グループ | 基準コード | 特徴 |
|---|---|---|
| `AUTOMATED_CODES` (16件) | `PVS1, PS1, PM1, PM2, PM4, PM5, PP2, PP3, PP5, BA1, BS1, BP1, BP3, BP4, BP6, BP7` | 決定的ルール判定。LLM不要・高速 |
| `LITERATURE_CODES` (3件) | `PS3, BS3, PS4` | PubMed全文検索 + LLM判定。低速(数十秒〜数分/variant) |
| `STUB_CODES` (9件) | `PS2, PM3, PM6, PP1, PP4, BS2, BS4, BP2, BP5` | 未実装。常に`NOT_EVALUATED`を返す |

---

## `GET /health`

死活監視用。vLLM等の設定に関わらず常に200を返す。

```json
{"status": "ok"}
```

---

## `POST /v1/classify_criteria`

全28基準を評価し、ACMG分類(Pathogenic/Likely Pathogenic/VUS/Likely
Benign/Benign)まで出す。`LITERATURE_CODES`を必ず評価するため**常に非同期**
(job_id + polling)。

**リクエスト**

```json
{"vcf": "##fileformat=VCFv4.2\n...", "clinical_note": "..."}
```

| フィールド | 型 | 必須 | 説明 |
|---|---|---|---|
| `vcf` | string | ✓ | 1バリアント分のVCFテキスト(ヘッダ込み) |
| `clinical_note` | string | – | 臨床所見の自然文(省略時は空文字) |

**レスポンス — 202 Accepted**

```json
{"job_id": "3e1b7c2a-...", "status": "queued", "poll_url": "/v1/classify_criteria/3e1b7c2a-..."}
```

**エラー**

| status | 条件 |
|---|---|
| 400 | `vcf` のパースに失敗(不正なVCF) |

## `GET /v1/classify_criteria/{job_id}`

ジョブの状態・結果をポーリングする。

```json
{
  "job_id": "3e1b7c2a-...",
  "status": "queued | running | succeeded | failed",
  "created_at": "2026-09-17T02:40:13.99Z",
  "variant": {"gene": "MYBPC3", "hgvsc": "c.1144C>T", "hgvsp": "p.(Arg382Trp)"},
  "started_at": "2026-09-17T02:40:13.99Z",
  "completed_at": "2026-09-17T02:42:38.95Z",
  "va_spec": { "...": "status=succeededのときのみ。ACMG分類+28基準のVA-Spec文書" },
  "error": { "type": "ExceptionGroup", "message": "..." }
}
```

`va_spec`は`status: "succeeded"`のときのみ、`error`は`status: "failed"`のときのみ値を持つ。
404は`job_id`が存在しない場合。

## `DELETE /v1/classify_criteria/{job_id}`

ジョブと結果を破棄する。成功時204、`job_id`が存在しなければ404。

---

## `POST /v1/get_evidence_line_by_target_criteria`

指定した基準だけを評価し、VA-Spec `EvidenceLine`を基準ごとに返す(ACMG分類は行わない)。
`LITERATURE_CODES`を1つも含まなければ**同期**、1つでも含めば**非同期**(job_id + polling)。

**リクエスト**

```json
{"vcf": "##fileformat=VCFv4.2\n...", "clinical_note": "...", "criteria": ["PM2", "BA1"]}
```

| フィールド | 型 | 必須 | 説明 |
|---|---|---|---|
| `vcf` | string | ✓ | 1バリアント分のVCFテキスト(ヘッダ込み) |
| `clinical_note` | string | – | 臨床所見の自然文(省略時は空文字) |
| `criteria` | string[] | ✓ | 評価したいACMGコード(例: `["PM2", "PS3"]`)。1件以上、既知コードのみ |

**レスポンス — 200 OK**(`criteria`が`AUTOMATED_CODES`/`STUB_CODES`のみの場合)

```json
{
  "status": "succeeded",
  "criteria": ["BA1", "PM2"],
  "evidence_lines": {"PM2": { "...": "VA-Spec EvidenceLine" }, "BA1": { "...": "VA-Spec EvidenceLine" }}
}
```

**レスポンス — 202 Accepted**(`criteria`に`PS3`/`BS3`/`PS4`のいずれかを含む場合)

```json
{"job_id": "c25e7233-...", "status": "queued", "poll_url": "/v1/get_evidence_line_by_target_criteria/c25e7233-..."}
```

**エラー**

| status | 条件 |
|---|---|
| 400 | `vcf` のパースに失敗 / `criteria`が空 / 未知のACMGコードを含む |

## `GET /v1/get_evidence_line_by_target_criteria/{job_id}`

`GET /v1/classify_criteria/{job_id}` と同じ形のジョブ状態を返すが、`va_spec`フィールドの中身は
`{code: EvidenceLine}`の辞書(ACMG分類は含まない)。404は`job_id`が存在しない場合。

## `DELETE /v1/get_evidence_line_by_target_criteria/{job_id}`

ジョブと結果を破棄する。成功時204、`job_id`が存在しなければ404。
