# 提案: キュレーター向け注意喚起(curator hints)フォーマットの統一

**対象**: 統合VA-Spec出力(`acmg_pipeline.pipeline.evaluate_variant_evidence_lines()`が返す28行のEvidenceLine)における、文献判定エンジン(PS3/BS3/PS4)と自動判定エンジン(Layer1、16基準)の間の、キュレーター向け注意喚起フォーマットの不一致。

**このドキュメントの使い方**: このまま別の開発者(またはその開発者が使うAIコーディングエージェント)への実装指示として渡せるように書いています。「現状」「提案」「具体的な変更箇所」「テストへの影響」「未決定事項」の順に必要な情報を揃えています。

---

## 1. 現状の問題

同じ「キュレーターが確認すべき注意事項」という目的の情報が、文献側と自動判定側で**別々のフォーマット**で出力されている。

### 文献側(PS3/BS3/PS4) — `acmg_pipeline/export.py`

`_hints_extension()`(export.py 323行目付近)が、`curatorHints`という独立した拡張(extension)を作る:

```python
def _hints_extension(hints: list[CuratorHint]) -> Optional[Extension]:
    if not hints:
        return None
    return Extension(
        name="curatorHints",
        value=[{"severity": h.severity, "message": h.message} for h in hints],
    )
```

`CuratorHint`(`acmg_pipeline/common.py` 55行目)は `severity: str` (`"info"` / `"caution"` / `"warning"` の3値) と `message: str` を持つ。出力されるEvidenceLineの`extensions`配列に、`curatorHints`という名前で1つ入る。

### 自動判定側(Layer1、16基準) — `acmg_pipeline/automated_va_spec.py`

`assessment_details()`(automated_va_spec.py 226行目付近)が、`bh26AssessmentDetails`という別の拡張の**値の中に**、3つの別々のキーとして格納する:

```python
optional = {
    ...
    "reviewPoints": result.review_points,      # 要確認事項(文字列のリスト)
    "conflictFlags": result.conflict_flags,    # 矛盾フラグ(文字列のリスト)
    ...
    "warnings": result.warnings,               # 警告(文字列のリスト)
    ...
}
```

`severity`という概念自体が無く、`reviewPoints`/`warnings`/`conflictFlags`という**キー名の違いだけ**で重要度・種別を表現している。

### 何が問題か

- 統合ドキュメント(1変異体あたり28行のEvidenceLine)を見る側(人間キュレーター、あるいは今後作るUI/レポートツール)は、「このEvidenceLineに注意喚起があるか」を確認するために、そのEvidenceLineがどちらのエンジン由来かに応じて**見る場所を変える**必要がある。
- 文献側は`severity`でソート・フィルタできるが、自動判定側はできない(3つの別々のリストを個別にチェックする必要がある)。
- `bh26AssessmentDetails`自体は必須(status/strength等、GA4GH標準フィールドだけでは表現できない情報を運ぶ役目がある — 別途確認済み)なので、**この拡張自体をなくす提案ではない**。あくまで`reviewPoints`/`warnings`/`conflictFlags`という3フィールドの扱いについての提案。

---

## 2. 提案するフォーマット

`curatorHints`(文献側の形式)に統一する。ただし自動判定側の`conflictFlags`(自動判定同士の矛盾)が持つ「これは特に重大」という区別を失わないよう、`category`フィールドを追加する:

```json
{
  "name": "curatorHints",
  "value": [
    {"severity": "warning", "category": "conflict", "message": "..."},
    {"severity": "caution", "category": "review",   "message": "..."},
    {"severity": "warning", "category": "warning",  "message": "..."}
  ]
}
```

- `severity`: 既存の`"info"`/`"caution"`/`"warning"`の3値をそのまま使う。
- `category`(新規): 元のキー名を保持する形で `"review"`(旧`reviewPoints`) / `"warning"`(旧`warnings`) / `"conflict"`(旧`conflictFlags`)。文献側の既存hint(`CuratorHint`由来)には`category`を付けない(または`null`)。

**なぜこちらの形式に寄せるか**(理由は2つ):
1. 消費側が1箇所(`curatorHints`)だけ見れば済み、`severity`で横断的にソート・フィルタできる。
2. `curatorHints`を使っているのはexport.py 1ファイルだけなので、自動判定側(16基準のモジュールに影響)を新形式に合わせて出力させる方が、結果的に1つの統一フォーマットに収束させやすい(下記「変更範囲」参照。ただし変更箇所自体は自動判定側の方が多い)。

---

## 3. 具体的な変更箇所

### 3-1. `acmg_pipeline/automated_va_spec.py`

`assessment_details()`(226行目付近)と`to_evidence_line()`(338行目付近)を変更する。

- `assessment_details()`から`reviewPoints`/`warnings`/`conflictFlags`を`bh26AssessmentDetails`の値に含めるのをやめる(あるいは、後方互換のため残すかは要相談 — 下記「未決定事項」参照)。
- `to_evidence_line()`内で、`result.review_points`/`result.warnings`/`result.conflict_flags`から`curatorHints`形式のリストを組み立て、`extensions`に`{"name": "curatorHints", "value": [...]}`として追加する:

```python
def _curator_hints_from_result(result) -> list[dict]:
    hints = []
    for message in result.conflict_flags:
        hints.append({"severity": "warning", "category": "conflict", "message": message})
    for message in result.review_points:
        hints.append({"severity": "caution", "category": "review", "message": message})
    for message in result.warnings:
        hints.append({"severity": "warning", "category": "warning", "message": message})
    return hints
```

(severityの割り当ては仮 — 実際にどのseverityにするかはチームで要確認)

- `export_record()`(400行目付近)も、`assessments`リストに積む`details`(=`assessment_details()`の戻り値)から該当キーが消えることを踏まえて確認が必要。

### 3-2. `acmg_pipeline/export.py`

変更不要(文献側は既に`curatorHints`形式)。ただし`_hints_extension()`が今後「両エンジン共通のヘルパー」になるなら、`acmg_pipeline/common.py`などの共有モジュールに移動する方が自然かもしれない(任意)。

### 3-3. 影響を受けるテスト

- `tests/test_va_spec.py` 191行目: `assessment["reviewPoints"]`を直接assertしている箇所を、新しい`curatorHints`の中身を見るように書き換える必要がある。
- `tests/test_curated.py` / `tests/test_population.py` / `tests/test_demo_pipeline.py`: `value.review_points`(`CriterionResult`自体の属性)をassertしている箇所は**変更不要**(このリファクタは「VA-Spec出力への変換方法」だけを変えるもので、`CriterionResult.review_points`等の内部属性自体は変えない)。ただし`test_demo_pipeline.py`150行目の`result["review_points"]`は`results.json`(内部フォーマット)を見ているだけなので、これも影響なし。

---

## 4. 未決定事項(チームで合意が必要)

1. **`bh26AssessmentDetails`から`reviewPoints`/`warnings`/`conflictFlags`を完全に削除するか、`curatorHints`と併記(重複)させるか。** 完全削除の方がフォーマットは綺麗になるが、既存の`bh26AssessmentDetails`だけを見ている他のコード(例: `pipeline_interface._evidence_from_line()`は`status`/`strength`しか見ていないので影響なし、と思われるが要確認)がないか、事前に洗い出す必要がある。
2. **severityの割り当てルール**: `conflictFlags`→`warning`、`review_points`→`caution`、`warnings`→`warning`、という割り当ては本ドキュメントの仮案。実際の重大度感覚に合わせて調整が必要。
3. **この統一を自動判定エンジン側(Layer1)のコードで行うか、統合レイヤー(`evaluate_variant_evidence_lines()`やその周辺)で後付け変換するか。** 前者(本ドキュメントの提案)は自動判定エンジン自身の出力が最初から正しい形になるが、Layer1担当チームの複数ファイルに影響する。後者は影響範囲を統合レイヤーだけに限定できるが、「なぜ自動判定エンジン自身のVA-Spec出力とpipeline.py経由の出力が違うのか」という新たな不一致を生む。

---

## 5. 背景(参考)

`bh26AssessmentDetails`自体が必要な理由: GA4GH VA-Specの標準フィールド(`directionOfEvidenceProvided`)は`supports`/`neutral`/`disputes`の3値しかなく、「NOT_MET」と「UNKNOWN(未評価)」がどちらも`neutral`に潰れてしまう。`bh26AssessmentDetails.status`だけが明示的に`met`/`not_met`/`unknown`の3値を保持しており、`pipeline_interface._evidence_from_line()`(28行のEvidenceLineを`classification.classify()`用の`CriterionEvidence`に変換する関数)はこの値だけを読んでいる。この部分は今回の提案の対象外で、変更しない。
