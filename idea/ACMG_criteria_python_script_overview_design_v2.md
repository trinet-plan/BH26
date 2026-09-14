# ACMG / ClinGen General Variant Interpretation
## Criterion単位 Pythonスクリプト 実装概要設計 v2

## 1. 目的

本資料は、ACMG/AMP 2015 + ClinGen General Guidance に基づく variant interpretation を、
**ACMG criterion ごとに独立した Python スクリプトとして実装するための概要設計**を定義する。

対象は ACMG/AMP 2015 の28 criteria とする。

```text
PVS1
PS1 PS2 PS3 PS4
PM1 PM2 PM3 PM4 PM5 PM6
PP1 PP2 PP3 PP4 PP5
BA1
BS1 BS2 BS3 BS4
BP1 BP2 BP3 BP4 BP5 BP6 BP7
```

本設計では、以下を重要な前提とする。

- 1 criterion = 1 Python スクリプト
- 1回の実行で、1 variant に対する1 criterionの評価結果を返す
- Variant を必須入力とする
- Disease / Phenotype / Family 等は optional input とする
- DB/APIアクセスはcriterionスクリプトから分離する
- 外部データソースは将来的に差し替え・追加可能とする
- 特にPopulation DBは国・地域ごとのデータソースを追加可能とする
- criterionロジックは特定DB名に依存させない
- Evidence取得結果とcriterion判定結果を分離する
- `MET` と `NOT_MET` と `NOT_EVALUATED` を明確に区別する
- 最終Evidence生成には **VA-SPEC-Python** を使用する
- criterionスクリプトからVA-SPEC-Pythonを直接呼ばず、Mapper層で変換する
- CSpec / VCEP-specific rule は本実装の基本スコープ外とする

---

# 2. システム全体像

```text
                         ┌──────────────────────┐
                         │ Variant Input        │
                         │ VCF / JSON           │
                         └──────────┬───────────┘
                                    │
                                    v
                         ┌──────────────────────┐
                         │ Input Normalization  │
                         │ split / normalize    │
                         └──────────┬───────────┘
                                    │
                                    v
                       ┌─────────────────────────┐
                       │ Logical Services        │
                       │ population              │
                       │ annotation              │
                       │ variant classification  │
                       │ gene-disease            │
                       │ functional              │
                       │ computational           │
                       │ phenotype / family      │
                       └────────────┬────────────┘
                                    │
                                    v
                       ┌─────────────────────────┐
                       │ Provider Registry       │
                       │ population:             │
                       │  ├ gnomAD               │
                       │  ├ jMorp                │
                       │  ├ Country A DB         │
                       │  └ Country B DB         │
                       └────────────┬────────────┘
                                    │
                                    v
                       ┌─────────────────────────┐
                       │ Normalized Evidence     │
                       │ common data model       │
                       └────────────┬────────────┘
                                    │
                ┌───────────────────┼────────────────────┐
                │                   │                    │
                v                   v                    v
             pvs1.py             pm2.py               pp3.py
                │                   │                    │
                v                   v                    v
         CriterionResult     CriterionResult      CriterionResult
                │                   │                    │
                └───────────────────┼────────────────────┘
                                    │
                                    v
                       ┌─────────────────────────┐
                       │ VA-Spec Mapper          │
                       └────────────┬────────────┘
                                    │
                                    v
                       ┌─────────────────────────┐
                       │ VA-SPEC-Python          │
                       │ model / validation      │
                       └────────────┬────────────┘
                                    │
                                    v
                       ┌─────────────────────────┐
                       │ VA-Spec EvidenceLine    │
                       │ JSON serialization      │
                       └─────────────────────────┘
```

---

# 3. 評価単位

評価単位は、

> **1 Variant / 1 ALT allele / 1 Criterion**

とする。

VCF:

```vcf
#CHROM POS      REF ALT
1      123456   A   G
3      456789   G   A,C
```

内部では以下に分解する。

```text
1:123456 A>G
3:456789 G>A
3:456789 G>C
```

multi-allelic variant は評価前にsplit / normalizeする。

---

# 4. 共通入力

## 4.1 最小入力

```json
{
  "variant": {
    "assembly": "GRCh38",
    "chrom": "13",
    "pos": 32316461,
    "ref": "C",
    "alt": "T"
  }
}
```

## 4.2 拡張入力

```json
{
  "variant": {
    "assembly": "GRCh38",
    "chrom": "13",
    "pos": 32316461,
    "ref": "C",
    "alt": "T"
  },
  "condition": {
    "mondo_id": "MONDO:xxxxxxx",
    "label": "example disease"
  },
  "inheritance": "autosomal_dominant",
  "case": {
    "case_id": "CASE001",
    "hpo_terms": ["HP:0001250", "HP:0001263"],
    "age": 35,
    "sex": "female"
  },
  "family": {
    "proband_genotype": "0/1",
    "father_genotype": "0/0",
    "mother_genotype": "0/0",
    "maternity_confirmed": true,
    "paternity_confirmed": true
  }
}
```

---

# 5. Criterionスクリプトの基本責務

各criterionスクリプトは以下のみを担当する。

```text
1. 入力確認
2. criterion適用対象確認
3. 必要Evidenceの要求
4. Evidence品質確認
5. criterion固有ルール評価
6. strength決定
7. CriterionResult生成
```

以下はcriterionスクリプトの責務外とする。

```text
VCF parserの詳細
HTTP/API低レベル処理
各DB固有レスポンス解析
Provider選択
Population DBの優先順位解決
VA-Spec schema生成
VA-SPEC-Pythonの直接操作
JSON serialization詳細
ログ基盤
```

各criterionファイルには、

> **そのcriterion固有の判定ロジックだけを書く**

ことを基本原則とする。

---

# 6. Criterionスクリプト 共通インターフェース

```python
def evaluate(input_data, services, config):
    # Evaluate one ACMG criterion for one variant.
    # Returns CriterionResult.
    ...
```

例:

```python
from core.result import CriterionResult

CRITERION = "PM2"

def evaluate(input_data, services, config):
    variant = input_data.variant

    population = services.population.get_resolved_evidence(
        variant=variant,
        context=input_data
    )

    if population is None:
        return CriterionResult.not_evaluated(
            criterion=CRITERION,
            reason="Population evidence unavailable"
        )

    # PM2-specific rule
    ...

    return result
```

CLI実行も可能にする。

```bash
python criteria/pm2.py   --input input/variant.json   --output output/PM2.internal.json
```

---

# 7. 共通評価ステータス

| Status | 意味 |
|---|---|
| `MET` | criterion成立 |
| `NOT_MET` | 必要情報が揃った上でcriterion不成立 |
| `NOT_EVALUATED` | 情報不足・品質不足等により評価不能 |
| `NOT_APPLICABLE` | variant type等から適用対象外 |
| `MANUAL_REVIEW` | Curator確認が必要 |
| `DEPRECATED` | ClinGen General方針上使用しない |

重要:

```text
NOT_MET != NOT_EVALUATED
```

---

# 8. CriterionResult

各criterionスクリプトは、VA-Specそのものではなく、
まず内部共通形式 `CriterionResult` を返す。

例:

```json
{
  "criterion": "PM2",
  "status": "MET",
  "direction": "supports",
  "strength": "supporting",
  "evidence_outcome": "PM2_supporting",
  "variant": {
    "assembly": "GRCh38",
    "chrom": "13",
    "pos": 32316461,
    "ref": "C",
    "alt": "T"
  },
  "summary": "Variant is absent or sufficiently rare in the resolved population evidence.",
  "evidence": [
    {
      "type": "population_frequency",
      "source": "jMorp",
      "source_version": "configured-version",
      "population": "JPN",
      "AC": 1,
      "AN": 100000,
      "AF": 0.00001
    },
    {
      "type": "population_frequency",
      "source": "gnomAD",
      "source_version": "configured-version",
      "population": "global",
      "AC": 2,
      "AN": 1200000,
      "AF": 0.00000167
    }
  ],
  "missing_inputs": [],
  "review_points": [],
  "conflict_flags": [],
  "provenance": {
    "rule_version": "PM2-v1.0",
    "provider_profile": "japan",
    "implementation_spec": "ACMG/AMP 2015 + ClinGen General Guidance",
    "evaluated_at": "ISO-8601 timestamp"
  }
}
```

---

# 9. Provider Plugin Architecture

## 9.1 基本思想

criterionロジックから、具体的なDB名を分離する。

悪い例:

```python
gnomad = GnomadClient()
result = gnomad.get_frequency(variant)
```

推奨:

```python
result = services.population.get_resolved_evidence(
    variant=variant,
    context=input_data
)
```

criterionは、

```text
Population Evidenceが必要
```

ということだけを知り、

```text
gnomADを使う
jMorpを使う
Country X DBを使う
```

という判断はService / Provider Registry側へ委譲する。

---

# 10. 層構造

```text
Criterion
    ↓
Logical Service
    ↓
Provider Registry
    ↓
Provider Adapter
    ↓
External DB / Local DB / API / File
```

例:

```text
PM2
 ↓
PopulationService
 ↓
PopulationProviderRegistry
 ↓
├ GnomADProvider
├ JMorpProvider
└ CountrySpecificProvider
```

---

# 11. Population Provider

## 11.1 共通Interface

各Population DBは共通interfaceを実装する。

```python
from abc import ABC, abstractmethod

class PopulationProvider(ABC):

    @abstractmethod
    def get_frequency(self, variant, context=None):
        # Return PopulationEvidence or None.
        raise NotImplementedError
```

## 11.2 Provider固有処理

各DBの物理フィールド差はProvider内部で吸収する。

```text
gnomAD raw response
       ↓
GnomADProvider
       ↓
PopulationEvidence

jMorp raw response
       ↓
JMorpProvider
       ↓
PopulationEvidence
```

criterion側は両者の違いを意識しない。

---

# 12. PopulationEvidence 共通モデル

国・DBによって取得形式が異なっても、
criterion側へ渡す形式は共通化する。

```json
{
  "type": "population_frequency",
  "source": "jMorp",
  "source_version": "configured-version",
  "provider_version": "1.0.0",
  "population": "Japanese",
  "population_code": "JPN",
  "assembly": "GRCh38",
  "AC": 3,
  "AN": 100000,
  "AF": 0.00003,
  "FAF": null,
  "homozygote_count": 0,
  "quality_status": "PASS",
  "retrieved_at": "ISO-8601 timestamp"
}
```

最低限の共通項目:

```text
source
source_version
provider_version
population
population_code
assembly
AC
AN
AF
FAF
homozygote_count
quality_status
retrieved_at
```

取得できない項目は `null` とし、
Provider固有値をcriterionロジックへ直接漏らさない。

---

# 13. Provider Registry

Providerをコードへ固定せず、Registryで管理する。

```python
class ProviderRegistry:

    def register(self, category, name, provider):
        ...

    def get(self, category, name):
        ...

    def get_enabled(self, category):
        ...
```

例:

```python
registry.register(
    category="population",
    name="gnomad",
    provider=GnomADProvider(...)
)

registry.register(
    category="population",
    name="jmorp",
    provider=JMorpProvider(...)
)
```

将来は、

```python
registry.register(
    category="population",
    name="country_x",
    provider=CountryXPopulationProvider(...)
)
```

のように追加可能にする。

---

# 14. Provider設定

## 14.1 Japan profile

```yaml
profile: japan

population:
  target_population: JPN

  providers:
    - name: jmorp
      enabled: true
      priority: 100
      role: primary

    - name: gnomad
      enabled: true
      priority: 200
      role: supporting
```

## 14.2 Korea profile

```yaml
profile: korea

population:
  target_population: KOR

  providers:
    - name: korean_population
      enabled: true
      priority: 100
      role: primary

    - name: gnomad
      enabled: true
      priority: 200
      role: supporting
```

criterionコードは変更しない。

---

# 15. 「切り替え」ではなく「複数Provider利用」

Population DBは排他的切り替えだけでなく、
複数Providerを同時利用できる構造とする。

例:

```text
gnomAD
  AF = 0

jMorp
  AF = 0.00012
```

この場合、

```text
gnomADで absent
```

だけを理由にPM2を判定しない。

そのため、複数Evidenceを統合する
**Population Evidence Resolver** を設ける。

---

# 16. Population Evidence Resolver

```text
Population Provider[]
       ↓
PopulationEvidence[]
       ↓
PopulationEvidenceResolver
       ↓
ResolvedPopulationEvidence
       ↓
PM2 / BA1 / BS1 / BS2 / PP2
```

概念例:

```python
evidence_list = population_service.get_all(
    variant=variant,
    context=input_data
)

resolved = population_resolver.resolve(
    evidence_list=evidence_list,
    target_population="JPN"
)
```

Resolverでは以下を扱う。

```text
target population
primary / supporting source
provider priority
ancestry matching
最大AF
FAF利用可否
quality gate
複数DBの矛盾
missing data
```

---

# 17. ResolvedPopulationEvidence

例:

```json
{
  "target_population": "JPN",
  "primary": {
    "source": "jMorp",
    "AF": 0.00012
  },
  "supporting": [
    {
      "source": "gnomAD",
      "population": "global",
      "AF": 0.000001
    }
  ],
  "max_observed_af": 0.00012,
  "max_observed_af_source": "jMorp",
  "quality_status": "PASS",
  "resolution_policy": "population-resolution-v1"
}
```

---

# 18. Populationを利用する主なcriterion

```text
PM2
BA1
BS1
BS2
PP2
```

各criterionはProvider固有ロジックではなく、
`ResolvedPopulationEvidence` を利用する。

---

# 19. Population以外への拡張

Provider Plugin方式はPopulationだけに限定しない。

将来的には以下にも適用可能。

```text
variant classification
gene-disease
functional assay
computational prediction
phenotype
literature
```

例:

```text
VariantClassificationService
   ↓
Provider Registry
   ↓
├ ClinVarProvider
├ ClinGenEvidenceProvider
└ CountrySpecificVariantDBProvider
```

これにより、特定国・施設の独自knowledge baseを追加できる。

---

# 20. Logical Service

criterionスクリプトからは、
ProviderではなくLogical Serviceを利用する。

```text
services.population
services.annotation
services.variant_classification
services.gene_disease
services.functional
services.computational
services.phenotype
services.family
services.literature
```

例:

```python
population = services.population.get_resolved_evidence(...)
classification = services.variant_classification.find_comparators(...)
annotation = services.annotation.annotate(...)
```

---

# 21. 推奨ディレクトリ構成

```text
acmg_interpreter/
│
├── criteria/
│   ├── pvs1.py
│   ├── ps1.py
│   ├── ps2.py
│   ├── ...
│   └── bp7.py
│
├── core/
│   ├── models.py
│   ├── statuses.py
│   ├── result.py
│   ├── input_loader.py
│   ├── normalizer.py
│   ├── provenance.py
│   └── exceptions.py
│
├── evidence/
│   ├── population.py
│   ├── computational.py
│   ├── functional.py
│   ├── classification.py
│   ├── phenotype.py
│   └── family.py
│
├── services/
│   ├── population.py
│   ├── annotation.py
│   ├── variant_classification.py
│   ├── gene_disease.py
│   ├── functional.py
│   ├── computational.py
│   ├── phenotype.py
│   ├── family.py
│   └── literature.py
│
├── providers/
│   ├── registry.py
│   ├── population/
│   │   ├── base.py
│   │   ├── gnomad.py
│   │   ├── jmorp.py
│   │   ├── alfa.py
│   │   └── country_example.py
│   ├── annotation/
│   │   ├── base.py
│   │   └── vep.py
│   ├── classification/
│   │   ├── base.py
│   │   ├── clinvar.py
│   │   └── clingen_evidence.py
│   ├── functional/
│   │   ├── base.py
│   │   └── mavedb.py
│   ├── computational/
│   │   ├── base.py
│   │   ├── revel.py
│   │   └── spliceai.py
│   └── ontology/
│       ├── hpo.py
│       └── mondo.py
│
├── resolvers/
│   ├── population.py
│   ├── classification.py
│   └── evidence_conflict.py
│
├── config/
│   ├── profiles/
│   │   ├── default.yaml
│   │   ├── japan.yaml
│   │   └── korea.yaml
│   ├── thresholds.yaml
│   ├── sources.yaml
│   └── rules.yaml
│
├── schemas/
│   ├── input.schema.json
│   └── criterion_result.schema.json
│
├── va_spec/
│   ├── mapper.py
│   ├── builders.py
│   ├── validation.py
│   └── serializer.py
│
├── tests/
│   ├── criteria/
│   ├── providers/
│   ├── services/
│   ├── resolvers/
│   └── fixtures/
│
└── README.md
```

---

# 22. VA-SPEC-Pythonの位置付け

VA-SPEC-Pythonは、
**criterion判定ロジックのライブラリとしてではなく、
最終Evidenceオブジェクト生成・型validation・serializationのために使用する**。

```text
Criterion Script
       ↓
CriterionResult
       ↓
VA-Spec Mapper
       ↓
VA-SPEC-Python
       ↓
Variant Pathogenicity EvidenceLine
```

各criterionスクリプトから直接VA-SPEC-Pythonを呼ばない。

理由:

- ACMG rule engineとVA-Spec schemaを分離できる
- VA-SPEC-Pythonのversion変更を局所化できる
- criterion単体テストが容易
- 内部workflow statusをVA-Specへ無理に埋め込まなくてよい
- VA-Spec以外の出力形式を将来追加しやすい

---

# 23. VA-Spec Mapper

概念:

```python
def to_va_spec(result: CriterionResult):
    ...
```

内部で、

```text
criterion
status
direction
strength
evidence_outcome
evidence
provenance
```

をVA-SPEC-Python modelへ変換する。

例:

```text
CriterionResult

criterion = PM2
status = MET
direction = supports
strength = supporting
evidence_outcome = PM2_supporting
```

↓

```text
Variant Pathogenicity EvidenceLine

methodType = PM2
directionOfEvidenceProvided = supports
strengthOfEvidenceProvided = supporting
evidenceOutcome = PM2_supporting
```

---

# 24. VA-SPEC-Python version管理

以下を別々にProvenance管理する。

```text
ACMG/AMP guideline version
ClinGen guidance version
implementation rule version
VA-Spec specification version
VA-SPEC-Python package version
Provider version
Dataset version
```

---

# 25. Evidence Provenance

Providerから取得したEvidenceには最低限以下を保持する。

```text
source
source_type
source_version
provider_name
provider_version
dataset_version
population
population_code
reference_assembly
retrieval_date
variant representation
```

これにより、

> 同じPM2でも、どの国のどのPopulation DBを使って評価したか

を再現できる。

---

# 26. Criterion間依存・二重計上

criterionスクリプトは独立実行可能とし、
他criterionの結果を必須入力にしない。

ただし、以下は重複・競合候補となる。

```text
PVS1 ↔ PM4
PVS1 ↔ PP3(splicing)
PS1  ↔ PP3(splicing)
PM4  ↔ BP3
PP3  ↔ BP4
PP3  ↔ BP7
BP4  ↔ BP7
```

各criterionは必要に応じて `conflict_flags` を返す。

最終的なdouble-counting resolutionは後段で行う。

---

# 27. Threshold管理

固定閾値はコードへ直接埋め込まず、
configurationとして管理する。

```yaml
PP3:
  REVEL:
    supporting: 0.644
    moderate: 0.773
    strong: 0.932

BP4:
  REVEL:
    supporting: 0.290
    moderate: 0.183
    strong: 0.016
    very_strong: 0.003
```

固定閾値はACMG/ClinGen universal ruleと区別し、
`Implementation Rule` としてversion管理する。

---

# 28. キャッシュ

複数criterionが同じProviderへ問い合わせるため、
Evidence取得結果はcache可能にする。

```text
Variant
  ↓
Evidence Cache
  ├ population
  ├ annotation
  ├ classification
  ├ computational
  └ functional
```

cache key例:

```text
provider
provider_version
dataset_version
assembly
chrom
pos
ref
alt
```

---

# 29. テスト

## 29.1 Criterion test

```text
MET
NOT_MET
NOT_EVALUATED
NOT_APPLICABLE
MANUAL_REVIEW
```

をcriterion単位でテストする。

## 29.2 Provider contract test

新しいPopulation Providerは、
共通interfaceを満たすことをテストする。

```text
get_frequency() が PopulationEvidence を返す
assemblyが保持される
AC/AN/AFの型が正しい
source/versionが保持される
取得不可値はnullになる
```

## 29.3 Resolver test

```text
複数Provider
target population
provider priority
矛盾するAF
quality failure
missing data
```

をテストする。

## 29.4 VA-Spec test

CriterionResultをVA-SPEC-Pythonへ変換し、
model validationが成功することを確認する。

---

# 30. 初期実装優先順位

まず自動化しやすく、
Provider architectureも検証しやすい以下を実装する。

```text
PM2
BA1
PP3
BP4
BP7
PM4
BP3
```

特にPM2 / BA1で、

```text
gnomAD + jMorp
```

の複数Population Provider構成を最初に検証するとよい。

次段階:

```text
PVS1
PS1
PM1
PM5
PP2
BP1
```

その後:

```text
PS2
PS3
PS4
PM3
PM6
PP1
PP4
BS1
BS2
BS3
BS4
BP2
BP5
```

PP5 / BP6は `DEPRECATED` として扱う。

---

# 31. Criterion実装のDefinition of Done

各criterionについて以下を満たす。

- criterion専用Pythonファイルが存在する
- 共通inputを受け取れる
- 特定DB名へ直接依存しない
- Logical Service経由でEvidenceを取得する
- MET / NOT_MET / NOT_EVALUATED / NOT_APPLICABLE / MANUAL_REVIEWを区別できる
- strengthを返せる
- evidence_outcomeを返せる
- Evidenceを結果へ含める
- missing_inputsを返せる
- provenanceを返せる
- provider / dataset / rule versionを記録できる
- unit testが存在する
- CriterionResult schema validationが通る
- VA-Spec Mapperへ渡せる
- VA-SPEC-Python model validationが通る

---

# 32. 新規Population DB追加のDefinition of Done

各国・機関のPopulation DBを追加する場合は、
以下のみで追加できることを目標とする。

```text
1. PopulationProviderを実装
2. DB固有レスポンスをPopulationEvidenceへ変換
3. Provider Registryへ登録
4. profile YAMLへ追加
5. Provider contract testを追加
```

以下は変更不要とする。

```text
PM2.py
BA1.py
BS1.py
BS2.py
PP2.py
VA-Spec Mapper
```

これをProvider Plugin Architectureの主要要件とする。

---

# 33. 最終的な設計原則

本システムは、

```text
Criterion Logic
        ↓
Logical Service
        ↓
Provider / Plugin
        ↓
External Evidence Source
```

と、

```text
CriterionResult
        ↓
VA-Spec Mapper
        ↓
VA-SPEC-Python
        ↓
VA-Spec Evidence
```

を分離する。

この構造により、

- criterionごとの独立開発
- 国別Population DBの追加
- DB切り替え
- 複数DB併用
- provider priority変更
- ClinGen guidance更新
- VA-Spec更新
- VA-SPEC-Python更新

を、それぞれ独立して扱える。

特にPopulation evidenceについては、

> **gnomADを前提とする実装ではなく、Population Providerを差し替え・追加可能な実装**

とすることを基本設計とする。
