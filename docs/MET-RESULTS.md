# MET一覧と成立理由

[俯瞰](OVERVIEW.md)から分離した、MET 52件（demo-data 32・fixture 20）の内訳。
判定に使った数値・校正区間・除外理由は`results.json`の`provenance`と`evidence`に全件残る。
demo-data分は下記で再現でき、fixture分は`tests/test_clingen_positive.py`と
`tests/test_clingen_gates.py`が同じ評価を実行する。

```powershell
$env:PYTHONPATH = 'src'
.venv/Scripts/python.exe -m acmg prepare-demo-online --input-dir demo-data --cache-dir tests/fixtures/ensembl-cache --evidence-cache-dir tests/fixtures/external-cache --output-dir work/full/prepared --ensembl-release 116 --with-gnomad --with-clinvar --clinvar-release 2026-09-15 --with-dbnsfp --with-pm1-hotspot --rules config/demo-rules.json --offline
.venv/Scripts/python.exe -m acmg evaluate --input work/full/prepared/variants.json --evidence work/full/prepared/evidence.json --config config/demo-rules.json --context config/curated-context.json --criteria all --offline --output-dir work/full/evaluated
```

## demo-data（32件）

PM1（2件）: ClinVar missense density、route=mutational_hotspot、窓±5aa。閾値は
`config/demo-rules.json`の`PM1.hotspot`から入り、評価対象の変異自身は密度に数えない。

| record | gene | 変異 | 強度 | P/LP | benign |
| --- | --- | --- | --- | ---: | ---: |
| case1:26:1 | KCNJ5 | p.Arg155Gln | moderate | 6 | 0 |
| case3:24:1 | MYH7 | p.Arg719Trp | moderate | 9 | 0 |

PP3（6件）: REVEL（dbNSFP 4.8a、Pejaver et al. 2022校正）。6件とも蛋白機序で成立し、
splicing機序は区間に届かない。PP3はORなので最強の区間を採る。

| record | gene | 変異 | 強度 | REVEL | SpliceAI |
| --- | --- | --- | --- | ---: | ---: |
| case1:26:1 | KCNJ5 | p.Arg155Gln | strong | 0.934 | 0 |
| case3:24:1 | MYH7 | p.Arg719Trp | moderate | 0.814 | 0.02 |
| case1:30:1 / case4:30:1 | TNNI3 | p.Pro82Ser | supporting | 0.76 | 0 |
| case3:25:1 | MYBPC3 | p.Glu334Lys | supporting | 0.757 | 0.02 |
| case4:26:1 | DSG2 | p.Phe531Cys | supporting | 0.744 | 0 |

BP4（12件）: REVEL + SpliceAIのAND。synonymous 4件はdbNSFPがnonsynonymous SNV用で
REVELを持たないため、splicing機序のみで成立する。

| gene / 変異 | 件数 | 強度 | REVEL | SpliceAI |
| --- | ---: | --- | ---: | ---: |
| MYBPC3 p.Ser236Gly | 2 | moderate | 0.173 | 0.01 |
| MYBPC3 p.Arg382Trp | 4 | supporting | 0.265 | 0.07 |
| MYBPC3 p.Val158Met | 2 | supporting | 0.228 | 0.01 |
| LMNA p.Ala287= | 2 | supporting | - | 0 |
| TNNI3 p.Arg68= | 1 | supporting | - | 0.01 |
| TNNI3 p.Glu66= | 1 | supporting | - | 0.02 |

BA1（11件）: gnomAD 4.1.1。全件がClinGen SVI例外リスト9変異に非該当。AN不足の集団は
`rejected_observations`に`INSUFFICIENT_AN`として理由が残る。下2行はglobalが5%未満でも
特定集団で超える例。

| record | 変異 | exome:global AF | 最大集団 |
| --- | --- | ---: | --- |
| case1:29:1 / case3:27:1 | 11:47348490T>C | 0.1202 | nfe 0.1331 |
| case2:28:1 / case3:28:1 | 11:47350047C>T | 0.0895 | fin 0.1047 |
| case1:31:1 / case4:27:1 | 1:156135237T>C | 0.0814 | afr 0.4586 |
| case2:30:1 | 19:55156279C>A | 0.0501 | fin 0.0809 |
| case2:26:1 / case3:30:1 / case4:28:1 | 1:201361238G>C | 0.0088 | eas 0.1884 |
| case4:31:1 | 19:55156285C>T | 0.0031 | afr 0.1096 |

PM5（1件）: case3:24:1 MYH7 p.Arg719Trp、moderate。ClinVar residue検索で同一残基の
別アミノ酸変化にP/LP報告。`condition_assessment: NOT_EVALUATED`のため確認事項を併記する。

## fixture（20件）

clingen-positive 15件が14 criterionすべてを覆う。demo-dataで未評価のPVS1・PP2・BP1・BS1が
ここで成立するのは、fixtureがgene_disease・region・disease_frequency_thresholdを
最初から積んでいるため。demo側に足りないのはロジックではなくキュレーションであることを示す。

| criterion | 強度 | gene | 変異 | ERepo metCode | 成立理由 |
| --- | --- | --- | --- | --- | --- |
| PVS1 | very_strong | PAH | 12:102852850GA>G | PVS1 | relevant transcript/exonでNMD予測 |
| PS1 | strong | GCK | 7:44149809C>A | PS1 | exact検索で同一アミノ酸変化のP報告 |
| PM1 | moderate | PAX6 | 11:31802793C>T | PM1 | critical_functional_domain（curated） |
| PM2 | supporting | GUCY2D | 17:8009531T>C | PM2_Supporting | 全観測がrarity policyを満たす |
| PM4 | moderate | OTC | X:38369845GAAG>G | PM4 | in-frame deletion + reviewed region |
| PM5 | moderate | OTC | X:38367331C>T | PM5 | residue検索で同一残基の別変化にP報告 |
| PP2 | supporting | GCK | 7:44150049A>G | PP2 | gene_level機序 + variant spectrum |
| PP3 | strong | PAX6 | 11:31802793C>T | PP3_Moderate | REVEL 0.967 |
| BA1 | stand_alone | ITGB3 | 17:47283530T>C | BA1 | AF>5%、例外リスト非該当 |
| BS1 | strong | MYH7 | 14:23420189C>T | BS1 | max_credible_AF 0.0002（MONDO:0005045、AD） |
| BP1 | supporting | BRCA2 | 13:32333210G>C | BP1_Strong | gene_level機序 + variant spectrum |
| BP3 | supporting | FOXG1 | 14:28767509C>CGCCGCC | BP3 | 反復領域のin-frame挿入、機能的重要性なし |
| BP4 | moderate | SLC6A8 | X:153688650G>A | BP4 | REVEL 0.079 |
| BP7 | supporting | PAX6 | 11:31790720G>A | BP7 | 位置・splice予測・保存性が全て良性側 |
| PM2 | supporting | MYH7 | 14:23420189C>T | - | bs1-myh7で付随的に成立 |

PP3（strong / ERepo `PP3_Moderate`）とBP1（supporting / ERepo `BP1_Strong`）はVCEPが
criterion-specificに調整した強度との差で、CSpec自動適用を対象外としているため意図的に残す。
個々のUUIDと差の根拠は[正例一覧](CLINGEN-POSITIVE-REFERENCES.md)に固定している。

clingen-gate 2件（PP3 supporting: MYBPC3 c.405A>G、SpliceAI 0.27 / BP4 supporting:
MYH7 c.3036C>T、SpliceAI 0）はBP7ゲート確認用synonymousの副産物。REVELを持たないため
splicing機序のみで成立する。synthetic・population-met 3件（PM2 2・BA1 1）は
`policy_source: synthetic-test-only`の配線確認用で、臨床閾値ではない。
