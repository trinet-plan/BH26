# ClinGen正例リファレンス

更新: 2026-09-16。

## 目的と境界

PP5/BP6を除く実装対象14 criterionについて、ClinGen Evidence Repository（ERepo）の公開VCEP解釈に
最低1件の適用例があり、現在のGeneral Guidance evaluatorが独立EvidenceからMETを再現できることを
固定する。これは最終分類の再現でも、VCEP/CSpec固有ロジックの実装でもない。

ERepoの`metCodes`は`clingen-positive-expected.json`にだけ保存する。evaluatorが読む
`clingen-positive-prepared.json`と`clingen-positive-evidence.json`には期待ラベルを入れない。
したがって「ClinGenが適用したからMET」という循環評価にはならない。

## 対応表

| Criterion | ClinGen variant | ERepo UUID | ERepo code | Engine strength |
| --- | --- | --- | --- | --- |
| PVS1 | PAH c.806delT | `f5d8dc3f-dba0-4dc9-98cd-4a22b6a11a83` | PVS1 | very_strong |
| PS1 | GCK c.630G>T (p.Met210Ile) | `d17187f9-d65a-4117-af07-cfcb9e907236` | PS1 | strong |
| PM1 | PAX6 c.52G>A (p.Gly18Arg) | `94abbf17-09a7-47b8-8d26-63bc6f30280f` | PM1 | moderate |
| PM2 | GUCY2D c.1694T>C (p.Phe565Ser) | `e4081a29-cfa7-40d5-b266-5ca0950eaed2` | PM2_Supporting | supporting |
| PM4 | OTC p.Arg89del | `1013949b-d83b-492b-adf2-b808a41e801c` | PM4 | moderate |
| PM5 | OTC p.Arg40Cys | `07236833-88e3-4b3a-9a4d-35a6c2230597` | PM5 | moderate |
| PP2 | GCK p.Trp167Arg | `8baf769d-4142-48b1-9144-f9044b76b177` | PP2 | supporting |
| PP3 | PAX6 c.52G>A (p.Gly18Arg) | `94abbf17-09a7-47b8-8d26-63bc6f30280f` | PP3_Moderate | strong |
| BA1 | ITGB3 c.342T>C (p.Ile114=) | `9b9e1625-5f63-4b04-aef5-9fa3f6098f24` | BA1 | stand_alone |
| BS1 | MYH7 c.3382G>A (p.Ala1128Thr) | `dc6a24be-6cd7-45ac-a097-1d5031608aa7` | BS1 | strong |
| BP1 | BRCA2 p.Gly578Arg | `60281fd2-8f05-4c28-abf1-89a1a3eded05` | BP1_Strong | supporting |
| BP3 | FOXG1 p.Pro79_Pro80dup | `6e2dafeb-f2f0-4001-98b4-3b72f6ccff6c` | BP3 | supporting |
| BP4 | SLC6A8 p.Gly26Arg | `49ff12e7-0969-4d36-a6ca-664d4a65e974` | BP4 | moderate |
| BP7 | PAX6 c.1215C>T (p.Thr405=) | `f6819b13-83d-4f10-864c-e25c671697bc` | BP7 | supporting |

PP3とBP1ではERepoのVCEP固有強度とengine強度が異なる。engineは`clingen-positive-rules.json`の
General Guidance校正を一貫適用し、ERepo codeを強度入力として使わないためである。

## データ監査

- ERepo summary API 2.5.6を2026-09-16に取得し、13,265件を走査した。
- MET件数はPVS1 2,894、PS1 197、PM1 1,438、PM2 9,308、PM4 245、PM5 1,386、
  PP2 757、PP3 4,346、BA1 1,138、BS1 747、BP1 52、BP3 23、BP4 2,370、BP7 1,103。
- GRCh38座標はERepoのHGVSに合わせた。FOXG1とOTC indelのVCF anchor/REFはEnsembl参照配列で
  確認した。
- PM2はGUCY2Dの1/1,613,704とVCEP閾値0.0004、BA1はITGB3のAfrican
  1,342/24,024（約5.586%）、BS1はMYH7の10/34,232とVCEP閾値0.0002を使用する。
  未公開の分母や便宜的な頻度は作らない。

## 再現

```powershell
$env:PYTHONPATH = 'src'
.venv/Scripts/python.exe -m unittest tests.test_clingen_positive -v
```

テストは14 criterionのstatus/strength、期待ラベルの入力分離、ERepo provenance、PVS1 trace、
VA-Spec 1.0.1 Evidence Lineの生成とschema検証を確認する。
