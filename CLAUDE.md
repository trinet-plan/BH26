# AGENTS.md

このファイルはAIコーディングエージェント向けの作業ガイドです。人間向けのセットアップ手順・アーキテクチャ概要は [`README.md`](README.md) を、設計判断の経緯は [`ps3_bs3_ps4_implementation_v10.md`](ps3_bs3_ps4_implementation_v10.md) を参照してください。

## プロジェクト概要

- LLM(vLLM上のgemma-4)とPubMed MCPを組み合わせ、ACMG/AMP 2015の28基準について、人間キュレーターの一次スクリーニングを高速化する下書き判定を生成するパイプライン
- 最終的にはdocker化してホストし、APIとして公開する予定

## 初期段階の構想

| 層 | 内容 | 対応するACMGコード | 実装 |
|---|---|---|---|
| 1. Automated Evidence | 構造化データから機械的に判定 | PVS1, PS1, PM1, PM2, PM4, PM5, PP2, PP3, PP5, BA1, BS1, BP1, BP3, BP4, BP6, BP7 | InterVar/AutoPVS1等既存の決定的ツール |
| 2. Semi-automated Evidence | 表現型を使う判定 | PP4 | PubCaseFinder(入力は自然文の臨床記述) |
| 3. Evidence Gap / Manual Review Required | 現行データでは判断できない点を明示 | PS2, PS3, PS4, PM3, PM6, PP1, BS2, BS3, BS4, BP2, BP5 | ルール+LLM(Evidence Gap Analysis) |
| 4. Documentation | Evidence・判断結果の自然文整理 | — | ローカルLLM |