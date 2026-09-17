# BH26 Demo Case 1-4:実症例概要・ClinGen登録状況・Ground Truth 統合版

**version: v6.0**(v5.0から、11節の重複・解消済み記述を整理。Case4臨床記述修正を10節の表に追加)

作成:トライネット高橋 / 2026-09-08(2026-09-10更新)
対象:`case1_variants.vcf`〜`case4_variants.vcf` + 対応する`case{1-4}_clinical_note.txt`(demo-cases-v3-real)

---

## 1. 信頼度階層の定義

比較ツールが参照する正解データは、単純に「これが正解」と隠すのではなく、**その正解自体がどれだけ確からしいか**を明示する。

| Tier | 定義 |
|---|---|
| **A** | ClinGen Variant Curation Expert Panel(VCEP)が基準ごとに正式検証済み(3-star) |
| **B** | ClinGen正式検証はないが、論文著者が根拠(ACMGコード等)を明記している |
| **C** | 実在する不一致(conflicting interpretations等)そのものが正解。単一の確定した分類はない |

## 2. ClinGen登録状況・全体サマリー

| Case | 変異 | ClinVar登録 | ClinGen VCEP登録 | Tier |
|---|---|---|---|---|
| 1 | MYBPC3 c.278delA(p.Lys93ArgfsTer3) | 未登録(新規) | ❌ | B |
| 1 | KCNJ5 c.464G>A(p.Arg155Gln) | 未登録(新規) | ❌(対応VCEP無し) | B |
| 2 | MYBPC3 c.2905+1G>A | ✅ VCV000042666 | ❌ | B |
| 2 | MYBPC3 c.836del(p.Gly279Valfs*21) | 未登録(新規) | ❌ | B |
| **3** | **MYH7 c.2155C>T(p.Arg719Trp)** | ✅ VCV000014104 | **✅ ClinGen Inherited Cardiomyopathy Expert Panel** | **A** |
| 3 | MYBPC3 c.1000G>A(p.Glu334Lys) | ✅(conflicting) | ❌ | C |
| 3 | MYH7 c.3382G>A(p.Ala1128Thr、ノイズ) | ✅(Likely benign) | ✅ ClinGen Inherited Cardiomyopathy Expert Panel(BS1) | **A** |
| 4 | DSG2 c.1592T>G(p.Phe531Cys) | ✅(conflicting) | ❌ | C |

**重要な発見**:ClinGen Cardiomyopathy VCEPは現状MYH7がPhase 1(完了)、MYBPC3・TNNT2・TNNI3・ACTC1等はPhase 2(適用可否を評価中)。**4症例中、正式なClinGen 3-starレビューを受けているのはCase 3のMYH7変異1つのみ**。

**なぜClinGenのみに絞らなかったか**:ClinGenがカバーする遺伝子・変異はごく一部であり、「1症例に複数の候補バリアント」という現実的な状況を、ClinGenレビュー済み変異だけで再現するのはほぼ不可能。また、ClinGenの証拠基盤自体が欧米中心データに依存しがちという偏りを抱えており、Pan-Asian実演の趣旨とも相性が悪い。**正解の確からしさに階層があること自体を明示する方が、藤原さんの論文が重視する「provenance and uncertaintyの保持」という原則により忠実**と判断した。

---

## 3. Case 1:MYBPC3 + KCNJ5(表現型ミスマッチ・パターン)

**出典**:Han S, Zhang Y-Y, Geng J (2026). Case Report: A novel MYBPC3 gene variant in a Chinese patient with hypertrophic cardiomyopathy and apical ventricular aneurysm, with a concurrent novel KCNJ5 gene variant. *Front Cardiovasc Med* 13:1841777. doi:10.3389/fcvm.2026.1841777. Open access, CC BY. **人間キュレーション(症例報告.pptx)により2026-09-08最終確認済み。**

**概要**:中国人HCM患者、心尖部瘤合併。新規のMYBPC3変異に加え、偶発的に新規KCNJ5変異(家族性高アルドステロン症III型・QT延長症候群13型の原因遺伝子、HCMとは無関係)も発見された。臨床医自身がFabry病(GLA)・心アミロイドーシスを実際の検査で除外した記述も含まれる。

**使いどころ**:表現型に一致しない偶発的所見(KCNJ5)に惑わされず、真の原因候補(MYBPC3)を正しく優先できるかを試すケース。**実際に`vep_parser.py`の課題(単独/新規報告のPathogenic変異の取りこぼし)がこの実データで再現される**(9節参照)。

### Ground Truth

| Variant | Criterion | Strength | Accept/Reject | Tier | 根拠 |
|---|---|---|---|---|---|
| case1-var1 MYBPC3 c.278delA | PVS1 | Very Strong | Accept | B | フレームシフト+早期終止コドン、NMD予測、MYBPC3はLOFが確立した病因(ClinGen SVI) |
| case1-var1 MYBPC3 c.278delA | PM2 | Moderate | Accept | B | 1000G/ESP6500/ExAC/地域コホートいずれにも未登録 |
| case1-var2 KCNJ5 c.464G>A | PM2 | Moderate | Accept | B | ExAC 8.236e-06のみ、1000G/ESP6500では未検出 |
| case1-var2 KCNJ5 c.464G>A | PP3 | Supporting | Accept | B | SIFT/PolyPhen2/MutationTaster全てD、REVEL=0.934 |
| case1-var2 KCNJ5 c.464G>A | PP4 | — | **Reject** | B | KCNJ5はHCM表現型と無関係。息子はKCNJ5のみ保有し心表現型なしと論文が報告 |
| case1-var3 MYH7 p.Ser1491Cys(ノイズ) | BS1, BP4 | Strong, Supporting | Accept | B | ClinVar VCV000043020、複数施設一致のBenign |

**ACMGガイドライン表(Table 5)とTavtigian点数法の判定差**:MYBPC3(PVS1_VeryStrong+PM2_Moderate)は、ACMGガイドラインのTable 5では**Likely Pathogenic**(Very Strong 1つ+Moderate 1つ)だが、Tavtigian点数法では**+10点でPathogenicの閾値に達する**。両手法で結論が変わりうる実例。KCNJ5(PM2_Moderate+PP3_Supporting)はTable 5のいずれの区分にも該当せず**VUS**、Tavtigian法でも+3点で同じくuncertain。

**最終Classification**:case1-var1(MYBPC3)= **Likely Pathogenic**。case1-var2(KCNJ5)= **VUS**、表現型不一致によりPP4はReject。**原因バリアント = case1-var1(MYBPC3)**

---

## 4. Case 2:MYBPC3 複合ヘテロ接合(二重正解パターン)

**出典**:Wang J, Hong L, Li Y, Mao Z, Zhu Y, Qi M, Zhou R, Hong X (2025). Case Report: Lethal neonatal hypertrophic cardiomyopathy from compound heterozygous MYBPC3 variants. *Front Cardiovasc Med*. doi:10.3389/fcvm.2025.1726463. CC BY. **人間キュレーション(症例報告.pptx)により2026-09-08最終確認済み。**

**概要**:新生児期致死性HCMの分子剖検症例。父方由来のスプライス部位変異(c.2905+1G>A、exon27、既知、ClinVar VCV000042666)と、母方由来の新規フレームシフト変異(c.836del、exon8、新規)——両方ともロスオブファンクション、複合ヘテロ接合。trio-WESで確認、ウエスタンブロットでタンパク確認失敗(nullフェノタイプ)。

**使いどころ**:「1つの正解に絞り込む」型ではなく、**「2つの変異が両方とも原因として必要」**という、Case 1・3とは異なるパターン。

### Ground Truth

| Variant | Criterion | Strength | Accept/Reject | Tier | 根拠 |
|---|---|---|---|---|---|
| case2-var1 MYBPC3 c.2905+1G>A(父方、exon27) | PVS1 | Very Strong | Accept | B | 標準スプライスドナー部位破壊、LOF確立病因 |
| case2-var1 | PS4 | Strong | Accept | B | ClinVar(VCV000042666)がケースコントロール比較論文(PMID:27532257)を引用。**7つの臨床検査機関が独立に評価し全て一致(Pathogenic/Likely Pathogenic)** |
| case2-var1 | PM2 | Moderate | Accept | B | gnomAD 0.0014%、UK Biobank保因者2人、1000genome未検出 |
| case2-var2 MYBPC3 c.836del(母方、exon8) | PVS1 | Very Strong | Accept | B | 早期終止コドン、NMD予測、ウエスタンブロットでタンパク確認失敗により機能喪失を直接裏付け |
| case2-var2 | PM2 | Moderate | Accept | B | gnomAD/UK Biobank/1000genomeいずれも未検出、ClinVar/dbSNP未登録 |
| case2-var3 TNNT2(ノイズ) | BS1, BP4 | Strong, Supporting | Accept | B | ClinVar VCV000188689 |

**ACMGガイドライン表とTavtigian点数法**:case2-var1(PVS1_VeryStrong+PS4_Strong+PM2_Moderate)=Table 5で**Pathogenic**、Tavtigian法+14点でもPathogenic(両手法一致)。case2-var2(PVS1_VeryStrong+PM2_Moderate)=Table 5で**Likely Pathogenic**、Tavtigian法+10点では**Pathogenicの閾値**(Case1のMYBPC3と同型の手法間差異)。

**最終Classification**:case2-var1 = **Pathogenic**、case2-var2 = **Likely Pathogenic**(単純に両方「Pathogenic」ではない)。**原因バリアント = case2-var1とcase2-var2の両方**(単一の勝者を選ぶ設計ではない)

---

## 5. Case 3:MYH7 + MYBPC3 p.Glu334Lys(Pan-Asian頻度パターン)★ClinGen正式レビュー済み変異を含む

**出典**:MYH7=ClinVar VCV000014104、ClinGen Inherited Cardiomyopathy Expert Panel(2016-12-15、3-star)。MYBPC3=The Need for Inclusive Genomic Research, *Circulation: Genomic and Precision Medicine*, DOI:10.1161/CIRCGEN.122.003736(Tomar et al.のSG10K解析を引用)、裏付け(香港コホート):PMC12123433。**患者の臨床記述は創作**(変異の同定情報は実在)。

**概要**:出典はHong Kong Chinese HCMコホート研究(53症例、P/LP 13変異・VUS 21変異)。MYH7 p.Arg719Trp(ClinGen承認済みPathogenic)を正解として含みつつ、以前Pathogenicと分類されていたMYBPC3 p.Glu334Lysが、東アジア集団データにより現在ClinVar上VUS(conflicting)という実例を組み合わせた。**人間キュレーション(症例報告_1.pptx)により2026-09-08最終確認済み。**

**使いどころ**:Pan-Asian核心変異。ClinGen承認済みの確立した正解(MYH7)と、まだ標準化されていない候補(MYBPC3)が同じ症例に混在する点で、①の自動化レベルの区分の教材としても使える。

### Ground Truth

| Variant | Criterion | Strength | Accept/Reject | Tier | 根拠 |
|---|---|---|---|---|---|
| case3-var1 MYH7 c.2155C>T | PS2 | Strong | Accept | **A** | ClinGen Evidence Repository。トリオ解析でde novo確認(PMID:10957787) |
| case3-var1 | PS3 | Strong | Accept | **A(確信度に留保あり)** | 血行動態低下による有害な心筋リモデリング(PMID:24829265)。**人間キュレーター自身が「ノックイン・ノックアウトマウスの知識があまりなく、確からしさの判定に自信がない」と明記** |
| case3-var1 | PS4 | Strong | Accept | **A** | 中国人非血縁6家系で罹患者8人に検出、非罹患者では非保因者のみ(PMID:19645038) |
| case3-var1 | PM1, PM2, PM5, PP1, PP3 | Moderate/Moderate/Moderate/Supporting/Supporting | Accept | **A** | ClinGen Evidence Repositoryが正式適用。PM5は同一残基の別アミノ酸変化(p.Arg719Gln、Variation ID 14107)がExpert PanelでPathogenic |
| case3-var2 MYBPC3 c.1000G>A | PS3 | Strong | Accept | **C(確信度に留保あり)** | E334K cMyBPCによるUPS障害がイオンチャネル・Ca2+ハンドリングタンパクレベルを変化。**同じくノックイン・ノックアウトマウス知見への自信のなさが明記されている** |
| case3-var2 | BS1 | — | Accept(地域データ使用時) | **C** | gnomAD全体AF=0.02368%、East Asian AF=0.3338%(約14倍)。患者が中国人であることを考慮しBenign方向 |
| case3-var2 | (総合) | — | — | **C** | PS3(病原性寄り)とBS1(良性寄り)が競合しVUSと判定。**gnomAD全体のAFのみ知られていた時期はPathogenicと判定されていた**——地域データの追加が判定を変えた核心的な実例 |
| case3-noise群(1-5) | 各種 | — | Reject | B | 実在の高頻度Benign変異 |
| case3-noise6 MYH7 c.3382G>A(p.Ala1128Thr) | BS1 | Strong | Reject | **A** | ClinGen Inherited Cardiomyopathy Expert Panel Evidence Repository(Kelly et al. 2018)。gnomAD Latino 0.02%・East Asian 0.016%でBS1適用。**両集団で同程度に高頻度であり、3-4節Bの逆方向パターンの実例ではない** |

**訂正(2026-09-08、修正版症例報告v2で解消)**:case3-var2のClinVar内訳(VCV000177902.53、Uncertain significance 9件・Benign 2件・Likely benign 3件)が、当初Case4のDSG2と完全に一致しており、コピー&ペーストミスが疑われていたが、**キュレーターが修正版で確認し、Case4側の内訳は実際には異なる(Likely pathogenic 3件・Uncertain significance 7件・Likely benign 1件)ことが判明した**。疑いは解消された。

**最終Classification**:case3-var1(MYH7)= **Pathogenic**(Tier A、PS3を除いてもPS2+PS4のStrong 2個で成立するため確信度は揺るがない)。case3-var2(MYBPC3)= **VUS**(Tier C、PS3の確信度に留保がある分、実質的にBS1寄りに傾く可能性もある)。**原因バリアント = case3-var1(MYH7)**

---

## 6. Case 4:DSG2 p.Phe531Cys(日本/ToMMo、Pan-Asian頻度パターン)

**出典**:Arrhythmogenic right ventricular cardiomyopathy in a Japanese patient with a homozygous founder variant of DSG2 in the East Asian population. *Human Genome Variation*. PMC9360431. **人間キュレーション(症例報告_1.pptx)により2026-09-08最終確認済み。**

**概要**:実際の日本人ARVC患者。ホモ接合のDSG2 p.Phe531Cys変異。gnomAD東アジア集団・ToMMo・HGVDのいずれでも、gnomAD全体より大幅に高頻度(最大約14倍)。東アジア系founder variant。

**使いどころ**:日本側のPan-Asianケース。TogoMCP経由でToMMo/NCBN/GEM-Jへのライブ照会と直接接続できる。

### Ground Truth

| Variant | Criterion | Strength | Accept/Reject | Tier | 根拠 |
|---|---|---|---|---|---|
| case4-var1 DSG2 c.1592T>G(ホモ接合) | PS3 | Strong | Accept(キュレーター独自の再導出) | **C(確信度に留保あり)** | マウスDsg2 p.Phe536Cysノックインモデル(ヒトp.Phe531Cysに対応)、ホモ接合で両心室拡大・LVEF低下・伝導異常・線維化を確認。**人間キュレーター自身が「ノックイン・ノックアウトマウスの知識があまりなく、確からしさの判定に自信がない」と明記** |
| case4-var1 | PS4 | Strong | Accept | **C(キュレーター独自の再導出)** | 家系スクリーニングでホモ接合保因者のみ完全浸透率でARVC表現型、ヘテロ接合保因者は非罹患 |
| case4-var1 | (総合) | — | — | **C** | キュレーターの再導出(PS3+PS4のStrong 2個)ではPathogenicとなる。**実際のClinVar登録(VCV000044283.23)はLikely pathogenic(3)/Uncertain significance(7)/Likely benign(1)**——訂正前に疑っていたほどの乖離はなく、11件中3件は既にLikely pathogenicと判定しているが、多数派(7件)はUncertainに留まる |
| case4-noise群 | 各種 | — | Reject | B | 実在のBenign変異(一部TogoMCPライブ照会値) |

**最終Classification**:case4-var1(DSG2)= **単一の正解なし、ただしPathogenic寄りに傾きつつある**(Tier C)。**キュレーター自身の再導出(PS3+PS4→Pathogenic)は、実際のClinVar登録の一部(11件中3件がLikely pathogenic)とは方向性が一致するが、多数派(7件)は依然Uncertainに留まる。PS3の確信度が低い(マウスモデル解釈への留保あり)ことを踏まえると、この変異は「Pathogenicへ向かいつつあるが未確定」という記述が最も正確。原因バリアント = case4-var1(DSG2)、ただし確信度は限定的である旨を明示して報告する**

---

## 7. Tier分布まとめ

| Case | Tier A | Tier B | Tier C |
|---|---|---|---|
| 1 | — | var1(MYBPC3)、var3(ノイズ) | — |
| 2 | — | var1・var2(MYBPC3×2)、var3(ノイズ) | — |
| 3 | var1(MYH7)、noise6(MYH7 BS1) | ノイズ群 | var2(MYBPC3) |
| 4 | — | ノイズ群 | var1(DSG2) |

**Tier Aに該当するのはCase 3のMYH7関連2変異のみ**。比較ツールでの結果報告時は、単純な一致率だけでなく、**どのTierの正解との一致/不一致だったか**を併記する設計とする(企画書3-1節参照)。

**確信度のさらなる内訳(重要)**:Tier A・Cの中にも、**PS3(機能実験、特にノックイン・ノックアウトマウスモデルによる評価)については、人間キュレーター自身が「確からしさの判定に自信がない」と明記した箇所が複数ある**(Case3のMYH7・MYBPC3両方、Case4のDSG2)。Tier区分だけでは表現しきれないため、該当するEvidenceには個別に「確信度に留保あり」と注記している。Case4は特に、PS3の確信度が低いことに加え、キュレーター独自の再導出(PS3+PS4→Pathogenic)と、実際のClinVar登録の合意(Uncertain寄り)が一致しないという、二重の不確実性を持つ。

---

## 8. デモデータの概要(ファイル構成・トリアージ検証結果)

| ファイル | 内容 |
|---|---|
| `case{1-4}_clinical_note.txt` | 各症例の臨床記述(英語)。#1・#2・#4は実際の症例報告から採録・要約、#3は変異は実在だが患者の物語は創作 |
| `case{1-4}_variants.vcf` | 各症例のバリアントリスト。標的変異に加え、TogoMCPライブ照会・ClinVarから収集した実在の良性ノイズを追加し、①トリアージ層(`vep_parser.py`)のテストにも使える規模(6〜7変異/症例)に拡張 |

**トリアージ実測結果**(`vep_parser.py`、mainブランチ未修正版):

| Case | 変異数(標的+ノイズ) | 実測結果 |
|---|---|---|
| 1 | 2標的+5ノイズ=7 | 人間キュレーション(症例報告.pptx)で確定した正しいデータ(MYBPC3=PVS1+PM2/Likely Pathogenic、KCNJ5=PM2+PP3/VUS)で最終検証。**vep_parser.pyの課題が実際に発現**:KCNJ5(VUS、Rank5)が真の原因候補MYBPC3(Likely Pathogenic、新規のためreview status未確立、Rank100)より上位に来る |
| 2 | 2標的+5ノイズ=7 | 人間キュレーションで確定したCLNREVSTATを反映して最終検証。**同じ複合ヘテロ接合ペアの中でも差が出た**:case2-var1(既知、7施設一致)はRank2で正しく最上位、case2-var2(新規、review status未確立)は課題によりRank100に留まる。ただし両方ともノイズ(200)よりは正しく上位 |
| 3 | 2標的+6ノイズ=8 | MYH7(正解)が正しく最上位(Rank2)、MYBPC3(Pan-Asian候補)も正しくノイズより上位。追加したMYH7 BS1変異(Tier A)も正しくノイズ最下位(Rank200) |
| 4 | 1標的+5ノイズ=6 | DSG2が正しくノイズより上位(Rank9 vs 200) |

**課題への対応の検証について(最終確定、2026-09-08)**:本項目は当日中に複数回の訂正を経ている。①当初、人間キュレーション入手前の推測でMYBPC3をPVS1+PM2/Likely Pathogenicと記録し、課題の再現を確認→②検索スニペットの誤読(KCNJ5の記述をMYBPC3のものと取り違え)により、誤ってPM2+PP3/VUSに訂正し、「課題は発現しない」と誤って報告→③人間キュレーション(症例報告.pptx、論文全文を読んだ結果)により、①の内容が正しかったことが確定的に確認された。**最終結論:MYBPC3=Likely Pathogenic(PVS1+PM2)、KCNJ5=VUS(PM2+PP3)という正しいデータで検証すると、vep_parser.pyの課題は実際に発現する。** Case 2でも、同一の複合ヘテロ接合ペア内で新規変異側だけがこの課題の影響を受けることを確認。課題への対応版(`vep_parser_bugfixed.py`)は当日Francisさんへの提案候補として参考資料に保管。

---

## 9. AlphaMissense・AlphaGenomeアノテーション(2026-09-10追加)

ユーザーが独自にAlphaMissense(ミスセンス病原性予測)とAlphaGenome(スプライシング・調節領域への影響予測、raw score・PHRED score)を全変異にアノテーションした(`annotation_alphamissense_alphagenome_v1.xlsx`)。この過程で、**4症例・28変異全てのゲノム座標(POS/REF/ALT)が実際に検証され、うち少なくとも5変異(case2-var1, case3-var1, case3-var2, case4-var1, case3-noise6)はClinVar記載と異なっていたことが判明し修正された**。**以前私が「デモ用の近似値、厳密な再マッピングをしていない」と正直に注記していた懸念が、まさに的中していた形になる。この座標検証は結果として全変異に及んでおり、11節の該当する未着手事項は解消されたとみなす。**

修正後、全ての変異でアノテーションが正常に取得できた。取得できなかった変異(case1-var1のMYBPC3欠失、case2-var2のMYBPC3欠失)は、フレームシフト/欠失であるため、ミスセンス変異を対象とするAlphaMissenseの対象外という、想定通りの理由だった。

**注目すべき発見**:case1-var2(KCNJ5 c.464G>A、表現型ミスマッチの罠として設計)について、AlphaMissenseが**am_pathogenicity=0.9485(likely_pathogenic)**という高い病原性スコアを独立に算出した。ACMG/ClinVarベースの評価ではPM2+PP3でVUSに留まるが、AlphaMissense単体で見ると「病原性が高そうな変異」に見える。**これはこの症例の教育的価値をむしろ強化する**:「変異自体の病原性スコアが高い」ことと、「その患者の表現型(HCM)の説明として正しい遺伝子か」は別の問題であり、KCNJ5はいくらAlphaMissenseスコアが高くても、HCMとは無関係な遺伝子であるため、原因候補として採用すべきではない。この対比は、決定論 vs LLM比較実証(企画書3-4節)やEvidence Copilotの設計思想(表現型との整合性を独立に評価する必要性)を補強する具体例になる。

再検証(`vep_parser.py`)の結果、座標修正はRank出力に影響を与えないことを確認した(vep_parser.pyはCLNSIG/CLNREVSTAT/頻度/IMPACTのみを使い、ゲノム座標自体は判定に使わないため)。

---

## 10. デモデータ作成メモ

本ドキュメント作成の過程で、**ACMG評価(Criterion/Strength/Classification)の数値だけでなく、症例の臨床記述(clinical note)にも同水準の検証が必要**であることが分かった。具体的に見つかった誤りは以下の通り。

| 箇所 | 誤りの内容 | 原因 | 発見のきっかけ |
|---|---|---|---|
| Case1 MYBPC3のACMG評価 | PVS1+PM2/Likely PathogenicをPM2+PP3/VUSに誤訂正 | 検索スニペットの誤読(KCNJ5に関する一文をMYBPC3のものと取り違え) | 人間キュレーション(症例報告.pptx)との突合 |
| Case2 MYBPC3(母方)のACMG評価 | PM3・PM6を誤って付与 | 論文に明記のない基準を「標準的拡張」として独自に補完 | 人間キュレーションとの突合 |
| Case2 MYBPC3(父方)のACMG評価 | PS4の見落とし | 検索スニペットの読み込み不足 | 人間キュレーションとの突合 |
| **Case1 clinical note** | **患者の年齢・臨床経過を誤って創作(34歳・労作時呼吸困難6か月→実際は63歳・動悸7日間)** | 検証せずに"それらしい"臨床像を補完 | 論文原文の直接検索 |
| **Case2 clinical note** | **全く別の論文(Alsters et al. 2019、オランダの別症例)の内容を混入**(母体糖尿病既往、家族歴陰性) | 複数の類似論文を検索する過程での取り違え | 論文原文の直接検索 |
| Case3 MYH7・MYBPC3、Case4 DSG2のACMG評価 | PS-series(PS2/PS3/PS4/PM5等)の詳細な内訳が不足していた(自分の検索だけではClinVarのPathogenic/Conflicting表示までしか分からなかった) | 検索スニペットではACMG評価の詳細な内訳(どの基準がなぜ適用されたか)まで拾いきれなかった | 人間キュレーション(症例報告_1.pptx)との突合 |
| Case3 MYBPC3・Case4 DSG2のgnomAD数値 | 概算値(illustrative)を使っていたが、実際の精密な数値(gnomAD global 0.02368%、East Asian 0.3338%等)が判明 | 検索で正確な数値まで確認できていなかった | 人間キュレーションとの突合 |
| **Case4 clinical note** | **患者の年齢・性別(58歳男性)と、両親の血族結婚(ホモ接合の理由そのもの)が完全に欠落** | 検索スニペットからは要約レベルの情報しか拾えておらず、フルテキストでの裏取りをしていなかった | 論文原文(フルテキスト)の直接検索 |

**追加の教訓(Case3・4、2026-09-08)**:
5. **機能実験(PS3)の評価、特にノックイン・ノックアウトマウスモデルを用いた研究の妥当性判断は、専門知識がないと自信を持って行えない**。人間キュレーター自身が「ノックイン・ノックアウトマウスの知識があまりなく、確からしさの判定に自信がない」とCase3(MYH7・MYBPC3両方)・Case4(DSG2)で明記しており、これはTierだけでは表現しきれない「基準ごとの確信度のばらつき」を示している。Evidence Recordには、Tier(出所の格)とは別に、**特定の基準タイプ(機能実験など)自体の評価の難しさ**も記録できる余地を残すべきかもしれない
6. **同じ内訳の数字が複数の異なる変異に登場した場合は、コピー&ペーストミスを疑う**。Case3のMYBPC3とCase4のDSG2で、当初ClinVarの提出者内訳が完全に一致しており、コピー&ペーストミスを疑っていた。**キュレーターが修正版(症例報告_v2.pptx)で確認し、実際には異なる内訳(DSG2:Likely pathogenic 3件・Uncertain significance 7件・Likely benign 1件)であることが判明し、疑いは解消された。** 疑わしい一致に気づいた時点で確認を依頼したことが、正しいデータへの修正につながった
7. **人間キュレーターの独自の再導出(ACMG基準からの計算結果)と、実際にデータベースに登録されている専門家集団の合意は、一致するとは限らない**(Case4で顕在化)。どちらか一方を「正解」として採用するのではなく、両方を記録し、その乖離自体を有用な情報として扱う

**教訓**:
1. **ACMG評価の数値は「際立って重要」に見えるため慎重に扱うが、物語的な臨床記述は"それらしければ良い"と気が緩みやすい**。しかし正解データとして使う以上、両者は同じ水準の検証が必要
2. **検索スニペットは、それが「どの変異・どの論文についての記述か」を、周辺の文脈まで確認してから採用する**。断片的な一致だけで判断すると、隣接する別の対象の記述を取り違えるリスクがある
3. **同じ疾患・同じ遺伝子(MYBPC3)を扱う論文は複数存在するため、複数の類似症例が検索結果に混在する際は、著者名・出版年・DOIを都度照合し、混同を防ぐ**
4. 人間による査読(今回は症例報告.pptx)との突合が、これらの誤りを発見する上で決定的に有効だった。**Day1でも、自動生成された内容を鵜呑みにせず、Ruthさん・Francisさんのような実務者による目視確認の機会を設ける価値がある**

---

## 11. 未着手・今後の課題

- ノイズ変異(良性バリアント)についてはClinGen登録状況を個別確認していない(標的変異のみ確認済み)
- Demo Case全体の疾患軸(HCMとARVCが混在)をどう扱うかは未決定
- 比較ツールの出力仕様(Tier別の一致率表示等)は設計のみで未実装
