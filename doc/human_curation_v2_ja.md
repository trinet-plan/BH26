# 人間キュレーション結果:Case 1〜4 ACMG評価

**version: v2**(元pptx: human_curation_v2.pptx を整形したもの)

---

## Case 1: A novel MYBPC3 gene variant in a Chinese patient with hypertrophic cardiomyopathy and apical ventricular aneurysm, with a concurrent novel KCNJ5 gene variant

症例報告:肥厚性心筋症および頂端心室動脈瘤を患う中国人患者における新規MYBPC3遺伝子変異と、同時に新規のKCNJ5遺伝子変異

### 概要

中国人HCM患者、心尖部瘤合併。新規のMYBPC3変異に加え、偶発的に新規KCNJ5変異(家族性高アルドステロン症III型・QT延長症候群13型の原因遺伝子、HCMとは無関係)も発見。

表現型に一致しない偶発的所見(KCNJ5)に惑わされず、真の原因候補(MYBPC3)を正しく判断できるか。

### 症状

- 肥大型心筋症(HCM)を発症
- MYBPC3遺伝子とMYH7遺伝子はHCM患者で最も一般的な原因遺伝子
- HCM患者では、冠動脈疾患がない場合に左心室動脈瘤が報告された症例は5%未満
- KCNJ5遺伝子の変異は、家族性高アルドステロン症タイプIIIや長QT症候群タイプ13を引き起こすことが多い
- 主に、他の心臓、全身、代謝的疾患が存在しない状態で、左心室肥大(LVH)が特徴的
- MYH7とMYBPC3が最も関与しており、変異陽性症例の大多数を占める
- 他の遺伝子(TNNI3、TNNT2、TPM1、MYL2、MYL3、ACTC1など)は、それぞれ患者の1〜5%に寄与

### 検知された2つの変異

- MYBPC3遺伝子エクソン2に新規ヘテロ接合バリアント NM_000256.3(MYBPC3):c.278delA (p.Lys93ArgfsTer3)
- KCNJ5遺伝子エクソン2に新規ヘテロ接合バリアント NM_000256.3(KCNJ5):c.464G>A (p.Arg155Gln)

### MYBPC3 c.278delA (p.Lys93ArgfsTer3) とは

- 集団頻度データベースの調査により、この変異は稀であり、1,000ゲノムプロジェクト、ESP6500、ExACの各データベースでは報告されていない
- 心筋症患者と対照群の両方を含む地域集団データベースにも含まれていなかった
- この単塩基欠失により、位置93の正に帯電した極性残基リジンが置換され、翻訳中にフレームシフトが生じる。2番目の次のアミノ酸位置に早期停止コドンが導入され、NMD(Nonsense-Mediated mRNA Decay)が予測される
- 同じ遺伝子座における下流のフレームシフトまたはナンセンス変異はClinVarで一貫して報告されており、肥厚性心筋症および関連する心臓疾患の病原性と分類されている
- ACMGガイドラインによれば、この変異はHCMに対して「可能性が高い病原性」と分類され、PVS1・PM2の基準を満たす

**ACMG評価:**

| 基準 | 根拠 | ACMG評価 | Tavtigian方式 |
|---|---|---|---|
| PVS1(機能喪失型変異の定義) | 1塩基欠失によるフレームシフト+早期終止コドン→NMD。ClinGen MYBPC3基準では機能喪失が確立した病因のためPVS1適用可 | Very Strong | +8点 |
| PM2(集団対照群データベースに認められない) | 1,000ゲノム・ESP6500・ExACいずれも未報告 | Moderate | +2点 |

ACMGガイドラインtable5ではVery Strong 1個+Moderate 1個でLikely Pathogenic。Tavtigian方式では+10点でPathogenic。(ClinGenはPM2をModerateではなくSupportingにすることを推奨)

参照:ACMGガイドライン https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/ / ClinGen MYBPC3基準 https://zenodo.org/records/21434332

### KCNJ5 c.464G>A (p.Arg155Gln) とは

- 集団頻度データベースの検索により稀な変異であることを確認(1,000ゲノム:欠如、ESP6500:不在、ExAC:8.236e-06)
- 複数のバイオインフォマティクスツール(SIFT、Polyphen-2等)によるクロス予測は一貫して有害な影響を示した(SIFT:D、Polyphen2:D、MutationTaster_pred:D、VEST4スコア:0.953、REVELスコア:0.934、その他ツール:7D/1H)
- ClinVarやHGMDでは見つからない
- 近隣のミスセンス変異(c.451G>A p.Gly151Arg、c.452G>A p.Gly151Glu、c.470T>G p.Ile157Ser、c.473C>G p.Thr158Arg等)は、ClinVarで長QT症候群やアルドステロン症に関連する病原性変異として繰り返し記録されている

**ACMG評価:**

| 基準 | 根拠 | ACMG評価 | Tavtigian方式 |
|---|---|---|---|
| PM2(集団対照群に含まれていない) | ExAC:8.236e-06でのみ報告 | Moderate | +2点 |
| PP3(複数の計算的証拠が悪影響を支持) | SIFT/Polyphen2/MutationTaster全てD、VEST4=0.953、REVEL=0.934、他7D/1H→蛋白機能への悪影響と判断 | Supporting | +1点 |

ACMGガイドラインtable5ではModerate 1個+Supporting 1個は該当区分がないためVUS。Tavtigian方式では+3点でuncertain(VUS)。(ExACでは報告があるが、論文著者は希少性を根拠にPM2を採用)

### 結論

**Likely Pathogenic**であるMYBPC3遺伝子エクソン2:c.278delA (p.Lys93ArgfsTer3)を採用。
**VUS**であるKCNJ5遺伝子エクソン2:c.464G>A (p.Arg155Gln)を棄却。

---

## Case 2: Lethal neonatal hypertrophic cardiomyopathy from compound heterozygous MYBPC3 variants

症例報告:複合ヘテロ接合体MYBPC3変異株による致死性新生児肥厚性心筋症

### 概要

新生児期致死性HCMの分子剖検症例。父方由来のスプライス部位変異(c.2905+1G>A)と、母方由来の新規フレームシフト変異(c.836del; p.Gly279Valfs*21)——両方ともロスオブファンクション、複合ヘテロ接合。trio-WESで確認、ACMG基準で解釈済み。両親は無症候性キャリア。

「1つの正解に絞り込む」ではなく、「2つの変異が両方とも原因として必要」ということを判断できるか。

### 症状

死後分子診断の臨床的有用性を強調する症例報告。生後2か月の乳児が突然発症した急性心不全で死亡。法医学解剖、プロバンドおよび親に対してWESを実施。解剖の結果、重度のHCM、心房中隔欠損(ASD)、広範な心筋壊死および線維化が認められた。WESの結果、MYBPC3において既知の父系スプライス部位変異(c.2905+1G>A)と新規の母体切断フレームシフト変異(c.836del; p.Gly279Valfs*21)を特定。両バリアントともタンパク質機能の完全な喪失をもたらすと予測。

### 検知された2つの変異

- MYBPC3遺伝子エクソン27、既知、父方、スプライス部位変異 c.2905+1G>A
- MYBPC3遺伝子エクソン8、新規、母方、フレームシフト欠失 c.836del (p.Gly279Valfs*21)

### MYBPC3 c.2905+1G>A とは

機能の完全な喪失をもたらすと予測。標準的なドナー部位を破壊し、エクソンスキップやイントロン保持を引き起こす。LOFはMYBPC3の既知の疾患メカニズムである。

**ACMG評価:**

| 基準 | 根拠 | ACMG評価 | Tavtigian方式 |
|---|---|---|---|
| PVS1(機能喪失型変異) | 標準スプライス供与部位を破壊。MYBPC3で機能喪失は既知の病因 | Very Strong | +8点 |
| PS4(対照群と比較した有病率の有意な増加) | ClinVar(VCV000042666, https://www.ncbi.nlm.nih.gov/clinvar/variation/42666/)を確認すると心筋症症例と集団対照の比較論文(PMID:27532257, https://pubmed.ncbi.nlm.nih.gov/27532257/)が記載されている | Strong | +4点 |
| PM2(集団対照群データベースに認められない) | gnomAD:0.0014%、UK Biobank:保因者2人、1000genome:未検出。ClinVar/dbSNPにID附番あり(VCV000042666, rs397515991) | Moderate | +2点 |

ACMGガイドラインtable5ではVery Strong 1個+Strong 1個+Moderate 1個でPathogenic。Tavtigian方式では+14点でPathogenic。(ClinGenはPM2をSupportingへの格下げを推奨。報告はあるが著者は希少性を根拠にPM2採用)

### MYBPC3 c.836del (p.Gly279Valfs*21) とは

機能の完全な喪失をもたらすと予測。初期コード配列に早期終了コドン(PTC)を導入。このようなPTCを含む転写産物は通常NMDによる分解対象となる。この複合ヘテロ接合状態は機能的cMyBP-Cタンパク質がほぼ全く欠如したnull表現型をもたらす可能性が高い。ウエスタンでのタンパク確認は失敗。

**ACMG評価:**

| 基準 | 根拠 | ACMG評価 | Tavtigian方式 |
|---|---|---|---|
| PVS1(機能喪失型変異) | 標準スプライス供与部位を破壊(訳注:原文ママ、フレームシフトについての記載)。MYBPC3で機能喪失は既知の病因 | Very Strong | +8点 |
| PM2(集団対照群データベースに認められない) | gnomAD:未検出、UK Biobank:未検出、1000genome:未検出。ClinVar/dbSNPにID附番なし | Moderate | +2点 |

ACMGガイドラインtable5ではVery Strong 1個+Moderate 1個でLikely Pathogenic。Tavtigian方式では+10点でPathogenic。(ClinGenはPM2のSupportingへの格下げを推奨)

### 結論

両者の変異を採用。
**Pathogenic**であるMYBPC3遺伝子エクソン27:c.2905+1G>A
**Likely Pathogenic**であるMYBPC3遺伝子エクソン8:c.836del (p.Gly279Valfs*21)

---

## Case 3: Genetic landscape of hypertrophic cardiomyopathy in Hong Kong Chinese population

香港中国人集団における肥大性心筋症の遺伝的状況

### 概要

以前Pathogenicと分類されていたMYBPC3 p.Glu334Lysが、SG10Kコホート(4810人)中13人(0.03%)に見つかり、ClinVar上は現在conflicting interpretations。祖先集団に合わせた参照データによって、変異の真の頻度が明らかになった実例。同じ変異は香港の中国系コホート研究でもVUSとして3例報告されており、複数の東アジア系集団で裏付けが取れている。

「gnomAD全体では稀に見えるが、アジア系集団では実は無視できない頻度」の発見。

(症状のスライドはCase1とほぼ同様のためスキップ)

### 対象の変異について

この論文はHCMの患者53人についてシーケンスを実施し、検出された変異をACMGガイドラインに沿って判定した。P/LPが13変異、VUSが21変異発見された。その内、P/LPのMYH7 p.Arg719Trp(ClinGen承認済みPathogenic)とVUSのMYBPC3 p.Glu334LysについてACMGガイドラインに沿って判定を行う。

### MYH7 p.Arg719Trp

ClinGenを参照してガイドラインに沿って判定(https://erepo.clinicalgenome.org/evrepo/ui/classification/7b17a8bd-b169-46a1-8efa-b7388c4ecdef?version=1.0)

**ACMG評価:**

| 基準 | 根拠 | ACMG評価 |
|---|---|---|
| PS2(家族歴のないde novo変異) | トリオ解析の結果、MYH7 p.Arg719Trpを有するde novo変異が検出(PMID:10957787, https://pubmed.ncbi.nlm.nih.gov/10957787/) | Strong |
| PS3(遺伝子産物への有害な影響、in vitro/in vivo研究) | 血行動態の低下は、交感神経系およびレニン・アンジオテンシン・アルドステロン系の刺激の両方により有害な心筋リモデリングを誘発(PMID:24829265, https://pubmed.ncbi.nlm.nih.gov/24829265/)。**⚠️ノックイン・ノックアウトマウスの知識があまりないため確からしさの判定にあまり自信がない** | Strong |
| PS4(対照群と比較した有病率の有意な増加) | 常染色体顕性HCMの中国人非血縁6家系を調査し8人で当該変異を検出。非罹患者では保因者以外検出されず(PMID:19645038, https://pubmed.ncbi.nlm.nih.gov/19645038/) | Strong |

その他、PP3(computational tools predict damaging)、PM2(ExACに不在)、PM5(c.2156G>A p.Arg719Gln - Variation ID 14107 - Expert PanelによりPathogenic)、PP1、PM1にも該当。

ACMGガイドラインではPS1-4の中でStrong 2個以上のためPathogenic。

### MYBPC3 p.Glu334Lys

ClinVarでは「Conflicting classifications of pathogenicity」——Uncertain significance (9); Benign (2); Likely benign (3) (VCV000177902.53, https://www.ncbi.nlm.nih.gov/clinvar/variation/177902/)

**ACMG評価:**

| 基準 | 根拠 | ACMG評価 |
|---|---|---|
| PS3(遺伝子産物への有害な影響) | E334K cMyBPCによるUPSの障害が心臓イオンチャネルおよびCa2+ハンドリングタンパク質のレベルを変化させ、カルシウム過渡応答の振幅増大・電気生理学的機能障害を誘発(https://www.sciencedirect.com/science/article/abs/pii/S002228361100996X )。**⚠️ノックイン・ノックアウトマウスの知識があまりないため確からしさの判定にあまり自信がない** | Strong |
| BS1(対立遺伝子頻度が疾患の予想より高い) | gnomAD global AF=0.0002368、gnomAD East Asia AF=0.003338。East Asiaで比較的高頻度なため、罹患者が中国人であることを考慮しbenign方向とする | — |

ACMGガイドラインにおいてPathogenic側のPS3とBenign側のBS1が拮抗するためVUSと判定。**なお、gnomAD globalのAFのみ公開されていた際はPathogenicと判定されていた。**

---

## Case 4: Arrhythmogenic right ventricular cardiomyopathy in a Japanese patient with a homozygous founder variant of DSG2 in the East Asian population

東アジア集団におけるDSG2のホモ接合創始変異を持つ日本人患者における不整脈性右心室心筋症

### 概要

日本人ARVC(不整脈原性右室心筋症)患者の症例報告。ホモ接合のDSG2 p.Phe531Cys変異を保有。HGMDでは「disease-causing」とされていたが、実際の集団頻度は東アジア集団で明確に高い。TogoMCP経由でToMMo/NCBN/GEM-Jへの照会が行える。

### 症状

不整脈性右心室心筋症(ARVC)は、致命的な不整脈と心不全を引き起こす遺伝性心筋症。若年層やアスリートにおける不整脈性心停止の主な原因の一つで、有病率は1000人から5000人に1人と推定される。最も一般的な臨床症状は、思春期や若年成人における動悸または努力性失神であり、ECGの右前胸線でT波反転、心室不整脈、画像検査での右心室異常が見られる。東アジア人の創始変異と推定されるDSG2のホモ接合変異、c.1592T>G (p.Phe531Cys)が特定された。

### 検知された変異

DSG2のホモ接合変異 c.1592T>G (p.Phe531Cys)

### DSG2 p.Phe531Cys

ClinVarでは「Conflicting classifications of pathogenicity」——Likely pathogenic (3); Uncertain significance (7); Likely benign (1) (VCV000044283.23, https://www.ncbi.nlm.nih.gov/clinvar/variation/44283/)

**ACMG評価:**

| 基準 | 根拠 | ACMG評価 |
|---|---|---|
| PS3(遺伝子産物への有害な影響) | ヒトp.Phe531Cysに対応するマウスDsg2 p.Phe536Cysのノックインモデルを作製、野生型・ヘテロ接合・ホモ接合を比較。ホモ接合マウスでは両心室拡大・左室収縮機能低下・心電図上の伝導異常・心筋細胞脱落・両心室の線維化を確認(https://link.springer.com/article/10.1186/s12916-024-03593-8 )。**⚠️ノックイン・ノックアウトマウスの知識があまりないため確からしさの判定にあまり自信がない** | Strong |
| PS4(対照群と比較した有病率の有意な増加) | 家族スクリーニングでは、ホモ接合型変異保因者のみが100%の浸透率で明確なARVC表現型を示し、ヘテロ接合型変異保因者は影響を受けなかった(Chen L et al. 2019, https://doi.org/10.1016/j.ijcard.2018.06.105 ) | Strong |

ACMGガイドラインではPS1-4の中でStrong 2個以上のためPathogenic。

---

## 付録

### ACMGガイドライン原文(Richards et al. 2015)

**PVS1**(Very strong evidence of pathogenicity):Null variant(nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multi-exon deletion)in a gene where loss of function(LOF)is a known mechanism of disease

Caveats:
- Beware of genes where LOF is not a known disease mechanism(**例:GFAP, MYH7**)——**MYH7でPVS1を安易に使わないよう、原文自体が明記している点は要注意**
- Use caution interpreting LOF variants at the extreme 3' end of a gene
- Use caution with splice variants that are predicted to lead to exon skipping but leave the remainder of the protein intact
- Use caution in the presence of multiple transcripts

**PM2**:Absent from controls(or at extremely low frequency if recessive)in Exome Sequencing Project, 1000 Genomes, or ExAC

Caveat:Population data for indels may be poorly called by next generation sequencing

**PP3**:Multiple lines of computational evidence support a deleterious effect on the gene or gene product(conservation, evolutionary, splicing impact, etc)

Caveat:As many in silico algorithms use the same or very similar input for their predictions, each algorithm should not be counted as an independent criterion. PP3 can be used only once in any evaluation of a variant.

### Table 5: Rules for Combining Criteria to Classify Sequence Variants

**Pathogenic**
1. 1 Very Strong(PVS1)AND
   a. ≥1 Strong(PS1-PS4)、OR
   b. ≥2 Moderate(PM1-PM6)、OR
   c. 1 Moderate(PM1-PM6)and 1 Supporting(PP1-PP5)、OR
   d. ≥2 Supporting(PP1-PP5)
2. ≥2 Strong(PS1-PS4)、OR
3. 1 Strong(PS1-PS4)AND
   a. ≥3 Moderate(PM1-PM6)、OR
   b. 2 Moderate(PM1-PM6)AND ≥2 Supporting(PP1-PP5)、OR
   c. 1 Moderate(PM1-PM6)AND ≥4 Supporting(PP1-PP5)

**Likely Pathogenic**
1. 1 Very Strong(PVS1)AND 1 Moderate(PM1-PM6)、OR
2. 1 Strong(PS1-PS4)AND 1-2 Moderate(PM1-PM6)、OR
3. 1 Strong(PS1-PS4)AND ≥2 Supporting(PP1-PP5)、OR
4. ≥3 Moderate(PM1-PM6)、OR
5. 2 Moderate(PM1-PM6)AND ≥2 Supporting(PP1-PP5)、OR
6. 1 Moderate(PM1-PM6)AND ≥4 Supporting(PP1-PP5)

**Benign**
1. 1 Stand-Alone(BA1)、OR
2. ≥2 Strong(BS1-BS4)

**Likely Benign**
1. 1 Strong(BS1-BS4)and 1 Supporting(BP1-BP7)、OR
2. ≥2 Supporting(BP1-BP7)

### Tavtigian点数法(Bayesian point-based system)

**ACMG/AMP強度カテゴリごとの点数**

| Evidence Strength | Pathogenic | Benign |
|---|---|---|
| Indeterminate | 0 | 0 |
| Supporting | 1 | -1 |
| Moderate | 2 | -2 |
| Strong | 4 | -4 |
| Very Strong | 8 | -8 |

**点数に基づく分類区分**

| Category | Point ranges |
|---|---|
| Pathogenic | ≥10 |
| Likely Pathogenic | 6〜9 |
| Uncertain | 0〜5 |
| Likely Benign | -1〜-6 |
| Benign | ≤-7 |
