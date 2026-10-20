# 标记基因数据库

标记基因数据库（`data/db/D1_marker/`，随仓库分发的 **mini** 层）驱动默认的
`--mode marker` 鉴定：`blastn` 将查询 contigs 比对到 81 条精选参考序列，
再由 38 条物种规则把命中转化为判定。gside 是物种鉴定引擎：每条参考均
锚定可核查的 NCBI/RefSeq accession，并经 291 基因组面板 blast 验证。
种下分型是 gapit 的职责——霍乱弧菌 O1/O139 血清群库已移交该处
（`gapit screen --db vc_serogroup`）。

## marker 鉴定流程

1. `blastn` 将 contigs 比对到预建的 `markers_blastdb`（evalue 1e-10，词长 11）。
2. 命中需通过全局阈值：identity >= 85% 且 coverage >= 60%；逐基因最佳命中
   在**通过阈值的命中中**选取，短小的高一致性片段不会挤掉全长匹配。
3. 每个基因的最佳命中进入物种规则：规则内基因命中数达到 `min_hits` 且
   identity 不低于 `min_identity` 时规则触发；命中数最多的规则获胜。
4. 命中基因平均 identity >= 90% 判 **high** 置信度，否则 **medium**；
   临界命中的 marker 结果会建议进行 ANI 复核。

规则文件格式见[标记规则](../reference/marker-rules.zh.md)，marker 与 ANI 方法的仲裁关系见
[鉴定模式](../usage/modes.zh.md)。

## 基因角色

| 角色 | 含义 | 数量 |
|---|---|---|
| `primary` | Species-defining marker; a rule's `genes` draw from these | 46 |
| `confirm` | Corroborating evidence reported alongside the verdict | 22 |

## 基因清单（80 条序列）

| 基因 | 角色 | 参考文献 |
|---|---|---|
| `inva` | primary | M90846.1 invasion protein A [Salmonella enterica] |
| `rpoD` | confirm | M14427.1 rpoD operon region (dnaG-rpoD) [Salmonella enterica] |
| `iroB` | confirm | PV799927.1 glycosyltransferase iroB [Salmonella enterica] |
| `uida` | primary | NC_000913.3 beta-D-glucuronidase uidA [Escherichia coli] |
| `lacy` | primary | NC_000913.3 lactose permease lacY [Escherichia coli] |
| `gada` | confirm | NC_000913.3 glutamate decarboxylase alpha [Escherichia coli] |
| `stx1` | virulence | AF125520.1 Shiga toxin 1 subunit A [Escherichia coli] |
| `stx2` | virulence | PZ374831.1 Shiga toxin 2 subunit A (stx2A) [Escherichia coli] |
| `eae` | virulence | PZ374841.1 intimin (eae) [Escherichia coli] |
| `ipah` | primary | NC_004337.2 invasion plasmid antigen H ipaH [Shigella flexneri] |
| `toxr` | primary | NC_004603.1 toxR regulatory protein [Vibrio parahaemolyticus] |
| `tlh` | confirm | M36437.1 thermolabile hemolysin tlh [Vibrio parahaemolyticus] |
| `tdh` | virulence | M10069.1 thermostable direct hemolysin [Vibrio parahaemolyticus] |
| `vvha` | primary | PQ308339.1 hemolysin vvhA [Vibrio vulnificus] |
| `hly` | primary | M24199.1 listeriolysin O hly [Listeria monocytogenes] |
| `prs` | confirm | LC700425.1 phosphoribosyl pyrophosphate synthetase prs [Listeria monocytogenes] |
| `inlj` | confirm | AL592102.1 internalin J inlJ [Listeria monocytogenes] |
| `nuc` | primary | V01281.1 thermonuclease nuc [Staphylococcus aureus] |
| `fema` | confirm | X17688.1 methicillin resistance femA [Staphylococcus aureus] |
| `khe` | primary | PZ143400.1 hemolysin khe [Klebsiella pneumoniae] |
| `gyrb-kpn` | confirm | PZ127011.1 DNA gyrase subunit B gyrB [Klebsiella pneumoniae] |
| `oxa51` | primary | MZ351777.1 beta-lactamase blaOXA-51 [Acinetobacter baumannii] |
| `rpob-aba` | confirm | QB079536.1 RNA polymerase beta rpoB [Acinetobacter baumannii] |
| `ecfx` | primary | AE004091.2 sigma factor ECF ecfX [Pseudomonas aeruginosa] |
| `gyrb-pae` | confirm | PZ521782.1 DNA gyrase subunit B gyrB [Pseudomonas aeruginosa] |
| `ddl-ef` | primary | QB064932.1 D-alanine-D-alanine ligase ddl [Enterococcus faecalis] |
| `ddl-efm` | primary | PZ904265.1 D-alanine-D-alanine ligase ddl [Enterococcus faecium] |
| `lyta` | primary | PP971805.1 autolysin lytA [Streptococcus pneumoniae] |
| `psaa` | confirm | AF030366.1 manganese transport psaA [Streptococcus pneumoniae] |
| `cfb` | primary | PP475782.1 CAMP factor cfb [Streptococcus agalactiae] |
| `ctra` | primary | AY336258.1 outer membrane protein ctrA [Neisseria meningitidis] |
| `sodc` | confirm | PP471286.1 superoxide dismutase sodC [Neisseria meningitidis] |
| `pora` | primary | M21289.1 porin protein A porA [Neisseria gonorrhoeae] |
| `tcda` | primary | M30307.1 toxin A tcdA [Clostridioides difficile] |
| `tcdb` | confirm | PQ863731.1 toxin B tcdB [Clostridioides difficile] |
| `rpob-cs` | primary | MW380341.1 RNA polymerase beta rpoB [Cronobacter sakazakii] |
| `aera` | primary | PQ871546.1 aerolysin aerA [Aeromonas hydrophila] |
| `ahh1` | confirm | PZ301496.1 hemolysin ahh1 [Aeromonas hydrophila] |
| `ail` | primary | M29945.1 attachment invasion locus ail [Yersinia enterocolitica] |
| `yada` | confirm | X13881.1 Yersinia adhesin yadA [Yersinia enterocolitica] |
| `nheb` | virulence | Z11886.1 non-hemolytic enterotoxin B nheB [Bacillus cereus] |
| `ces` | virulence | AF110410.1 cereulide synthetase ces [Bacillus cereus] |
| `cpe` | virulence | L43549.1 enterotoxin cpe [Clostridium perfringens] |
| `mip` | primary | PX526202.1 macrophage infectivity potentiator mip [Legionella pneumophila] |
| `dota` | confirm | MW092166.1 dotA [Legionella pneumophila] |
| `cards` | confirm | AF390408.1 CARDS toxin [Mycoplasma pneumoniae] |
| `ompa-cp` | primary | Z31593.1 major outer membrane protein [Chlamydia pneumoniae] |
| `ompa-cps` | primary | PZ860298.1 major outer membrane protein ompA [Chlamydia psittaci] |
| `urea` | primary | M60398.1 urease subunit alpha [Helicobacter pylori] |
| `urec` | confirm | M60398.1 urease subunit gamma/ureC [Helicobacter pylori] |
| `caga` | virulence | AB015416.1 cytotoxin-associated gene A [Helicobacter pylori] |
| `ptxs1` | confirm | M13223.1 pertussis toxin S1 subunit [Bordetella pertussis] |
| `tox` | primary | KX703000.1 diphtheria toxin tox [Corynebacterium diphtheriae] |
| `hpd` | primary | U32723.1 hemoglobin-binding protein D [Haemophilus influenzae] |
| `gdh` | primary | MT542692.1 glutamate dehydrogenase gdh [Streptococcus suis] |
| `bimabp` | primary | AF505187.1 bacterial actin motility BimA(Bp) [Burkholderia pseudomallei] |
| `tts1` | confirm | AY089501.1 type III secretion system 1 [Burkholderia pseudomallei] |
| `lipl32` | primary | PZ599989.1 lipoprotein lipL32 [Leptospira interrogans] |
| `tpp47` | primary | AE000520.1:622266-623570 47 kDa membrane lipoprotein tpp47 [Treponema pallidum] |
| `gyrb-ps` | primary | AB015801.1 DNA gyrase subunit B [Plesiomonas shigelloides] |
| `eacf` | primary | CP099906.1:3631386-3634232 intimin eae [Escherichia albertii] |
| `toxrf` | primary | OR123876.1 transmembrane regulatory protein toxR [Vibrio fluvialis] |
| `bontA` | primary | M30171.1 botulinum neurotoxin type A [Clostridium botulinum] |
| `bontB` | typing | M81186.1 botulinum neurotoxin type B [Clostridium botulinum] |
| `bontE` | typing | X62089.1 botulinum neurotoxin type E [Clostridium botulinum] |
| `paga` | primary | M22504.1 protective antigen pagA [Bacillus anthracis] |
| `capb` | confirm | M64091.1 capsule biosynthesis capB [Bacillus anthracis] |
| `caf1` | primary | X13880.1 F1 capsule antigen caf1 [Yersinia pestis] |
| `pla` | confirm | PX130171.1 plasminogen activator pla pPCP1 [Yersinia pestis] |
| `trh` | virulence | AB663491.1 thermostable direct hemolysin-related trh [Vibrio parahaemolyticus] |
| `p1` | primary | U00089.2:180858-185741 P1 adhesin operon mgpA [Mycoplasma pneumoniae] |
| `is481` | primary | MK841515.1 insertion sequence IS481 transposase [Bordetella pertussis] |
| `cadf` | primary | GCF_000009085.1 cadf [Campylobacter jejuni] |
| `ceue` | primary | GCF_009730395.1 ceue [Campylobacter coli] |
| `hipo` | primary | GCF_000009085.1 hipo [Campylobacter jejuni] |
| `mapa` | primary | GCF_000009085.1 mapa [Campylobacter jejuni] |
| `speb` | primary | GCF_000006785.1 speb [Streptococcus pyogenes] |
| `ompw` | primary | X51948.1 outer membrane protein W [Vibrio cholerae] |
| `cpa` | primary | PV804807.1 phospholipase C alpha toxin (cpa/plc) [Clostridium perfringens] |
| `safC` | typing | NC_003197.2:343560-346076 saf atypical fimbrial operon [Salmonella enterica] |

## 物种覆盖（38 条规则）

| 物种 | 规则基因 | 最少命中 | 最低 identity |
|---|---|---|---|
| *Salmonella* | `inva` | 1 | 90%
| *V parahaemolyticus* | `tlh` | 1 | 90%
| *Listeria monocytogenes* | `hly` | 1 | 90%
| *Vibrio cholerae* | `ompw` | 1 | 90%
| *Vibrio vulnificus* | `vvha` | 1 | 90%
| *Staphylococcus aureus* | `nuc` | 1 | 90%
| *Klebsiella pneumoniae* | `khe` | 1 | 90%
| *Acinetobacter baumannii* | `oxa51` | 1 | 90%
| *Pseudomonas aeruginosa* | `ecfx` | 1 | 90%
| *Enterococcus faecalis* | `ddl-ef` | 1 | 90%
| *Enterococcus faecium* | `ddl-efm` | 1 | 90%
| *Streptococcus pneumoniae* | `lyta` | 1 | 90%
| *Streptococcus pyogenes* | `speb` | 1 | 90%
| *Streptococcus agalactiae* | `cfb` | 1 | 90%
| *Neisseria meningitidis* | `ctra` | 1 | 90%
| *Neisseria gonorrhoeae* | `pora` | 1 | 90%
| *Clostridioides difficile* | `tcda, tcdb` | 1 | 90%
| *Cronobacter sakazakii* | `rpob-cs` | 1 | 95% (*excludes* `speb`)
| *Legionella pneumophila* | `mip` | 1 | 90%
| *Mycoplasma pneumoniae* | `p1` | 1 | 90%
| *Helicobacter pylori* | `urea, urec` | 1 | 90%
| *Bordetella pertussis* | `is481, ptxs1` | 1 | 90%
| *Corynebacterium diphtheriae* | `tox` | 1 | 90%
| *Haemophilus influenzae* | `hpd` | 1 | 90%
| *Streptococcus suis* | `gdh` | 1 | 90%
| *Burkholderia pseudomallei* | `tts1, bimabp` | 1 | 90%
| *Leptospira interrogans* | `lipl32` | 1 | 90%
| *Treponema pallidum* | `tpp47` | 1 | 90%
| *Bacillus anthracis* | `paga, capb` | 2 | 90%
| *Yersinia pestis* | `caf1, pla` | 2 | 90%
| *Yersinia enterocolitica* | `ail, yada` | 1 | 90%
| *Aeromonas hydrophila* | `aera, ahh1` | 1 | 90%
| *Clostridium botulinum* | `bontA, bontB, bontE` | 1 | 90%
| *Clostridium perfringens* | `cpa` | 1 | 90%
| *DEC* | `uida, lacy, gada` | 2 | 95%
| *Shigella EIEC* | `ipah` | 1 | 95%
| *Campylobacter jejuni* | `mapa, hipo, cadf` | 2 | 90%
| *Campylobacter coli* | `ceue` | 1 | 90% (*excludes* `mapa, hipo, cadf`)

## 自定义

编辑 `data/db/D1_marker/marker_rules.yaml`（字段说明见
[标记规则](../reference/marker-rules.zh.md)）；修改 `markers.fasta` 后需重建 BLAST 库：

```bash
makeblastdb -in data/db/D1_marker/markers.fasta -dbtype nucl -out data/db/D1_marker/markers_blastdb
```
