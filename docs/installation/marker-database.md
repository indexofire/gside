# Marker Gene Database

The marker database (`data/db/D1_marker/`, the bundled **mini** tier) powers the
default `-m marker` identification: `blastn` maps query contigs against
188 curated reference sequences, and 48 species rules turn the hits into a
verdict. gside is a species-identification engine: every reference is anchored
to a verifiable NCBI/RefSeq accession and was validated by blast against its
reference genome and, where available, a 100-genome real-world panel.
Sub-species typing is gspec's job — the V. cholerae
O1/O139 serogroup database lives there (`gspec screen --db vc_serogroup`).
Shigella follows the same rule: gside reports the **genus** `Shigella` only
(uidA+ipaH, lacY-negative) and adds a note pointing at gspec for
species/serotype discrimination.

For a user-provided marker database, see
[Custom marker databases](../usage/custom-databases.md). The algorithm is the
same; only the FASTA/rules change.

## How marker identification works

**FASTA**: the path below. **FASTQ**: minimap2 maps reads to the same
`markers.fasta` (`-x sr` for short reads, `-x map-ont` for long), per-gene
breadth is the union of the covered reference intervals, and the identity and
rule thresholds are identical. Reads only need ~1x genome coverage; see
[Identification Modes](../usage/modes.md#l1-marker-target-gene-combination).

1. `blastn` aligns contigs against the pre-built `markers_blastdb`
   (evalue 1e-10, word size 11).
2. Hits pass the global thresholds: identity >= 85% and coverage >= 60%; the
   best hit per gene is selected **among passing hits**, so a short
   high-identity fragment cannot displace a full-length match.
3. The best hit per gene feeds the species rules: a rule fires when enough of
   its genes (`min_hits`) match at `min_identity` or better; the rule with the
   most matched genes wins.
4. Mean identity of matched genes >= 90% reports **high** confidence,
   otherwise **medium**; near-threshold marker calls suggest an ANI recheck.

See [Marker Rules](../reference/marker-rules.md) for the rule-file format and
[Identification Modes](../usage/modes.md) for how marker results arbitrate
against ANI methods.

## Gene roles

| Role | Meaning | Count |
|---|---|---|
| `primary` | Species-defining marker; a rule's `genes` draw from these | 139 |
| `confirm` | Corroborating evidence reported alongside the verdict | 26 |
| `typing` | Subtype/toxin-type reporting (e.g. BoNT serotypes) | 3 |
| `virulence` | Virulence/toxigenicity reporting; never gates the species call | 2 |
| `control` | Single-copy reads-mode depth controls; never rule genes | 18 |

## Gene inventory (188 sequences)

| Gene | Role | Reference |
|---|---|---|
| `ail` | primary | NC_008800.1:2007722-2008258 attachment invasion locus ail, pathogenic-bioserotype marker [Yersinia enterocolitica 8081] |
| `asst` | primary | NC_015761.1:34270-35988 aryl-sulfate sulfotransferase [Salmonella bongori NCTC 12419] |
| `ba5345` | primary | CP009981.1:4873800-4874400 lambdaBa03 prophage chromosomal marker BA_5345 [Bacillus anthracis Ames BA1004] |
| `bonta` | primary | NC_009495.1:909709-913599 botulinum neurotoxin type A1 [Clostridium botulinum ATCC 3502] |
| `bontc` | primary | D90210.1&rettype=fasta botulinum neurotoxin type C1 [Clostridium botulinum phage c-st] |
| `bontd` | primary | X54254.1:1-3870 botulinum neurotoxin type D [Clostridium botulinum] |
| `bontf` | primary | NZ_CP013710.1:158109-161942 botulinum neurotoxin subtype F5 [Clostridium botulinum F634] |
| `btub-yp` | primary | NC_003143.1:4391024-4392901 TonB-dependent vitamin B12 receptor BtuB, chromosomal marker [Yersinia pestis CO92] |
| `cadf` | primary | GCF_000009085.1 cadf [Campylobacter jejuni] |
| `caf1` | primary | X13880.1 F1 capsule antigen caf1 [Yersinia pestis] |
| `ceue` | primary | GCF_009730395.1 ceue [Campylobacter coli] |
| `cfb` | primary | PP475782.1 CAMP factor cfb [Streptococcus agalactiae] |
| `cpa` | primary | PV804807.1 phospholipase C alpha toxin (cpa/plc) [Clostridium perfringens] |
| `ctra` | primary | AY336258.1 outer membrane protein ctrA [Neisseria meningitidis] |
| `cupin-bp` | primary | NZ_CP025371.1:2754528-2755208 cupin domain-containing protein, chromosomal species marker [Bordetella pertussis] |
| `ddl-ef` | primary | QB064932.1 D-alanine-D-alanine ligase ddl [Enterococcus faecalis] |
| `ddl-efm` | primary | PZ904265.1 D-alanine-D-alanine ligase ddl [Enterococcus faecium] |
| `eacf` | primary | CP099906.1:3631386-3634232 intimin eae [Escherichia albertii] |
| `ecfx` | primary | NC_002516.2:169361-169906 ECF subfamily sigma-70 factor ecfX/PA0149 [Pseudomonas aeruginosa] |
| `gdh` | primary | MT542692.1 glutamate dehydrogenase gdh [Streptococcus suis] |
| `glgb-bp` | primary | NZ_AP028077.1:1548836-1551052 GlgB 1,4-alpha-glucan branching enzyme (species anchor); 99.6-100% identity in 100 clinical B. pseudomallei genomes; nearest cross B. thailandensis 73.3%, B. cepacia 82.5% |
| `glk-ye` | primary | NC_008800.1:412574-413737 glycerate kinase, chromosomal species anchor [Yersinia enterocolitica 8081] |
| `gyrb-cf` | primary | GCF_011064845.1:G4551_RS00020 ref locus |
| `gyrb-mp` | primary | NZ_LR214945.1:810714-812666 DNA topoisomerase subunit B, second chromosomal marker [Mycoplasma pneumoniae NCTC10119] |
| `gyrb-ps` | primary | AB015801.1 DNA gyrase subunit B [Plesiomonas shigelloides] |
| `hipo` | primary | GCF_000009085.1 hipo [Campylobacter jejuni] |
| `hly` | primary | M24199.1 listeriolysin O hly [Listeria monocytogenes] |
| `hpd` | primary | U32723.1 hemoglobin-binding protein D [Haemophilus influenzae] |
| `inva` | primary | M90846.1 invasion protein A [Salmonella enterica] |
| `ipah` | primary | NC_004851.1:182187-182294 invasion plasmid antigen H ipaH conserved amplicon [Shigella flexneri] |
| `is1081` | primary | PZ498641.1 IS1081 family transposase fragment (300bp); MTBC 96/100 at 95%, non-TB Mycobacterium 2/43 (MTBC members) |
| `is6110` | primary | cpin-build insertion sequence IS6110 [GCF_000195955.2] |
| `isph-cdf` | primary | NZ_CP076401.1:1943117-1943956 4-hydroxy-3-methylbut-2-enyl diphosphate reductase IspH, second chromosomal marker [Clostridioides difficile CD6] |
| `khe` | primary | PZ143400.1 hemolysin khe [Klebsiella pneumoniae] |
| `lipl32` | primary | PZ599989.1 lipoprotein lipL32 [Leptospira interrogans] |
| `lyta3` | primary | NC_003028.3:1840882-1840405 autolysin lytA 3' half (fragmentation-resistant split) [Streptococcus pneumoniae TIGR4] |
| `lyta5` | primary | NC_003028.3:1841361-1840883 autolysin lytA 5' half (fragmentation-resistant split) [Streptococcus pneumoniae TIGR4] |
| `mapa` | primary | GCF_000009085.1 mapa [Campylobacter jejuni] |
| `mena-cd` | primary | NZ_LN831026.1:333379-334341 o-succinylbenzoate synthase MenC, second chromosomal marker [Corynebacterium diphtheriae NCTC 11397] |
| `mip` | primary | PX526202.1 macrophage infectivity potentiator mip [Legionella pneumophila] |
| `muri-cdf` | primary | NZ_CP076401.1:3962663-3963469 glutamate racemase MurI, chromosomal species anchor [Clostridioides difficile CD6] |
| `nuc` | primary | V01281.1 thermonuclease nuc [Staphylococcus aureus] |
| `ompa-cp` | primary | Z31593.1 major outer membrane protein [Chlamydia pneumoniae] |
| `ompa-cpn2` | primary | cpin-build major outer membrane protein [C. pneumoniae J138] |
| `ompa-cps` | primary | PZ860298.1 major outer membrane protein ompA [Chlamydia psittaci] |
| `ompw` | primary | X51948.1 outer membrane protein W [Vibrio cholerae] |
| `oxa51` | primary | MZ351777.1 beta-lactamase blaOXA-51 [Acinetobacter baumannii] |
| `p1a` | primary | U00089.2:180858-185741[p1a] P1 adhesin quarter 1/4 (fragmentation-resistant split) [Mycoplasma pneumoniae] |
| `p1c` | primary | U00089.2:180858-185741[p1c] P1 adhesin quarter 3/4 (fragmentation-resistant split) [Mycoplasma pneumoniae] |
| `paga` | primary | NC_007322.2:143779-146073 protective antigen pagA on pXO1 [Bacillus anthracis Ames Ancestor] |
| `pora` | primary | M21289.1 porin protein A porA [Neisseria gonorrhoeae] |
| `porb1b` | primary | NZ_KI391931.1:186504-187487 trimeric porin PorB.IB porB1b allele family [Neisseria gonorrhoeae FA19] |
| `rela-cd` | primary | NZ_LN831026.1:1595158-1596177 GTP pyrophosphokinase RelA/SpoT, chromosomal species anchor [Corynebacterium diphtheriae NCTC 11397] |
| `rfbe` | primary | NC_004337.2:2113945-2115201 polysaccharide biosynthesis protein rfbE [Shigella flexneri] |
| `rfbi` | primary | NC_004337.2:2110581-2110955 glycosyl translocase rfbI [Shigella flexneri] |
| `rfc` | primary | NC_004337.2:2110959-2112107 O-antigen polymerase rfc [Shigella flexneri] |
| `rfpb` | primary | CP000640.1:2979-4112 glycosyltransferase rfpB [Shigella dysenteriae] |
| `ribf-hpi` | primary | NZ_CP065991.1:1099665-1100312 GTP cyclohydrolase II RibF, chromosomal species anchor [Haemophilus parainfluenzae FDAARGOS_1000] |
| `rpob-cs` | primary | MW380341.1 RNA polymerase beta rpoB [Cronobacter sakazakii] |
| `sb11_wzx` | primary | ShigaTyper ShigellaRef5 Sb11_wzx [Shigella boydii] |
| `sb11_wzy` | primary | ShigaTyper ShigellaRef5 Sb11_wzy [Shigella boydii] |
| `sb12_wzx` | primary | ShigaTyper ShigellaRef5 Sb12_wzx [Shigella boydii] |
| `sb12_wzy` | primary | ShigaTyper ShigellaRef5 Sb12_wzy [Shigella boydii] |
| `sb14_wzx` | primary | ShigaTyper ShigellaRef5 Sb14_wzx [Shigella boydii] |
| `sb14_wzy` | primary | ShigaTyper ShigellaRef5 Sb14_wzy [Shigella boydii] |
| `sb15_wzx` | primary | ShigaTyper ShigellaRef5 Sb15_wzx [Shigella boydii] |
| `sb15_wzy` | primary | ShigaTyper ShigellaRef5 Sb15_wzy [Shigella boydii] |
| `sb16_wzy` | primary | ShigaTyper ShigellaRef5 Sb16_wzy [Shigella boydii] |
| `sb17_wzx` | primary | ShigaTyper ShigellaRef5 Sb17_wzx [Shigella boydii] |
| `sb17_wzy` | primary | ShigaTyper ShigellaRef5 Sb17_wzy [Shigella boydii] |
| `sb18_wzx` | primary | ShigaTyper ShigellaRef5 Sb18_wzx [Shigella boydii] |
| `sb18_wzy` | primary | ShigaTyper ShigellaRef5 Sb18_wzy [Shigella boydii] |
| `sb19_wzx` | primary | ShigaTyper ShigellaRef5 Sb19_wzx [Shigella boydii] |
| `sb19_wzy` | primary | ShigaTyper ShigellaRef5 Sb19_wzy [Shigella boydii] |
| `sb1_wzx` | primary | ShigaTyper ShigellaRef5 Sb1_wzx [Shigella boydii] |
| `sb1_wzy` | primary | ShigaTyper ShigellaRef5 Sb1_wzy [Shigella boydii] |
| `sb2_wzx` | primary | ShigaTyper ShigellaRef5 Sb2_wzx [Shigella boydii] |
| `sb2_wzy` | primary | ShigaTyper ShigellaRef5 Sb2_wzy [Shigella boydii] |
| `sb3_wzx` | primary | ShigaTyper ShigellaRef5 Sb3_wzx [Shigella boydii] |
| `sb3_wzy` | primary | ShigaTyper ShigellaRef5 Sb3_wzy [Shigella boydii] |
| `sb4_wzx` | primary | ShigaTyper ShigellaRef5 Sb4_wzx [Shigella boydii] |
| `sb4_wzy` | primary | ShigaTyper ShigellaRef5 Sb4_wzy [Shigella boydii] |
| `sb5_wzx` | primary | ShigaTyper ShigellaRef5 Sb5_wzx [Shigella boydii] |
| `sb5_wzy` | primary | ShigaTyper ShigellaRef5 Sb5_wzy [Shigella boydii] |
| `sb6_wzx` | primary | ShigaTyper ShigellaRef5 Sb6_wzx [Shigella boydii] |
| `sb6_wzy` | primary | ShigaTyper ShigellaRef5 Sb6_wzy [Shigella boydii] |
| `sb7_wzx` | primary | ShigaTyper ShigellaRef5 Sb7_wzx [Shigella boydii] |
| `sb7_wzy` | primary | ShigaTyper ShigellaRef5 Sb7_wzy [Shigella boydii] |
| `sb8_wzx` | primary | ShigaTyper ShigellaRef5 Sb8_wzx [Shigella boydii] |
| `sb8_wzy` | primary | ShigaTyper ShigellaRef5 Sb8_wzy [Shigella boydii] |
| `sb9_wzx` | primary | ShigaTyper ShigellaRef5 Sb9_wzx [Shigella boydii] |
| `sb9_wzy` | primary | ShigaTyper ShigellaRef5 Sb9_wzy [Shigella boydii] |
| `sbprov_wzx` | primary | ShigaTyper ShigellaRef5 SbProv_wzx [Shigella boydii] |
| `sbprov_wzy` | primary | ShigaTyper ShigellaRef5 SbProv_wzy [Shigella boydii] |
| `sd10_wzx` | primary | ShigaTyper ShigellaRef5 Sd10_wzx [Shigella dysenteriae] |
| `sd10_wzy` | primary | ShigaTyper ShigellaRef5 Sd10_wzy [Shigella dysenteriae] |
| `sd12_wzx` | primary | ShigaTyper ShigellaRef5 Sd12_wzx [Shigella dysenteriae] |
| `sd12_wzy` | primary | ShigaTyper ShigellaRef5 Sd12_wzy [Shigella dysenteriae] |
| `sd13_wzx` | primary | ShigaTyper ShigellaRef5 Sd13_wzx [Shigella dysenteriae] |
| `sd13_wzy` | primary | ShigaTyper ShigellaRef5 Sd13_wzy [Shigella dysenteriae] |
| `sd14_wzx` | primary | ShigaTyper ShigellaRef5 Sd14_wzx [Shigella dysenteriae] |
| `sd14_wzy` | primary | ShigaTyper ShigellaRef5 Sd14_wzy [Shigella dysenteriae] |
| `sd15_wzx` | primary | ShigaTyper ShigellaRef5 Sd15_wzx [Shigella dysenteriae] |
| `sd15_wzy` | primary | ShigaTyper ShigellaRef5 Sd15_wzy [Shigella dysenteriae] |
| `sd1_rfp` | primary | ShigaTyper ShigellaRef5 Sd1_rfp [Shigella dysenteriae] |
| `sd1_wzx` | primary | ShigaTyper ShigellaRef5 Sd1_wzx [Shigella dysenteriae] |
| `sd1_wzy` | primary | ShigaTyper ShigellaRef5 Sd1_wzy [Shigella dysenteriae] |
| `sd2_wzx` | primary | ShigaTyper ShigellaRef5 Sd2_wzx [Shigella dysenteriae] |
| `sd2_wzy` | primary | ShigaTyper ShigellaRef5 Sd2_wzy [Shigella dysenteriae] |
| `sd3_wzx` | primary | ShigaTyper ShigellaRef5 Sd3_wzx [Shigella dysenteriae] |
| `sd3_wzy` | primary | ShigaTyper ShigellaRef5 Sd3_wzy [Shigella dysenteriae] |
| `sd4_wzx` | primary | ShigaTyper ShigellaRef5 Sd4_wzx [Shigella dysenteriae] |
| `sd4_wzy` | primary | ShigaTyper ShigellaRef5 Sd4_wzy [Shigella dysenteriae] |
| `sd6_wzx` | primary | ShigaTyper ShigellaRef5 Sd6_wzx [Shigella dysenteriae] |
| `sd6_wzy` | primary | ShigaTyper ShigellaRef5 Sd6_wzy [Shigella dysenteriae] |
| `sd8_wzx` | primary | ShigaTyper ShigellaRef5 Sd8_wzx [Shigella dysenteriae] |
| `sd8_wzy` | primary | ShigaTyper ShigellaRef5 Sd8_wzy [Shigella dysenteriae] |
| `sd9_wzx` | primary | ShigaTyper ShigellaRef5 Sd9_wzx [Shigella dysenteriae] |
| `sd9_wzy` | primary | ShigaTyper ShigellaRef5 Sd9_wzy [Shigella dysenteriae] |
| `sdprov_wzx` | primary | ShigaTyper ShigellaRef5 SdProv_wzx [Shigella dysenteriae] |
| `sdprov_wzy` | primary | ShigaTyper ShigellaRef5 SdProv_wzy [Shigella dysenteriae] |
| `sdprove_wzx` | primary | ShigaTyper ShigellaRef5 SdProvE_wzx [Shigella dysenteriae] |
| `sdprove_wzy` | primary | ShigaTyper ShigellaRef5 SdProvE_wzy [Shigella dysenteriae] |
| `sera-hpi` | primary | NZ_CP065991.1:1996188-1997420 phosphoglycerate dehydrogenase SerA, second chromosomal marker [Haemophilus parainfluenzae FDAARGOS_1000] |
| `sonnei_psuk` | primary | GCF_002950395.1:595851-596792 pseudouridine kinase [Shigella sonnei] |
| `speb` | primary | GCF_000006785.1 speb [Streptococcus pyogenes] |
| `spi-yp` | primary | NC_003143.1:1900484-1901470 sugar phosphate isomerase/epimerase, chromosomal species anchor [Yersinia pestis CO92] |
| `ssav` | primary | NC_003197.2:1493647-1495692 SPI-2 type III secretion system apparatus protein SsaV [Salmonella enterica Typhimurium LT2] |
| `ssec` | primary | NC_003197.2:1484988-1486442 SPI-2 type III secretion system translocon protein SseC [Salmonella enterica Typhimurium LT2] |
| `tcda` | primary | M30307.1 toxin A tcdA [Clostridioides difficile] |
| `tox` | primary | KX703000.1 diphtheria toxin tox [Corynebacterium diphtheriae] |
| `toxrf` | primary | OR123876.1 transmembrane regulatory protein toxR [Vibrio fluvialis] |
| `tpp47` | primary | AE000520.1:622266-623570 47 kDa membrane lipoprotein tpp47 [Treponema pallidum] |
| `tuf-mp` | primary | NZ_LR214945.1:25477-26661 elongation factor Tu, chromosomal species anchor [Mycoplasma pneumoniae NCTC10119] |
| `uida` | primary | NC_000913.3:1694260-1696071 beta-D-glucuronidase uidA [Escherichia coli] |
| `urea` | primary | M60398.1 urease subunit alpha [Helicobacter pylori] |
| `vvha` | primary | PQ308339.1 hemolysin vvhA [Vibrio vulnificus] |
| `wbgz` | primary | CP000039.1:204752-206671 putative epimerase/dehydratase wbgZ [Shigella sonnei] |
| `zadh` | primary | NC_015761.1:1659167-1660180 zinc-dependent alcohol dehydrogenase [Salmonella bongori NCTC 12419] |
| `adk` | confirm | CP009981.1:4359650-4360300 adenylate kinase adk, Bacillus cereus group common marker [Bacillus anthracis Ames BA1004] |
| `capb` | confirm | NC_007323.3:55599-56993 capsule biosynthesis capB on pXO2 [Bacillus anthracis Ames Ancestor] |
| `cards` | confirm | AF390408.1 CARDS toxin [Mycoplasma pneumoniae] |
| `cdsa-citrobacter` | confirm | GCF_033726075.1:U4W25_RS05535 Citrobacter cdsA phosphatidate cytidylyltransferase (858bp); covers C. amalonaticus/farmeri at 95-100%, Enterobacteriaceae max 68-86% |
| `dota` | confirm | MW092166.1 dotA [Legionella pneumophila] |
| `fema` | confirm | X17688.1 methicillin resistance femA [Staphylococcus aureus] |
| `gyrb-kpn` | confirm | PZ127011.1 DNA gyrase subunit B gyrB [Klebsiella pneumoniae] |
| `gyrb-pae` | confirm | PZ521782.1 DNA gyrase subunit B gyrB [Pseudomonas aeruginosa] |
| `inlj` | confirm | AL592102.1 internalin J inlJ [Listeria monocytogenes] |
| `irob` | confirm | PV799927.1 glycosyltransferase iroB [Salmonella enterica] |
| `lacy` | confirm | NC_000913.3:361918-362218 lactose permease lacY 5' region [Escherichia coli] |
| `leus-citrobacter` | confirm | GCF_011064845.1:G4551_RS07165 Citrobacter leuS leucine--tRNA ligase (2583bp); E. coli/Salmonella/Klebsiella/Cronobacter max 85-89% at 400 genomes |
| `pla` | confirm | PX130171.1 plasminogen activator pla pPCP1 [Yersinia pestis] |
| `prs` | confirm | LC700425.1 phosphoribosyl pyrophosphate synthetase prs [Listeria monocytogenes] |
| `psaa` | confirm | AF030366.1 manganese transport psaA [Streptococcus pneumoniae] |
| `ptxs1` | confirm | NZ_CP025371.1:84308-85117 pertussis toxin S1 subunit ptxA [Bordetella pertussis] |
| `rpob-aba` | confirm | QB079536.1 RNA polymerase beta rpoB [Acinetobacter baumannii] |
| `rpod` | confirm | M14427.1 rpoD operon region (dnaG-rpoD) [Salmonella enterica] |
| `sodc` | confirm | PP471286.1 superoxide dismutase sodC [Neisseria meningitidis] |
| `tcdb` | confirm | PQ863731.1 toxin B tcdB [Clostridioides difficile] |
| `tcfa3` | confirm | NZ_CP025371.1:1277887-1276858R tracheal colonization factor tcfA 3' half (fragmentation-resistant split) [Bordetella pertussis] |
| `tcfa5` | confirm | NZ_CP025371.1:1274915-1277886 tracheal colonization factor tcfA 5' half (fragmentation-resistant split) [Bordetella pertussis] |
| `tlh` | confirm | M36437.1 thermolabile hemolysin tlh [Vibrio parahaemolyticus] |
| `urec` | confirm | M60398.1 urease subunit gamma/ureC [Helicobacter pylori] |
| `yada` | confirm | NC_008791.1:47054-48322 Yersinia adhesin yadA on pYV plasmid, virulence reporting [Yersinia enterocolitica 8081] |
| `yddg-ye` | confirm | NC_008800.1:823440-824327 aromatic amino acid transporter YddG, second chromosomal marker [Yersinia enterocolitica 8081] |
| `bontb` | typing | M81186.1 botulinum neurotoxin type B [Clostridium botulinum] |
| `bonte` | typing | X62089.1 botulinum neurotoxin type E [Clostridium botulinum] |
| `safc` | typing | NC_003197.2:343560-346076 saf atypical fimbrial operon [Salmonella enterica] |
| `bima-bp` | virulence | NZ_AP028078.1:2020709-2022298 BimA(Bp) actin-nucleating factor (virulence reporter); present in 95/100 clinical genomes at 95%+ identity |
| `sctn-bp` | virulence | NZ_AP028078.1:1898238-1899587 T3SS ATPase SctN (virulence reporter); present in 99/100 clinical genomes at 95%+ identity |
| `lpxa` | control | AL111168.1:251030-251821 lpxA [Campylobacter jejuni] |
| `lpxa` | control | lpxA_2 lpxA [Campylobacter coli] |
| `rpob-bacillus_anthracis` | control | NC_007530.2:108391-111924 RNA polymerase beta subunit rpoB, single-copy chromosome control [Bacillus anthracis Ames Ancestor] |
| `rpob-bordetella_pertussis` | control | NC_002929.2:12806-13978 rpoB [Bordetella pertussis] |
| `rpob-chlamydia_pneumoniae` | control | cpin-build rpoB window [GCF_000008745.1; CPN_RS00410] |
| `rpob-clostridioides_difficile` | control | NC_009089.1:91896-92362 rpoB [Clostridioides difficile] |
| `rpob-clostridium_perfringens` | control | NZ_CP065681.1:591382-592221 rpoB [Clostridium perfringens] |
| `rpob-corynebacterium_diphtheriae` | control | NC_002935.2:410798-411818 rpoB [Corynebacterium diphtheriae] |
| `rpob-enterococcus_faecalis` | control | NZ_KB944666.1:168618-169154 rpoB [Enterococcus faecalis] |
| `rpob-escherichia_coli` | control | NC_002695.2:4991564-4992488 rpoB [Escherichia coli] |
| `rpob-haemophilus_influenzae` | control | NC_000907.1:536660-537984 rpoB [Haemophilus influenzae] |
| `rpob-legionella_pneumophila` | control | NC_002942.5:371299-372178 rpoB [Legionella pneumophila] |
| `rpob-mycobacterium_tuberculosis` | control | cpin-build rpoB window [GCF_000195955.2; Rv0667] |
| `rpob-neisseria_gonorrhoeae` | control | NC_002946.2:1821201-1822430 rpoB [Neisseria gonorrhoeae] |
| `rpob-shigella` | control | NC_004337.2:4202294-4203218 rpoB [Shigella] |
| `rpob-streptococcus_agalactiae` | control | NC_004116.1:180465-181082 rpoB [Streptococcus agalactiae] |
| `rpob-streptococcus_pyogenes` | control | NC_002737.1:94763-95380 rpoB [Streptococcus pyogenes] |
| `rpob-vibrio_parahaemolyticus` | control | NC_004603.1:3114701-3115624 rpoB [Vibrio parahaemolyticus] |

## Species coverage (48 rules)

| Species | Rule genes | min hits | min identity | Constraints |
|---|---|---|---|---|
| *Salmonella bongori* | `inva`, `asst`, `zadh` | 2 | 90% | require `inva`; priority over *Salmonella*, *Salmonella enterica* |
| *Salmonella enterica* | `inva`, `ssav`, `ssec` | 2 | 90% | require `inva`; exclude `asst`, `zadh`; priority over *Salmonella* |
| *Salmonella* | `inva` | 1 | 90% | — |
| *Vibrio parahaemolyticus* | `tlh` | 1 | 90% | — |
| *Haemophilus parainfluenzae* | `ribf-hpi`, `sera-hpi` | 1 | 90% | require `ribf-hpi` |
| *Listeria monocytogenes* | `hly` | 1 | 90% | — |
| *Vibrio cholerae* | `ompw` | 1 | 85% | — |
| *Vibrio vulnificus* | `vvha` | 1 | 90% | — |
| *Staphylococcus aureus* | `nuc` | 1 | 90% | — |
| *Klebsiella pneumoniae* | `khe` | 1 | 90% | — |
| *Acinetobacter baumannii* | `oxa51` | 1 | 90% | — |
| *Pseudomonas aeruginosa* | `ecfx` | 1 | 90% | — |
| *Enterococcus faecalis* | `ddl-ef` | 1 | 90% | — |
| *Enterococcus faecium* | `ddl-efm` | 1 | 90% | — |
| *Streptococcus pneumoniae* | `lyta5`, `lyta3` | 1 | 90% | — |
| *Streptococcus pyogenes* | `speb` | 1 | 90% | — |
| *Streptococcus agalactiae* | `cfb` | 1 | 90% | — |
| *Neisseria meningitidis* | `sodc`, `ctra` | 1 | 90% | require `sodc` |
| *Neisseria gonorrhoeae* | `pora`, `porb1b` | 1 | 90% | — |
| *Clostridioides difficile* | `muri-cdf`, `isph-cdf`, `tcda`, `tcdb` | 1 | 90% | require `muri-cdf` |
| *Cronobacter sakazakii* | `rpob-cs` | 1 | 95% | exclude `speb` |
| *Legionella pneumophila* | `mip` | 1 | 90% | — |
| *Mycoplasma pneumoniae* | `tuf-mp`, `gyrb-mp`, `p1a`, `p1c` | 1 | 90% | require `tuf-mp` |
| *Helicobacter pylori* | `urea`, `urec` | 1 | 90% | — |
| *Bordetella pertussis* | `cupin-bp`, `ptxs1`, `tcfa5`, `tcfa3` | 2 | 90% | require `cupin-bp` |
| *Corynebacterium diphtheriae* | `rela-cd`, `mena-cd`, `tox` | 1 | 90% | require `rela-cd` |
| *Haemophilus influenzae* | `hpd` | 1 | 90% | — |
| *Streptococcus suis* | `gdh` | 1 | 90% | — |
| *Burkholderia pseudomallei* | `glgb-bp`, `bima-bp`, `sctn-bp` | 1 | 95% | require `glgb-bp` |
| *Leptospira* | `lipl32` | 1 | 90% | — |
| *Treponema pallidum* | `tpp47` | 1 | 90% | — |
| *Bacillus anthracis* | `adk`, `ba5345`, `paga`, `capb` | 2 | 90% | require `adk`, `ba5345`; priority over *Bacillus cereus group* |
| *Bacillus cereus group* | `adk` | 1 | 90% | — |
| *Yersinia pestis* | `spi-yp`, `btub-yp`, `caf1`, `pla` | 1 | 90% | require `spi-yp` |
| *Yersinia enterocolitica* | `glk-ye`, `yddg-ye`, `ail`, `yada` | 1 | 90% | require `glk-ye` |
| *Clostridium botulinum* | `bonta`, `bontb`, `bonte`, `bontc`, `bontd`, `bontf` | 1 | 80% | — |
| *Clostridium perfringens* | `cpa` | 1 | 90% | — |
| *Escherichia coli* | `uida`, `lacy` | 2 | 95% | — |
| *Shigella flexneri* | `uida`, `ipah`, `rfc`, `rfbe`, `rfbi` | 3 | 95% | require `uida`, `ipah`; exclude `lacy`; priority over *Shigella*, *Escherichia coli* |
| *Shigella sonnei* | `uida`, `ipah`, `wbgz`, `sonnei_psuk` | 3 | 95% | require `uida`, `ipah`; exclude `lacy`; priority over *Shigella*, *Escherichia coli* |
| *Shigella dysenteriae* | `uida`, `ipah`, `rfpb`, `sd10_wzx`, `sd10_wzy`, `sd12_wzx`, `sd12_wzy`, `sd13_wzx`, `sd13_wzy`, `sd14_wzx`, `sd14_wzy`, `sd15_wzx`, `sd15_wzy`, `sd1_rfp`, `sd1_wzx`, `sd1_wzy`, `sd2_wzx`, `sd2_wzy`, `sd3_wzx`, `sd3_wzy`, `sd4_wzx`, `sd4_wzy`, `sd6_wzx`, `sd6_wzy`, `sd8_wzx`, `sd8_wzy`, `sd9_wzx`, `sd9_wzy`, `sdprove_wzx`, `sdprove_wzy`, `sdprov_wzx`, `sdprov_wzy` | 3 | 85% | require `uida`, `ipah`; exclude `lacy`; priority over *Shigella*, *Escherichia coli* |
| *Shigella boydii* | `uida`, `ipah`, `sb11_wzx`, `sb11_wzy`, `sb12_wzx`, `sb12_wzy`, `sb14_wzx`, `sb14_wzy`, `sb15_wzx`, `sb15_wzy`, `sb16_wzy`, `sb17_wzx`, `sb17_wzy`, `sb18_wzx`, `sb18_wzy`, `sb19_wzx`, `sb19_wzy`, `sb1_wzx`, `sb1_wzy`, `sb2_wzx`, `sb2_wzy`, `sb3_wzx`, `sb3_wzy`, `sb4_wzx`, `sb4_wzy`, `sb5_wzx`, `sb5_wzy`, `sb6_wzx`, `sb6_wzy`, `sb7_wzx`, `sb7_wzy`, `sb8_wzx`, `sb8_wzy`, `sb9_wzx`, `sb9_wzy`, `sbprov_wzx`, `sbprov_wzy` | 3 | 85% | require `uida`, `ipah`; exclude `lacy`; priority over *Shigella*, *Escherichia coli* |
| *Shigella* | `uida`, `ipah` | 2 | 95% | exclude `lacy`; priority over *Escherichia coli* |
| *Campylobacter jejuni* | `mapa`, `hipo`, `cadf` | 2 | 90% | — |
| *Campylobacter coli* | `ceue` | 1 | 85% | exclude `hipo`, `cadf` |
| *Mycobacterium tuberculosis* | `is6110`, `is1081` | 1 | 96% | — |
| *Citrobacter* | `gyrb-cf`, `leus-citrobacter`, `cdsa-citrobacter` | 1 | 91% | — |
| *Chlamydia pneumoniae* | `ompa-cp`, `ompa-cpn2` | 1 | 94% | — |

## Customizing

Edit `data/db/D1_marker/marker_rules.yaml` (fields documented in
[Marker Rules](../reference/marker-rules.md)), then rebuild the BLAST database
after changing `markers.fasta`:

```bash
makeblastdb -in data/db/D1_marker/markers.fasta -dbtype nucl -out data/db/D1_marker/markers_blastdb
```
