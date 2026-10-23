# Marker Gene Database

The marker database (`data/db/D1_marker/`, the bundled **mini** tier) powers the
default `-m marker` identification: `blastn` maps query contigs against
73 curated reference sequences, and 42 species rules turn the hits into a
verdict. gside is a species-identification engine: every reference is anchored
to a verifiable NCBI/RefSeq accession and was validated by blast against the
291-genome panel. Sub-species typing is gapit's job — the V. cholerae
O1/O139 serogroup database lives there (`gapit screen --db vc_serogroup`).

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
| `primary` | Species-defining marker; a rule's `genes` draw from these | 113 |
| `confirm` | Corroborating evidence reported alongside the verdict | 23 |

## Gene inventory (139 sequences)

| Gene | Role | Reference |
|---|---|---|
| `inva` | primary | M90846.1 invasion protein A [Salmonella enterica] role=primary |
| `rpoD` | confirm | M14427.1 rpoD operon region (dnaG-rpoD) [Salmonella enterica] role=confirm |
| `iroB` | confirm | PV799927.1 glycosyltransferase iroB [Salmonella enterica] role=confirm |
| `uida` | primary | NC_000913.3:1694260-1696071 beta-D-glucuronidase uidA [Escherichia coli] role=primary |
| `lacy` | confirm | NC_000913.3:361918-362218 lactose permease lacY 5' region [Escherichia coli] role=confirm |
| `ipah` | primary | NC_004851.1:182187-182294 invasion plasmid antigen H ipaH conserved amplicon [Shigella flexneri] role=primary |
| `tlh` | confirm | M36437.1 thermolabile hemolysin tlh [Vibrio parahaemolyticus] role=confirm |
| `vvha` | primary | PQ308339.1 hemolysin vvhA [Vibrio vulnificus] role=primary |
| `hly` | primary | M24199.1 listeriolysin O hly [Listeria monocytogenes] role=primary |
| `prs` | confirm | LC700425.1 phosphoribosyl pyrophosphate synthetase prs [Listeria monocytogenes] role=confirm |
| `inlj` | confirm | AL592102.1 internalin J inlJ [Listeria monocytogenes] role=confirm |
| `nuc` | primary | V01281.1 thermonuclease nuc [Staphylococcus aureus] role=primary |
| `fema` | confirm | X17688.1 methicillin resistance femA [Staphylococcus aureus] role=confirm |
| `khe` | primary | PZ143400.1 hemolysin khe [Klebsiella pneumoniae] role=primary |
| `gyrb-kpn` | confirm | PZ127011.1 DNA gyrase subunit B gyrB [Klebsiella pneumoniae] role=confirm |
| `oxa51` | primary | MZ351777.1 beta-lactamase blaOXA-51 [Acinetobacter baumannii] role=primary |
| `rpob-aba` | confirm | QB079536.1 RNA polymerase beta rpoB [Acinetobacter baumannii] role=confirm |
| `ecfx` | primary | NC_002516.2:169361-169906 ECF subfamily sigma-70 factor ecfX/PA0149 [Pseudomonas aeruginosa] role=primary |
| `gyrb-pae` | confirm | PZ521782.1 DNA gyrase subunit B gyrB [Pseudomonas aeruginosa] role=confirm |
| `ddl-ef` | primary | QB064932.1 D-alanine-D-alanine ligase ddl [Enterococcus faecalis] role=primary |
| `ddl-efm` | primary | PZ904265.1 D-alanine-D-alanine ligase ddl [Enterococcus faecium] role=primary |
| `lyta` | primary | PP971805.1 autolysin lytA [Streptococcus pneumoniae] role=primary |
| `psaa` | confirm | AF030366.1 manganese transport psaA [Streptococcus pneumoniae] role=confirm |
| `cfb` | primary | PP475782.1 CAMP factor cfb [Streptococcus agalactiae] role=primary |
| `ctra` | primary | AY336258.1 outer membrane protein ctrA [Neisseria meningitidis] role=primary |
| `sodc` | confirm | PP471286.1 superoxide dismutase sodC [Neisseria meningitidis] role=confirm |
| `pora` | primary | M21289.1 porin protein A porA [Neisseria gonorrhoeae] role=primary |
| `tcda` | primary | M30307.1 toxin A tcdA [Clostridioides difficile] role=primary |
| `tcdb` | confirm | PQ863731.1 toxin B tcdB [Clostridioides difficile] role=confirm |
| `rpob-cs` | primary | MW380341.1 RNA polymerase beta rpoB [Cronobacter sakazakii] role=primary |
| `aera` | primary | PQ871546.1 aerolysin aerA [Aeromonas hydrophila] role=primary |
| `ahh1` | confirm | PZ301496.1 hemolysin ahh1 [Aeromonas hydrophila] role=confirm |
| `ail` | primary | M29945.1 attachment invasion locus ail [Yersinia enterocolitica] role=primary |
| `yada` | confirm | X13881.1 Yersinia adhesin yadA [Yersinia enterocolitica] role=confirm |
| `mip` | primary | PX526202.1 macrophage infectivity potentiator mip [Legionella pneumophila] role=primary |
| `dota` | confirm | MW092166.1 dotA [Legionella pneumophila] role=confirm |
| `cards` | confirm | AF390408.1 CARDS toxin [Mycoplasma pneumoniae] role=confirm |
| `ompa-cp` | primary | Z31593.1 major outer membrane protein [Chlamydia pneumoniae] role=primary |
| `ompa-cps` | primary | PZ860298.1 major outer membrane protein ompA [Chlamydia psittaci] role=primary |
| `urea` | primary | M60398.1 urease subunit alpha [Helicobacter pylori] role=primary |
| `urec` | confirm | M60398.1 urease subunit gamma/ureC [Helicobacter pylori] role=confirm |
| `ptxs1` | confirm | NZ_CP025371.1:84308-85117 pertussis toxin S1 subunit ptxA [Bordetella pertussis] role=confirm |
| `tcfa` | confirm | NZ_CP025371.1:1274915-1276858 tracheal colonization factor tcfA [Bordetella pertussis] role=confirm |
| `tox` | primary | KX703000.1 diphtheria toxin tox [Corynebacterium diphtheriae] role=primary |
| `hpd` | primary | U32723.1 hemoglobin-binding protein D [Haemophilus influenzae] role=primary |
| `gdh` | primary | MT542692.1 glutamate dehydrogenase gdh [Streptococcus suis] role=primary |
| `bimabp` | primary | AF505187.1 bacterial actin motility BimA(Bp) [Burkholderia pseudomallei] role=primary |
| `tts1` | confirm | AY089501.1 type III secretion system 1 [Burkholderia pseudomallei] role=confirm |
| `lipl32` | primary | PZ599989.1 lipoprotein lipL32 [Leptospira interrogans] role=primary |
| `tpp47` | primary | AE000520.1:622266-623570 47 kDa membrane lipoprotein tpp47 [Treponema pallidum] role=primary |
| `gyrb-ps` | primary | AB015801.1 DNA gyrase subunit B [Plesiomonas shigelloides] role=primary |
| `eacf` | primary | CP099906.1:3631386-3634232 intimin eae [Escherichia albertii] role=primary |
| `toxrf` | primary | OR123876.1 transmembrane regulatory protein toxR [Vibrio fluvialis] role=primary |
| `bontA` | primary | M30171.1 botulinum neurotoxin type A [Clostridium botulinum] role=primary |
| `bontB` | typing | M81186.1 botulinum neurotoxin type B [Clostridium botulinum] role=typing |
| `bontE` | typing | X62089.1 botulinum neurotoxin type E [Clostridium botulinum] role=typing |
| `paga` | primary | M22504.1 protective antigen pagA [Bacillus anthracis] role=primary |
| `capb` | confirm | M64091.1 capsule biosynthesis capB [Bacillus anthracis] role=confirm |
| `caf1` | primary | X13880.1 F1 capsule antigen caf1 [Yersinia pestis] role=primary |
| `pla` | confirm | PX130171.1 plasminogen activator pla pPCP1 [Yersinia pestis] role=confirm |
| `p1` | primary | U00089.2:180858-185741 P1 adhesin operon mgpA [Mycoplasma pneumoniae] role=primary |
| `cadf` | primary | GCF_000009085.1 cadf [Campylobacter jejuni] role=primary |
| `ceue` | primary | GCF_009730395.1 ceue [Campylobacter coli] role=primary |
| `hipo` | primary | GCF_000009085.1 hipo [Campylobacter jejuni] role=primary |
| `mapa` | primary | GCF_000009085.1 mapa [Campylobacter jejuni] role=primary |
| `speb` | primary | GCF_000006785.1 speb [Streptococcus pyogenes] role=primary |
| `ompw` | primary | X51948.1 outer membrane protein W [Vibrio cholerae] role=primary |
| `cpa` | primary | PV804807.1 phospholipase C alpha toxin (cpa/plc) [Clostridium perfringens] role=primary |
| `safC` | typing | NC_003197.2:343560-346076 saf atypical fimbrial operon [Salmonella enterica] role=typing |
| `rfc` | primary | NC_004337.2:2110959-2112107 O-antigen polymerase rfc [Shigella flexneri] role=primary |
| `wbgz` | primary | CP000039.1:204752-206671 putative epimerase/dehydratase wbgZ [Shigella sonnei] role=primary |
| `rfpb` | primary | CP000640.1:2979-4112 glycosyltransferase rfpB [Shigella dysenteriae] role=primary |
| `sonnei_psuk` | primary | GCF_002950395.1:595851-596792 pseudouridine kinase [Shigella sonnei] role=primary |
| `rfbe` | primary | NC_004337.2:2113945-2115201 polysaccharide biosynthesis protein rfbE [Shigella flexneri] role=primary |
| `rfbi` | primary | NC_004337.2:2110581-2110955 glycosyl translocase rfbI [Shigella flexneri] role=primary |
| `sb11_wzx` | primary | ShigaTyper ShigellaRef5 Sb11_wzx [Shigella boydii] role=primary |
| `sb11_wzy` | primary | ShigaTyper ShigellaRef5 Sb11_wzy [Shigella boydii] role=primary |
| `sb12_wzx` | primary | ShigaTyper ShigellaRef5 Sb12_wzx [Shigella boydii] role=primary |
| `sb12_wzy` | primary | ShigaTyper ShigellaRef5 Sb12_wzy [Shigella boydii] role=primary |
| `sb14_wzx` | primary | ShigaTyper ShigellaRef5 Sb14_wzx [Shigella boydii] role=primary |
| `sb14_wzy` | primary | ShigaTyper ShigellaRef5 Sb14_wzy [Shigella boydii] role=primary |
| `sb15_wzx` | primary | ShigaTyper ShigellaRef5 Sb15_wzx [Shigella boydii] role=primary |
| `sb15_wzy` | primary | ShigaTyper ShigellaRef5 Sb15_wzy [Shigella boydii] role=primary |
| `sb16_wzy` | primary | ShigaTyper ShigellaRef5 Sb16_wzy [Shigella boydii] role=primary |
| `sb17_wzx` | primary | ShigaTyper ShigellaRef5 Sb17_wzx [Shigella boydii] role=primary |
| `sb17_wzy` | primary | ShigaTyper ShigellaRef5 Sb17_wzy [Shigella boydii] role=primary |
| `sb18_wzx` | primary | ShigaTyper ShigellaRef5 Sb18_wzx [Shigella boydii] role=primary |
| `sb18_wzy` | primary | ShigaTyper ShigellaRef5 Sb18_wzy [Shigella boydii] role=primary |
| `sb19_wzx` | primary | ShigaTyper ShigellaRef5 Sb19_wzx [Shigella boydii] role=primary |
| `sb19_wzy` | primary | ShigaTyper ShigellaRef5 Sb19_wzy [Shigella boydii] role=primary |
| `sb1_wzx` | primary | ShigaTyper ShigellaRef5 Sb1_wzx [Shigella boydii] role=primary |
| `sb1_wzy` | primary | ShigaTyper ShigellaRef5 Sb1_wzy [Shigella boydii] role=primary |
| `sb2_wzx` | primary | ShigaTyper ShigellaRef5 Sb2_wzx [Shigella boydii] role=primary |
| `sb2_wzy` | primary | ShigaTyper ShigellaRef5 Sb2_wzy [Shigella boydii] role=primary |
| `sb3_wzx` | primary | ShigaTyper ShigellaRef5 Sb3_wzx [Shigella boydii] role=primary |
| `sb3_wzy` | primary | ShigaTyper ShigellaRef5 Sb3_wzy [Shigella boydii] role=primary |
| `sb4_wzx` | primary | ShigaTyper ShigellaRef5 Sb4_wzx [Shigella boydii] role=primary |
| `sb4_wzy` | primary | ShigaTyper ShigellaRef5 Sb4_wzy [Shigella boydii] role=primary |
| `sb5_wzx` | primary | ShigaTyper ShigellaRef5 Sb5_wzx [Shigella boydii] role=primary |
| `sb5_wzy` | primary | ShigaTyper ShigellaRef5 Sb5_wzy [Shigella boydii] role=primary |
| `sb6_wzx` | primary | ShigaTyper ShigellaRef5 Sb6_wzx [Shigella boydii] role=primary |
| `sb6_wzy` | primary | ShigaTyper ShigellaRef5 Sb6_wzy [Shigella boydii] role=primary |
| `sb7_wzx` | primary | ShigaTyper ShigellaRef5 Sb7_wzx [Shigella boydii] role=primary |
| `sb7_wzy` | primary | ShigaTyper ShigellaRef5 Sb7_wzy [Shigella boydii] role=primary |
| `sb8_wzx` | primary | ShigaTyper ShigellaRef5 Sb8_wzx [Shigella boydii] role=primary |
| `sb8_wzy` | primary | ShigaTyper ShigellaRef5 Sb8_wzy [Shigella boydii] role=primary |
| `sb9_wzx` | primary | ShigaTyper ShigellaRef5 Sb9_wzx [Shigella boydii] role=primary |
| `sb9_wzy` | primary | ShigaTyper ShigellaRef5 Sb9_wzy [Shigella boydii] role=primary |
| `sbprov_wzx` | primary | ShigaTyper ShigellaRef5 SbProv_wzx [Shigella boydii] role=primary |
| `sbprov_wzy` | primary | ShigaTyper ShigellaRef5 SbProv_wzy [Shigella boydii] role=primary |
| `sd10_wzx` | primary | ShigaTyper ShigellaRef5 Sd10_wzx [Shigella dysenteriae] role=primary |
| `sd10_wzy` | primary | ShigaTyper ShigellaRef5 Sd10_wzy [Shigella dysenteriae] role=primary |
| `sd12_wzx` | primary | ShigaTyper ShigellaRef5 Sd12_wzx [Shigella dysenteriae] role=primary |
| `sd12_wzy` | primary | ShigaTyper ShigellaRef5 Sd12_wzy [Shigella dysenteriae] role=primary |
| `sd13_wzx` | primary | ShigaTyper ShigellaRef5 Sd13_wzx [Shigella dysenteriae] role=primary |
| `sd13_wzy` | primary | ShigaTyper ShigellaRef5 Sd13_wzy [Shigella dysenteriae] role=primary |
| `sd14_wzx` | primary | ShigaTyper ShigellaRef5 Sd14_wzx [Shigella dysenteriae] role=primary |
| `sd14_wzy` | primary | ShigaTyper ShigellaRef5 Sd14_wzy [Shigella dysenteriae] role=primary |
| `sd15_wzx` | primary | ShigaTyper ShigellaRef5 Sd15_wzx [Shigella dysenteriae] role=primary |
| `sd15_wzy` | primary | ShigaTyper ShigellaRef5 Sd15_wzy [Shigella dysenteriae] role=primary |
| `sd1_rfp` | primary | ShigaTyper ShigellaRef5 Sd1_rfp [Shigella dysenteriae] role=primary |
| `sd1_wzx` | primary | ShigaTyper ShigellaRef5 Sd1_wzx [Shigella dysenteriae] role=primary |
| `sd1_wzy` | primary | ShigaTyper ShigellaRef5 Sd1_wzy [Shigella dysenteriae] role=primary |
| `sd2_wzx` | primary | ShigaTyper ShigellaRef5 Sd2_wzx [Shigella dysenteriae] role=primary |
| `sd2_wzy` | primary | ShigaTyper ShigellaRef5 Sd2_wzy [Shigella dysenteriae] role=primary |
| `sd3_wzx` | primary | ShigaTyper ShigellaRef5 Sd3_wzx [Shigella dysenteriae] role=primary |
| `sd3_wzy` | primary | ShigaTyper ShigellaRef5 Sd3_wzy [Shigella dysenteriae] role=primary |
| `sd4_wzx` | primary | ShigaTyper ShigellaRef5 Sd4_wzx [Shigella dysenteriae] role=primary |
| `sd4_wzy` | primary | ShigaTyper ShigellaRef5 Sd4_wzy [Shigella dysenteriae] role=primary |
| `sd6_wzx` | primary | ShigaTyper ShigellaRef5 Sd6_wzx [Shigella dysenteriae] role=primary |
| `sd6_wzy` | primary | ShigaTyper ShigellaRef5 Sd6_wzy [Shigella dysenteriae] role=primary |
| `sd8_wzx` | primary | ShigaTyper ShigellaRef5 Sd8_wzx [Shigella dysenteriae] role=primary |
| `sd8_wzy` | primary | ShigaTyper ShigellaRef5 Sd8_wzy [Shigella dysenteriae] role=primary |
| `sd9_wzx` | primary | ShigaTyper ShigellaRef5 Sd9_wzx [Shigella dysenteriae] role=primary |
| `sd9_wzy` | primary | ShigaTyper ShigellaRef5 Sd9_wzy [Shigella dysenteriae] role=primary |
| `sdprove_wzx` | primary | ShigaTyper ShigellaRef5 SdProvE_wzx [Shigella dysenteriae] role=primary |
| `sdprove_wzy` | primary | ShigaTyper ShigellaRef5 SdProvE_wzy [Shigella dysenteriae] role=primary |
| `sdprov_wzx` | primary | ShigaTyper ShigellaRef5 SdProv_wzx [Shigella dysenteriae] role=primary |
| `sdprov_wzy` | primary | ShigaTyper ShigellaRef5 SdProv_wzy [Shigella dysenteriae] role=primary |

## Species coverage (42 rules)

| Species | Rule genes | min hits | min identity |
|---|---|---|---|
| *Salmonella enterica* | `inva` | 1 | 90% |
| *Vibrio parahaemolyticus* | `tlh` | 1 | 90% |
| *Listeria monocytogenes* | `hly` | 1 | 90% |
| *Vibrio cholerae* | `ompw` | 1 | 90% |
| *Vibrio vulnificus* | `vvha` | 1 | 90% |
| *Staphylococcus aureus* | `nuc` | 1 | 90% |
| *Klebsiella pneumoniae* | `khe` | 1 | 90% |
| *Acinetobacter baumannii* | `oxa51` | 1 | 90% |
| *Pseudomonas aeruginosa* | `ecfx` | 1 | 90% |
| *Enterococcus faecalis* | `ddl-ef` | 1 | 90% |
| *Enterococcus faecium* | `ddl-efm` | 1 | 90% |
| *Streptococcus pneumoniae* | `lyta` | 1 | 90% |
| *Streptococcus pyogenes* | `speb` | 1 | 90% |
| *Streptococcus agalactiae* | `cfb` | 1 | 90% |
| *Neisseria meningitidis* | `ctra` | 1 | 90% |
| *Neisseria gonorrhoeae* | `pora` | 1 | 90% |
| *Clostridioides difficile* | `tcda`, `tcdb` | 1 | 90% |
| *Cronobacter sakazakii* | `rpob-cs` | 1 | 95% (*excludes* `speb`) |
| *Legionella pneumophila* | `mip` | 1 | 90% |
| *Mycoplasma pneumoniae* | `p1` | 1 | 90% |
| *Helicobacter pylori* | `urea`, `urec` | 1 | 90% |
| *Bordetella pertussis* | `ptxs1`, `tcfa` | 2 | 90% |
| *Corynebacterium diphtheriae* | `tox` | 1 | 90% |
| *Haemophilus influenzae* | `hpd` | 1 | 90% |
| *Streptococcus suis* | `gdh` | 1 | 90% |
| *Burkholderia pseudomallei* | `tts1`, `bimabp` | 1 | 90% |
| *Leptospira interrogans* | `lipl32` | 1 | 90% |
| *Treponema pallidum* | `tpp47` | 1 | 90% |
| *Bacillus anthracis* | `paga`, `capb` | 2 | 90% |
| *Yersinia pestis* | `caf1`, `pla` | 2 | 90% |
| *Yersinia enterocolitica* | `ail`, `yada` | 1 | 90% |
| *Aeromonas hydrophila* | `aera`, `ahh1` | 1 | 90% |
| *Clostridium botulinum* | `bontA`, `bontB`, `bontE` | 1 | 90% |
| *Clostridium perfringens* | `cpa` | 1 | 90% |
| *Escherichia coli* | `uida`, `lacy` | 2 | 95% |
| *Shigella flexneri* | `uida`, `ipah`, `rfc`, `rfbe`, `rfbi` | 3 | 95% (*excludes* `lacy`) |
| *Shigella sonnei* | `uida`, `ipah`, `wbgz`, `sonnei_psuk` | 3 | 95% (*excludes* `lacy`) |
| *Shigella dysenteriae* | `uida`, `ipah`, `rfpb`, `sd10_wzx`, `sd10_wzy`, `sd12_wzx`, … (+26 more) | 3 | 85% (*excludes* `lacy`) |
| *Shigella boydii* | `uida`, `ipah`, `sb11_wzx`, `sb11_wzy`, `sb12_wzx`, `sb12_wzy`, … (+31 more) | 3 | 85% (*excludes* `lacy`) |
| *Shigella* | `uida`, `ipah` | 2 | 95% (*excludes* `lacy`) |
| *Campylobacter jejuni* | `mapa`, `hipo`, `cadf` | 2 | 90% |
| *Campylobacter coli* | `ceue` | 1 | 90% (*excludes* `mapa`, `hipo`, `cadf`) |

## Customizing

Edit `data/db/D1_marker/marker_rules.yaml` (fields documented in
[Marker Rules](../reference/marker-rules.md)), then rebuild the BLAST database
after changing `markers.fasta`:

```bash
makeblastdb -in data/db/D1_marker/markers.fasta -dbtype nucl -out data/db/D1_marker/markers_blastdb
```
