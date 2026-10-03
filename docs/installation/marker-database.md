# Marker Gene Database

The marker database (`data/db/D1_marker/`, the bundled **mini** tier) powers the
default `--mode marker` identification: `blastn` maps query contigs against
81 curated reference sequences, and 38 species rules turn the hits into a
verdict. It ships with the repository — no `db setup` needed — and its
fingerprint (`database.version`, a SHA-256 prefix of `markers.fasta`) is
reported in every marker result.

## How marker identification works

1. `blastn` aligns contigs against the pre-built `markers_blastdb`
   (evalue 1e-10, word size 11).
2. Hits pass the global thresholds: identity >= 85% and coverage >= 60%.
3. The best hit per gene feeds the species rules: a rule fires when enough of
   its genes (`min_hits`) match at `min_identity` or better.
4. Mean identity of matched genes >= 90% reports **high** confidence,
   otherwise **medium**; near-threshold marker calls suggest an ANI recheck.

See [Marker Rules](../reference/marker-rules.md) for the rule-file format and
[Identification Modes](../usage/modes.md) for how marker results arbitrate
against ANI methods.

## Gene roles

| Role | Meaning | Count |
|---|---|---|
| `primary` | Species-defining marker; a rule's `genes` draw from these | 46 |
| `confirm` | Corroborating evidence reported alongside the verdict | 22 |
| `typing` | Sub-typing / serovar context, not species-defining | 3 |
| `virulence` | Toxin / virulence factors, reported as findings | 10 |

## Gene inventory (81 sequences)

| Gene | Role | Reference |
|---|---|---|
| `inva` | primary | M90846.1 invasion protein A [Salmonella enterica] |
| `rpoD` | confirm | CP026417.1 RNA polymerase sigma 70 [Salmonella enterica] |
| `iroB` | confirm | AF080424.1 iroB C-glycosyltransferase [Salmonella enterica] |
| `safC` | typing | FM201984.1 Salmonella atypical fimbriae usher [Salmonella enterica] |
| `uida` | primary | NC_000913.3 beta-D-glucuronidase uidA [Escherichia coli] |
| `lacy` | primary | NC_000913.3 lactose permease lacY [Escherichia coli] |
| `gada` | confirm | NC_000913.3 glutamate decarboxylase alpha [Escherichia coli] |
| `stx1` | virulence | AF125520.1 Shiga toxin 1 subunit A [Escherichia coli] |
| `stx2` | virulence | AF125522.1 Shiga toxin 2 subunit A [Escherichia coli] |
| `eae` | virulence | AF022231.1 intimin eae [Escherichia coli] |
| `ipah` | primary | NC_004337.2 invasion plasmid antigen H ipaH [Shigella flexneri] |
| `toxr` | primary | NC_004603.1 toxR regulatory protein [Vibrio parahaemolyticus] |
| `tlh` | confirm | M36437.1 thermolabile hemolysin tlh [Vibrio parahaemolyticus] |
| `tdh` | virulence | M10069.1 thermostable direct hemolysin [Vibrio parahaemolyticus] |
| `ompw` | primary | AF055890.1 outer membrane protein W [Vibrio cholerae] |
| `ctxa` | virulence | X00171.1 cholera enterotoxin A subunit [Vibrio cholerae] |
| `vvha` | primary | M34462.1 hemolysin vvhA [Vibrio vulnificus] |
| `hly` | primary | M24199.1 listeriolysin O hly [Listeria monocytogenes] |
| `prs` | confirm | AF261779.1 putrescine transport prs [Listeria monocytogenes] |
| `inlj` | confirm | AL592102.1 internalin J inlJ [Listeria monocytogenes] |
| `nuc` | primary | V01281.1 thermonuclease nuc [Staphylococcus aureus] |
| `fema` | confirm | X17688.1 methicillin resistance femA [Staphylococcus aureus] |
| `khe` | primary | AF072244.1 hemolysin gene khe [Klebsiella pneumoniae] |
| `gyrb-kpn` | confirm | AF318700.1 DNA gyrase subunit B [Klebsiella pneumoniae] |
| `oxa51` | primary | AF300835.1 OXA-51-like carbapenemase [Acinetobacter baumannii] |
| `rpob-aba` | confirm | X82164.1 RNA polymerase beta rpoB [Acinetobacter baumannii] |
| `ecfx` | primary | AE004091.2 sigma factor ECF ecfX [Pseudomonas aeruginosa] |
| `gyrb-pae` | confirm | L05639.1 DNA gyrase subunit B [Pseudomonas aeruginosa] |
| `ddl-ef` | primary | AF181880.1 D-ala-D-ala ligase ddl [Enterococcus faecalis] |
| `ddl-efm` | primary | AF181882.1 D-ala-D-ala ligase ddl [Enterococcus faecium] |
| `lyta` | primary | M81227.1 autolysin lytA [Streptococcus pneumoniae] |
| `psaa` | confirm | AF030366.1 manganese transport psaA [Streptococcus pneumoniae] |
| `cfb` | primary | M33317.1 CAMP factor cfb [Streptococcus agalactiae] |
| `ctra` | primary | M57681.1 capsule transport ctrA [Neisseria meningitidis] |
| `sodc` | confirm | AF322864.1 superoxide dismutase C [Neisseria meningitidis] |
| `pora` | primary | M21289.1 porin protein A porA [Neisseria gonorrhoeae] |
| `tcda` | primary | M30307.1 toxin A tcdA [Clostridioides difficile] |
| `tcdb` | confirm | M19030.1 toxin B tcdB [Clostridioides difficile] |
| `rpob-cs` | primary | AF064441.1 RNA polymerase beta rpoB [Cronobacter sakazakii] |
| `aera` | primary | M64716.1 aerolysin aerA [Aeromonas hydrophila] |
| `ahh1` | confirm | S57479.1 hemolysin ahh1 [Aeromonas hydrophila] |
| `ail` | primary | M29945.1 attachment invasion locus ail [Yersinia enterocolitica] |
| `yada` | confirm | X13881.1 Yersinia adhesin yadA [Yersinia enterocolitica] |
| `nheb` | virulence | Z11886.1 non-hemolytic enterotoxin B nheB [Bacillus cereus] |
| `ces` | virulence | AF110410.1 cereulide synthetase ces [Bacillus cereus] |
| `cpa` | primary | M24905.1 alpha toxin cpa/plc [Clostridium perfringens] |
| `cpe` | virulence | L43549.1 enterotoxin cpe [Clostridium perfringens] |
| `mip` | primary | M32024.1 macrophage infectivity potentiator [Legionella pneumophila] |
| `dota` | confirm | U91654.1 DotA type IV secretion [Legionella pneumophila] |
| `cards` | confirm | AF390408.1 CARDS toxin [Mycoplasma pneumoniae] |
| `ompa-cp` | primary | Z31593.1 major outer membrane protein [Chlamydia pneumoniae] |
| `ompa-cps` | primary | AY006345.1 major outer membrane protein [Chlamydia psittaci] |
| `urea` | primary | M60398.1 urease subunit alpha [Helicobacter pylori] |
| `urec` | confirm | M60398.1 urease subunit gamma/ureC [Helicobacter pylori] |
| `caga` | virulence | AB015416.1 cytotoxin-associated gene A [Helicobacter pylori] |
| `ptxs1` | confirm | M13223.1 pertussis toxin S1 subunit [Bordetella pertussis] |
| `tox` | primary | K01722.1 diphtheria toxin [Corynebacterium diphtheriae] |
| `hpd` | primary | U32723.1 hemoglobin-binding protein D [Haemophilus influenzae] |
| `gdh` | primary | AF363735.1 glutamate dehydrogenase gdh [Streptococcus suis] |
| `bimabp` | primary | AF505187.1 bacterial actin motility BimA(Bp) [Burkholderia pseudomallei] |
| `tts1` | confirm | AY089501.1 type III secretion system 1 [Burkholderia pseudomallei] |
| `lipl32` | primary | AF037900.1 outer membrane protein LipL32 [Leptospira interrogans] |
| `tpp47` | primary | M88726.1 47kDa membrane lipoprotein [Treponema pallidum] |
| `gyrb-ps` | primary | AB015801.1 DNA gyrase subunit B [Plesiomonas shigelloides] |
| `eacf` | primary | AB074997.1 E. albertii adherence factor [Escherichia albertii] |
| `toxrf` | primary | AF449565.1 toxR transcriptional regulator [Vibrio fluvialis] |
| `bontA` | primary | M30171.1 botulinum neurotoxin type A [Clostridium botulinum] |
| `bontB` | typing | M81186.1 botulinum neurotoxin type B [Clostridium botulinum] |
| `bontE` | typing | X62089.1 botulinum neurotoxin type E [Clostridium botulinum] |
| `paga` | primary | M22504.1 protective antigen pagA [Bacillus anthracis] |
| `capb` | confirm | M64091.1 capsule biosynthesis capB [Bacillus anthracis] |
| `caf1` | primary | X13880.1 F1 capsule antigen caf1 [Yersinia pestis] |
| `pla` | confirm | M77367.1 plasminogen activator pla [Yersinia pestis] |
| `trh` | virulence | AB025286.1 TDH-related hemolysin trh [Vibrio parahaemolyticus] |
| `p1` | primary | U00089.2 P1 adhesin protein [Mycoplasma pneumoniae] |
| `is481` | primary | U15171.1 insertion sequence IS481 [Bordetella pertussis] |
| `cadf` | primary | GCF_000009085.1 cadf [Campylobacter jejuni] |
| `ceue` | primary | GCF_009730395.1 ceue [Campylobacter coli] |
| `hipo` | primary | GCF_000009085.1 hipo [Campylobacter jejuni] |
| `mapa` | primary | GCF_000009085.1 mapa [Campylobacter jejuni] |
| `speb` | primary | GCF_000006785.1 speb [Streptococcus pyogenes] |

## Species coverage (38 rules)

| Species | Rule genes | min hits | min identity |
|---|---|---|---|
| *Salmonella* | `inva` | 1 | 90% |
| *V parahaemolyticus* | `tlh` | 1 | 90% |
| *Listeria monocytogenes* | `hly` | 1 | 90% |
| *Vibrio cholerae* | `ompw, ctxa` | 1 | 90% |
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
| *Clostridioides difficile* | `tcda, tcdb` | 1 | 90% |
| *Cronobacter sakazakii* | `rpob-cs` | 1 | 90% (*excludes* `speb`) |
| *Legionella pneumophila* | `mip` | 1 | 90% |
| *Mycoplasma pneumoniae* | `p1` | 1 | 90% |
| *Helicobacter pylori* | `urea, urec` | 1 | 90% |
| *Bordetella pertussis* | `is481, ptxs1` | 1 | 90% |
| *Corynebacterium diphtheriae* | `tox` | 1 | 90% |
| *Haemophilus influenzae* | `hpd` | 1 | 90% |
| *Streptococcus suis* | `gdh` | 1 | 90% |
| *Burkholderia pseudomallei* | `tts1, bimabp` | 1 | 90% |
| *Leptospira interrogans* | `lipl32` | 1 | 90% |
| *Treponema pallidum* | `tpp47` | 1 | 90% |
| *Bacillus anthracis* | `paga, capb` | 2 | 90% |
| *Yersinia pestis* | `caf1, pla` | 2 | 90% |
| *Yersinia enterocolitica* | `ail, yada` | 1 | 90% |
| *Aeromonas hydrophila* | `aera, ahh1` | 1 | 90% |
| *Clostridium botulinum* | `bontA, bontB, bontE` | 1 | 90% |
| *Clostridium perfringens* | `cpa` | 1 | 90% |
| *DEC* | `uida, lacy, gada` | 2 | 95% |
| *Shigella EIEC* | `ipah` | 1 | 95% |
| *Campylobacter jejuni* | `mapa, hipo, cadf` | 2 | 90% |
| *Campylobacter coli* | `ceue` | 1 | 90% (*excludes* `mapa, hipo, cadf`) |

## Customizing

Edit `data/db/D1_marker/marker_rules.yaml` (fields documented in
[Marker Rules](../reference/marker-rules.md)), then rebuild the BLAST database
after changing `markers.fasta`:

```bash
makeblastdb -in data/db/D1_marker/markers.fasta -dbtype nucl -out data/db/D1_marker/markers_blastdb
```
