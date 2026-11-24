# Marker Rules Format

`marker_rules.yaml` defines target gene combination identification rules —
adding new species requires only a rule entry, no code changes.

## File Location

```
data/db/D1_marker/marker_rules.yaml
```

## Rule Structure

```yaml
rules:
  - species: <species_name>          # required: registry key (underscored)
    genes:                            # required: candidate marker gene list
      - <gene1>
      - <gene2>
    min_hits: <int>                   # optional: minimum hits (default 1)
    min_identity: <float>             # optional: minimum identity % (default 90)
    exclude_genes:                    # optional: exclusion guard
      - <cross_reactive_gene>
    min_coverage: <float>             # optional: minimum coverage % (default 60)
    note: <str>                       # optional: documentation
    co_detect_note: <str>             # optional: documentation for
                                      # co-detection conflicts
```

## Field Reference

### species

Species name using registry key format (spaces → underscores).
Examples: `Campylobacter_jejuni`, `Vibrio_parahaemolyticus`, `Escherichia_coli`.

### genes

Candidate marker gene names (lowercase). These genes must exist in
`markers.fasta`.

### min_hits

Minimum number of genes that must hit (default 1). Setting 2 means at least
2 genes from the `genes` list must reach thresholds for the species to be
called.

Multi-gene combinations significantly reduce cross-reactions (e.g.,
C.jejuni's 2-of-3 design).

### require

Genes that MUST be among the matched hits, in addition to the `min_hits`
quorum (default: none). Use it when a rule mixes a genus signal with
species-specific markers — e.g. `Salmonella_enterica` = `require: [inva]`
plus at least one SPI-2 marker — so the specific markers alone can never
fire the species call without the genus anchor. Every `require` entry must
also appear in `genes` (validated by `gside db validate`).

### min_identity

Minimum identity percentage per gene (blastn pident). Default 90.

Lower to 85 for short sequences or conserved genes (but increases false
positive risk).

### exclude_genes

**Exclusion guard**: if these genes are detected in the sample, the species
call is suppressed.

Typical usage — Campylobacter species discrimination:

```yaml
- species: Campylobacter_coli
  genes: [ceue]
  exclude_genes: [mapa, hipo, cadf]    # jejuni markers present → prefer jejuni
- species: Campylobacter_jejuni
  genes: [mapa, hipo, cadf]
  min_hits: 2
  # no exclusion: jejuni genomes naturally contain ceuE homologs
```

### min_coverage

Minimum query coverage percentage for a gene to count toward this rule
(default 60, matching the global gate). The check runs per rule against
every raw hit of the gene, so a lower-identity hit with better coverage can
still satisfy the rule.

Set it lower (e.g. 50) only for rules where partial gene matches are
accepted; those two *Campylobacter* rules do exactly that.

### note / co_detect_note

Free-text documentation fields. They carry no runtime effect. Use `note` for
general context and `co_detect_note` to record how to interpret the rule
when its markers co-occur with another species' markers (e.g., Shigella vs
Escherichia_coli).

## Global Thresholds

The engine layer defines global gates (in `multigene_identifier.py`):

| Parameter | Value | Description |
|-----------|-------|-------------|
| `_MIN_IDENTITY` | 85.0 | Global minimum identity |
| `_MIN_COVERAGE` | 60.0 | Global minimum coverage (prevents short-HSP false fires) |
| `_HIGH_CONF` | 90.0 | High confidence threshold (average identity of matched genes) |

Rule-level `min_identity` / `min_coverage` apply to each gene's raw hits;
rules that omit them use 90 / 60. The global `_MIN_IDENTITY` (85) and
`_MIN_COVERAGE` (60) still gate `detected_markers` reporting and the
`exclude_genes` guard.

`priority_over` is supported by the engine (a matching priority species
suppresses non-priority ones). The four Shigella species rules
(`Shigella_flexneri`, `Shigella_sonnei`, `Shigella_dysenteriae`,
`Shigella_boydii`) use it with `priority_over: [Shigella, Escherichia_coli]`
and `require: [uida, ipah]`: a species determinant (e.g. `rfc`, `wbgZ`,
`sb1_wzx`) upgrades the genus signature — and the O-antigen determinants
alone can never fire the species call without the uidA+ipaH anchors — while
the plain `uidA+ipaH` signature stays at the `Shigella` genus and every
Shigella rule is lacY-negative (`exclude_genes: [lacy]`) so EIEC never fires
a Shigella call.

## CPIN coverage block

The rules file carries a marked block covering the `data/CPIN.txt` species
list. Its entries follow two patterns:

- **Classical targets** extracted from verified RefSeq reference genomes by
  `locus_tag` (IS6110, TUL4/FopA, invasin+VirF, BimA(Bm), ompA, hugA,
  BCSP31-family genes, …). Rule notes record the reference strain and the
  measured nearest cross-identity; thresholds sit above the measured cross.
- **Barcodes** (`gyrB` for `Citrobacter`/`Escherichia_albertii`,
  `rpoB+gyrB` for `Vibrio_fluvialis`/`Burkholderia_gladioli`,
  `gltA+rpoB` for `Bartonella`) where classical targets do not separate
  neighbours: Enterobacteriaceae `rpoB` cross-identity reaches 97%, so
  `gyrB` (≤89.5% cross) carries those rules instead.
- **Fungal ITS barcodes** (`its-*`, ITS1/ITS1-5.8S-ITS2 from type or
  clinical-strain records) provide genus-level calls at 92% identity — the
  Onygenales genera cross at up to 90.5%.

Every rule was validated end-to-end against its reference genome plus
near-neighbour adversarial genomes (`Y. pestis` vs `Y. pseudotuberculosis`,
`B. mallei` vs `B. pseudomallei`); see
`tests/integration/test_cpin_markers.py`.

The `B. mallei` BimA(Bm) rule sits at 98% identity: 3/100 clinical
*B. pseudomallei* genomes carry a Bm-like allele at 96.2–96.4%, while true
*B. mallei* measures 99.6–100%. At the old 95% threshold those three
pseudomallei genomes were suppressed by `priority_over` and reported as
`Burkholderia_mallei`. The two taxa are one genomic species (ANI ~99.5%);
see the *Burkholderia pseudomallei* entry below.

### Genus/family/MTBC barcode validation (100 random genomes)

The CPIN genus/family rules were re-tested against 100 random genomes per
taxon (RefSeq/GenBank, one 1-2 gene rule each). The original thresholds were
tuned on a single reference strain and missed most of the genus:

| Rule | Old design | Sensitivity | New design | Sensitivity |
|---|---|---|---|---|
| `Bartonella` | gltA + rpoB, both ≥92 | 6/97 (6%) | both ≥85 | 92/97 (95%) |
| `Nocardia` | hsp65 + gyrB, both ≥90 | 29/98 (30%) | hsp65 or gyrB, ≥91 | 87/98 (89%) |
| `Rickettsiaceae` | gltA + OmpB, both ≥90 | 48/100 (48%) | + 16S rRNA, any ≥90 | 85/100 (85%) |
| `Anaplasmataceae` | GroEL ≥90 | 1/99 (1%) | + 16S rRNA, any ≥85 | 96/99 (97%) |

The added 16S barcodes are a 1427bp *Rickettsia rickettsii* fragment
(PV716478.1) and a 1435bp *Ehrlichia chaffeensis* fragment (AF147752.2); 16S
is the standard family-level marker when genus-level genes (`gltA`, `ompB`,
`groEL`) are absent or divergent. The old gltA/OmpB/GroEL markers remain as
specific confirmatory reporters.

Specificity was checked against the other random panels: the 16S barcodes
gave zero hits in the reciprocal family and in *Brucella* (99), *Bartonella*
(97), *Citrobacter* (97) and *Nocardia* (98) at their thresholds. The
*Nocardia* 91% floor sits above the highest full-length *Mycobacterium*
cross (90.5%) and *Rhodococcus* cross (89.4%). *Bartonella* gltA/rpoB give
no hits in *Brucella*/Ochrobactrum at 85%.

Residual limits:

- *Brucella*: bp26+omp25 calls 79/99; all misses are the former
  *Ochrobactrum* species (*B. anthropi*, *B. intermedia*, *B. pituitosa*,
  *B. pseudogrignonensis*), while the classical species (*B. abortus*,
  *B. melitensis*, *B. suis*, *B. canis*, *B. ovis*, *B. ceti*) are 100%.
- *Citrobacter*: gyrB alone calls 90/100; leuS (leucine--tRNA ligase) lifts
  the union to 94/100, and cdsA (phosphatidate cytidylyltransferase) covers
  the divergent *C. amalonaticus*/*C. farmeri* clade at 95-100%, giving
  **100/100**. Negative controls stay at 68-86% across *E. coli*,
  *Salmonella*, *Klebsiella*, *Cronobacter*, *Enterobacter*, *Serratia*,
  *Proteus* and *Yersinia*.
- *M. tuberculosis*: IS6110 calls 95/100 random MTBC genomes and IS1081
  96/100; the union is **99/100**. The 43-genome non-TB *Mycobacterium*
  panel is negative except *M. orygis* and *M. canetti* (both MTBC members).
  One fragmented *M. tuberculosis* draft lacks both IS elements.
- *Rickettsiaceae*: the misses are mostly uncultured/endosymbiont genomes and
  *Orientia*; clinical *Rickettsia* (spotted-fever and typhus groups) are
  covered. *Anaplasmataceae* misses are *Wolbachia* and *Neorickettsia*
  outliers, not the clinical *Anaplasma*/*Ehrlichia* species.

### Fungal ITS barcode survey

Whole-genome FTP downloads for the fungal genera were too slow for a
100-genome panel (NCBI throttled the ~20-40Mb assemblies). The existing ITS
markers were therefore tested against 100 NCBI Nucleotide ITS records per
genus (full-length subset ≥500bp, marker identity ≥92% with ≥60% marker
coverage):

| Marker | Full-length records | ≥92% | Nearest cross at 92% |
|---|---|---|---|
| *Histoplasma* | 87 | 84 | — |
| *Blastomyces* | 66 | 66 | — |
| *Coccidioides* | 54 | 54 | — |
| *Paracoccidioides* | 26 | 25 | — |
| *Talaromyces marneffei* | 6 | 6 | — |
| *Rhinocladiella mackenziei* | 19 | 19 | — |
| *Cladophialophora bantiana* | 95 | 90 | — |
| *Sporothrix* | 96 | 94 | *Ophiostoma* 0/86 |
| *Cryptococcus* | 53 | 36 | *Rhodotorula* 0/90 |
| *Candida* | 33 | 17 | *Saccharomyces/Pichia* 0 |
| *Aspergillus* | 92 | 32 | *Penicillium* 4/87 at 92, 85/87 at 85 |
| *Fusarium* | 71 | 25 | *Trichoderma* 0/94 at 88, 47/94 at 85 |
| *Pneumocystis* | 59 | 9 | — |

The single-gene ITS rule is already adequate for the endemic dimorphic fungi
and for *Cladophialophora*. *Sporothrix* already works at 92% (94/96);
*Cryptococcus* moved to 88% (36/53 -> 46/53) without crossing *Rhodotorula*;
*Candida* moved to 85% in the tested neighbour set (17/33 -> 29/33, zero
*Saccharomyces*/*Pichia*). *Aspergillus* and *Fusarium* cannot lower their
ITS threshold: at 85% they cross *Penicillium* (85/87) and *Trichoderma*
(47/94). They now carry a second locus:

- *Aspergillus* + `bena-aspergillus` (β-tubulin, 506bp; PZ612096.1): 0/199
  *Penicillium*+*Talaromyces* records at 85-95%.
- *Fusarium* + `tef1-fusarium` (TEF1-α, 683bp; QB150418.1): 35/92 *Fusarium*
  records at 92% with 0/87 *Trichoderma* at 85-95%.
- *Pneumocystis* + `mtlsu-pneumocystis` (mtLSU rRNA, 255bp; PZ378845.1):
  42/100 records, retained alongside the weak ITS marker.

These are Nucleotide-record estimates, not whole-genome panels; the
Cryptococcus/Candida threshold changes and the new second markers should
still be confirmed on assembled genomes.

## Reads-Mode Depth Gate (FASTQ input)

On FASTQ input the marker engine applies a gspec-style abundance gate on top
of the identity/coverage rules. Per gene the screening reports `mean_depth`
(summed aligned target-interval length over gene length — not a read count),
`reads`, and `mapq_mean`; alignments below `min_mapq` (default 20) are
dropped before aggregation.

A marker passes when:

```
mean_depth >= min_depth  and  mean_depth / control_depth >= min_relative_depth
```

Defaults (in `metadata.yaml`, overridable per call via `min_depth` /
`min_relative_depth` / `min_mapq`):

| Parameter | Default | Description |
|---|---|---|
| `min_depth` | 5.0 | absolute mean-depth floor |
| `min_relative_depth` | 0.10 | fraction of the control depth |
| `min_mapq` | 20 | minimap2 MAPQ floor applied before aggregation |

`control_depth` is the median mean depth of the species' configured
single-copy control genes (`metadata.yaml: control_genes`, e.g. `rpoD` for
*Salmonella*, two `lpxA` alleles for *Campylobacter*). The deepest allele per
gene is used so a few contaminant reads cannot cap the control. Control genes
carry `role=control` in `markers.fasta` and are never shown in
`detected_markers`. Species without a dedicated control fall back to the
sample median depth of detected markers (`control_source: "fallback"`);
if no control is available at all, the relative gate is disabled and only the
absolute floor applies. Assembly (BLAST) input has no depth and keeps the
legacy behaviour.

## Cross-species co-detection and arbitration {: #cross-species-arbitration }

A single genome can produce marker hits from more than one species even when
it is a pure isolate. The usual causes are shared ancestry (conserved
housekeeping genes), horizontal transfer (virulence plasmids, pathogenicity
islands, prophages), multi-copy paralogs, and mixed or contaminated
assemblies. Marker hits are evidence, not a conclusion; the engine arbitrates
in a fixed order:

1. **Hit gates** — a hit is discarded below the global floor (85% identity,
   60% coverage). Each rule then applies its own `min_identity` (default 90)
   and `min_coverage` (default 60) to every raw hit of its genes.
2. **Rule quorum** — `min_hits` sets how many of the rule's genes must pass;
   `require` genes must all be among them (hard AND). A rule that cannot
   reach its quorum is not a candidate.
3. **Exclusion guards** — `exclude_genes` suppresses a rule when a
   cross-reactive gene is detected at the global floor, even if the rule's
   own markers pass (e.g. *Campylobacter coli* is suppressed by hipo/cadf;
   every Shigella species rule is lacY-negative, so EIEC stays *E. coli*).
4. **Priority** — when a candidate rule declares `priority_over`, it wins
   over the listed species regardless of score (*B. mallei* BimA(Bm) over
   *B. pseudomallei*; Shigella species over Shigella/*E. coli*;
   *S. bongori* over *S. enterica*). Priority is the documented way to
   resolve shared-marker competition.
5. **Score** — among the remaining candidates, the rule with the most
   matched genes wins; a tie prefers high confidence (average identity
   ≥90%), and a remaining tie is broken by rule order in the YAML.

The output keeps the evidence visible: `detected_markers` lists every
significant gene, `matched_rule` names the winning gene combination, and
`notes` / `co_detect_note` explain mixed profiles. A profile with, say,
*B. anthracis* plasmid genes but no chromosomal marker reports the cereus
group plus an "Anthrax-like" note rather than forcing the species call.

When marker, ANI and GTDB-Tk results disagree (`all` mode), arbitration is
layer-first: GTDB-Tk (L3) > ANI/mash/sourmash (L2) > marker (L1); within a
layer, high confidence beats medium, and `verdict.supporting` /
`verdict.agreement` record cross-method corroboration. A marker hit never
overrides a whole-genome ANI call, and a whole-genome ANI call never hides a
marker profile — inspect both when they disagree.

Interpretation guidance:

- A single off-target marker with no rule quorum is background, not a mixed
  sample.
- Multiple virulence markers from a foreign species suggest horizontal
  transfer or a pathotype hybrid; confirm with L2.
- Broad, high-identity hits across unrelated genera suggest contamination or
  a mixed assembly; re-check the input.

## Species-level Salmonella rules

The bundled database distinguishes the two *Salmonella* species with three
layered rules (empirical basis: 26 *Salmonella* genomes across all six
*S. enterica* subspecies plus 13 *S. bongori* strains, and five outgroup
genomes; background in Fookes 2011, PLoS Pathog e1002191):

| Rule | Genes | Threshold | Meaning |
|---|---|---|---|
| `Salmonella_enterica` | `inva` (required) + (`ssav` \| `ssec`), exclude `asst`/`zadh` | 2 of 3 genes, ≥ 90 | genus anchor AND a SPI-2 marker (T3SS-2 apparatus `ssav`/translocon `ssec`, core in all 11 *S. enterica* panel genomes across 6 subspecies, absent from every *S. bongori* strain and outgroup); genus-level threshold keeps divergent enterica lineages species-called |
| `Salmonella_bongori` | `inva` (required) + (`asst` \| `zadh`) | 2 of 3 genes, ≥ 90 | genus anchor AND a bongori-specific marker (`asst` = aryl-sulfate sulfotransferase, `zadh` = zinc-dependent alcohol dehydrogenase) |
| `Salmonella` | `inva` | identity ≥ 90 | genus fallback — invA-positive genomes where no species rule fires (bongori drafts missing both markers; SPI-2-negative enterica; markers in the 85–90 detection band) |

Both species rules `require: inva` (the genus signal must be among the
matched hits) and `Salmonella_bongori` carries
`priority_over: [Salmonella, Salmonella_enterica]`, so specific markers
alone never fire a species call, and a marker-carrying bongori genome always
resolves correctly. Known trade-off (deliberate): *S. enterica* isolates
lacking SPI-2 (rare, historic subsp. arizonae) fall back to the genus label.
Serovar/subspecies typing stays with gspec.

## Bordetella pertussis: IS fragmentation-resistant design

B. pertussis carries hundreds of IS481 copies, so draft assemblies fragment
badly (700+ contigs are routine; 30% of RefSeq drafts in a 20-genome panel
split tcfA across two contigs). The rule design absorbs that:

| Gene | Role |
|---|---|
| `cupin-bp` (681bp, required) | chromosomal species anchor — 100% conserved across 20 diverse strains (incl. ptxP3 era), max 71% identity in near relatives (ptxS1 cross-reacts at 96–99% with B. parapertussis/holmesii/bronchiseptica and IS481 is shared with B. holmesii, so neither can anchor) |
| `tcfa5` / `tcfa3` (~972bp each) | the 1944bp tcfA split in halves — each half survives on a single contig when an IS breakpoint fragments the gene |
| `ptxs1` (810bp) | pertussis toxin S1 — virulence-locus reporting |

Virulence in B. pertussis is chromosomal (no plasmid issue as in B. anthracis):
Bvg⁻ phase variants and detoxified vaccine derivatives remain detectable —
the historical misses were IS-fragmentation artifacts, not avirulence.

## Yersinia pestis: chromosomal anchor + plasmid virulence reporting

Y. pestis is a recently emerged clone of Y. pseudotuberculosis (chromosomal
identity >97%), making chromosomal discrimination the hardest in the genus.
A 50-strain real-world panel (spanning all three biovars) against 5
near-relative genomes found only 4 chromosomal genes specific to Y. pestis —
the species is extremely clonal (all candidates 100% conserved):

| Gene | Role |
|---|---|
| `spi-yp` (987bp, required) | chromosomal sugar phosphate isomerase — species anchor, 100% conserved across 50 strains, zero homology in near relatives |
| `btub-yp` (1878bp) | TonB-dependent BtuB — second chromosomal marker (literature precedent) |
| `caf1` (pMT1) | F1 capsule antigen — virulence plasmid reporting |
| `pla` (pPCP1) | plasminogen activator — virulence reporting; 13% of real-world strains have lost the small pPCP1 plasmid |

The original `min_hits: 2` on `[caf1, pla]` made every pPCP1-negative strain
Unknown. With the chromosomal anchor, plasmid loss no longer gates the
species call.

## Yersinia enterocolitica: species anchor + pathogenicity reporting

The original `ail`/`yada` entries were corrupted sequences matching no genome
(the rule never fired — 76/100 real-world genomes reported Unknown). The
redesigned rule separates species identity from pathogenicity:

| Gene | Role |
|---|---|
| `glk-ye` (1164bp, required) | glycerate kinase chromosomal anchor — 100% present at ≥96% across a 50-strain real-world panel (spanning subsp. enterocolitica/palearctica and biotype 1A), zero homology in Y. pseudotuberculosis/frederiksenii/intermedia/kristensenii/pestis and E. coli |
| `yddg-ye` (888bp) | second chromosomal marker (redundancy rider) |
| `ail` (537bp) | real AIL CDS — present only in pathogenic bioserotypes (1B/2–5); absent from non-pathogenic biotype 1A |
| `yada` (1269bp) | real pYV-plasmid YadA — virulence reporting; the plasmid is easily cured |

Non-pathogenic biotype 1A and pYV-cured strains therefore still report
`Yersinia_enterocolitica` (anchor-only) with a note flagging the absence of
pathogenicity markers.

## Bacillus cereus group: species vs complex

The cereus group (B. anthracis / B. cereus / B. thuringiensis / B. mycoides /
B. weihenstephanensis / ...) shares >99% chromosomal identity; only the
B. anthracis clade is chromosomally separable, via the lambdaBa03 prophage.
The bundled rules mirror that biology (empirical basis: 404-genome panel —
150 B. anthracis + 246 non-anthracis cereus group + 8 outgroups; ring-trial
literature: Wielinga 2013, Bentahir 2018):

| Rule | Genes | Threshold | Meaning |
|---|---|---|---|
| `Bacillus_anthracis` | `adk` (required) + `ba5345` (required), `paga`/`capb` ride along | 2 of 4, ≥ 90 | group anchor AND the lambdaBa03 chromosomal marker; the plasmid genes report virulence status but never gate the species call — plasmid-cured/vaccine strains (Sterne) stay anthracis |
| `Bacillus_cereus_group` | `adk` | ≥ 90 | complex-level fallback: B. cereus, B. thuringiensis, other group members, and lambdaBa03-negative B. anthracis drafts (the prophage assembles poorly from short reads — 27% of RefSeq drafts) |

When `paga`/`capb` fire without `ba5345` (the B. cereus biovar anthracis /
cross-over profile, e.g. Bcbva CI, BcFL2013, 03BB87), the result carries an
"Anthrax-like" note directing confirmation at a reference laboratory — the
plasmid pair alone never fires `Bacillus_anthracis`. `verdict.scope` for the
group label is `complex`.

## Species Marker Reference (with pathogen biology)

Each rule's marker set reflects the underlying pathogen biology. The tables
below explain what each gene measures and why it was chosen, so users can
interpret `detected_markers` and `notes` in the output with full context.

### Salmonella (S. enterica vs S. bongori)

| Gene | Role | Biology |
|---|---|---|
| `inva` | Genus anchor (required) | Invasion protein A; conserved across the genus (invA⁻ mutants are non-invasive lab constructs). Identity band separates species: ≥93% = S. enterica, 90–93% = S. bongori (verified 26-genome panel). |
| `ssav` | S. enterica specific (required pair) | SPI-2 type III secretion system apparatus. SPI-2 is the intracellular survival island, absent from every S. bongori strain and outgroup screened. |
| `ssec` | S. enterica specific (required pair) | SPI-2 translocon protein. Same island, independent gene for redundancy. |
| `asst` | S. bongori specific | Aryl-sulfate sulfotransferase, chromosomal. Zero full-length homologs in S. enterica or outgroups. |
| `zadh` | S. bongori specific | Zinc-dependent alcohol dehydrogenase. Same species specificity as asst. |

**Interpretation**: `Salmonella_enterica` = genus signal + species-specific
SPI-2 marker. `Salmonella_bongori` = genus signal + bongori metabolic marker.
Bare `Salmonella` = invA-positive but neither species rule fired (partial
assembly, divergent lineage). Serovar typing (Typhi, Enteritidis, etc.)
belongs to **gspec**, not gside.

### Bacillus cereus group

| Gene | Role | Biology |
|---|---|---|
| `adk` | Group anchor (required for B. anthracis) | Adenylate kinase, conserved across the entire cereus group (>99% identity). Distinguishes the group from other Bacillus and outgroups. |
| `ba5345` | B. anthracis specific (required) | λBa03 prophage chromosomal marker; the only validated chromosomal separator of B. anthracis from the cereus group (Wielinga 2013 ring trial, 134-genome in silico). |
| `paga` | Virulence reporting (pXO1 plasmid) | Protective antigen, anthrax toxin component. Reports toxigenicity; plasmid-cured/vaccine strains (Sterne) remain B. anthracis without it. |
| `capb` | Virulence reporting (pXO2 plasmid) | Capsule biosynthesis. Full virulence requires both plasmids. |

**Interpretation**: `Bacillus_anthracis` = chromosomal identity confirmed.
`Bacillus_cereus_group` = cereus-group member without the anthracis
chromosome — if paga/capb are present, an "Anthrax-like" note is emitted
(B. cereus biovar anthracis / cross-over strain, or λBa03-negative draft).
The **scope** is `complex` (chromosomally inseparable group).

### Bordetella pertussis

| Gene | Role | Biology |
|---|---|---|
| `cupin-bp` | Species anchor (required) | Cupin-domain chromosomal protein; 100% conserved across diverse strains. S. pseudotuberculosis and other near relatives share ≤71% identity. |
| `ptxs1` | Virulence reporting | Pertussis toxin S1 subunit (chromosomal, Bvg-regulated). Cross-reacts 96–99% with B. parapertussis/holmesii — reporting only, never gating. |
| `tcfa5`/`tcfa3` | Virulence reporting (split) | Tracheal colonization factor, 1944bp split in halves. IS481 (hundreds of copies) fragments assemblies; each half survives on a single contig. |

**Interpretation**: All B. pertussis (including Bvg⁻ phase variants and
detoxified vaccine derivatives) carry the chromosomal anchor. ptxs1/tcfa
report toxin locus presence. IS481-based assays were removed (shared with
B. holmesii).

### Yersinia pestis

| Gene | Role | Biology |
|---|---|---|
| `spi-yp` | Species anchor (required) | Chromosomal sugar phosphate isomerase. Y. pestis is a recently emerged clone of Y. pseudotuberculosis (>97% chromosomal identity) — only 4 chromosomal genes distinguish them. |
| `btub-yp` | Second chromosomal marker | Vitamin B12 receptor BtuB (literature precedent). |
| `caf1` | Virulence reporting (pMT1 plasmid) | F1 capsule antigen. 13% of real-world strains have lost pMT1. |
| `pla` | Virulence reporting (pPCP1 plasmid) | Plasminogen activator. The small 9.6kb pPCP1 is easily cured during lab passage. |

**Interpretation**: Species identity from the chromosome, toxigenicity from
plasmid status. A pPCP1-negative strain (common lab derivative) still
reports as Y. pestis.

### Yersinia enterocolitica

| Gene | Role | Biology |
|---|---|---|
| `glk-ye` | Species anchor (required) | Glycerate kinase, chromosomal. 100% conserved across a 50-strain real-world panel (biotypes 1A–4). Zero homology in Y. pseudotuberculosis/frederiksenii/intermedia/kristensenii/pestis. |
| `yddg-ye` | Second chromosomal marker | Aromatic amino acid transporter. Redundancy rider. |
| `ail` | Pathogenicity reporting | Attachment-invasion locus, chromosomal. Present only in pathogenic bioserotypes (1B, 2–5); absent from non-pathogenic biotype 1A. |
| `yada` | Virulence reporting (pYV plasmid) | Yersinia adhesin A. The 70kb pYV plasmid is easily cured. |

**Interpretation**: Non-pathogenic biotype 1A strains (common in food
surveys) still report as Y. enterocolitica with a note: "No pathogenicity
markers detected". Pathogenic bioserotypes additionally show ail/yada.

### Neisseria meningitidis

| Gene | Role | Biology |
|---|---|---|
| `sodc` | Species anchor (required) | Periplasmic Cu/Zn superoxide dismutase. The classic N. meningitidis vs N. gonorrhoeae discriminator (Corless 2001). Present in all meningococci including capsule-null carriage isolates. |
| `ctra` | Capsule reporting | Capsule transport protein. Non-groupable (NG) and capsule-null strains lack this; the species call remains N. meningitidis. |

### Neisseria gonorrhoeae

| Gene | Role | Biology |
|---|---|---|
| `pora` | Allele family 1 (porB1a) | Major porin, originally named "porA". Covers the porB1a/PIA allele family. |
| `porb1b` | Allele family 2 (porB1b) | The second antigenically divergent porin allele family (~78% inter-family identity). Either family firing is sufficient. N. meningitidis porB tops at 87.8% (below threshold). |

### Streptococcus pneumoniae

| Gene | Role | Biology |
|---|---|---|
| `lyta5`/`lyta3` | Species markers (split) | Autolysin (N-acetylmuramoyl-L-alanine amidase), the gold-standard pneumococcal marker. 957bp split in halves for fragmentation resistance (S. pneumoniae's high recombination rate breaks the gene in 6/100 draft assemblies). |

### Haemophilus influenzae

| Gene | Role | Biology |
|---|---|---|
| `hpd` | Species marker | Hemoglobin-binding protein D. Specific to H. influenzae; H. parainfluenzae/haemolyticus/pittmaniae lack full-length homologs. |

### Haemophilus parainfluenzae

| Gene | Role | Biology |
|---|---|---|
| `ribf-hpi` | Species anchor (required) | GTP cyclohydrolase II (riboflavin biosynthesis). 98.6% in confirmed strain; max 83.4% in near relatives. H. parainfluenzae is commonly misidentified as H. influenzae in clinical labs. |
| `sera-hpi` | Second marker | Phosphoglycerate dehydrogenase. |

### Mycoplasma pneumoniae

| Gene | Role | Biology |
|---|---|---|
| `tuf-mp` | Species anchor (required) | Elongation factor Tu, 1185bp. 100% conserved across 30 diverse strains; zero homology in M. genitalium/hominis. |
| `gyrb-mp` | Second chromosomal marker | DNA topoisomerase B. |
| `p1a`/`p1c` | Cytadherence reporting (2 of 4 quarters) | P1 adhesin (4884bp), split in quarters. The full gene requires 0.6% of the entire 810kb genome on one contig; quarters survive fragmentation. Virulence gene, not a species anchor. |

### Corynebacterium diphtheriae

| Gene | Role | Biology |
|---|---|---|
| `rela-cd` | Species anchor (required) | GTP pyrophosphokinase (stringent response). 98.3–100% across diverse strains; max 70.5% in C. sporogenes/tetani. |
| `mena-cd` | Second chromosomal marker | o-Succinylbenzoate synthase (vitamin K2 biosynthesis). |
| `tox` | Toxigenicity reporting | Diphtheria toxin, carried on β-corynephage. 37% of real-world strains are non-toxigenic (phage absent). Non-toxigenic C. diphtheriae causes invasive disease (reported since 1990s). |

### Clostridium botulinum

| Gene | Role | Biology |
|---|---|---|
| `bonta` | Neurotoxin type A | 3891bp from ATCC 3502 (type A1 reference). BoNT serotypes share only 30–40% inter-type identity — each type needs its own marker. |
| `bontb` | Neurotoxin type B | Covers type B strains (some subtypes). |
| `bonte` | Neurotoxin type E | Covers type E strains (foodborne, fish-associated). |
| `bontc` | Neurotoxin type C | From phage c-st (D90210). Covers animal-associated type C strains. |
| `bontd` | Neurotoxin type D | From X54254. Covers type D and C/D mosaic strains. |
| `bontf` | Neurotoxin type F | From F634 (subtype F5). Covers rare but clinically important type F. |

**Interpretation**: C. botulinum is a **toxotype**, not a genomic species —
four genomic groups (I–IV) unified only by BoNT production. Species
identity is toxin-based; the chromosomal background (C. botulinum vs
C. sporogenes, ANI 93.4%) requires L2 ANI methods (`-m marker+mash`).
Non-toxigenic Clostridium correctly reports Unknown.

### Vibrio cholerae

| Gene | Role | Biology |
|---|---|---|
| `ompw` | Species marker (identity 85%) | Outer membrane protein W. Non-O1/non-O139 environmental strains carry divergent alleles at 86–90% identity (the marker derives from an O1 strain). Nearest non-cholerae: V. mimicus at 80.8%. |

**Interpretation**: Threshold lowered to 85% to accommodate non-O1
environmental strain diversity. O1/O139 pandemic lineages hit at ~99–100%.

### Burkholderia pseudomallei

| Gene | Role | Biology |
|---|---|---|
| `glgb-bp` | Species anchor (required) | GlgB 1,4-alpha-glucan branching enzyme, 2217bp chromosomal. 99.6–100% identity and 100% coverage across 100 clinical melioidosis genomes; nearest cross *B. thailandensis* 73.3%, *B. cepacia* 82.5%. |
| `bima-bp` | Virulence reporting | BimA(Bp) actin-nucleating factor; mediates actin-based motility and cell-to-cell spread. Present in 95/100 panel genomes (some strains carry Bm-like or divergent alleles); never gates the call. |
| `sctn-bp` | Virulence reporting | T3SS ATPase SctN (TTS1 cluster). Present in 99/100 panel genomes at ≥95% (one TTS-cluster deletion); never gates the call. |

**Interpretation**: *B. pseudomallei* (melioidosis) and *B. mallei*
(glanders) are one genomic species (ANI ~99.5%) — *B. mallei* is a clone that
lost ~1.5Mb and became host-restricted. The marker layer cannot separate
them. The bundled `cpin` database adds a BimA(Bm) allele rule with
`priority_over` for glanders, thresholded at 98% because 3/100 clinical
pseudomallei genomes carry a Bm-like allele at 96.2–96.4%. GTDB rs226 labels
the whole clade `s__Burkholderia mallei`; NCBI and clinical microbiology use
*B. pseudomallei*, which this rule reports. A marker/ANI disagreement
(pseudomallei vs mallei) reflects the naming convention, not a
misidentification — confirm with L2 ANI and report both.

## Adding a New Species

```yaml
# Add to the rules list:
- species: Mycobacterium_tuberculosis
  genes:
  - mpb64
  - is6110
  - rpoB_mtbc
  min_hits: 2
  min_identity: 90
  exclude_genes: [rpob_maa]    # exclude M.africanum
  note: MTBC complex identification
```

Then:
1. Append corresponding gene sequences to `markers.fasta`
   (format: `>gside|db=<database>|gene=<gene>|acc=<accession>|func=<role>|note=<text>`;
   `func` is `primary`/`confirm`/`typing`/`virulence`/`control`)
2. Rebuild BLAST database:
   `makeblastdb -in markers.fasta -dbtype nucl -out markers_blastdb`
3. Run `gside scan test.fna -m marker` to verify
