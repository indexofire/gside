# Marker DB control genes (reads-mode relative-depth gate)

The reads path uses a gspec-style abundance gate:

```
mean_depth = sum(aligned target interval lengths) / target length
relative_depth = mean_depth / control_depth
keep if mean_depth >= min_depth and relative_depth >= min_relative_depth
```

Defaults (also in `metadata.yaml`):

```
min_depth: 5.0
min_relative_depth: 0.10
min_mapq: 20
```

`mean_depth` is not a read count, so genes of different lengths are directly
comparable.  MAPQ filtering happens before aggregation.

## Dedicated controls (28 species)

| species | control gene | source |
|---|---|---|
| **Campylobacter** | | |
| Campylobacter_jejuni | lpxa | C. jejuni allele + C. coli lpxA_2 allele |
| Campylobacter_coli | lpxa | shared with C. jejuni |
| **Confirm markers (pre-existing)** | | |
| Salmonella | rpoD | confirm marker |
| Listeria_monocytogenes | prs | confirm marker |
| Staphylococcus_aureus | femA | confirm marker |
| Klebsiella_pneumoniae | gyrb-kpn | confirm marker |
| Acinetobacter_baumannii | rpob-aba | confirm marker |
| Pseudomonas_aeruginosa | gyrb-pae | confirm marker |
| Streptococcus_pneumoniae | psaa | confirm marker |
| Neisseria_meningitidis | sodc | confirm marker |
| **New rpoB controls (RefSeq reference genomes)** | | |
| Escherichia_coli | rpob-escherichia_coli | E. coli O157:H7 Sakai |
| Shigella (all species) | rpob-escherichia_coli | shared with E. coli |
| Vibrio_parahaemolyticus | rpob-vibrio_parahaemolyticus | RIMD 2210633 |
| Bacillus_anthracis | rpob-bacillus_anthracis | B. cereus ATCC 14579 (close relative) |
| Bordetella_pertussis | rpob-bordetella_pertussis | Tohama I |
| Clostridioides_difficile | rpob-clostridioides_difficile | strain 630 |
| Clostridium_perfringens | rpob-clostridium_perfringens | FDAARGOS_903 |
| Corynebacterium_diphtheriae | rpob-corynebacterium_diphtheriae | NCTC 13129 |
| Enterococcus_faecalis | rpob-enterococcus_faecalis | EnGen0336 |
| Haemophilus_influenzae | rpob-haemophilus_influenzae | Rd KW20 |
| Legionella_pneumophila | rpob-legionella_pneumophila | Philadelphia 1 |
| Neisseria_gonorrhoeae | rpob-neisseria_gonorrhoeae | FA 1090 |
| Streptococcus_agalactiae | rpob-streptococcus_agalactiae | 2603 V/R |
| Streptococcus_pyogenes | rpob-streptococcus_pyogenes | SF370 |

## Species on sample_median fallback (14)

These species use `control_fallback: sample_median` (median mean depth of
detected markers).  The absolute `min_depth` + `min_mapq` still apply, but
there is no independent single-copy control.  Add a dedicated rpoB/gyrB
sequence when a verified reference genome is available.

| species | recommended gene |
|---|---|
| Vibrio_cholerae | rpoB / gyrB |
| Vibrio_vulnificus | rpoB / gyrB |
| Enterococcus_faecium | rpoB / gyrB |
| Cronobacter_sakazakii | gyrB (rpoB is the rule gene) |
| Mycoplasma_pneumoniae | rpoB / gyrB / tufA |
| Helicobacter_pylori | rpoB / gyrB |
| Streptococcus_suis | rpoB / gyrB (gdh is the rule gene) |
| Burkholderia_pseudomallei | rpoB / gyrB |
| Leptospira_interrogans | rpoB / gyrB |
| Treponema_pallidum | rpoB / gyrB |
| Yersinia_pestis | rpoB / gyrB |
| Yersinia_enterocolitica | rpoB / gyrB |
| Aeromonas_hydrophila | rpoB / gyrB |
| Clostridium_botulinum | rpoB / gyrB |

## Implementation notes

- Control genes are added to `markers.fasta` with `role=control` and are not
  referenced by any species rule in `marker_rules.yaml`.
- `metadata.yaml` maps species → control gene names; the code resolves these
  at runtime.
- E. coli rpoB serves all Shigella species rules (E. coli and Shigella are
  genomically the same species).
- B. cereus rpoB serves B. anthracis (very close relatives, >95% genome identity).
- For Campylobacter, two `lpxa` alleles are present.  The depth control uses
  the **deepest** allele per gene, not the highest-identity allele, to avoid
  a few contaminant reads overriding the sample's real control depth.
- The control identity floor is 80% (vs 90% for species markers) so divergent
  control alleles still contribute depth.

## How to add a new control gene

1. Obtain a verified reference genome for the species (verify the contig
   descriptions match the expected organism).
2. BLAST an existing rpoB query against the genome to locate the gene.
3. Extract ~1 kb around the best hit.
4. Add to `markers.fasta`:
   ```
   >markers~~~rpob-<species>~~~<accession> rpoB [<species name>] role=control
   ```
5. Add to `metadata.yaml`:
   ```yaml
   control_genes:
     <Species_name>: [rpob-<species>]
   ```
6. Rebuild the BLAST database:
   ```bash
   makeblastdb -in markers.fasta -dbtype nucl -out markers_blastdb -parse_seqids
   ```
7. Validate: `gside db validate <D1_marker path>`
