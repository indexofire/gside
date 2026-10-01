# Marker Rules Format

`marker_rules.yaml` defines target gene combination identification rules —
adding new species requires only a rule entry, no code changes.

## File Location

```
data/reference/species/marker_rules.yaml
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
    min_coverage: <float>             # optional: minimum coverage % (default 50)
    exclude_genes:                    # optional: exclusion guard
      - <cross_reactive_gene>
    priority_over:                    # optional: priority declaration
      - <other_species>
    note: <str>                       # optional: documentation
```

## Field Reference

### species

Species name using registry key format (spaces → underscores).
Examples: `Campylobacter_jejuni`, `V_parahaemolyticus`.

### genes

Candidate marker gene names (lowercase). These genes must exist in
`markers_v2.fasta`.

### min_hits

Minimum number of genes that must hit. Setting 2 means at least 2 genes from
the `genes` list must reach thresholds for the species to be called.

Multi-gene combinations significantly reduce cross-reactions (e.g.,
C.jejuni's 2-of-3 design).

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

### priority_over

Declares this species takes priority over others (e.g., DEC over
Shigella_EIEC).

## Global Thresholds

The engine layer defines global gates (in `multigene_identifier.py`):

| Parameter | Value | Description |
|-----------|-------|-------------|
| `_MIN_IDENTITY` | 85.0 | Global minimum identity |
| `_MIN_COVERAGE` | 60.0 | Global minimum coverage (prevents short-HSP false fires) |
| `_HIGH_CONF` | 90.0 | High confidence threshold |

Rule-level `min_identity` can override the global value (but not below 85).

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
1. Append corresponding gene sequences to `markers_v2.fasta`
   (format: `>markers_v2~~~mpb64~~~ACCESSION description [species] role=primary`)
2. Rebuild BLAST database:
   `makeblastdb -in markers_v2.fasta -dbtype nucl -out markers_v2_blastdb`
3. Run `gside species test.fna --mode marker` to verify
