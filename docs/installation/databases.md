# Database Management

gside reference databases are managed via the `gside db` subcommand with tiered installation.

## Database Tiers

| Tier | Contents | Size | Purpose |
|------|----------|------|---------|
| **mini** | marker_rules.yaml + markers_v2.fasta + BLAST db | ~1MB | Target gene identification (L1) |
| **panel** | refseq_panel (291 genomes, skani sketch) | ~88MB | ANI identification (L2) |
| **mash** | mash_refseq (RefSeq MinHash sketch) | ~179MB | Distance identification (L2) |
| **all** | panel + mash | ~267MB | Full capability |

## Commands

### Check status

```bash
gside db status
```

```
gside database status
──────────────────────────────────────────
  ✅ markers_v2       tier=mini   path=data/reference/species
  ✅ refseq_panel     tier=panel  path=data/db/refseq_panel
  ✅ mash_refseq      tier=mash   path=data/db/mash_refseq
```

### Install / update

```bash
gside db setup --tier panel          # Panel only (ANI)
gside db setup --tier mash           # MinHash only
gside db setup --tier all            # Everything
```

### Copy from bacmap

Users with an existing hermes-bacmap installation can copy databases directly:

```bash
gside db setup --tier all --source ~/repos/github/hermes-bacmap/data/db
```

### List available tiers

```bash
gside db list
```

## Database Paths

Default paths can be overridden via environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `GSIDE_DATA_DIR` | `<repo>/data` | Small reference data (marker rules) |
| `GSIDE_DB_DIR` | `<data>/db` | Large databases (panel/mash) |
| `GSIDE_PIXI_BIN` | Auto-detected | Bioinformatics binary directory |

## Custom Marker Rules

Edit `data/reference/species/marker_rules.yaml`:

```yaml
- species: Mycobacterium_tuberculosis
  genes:
  - mpb64
  - is6110
  min_hits: 2
  min_identity: 90
  exclude_genes: [some_cross_reactive_gene]
```

Then rebuild the BLAST database:

```bash
cd data/reference/species
makeblastdb -in markers_v2.fasta -dbtype nucl -out markers_v2_blastdb
```

No code changes required — rules are data.
