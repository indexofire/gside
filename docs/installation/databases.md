# Database Management

gside reference databases are managed via the `gside db` subcommand with tiered installation.

## Database Tiers

| Tier | Contents | Size | Purpose |
|------|----------|------|---------|
| **mini** | marker_rules.yaml + markers_v2.fasta + BLAST db | ~1MB | Target gene identification (L1) |
| **panel** | refseq_panel (291 genomes, skani sketch) | ~130MB | ANI identification (L2) |
| **mash** | mash_refseq (RefSeq MinHash sketch) | ~179MB | Distance identification (L2) |
| **all** | panel + mash | ~310MB | Full capability |

The `sourmash` mode uses a separate GTDB database (`sourmash_gtdb`) that is
not covered by the tiers above.

## Commands

### Check status

```bash
gside db status
```

```
gside database status
───────────────────────────────────────────────
  ✅ markers_v2       tier=mini   path=data/reference/species
  ✅ refseq_panel     tier=panel  path=data/db/refseq_panel
  ✅ mash_refseq      tier=mash   path=data/db/mash_refseq

  Run 'gside db setup --tier <tier>' to provision
```

### Install / update

```bash
gside db setup --tier panel          # Panel only (ANI), the default tier
gside db setup --tier mash           # MinHash only
gside db setup --tier all            # Everything
```

With no `--tier` given, setup defaults to `panel`. Panel tries a pre-built
sketch from the GitHub release first and falls back to building from the
bundled manifest (downloading genomes from NCBI). The mash sketch is
downloaded from Zenodo and MD5-verified. Setup returns a non-zero exit code
if any step reports an error.

### Copy from an existing local database directory

If you already have the databases on disk (for example from another machine
or a previous checkout), skip the download and copy them directly:

```bash
gside db setup --tier all --source /path/to/existing/data/db
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
