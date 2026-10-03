# Database Management

gside reference databases are managed via the `gside db` subcommand with tiered installation.

## Database Tiers

| Tier | Contents | Size | Purpose |
|------|----------|------|---------|
| **mini** | marker_rules.yaml + markers.fasta + BLAST db | ~1MB | Target gene identification (L1) |
| **panel** | D2_ani (291 genomes, skani sketch) | ~130MB | ANI identification (L2) |
| **mash** | D3_mash (RefSeq MinHash sketch) | ~179MB | Distance identification (L2) |
| **sourmash** | D4_sourmash (GTDB reps k=31 + lineages) | ~3.9GB | Gather classification (sourmash) |
| **all** | panel + mash | ~310MB | Full capability (excludes sourmash) |

### Directory naming

`data/db/` directories are named per identification layer:

| Old name | New name | Method |
|----------|----------|--------|
| `refseq_panel` | `D2_ani` | panel (skani ANI) |
| `mash_refseq` | `D3_mash` | mash_refseq (MinHash) |
| `sourmash_gtdb` | `D4_sourmash` | sourmash (GTDB gather) |

Migrate an existing checkout with:

```bash
mv data/db/refseq_panel data/db/D2_ani
mv data/db/mash_refseq data/db/D3_mash
mkdir -p data/db/D4_sourmash
```

Note: the prebuilt GitHub Release archive still contains the inner
`panel.sketch/` directory unchanged, so download-and-extract works as before.

## Commands

### Check status

```bash
gside db status
```

```
gside database status
───────────────────────────────────────────────
  ✅ markers          tier=mini   path=data/db/D1_marker
  ✅ D2_ani        tier=panel  path=data/db/D2_ani
  ✅ D3_mash         tier=mash   path=data/db/D3_mash
  ✅ D4_sourmash     tier=sourmash path=data/db/D4_sourmash

  Run 'gside db setup --tier <tier>' to provision
```

### Install / update

```bash
gside db setup --tier panel          # Panel only (ANI), the default tier
gside db setup --tier mash           # MinHash only
gside db setup --tier sourmash       # GTDB gather db (~3.9GB)
gside db setup --tier all            # panel + mash (excludes sourmash)
```

With no `--tier` given, setup defaults to `panel`. Panel tries a pre-built
sketch from the GitHub release first and falls back to building from the
bundled manifest (downloading genomes from NCBI). The mash sketch is
downloaded from Zenodo and MD5-verified. Setup returns a non-zero exit code
if any step reports an error.

`--tier all` installs panel + mash only. The sourmash GTDB database
(~3.9GB) is never included in `all`; run `gside db setup --tier sourmash`
separately to provision it.

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

Edit `data/db/D1_marker/marker_rules.yaml`:

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
cd data/db/D1_marker
makeblastdb -in markers.fasta -dbtype nucl -out markers_blastdb
```

No code changes required — rules are data.
