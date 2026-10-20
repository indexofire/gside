# Database Management

gside reference databases are managed via the `gside db` subcommand with tiered installation.

## Database Tiers

| Tier | Contents | Size | Purpose |
|------|----------|------|---------|
| **mini** | marker_rules.yaml + markers.fasta + BLAST db | ~1MB | Target gene identification (L1) |
| **panel** | D2_ani (291 genomes, skani sketch) | ~130MB | ANI identification (L2) |
| **mash** | D3_mash (RefSeq MinHash sketch) | ~179MB | Distance identification (L2) |
| **sourmash** | D4_sourmash (GTDB reps k=31 + lineages) | ~3.9GB | Gather classification (sourmash) |
| **meta** | D5_kraken2 (Kraken2 Standard-8 index) | ~5.5GB | Metagenomic profiling (`gside meta`) |
| **meta_sylph** | D6_sylph (sylph GTDB database) | ~14GB | ANI-based metagenomic profiling (`gside meta`) |
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
gside db setup --tier meta           # Kraken2 index (manual, ~5.5GB)
gside db setup --tier meta_sylph     # sylph GTDB db (manual, ~14GB)
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

The metagenomic tiers are **manual setup** — there is no auto-download.
`gside db setup --tier meta` / `--tier meta_sylph` print an actionable
pointer and exit non-zero until the database is placed by hand:

```bash
# D5_kraken2: download a prebuilt index (e.g. Standard-8) from
# https://benlangmead.github.io/aws-indexes/k2 and extract so that
# hash.k2d / opts.k2d / taxo.k2d sit in data/db/D5_kraken2/
# (Bracken additionally needs database100mer.kmer_distrib there)

# D6_sylph: download a prebuilt *.syl2db (e.g. GTDB) from
# https://sylph-docs.github.io/ and place it under data/db/D6_sylph/
```

`--source` also works for these tiers, copying an existing local database
directory into place.

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
