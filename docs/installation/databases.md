# Database Management

gside reference databases are managed via the `gside db` subcommand with tiered installation.

## Database Tiers

| Tier | Contents | Size | Purpose |
|------|----------|------|---------|
| **mini** | marker_rules.yaml + markers.fasta + BLAST db | ~1MB | Target gene identification (L1) |
| **panel** | D2_ani (291 genomes, skani sketch) | ~130MB | ANI identification (L2) |
| **mash** | D3_mash (RefSeq MinHash sketch) | ~179MB | Distance identification (L2) |
| **sourmash** | D4_sourmash (GTDB reps k=31 + lineages) | ~3.9GB | Gather classification (sourmash) |
| **meta** | D5_kraken2 (Kraken2 Standard-8 index) | ~5.5GB | Metagenomic profiling (`gside scan --meta`) |
| **meta_sylph** | D6_sylph (sylph GTDB database) | ~14GB | ANI-based metagenomic profiling (`gside scan --meta`) |
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
  ✅ D1_marker          tier=mini   path=data/db/D1_marker
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
gside db setup --tier meta           # Kraken2 + Bracken (~5.5GB)
gside db setup --tier meta_sylph     # sylph GTDB db (~14GB)
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

The metagenomic tiers download automatically from upstream:

```bash
gside db setup --tier meta         # Kraken2 Standard-8 + Bracken (~5.5GB)
gside db setup --tier meta_sylph   # sylph GTDB r226, .syl2db (~14GB)
```

- **meta** pulls Ben Langmead's Standard-8 archive, which bundles the Kraken2
  index *and* the Bracken 50-300-mer distributions. Extracted files are
  verified against the published `.md5` manifest; readiness requires
  `database150mers.kmer_distrib` (Bracken's read length in gside).
- **meta_sylph** pulls the two-stage `gtdb-r226-c200-dbv2.syl2db` from the
  HTTPS mirror and verifies the GCS-published MD5.

Both are large, one-time downloads. `gside scan --meta --auto-setup` installs the
missing tier on first use (progress on stderr) and retries the method once.
`--source` still copies an existing local database instead of downloading.

### Integrity checks

`gside db setup` verifies digests for everything it can pin:

| Tier | Digest | Notes |
|---|---|---|
| panel | SHA256 (`PANEL_SHA256`) | GitHub release asset |
| mash | MD5 (`MASH_MD5`) | Zenodo publishes MD5 only |
| sourmash | SHA256 for lineages + signature (`SOURMASH_SIG_SHA256`) | signature pinned from a verified download |

The sourmash signature zip (~3.9GB) has no upstream digest. gside pins one
recorded after a locally verified download (complete size, zip integrity,
143,385 k=31 DNA signatures, CC0):

```
SOURMASH_SIG_SHA256 = 6f3cc74f0bdcf6f84466ba91f87d6f2bf30081d5814bdf006cafd1e9cbbeada9
```

Override it at setup time if you mirror a different build:

```bash
export GSIDE_SOURMASH_SIG_SHA256=<sha256-of-your-zip>
gside db setup --tier sourmash
```

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

## Genome quality databases

`gside validate` needs two external databases that are **not** `db setup`
tiers: CheckM2 (completeness/contamination) and GTDB-Tk (classification).
They are large, version-coupled and GTDB-Tk is memory-hungry, so gside
deliberately leaves them to the official installers.

| Database | Download | On disk | Runtime | Env var |
|---|---|---|---|---|
| CheckM2 (v1.1.0) | ~1.7GB | ~3GB | modest (DIAMOND) | `CHECKM2DB` → `.../uniref100.KO.1.dmnd` |
| GTDB-Tk R232 | ~57GB | ~98GB | ~140GB RAM (bacteria) | `GTDBTK_DATA_PATH` → DB directory |

### CheckM2

```bash
pixi add checkm2          # or: conda install -c bioconda checkm2
checkm2 database --download --path /path/to/D7_checkm2
export CHECKM2DB=$(find /path/to/D7_checkm2 -name uniref100.KO.1.dmnd | head -1)
```

Manual fallback (no CheckM2 binary needed):

```bash
mkdir -p /path/to/D7_checkm2
wget -c https://zenodo.org/records/14897628/files/checkm2_database.tar.gz
md5sum checkm2_database.tar.gz    # 07c10655620843b517d0df0c160d911f
tar xzf checkm2_database.tar.gz -C /path/to/D7_checkm2
export CHECKM2DB=$(find /path/to/D7_checkm2 -name uniref100.KO.1.dmnd | head -1)
```

The database must be compatible with the installed CheckM2 version; prefer
`checkm2 database --download`, which selects the matching release.

### GTDB-Tk

```bash
pixi add gtdbtk
download-db.sh /path/to/D8_gtdbtk    # official script; target dir must be empty
export GTDBTK_DATA_PATH=/path/to/D8_gtdbtk
```

Manual fallback:

```bash
mkdir -p /path/to/D8_gtdbtk
wget -c https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
md5sum gtdbtk_data.tar.gz            # R232: 25a59e0352b1fd150c589f56559767d4
tar xvzf gtdbtk_data.tar.gz -C /path/to/D8_gtdbtk --strip 1
export GTDBTK_DATA_PATH=/path/to/D8_gtdbtk
```

| Release | GTDB-Tk version | MD5 |
|---|---|---|
| R232 | ≥ 2.7.0 | `25a59e0352b1fd150c589f56559767d4` |
| R226 | 2.4.1 – 2.6.1 | `24b476ea5a4ef30519d461e56cc4a27f` |

Hardware: ~140GB RAM for bacteria (~100GB archaea) and ~98GB disk for R232.
On smaller machines `validate` records GTDB-Tk as `missing`, writes the
CheckM2 results anyway and reports `status: partial` (exit 1).

### What gside reports

`gside db status` covers the `db setup` tiers only. `gside validate` reads
the two environment variables above and records per-tool status in
`validation.json` (`tools.checkm2`, `tools.gtdbtk`, `errors`). Exit code is
0 only when both tools ran.

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
