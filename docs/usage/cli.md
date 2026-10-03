# CLI Reference

## gside species

Main species identification command.

### Usage

```bash
gside species <contigs.fasta> [options]
```

### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `contigs` | positional | — | Assembled contigs FASTA file |
| `--mode` | option | `marker` | Identification mode (see below) |
| `--db-dir` | option | auto | Database root (`$GSIDE_DB_DIR` or `data/db`) |

### Identification Modes

Five modes are available. `marker` is the default.

| Mode | Description | Database Required |
|------|-------------|-------------------|
| `marker` | Target gene combination rules (blastn) | mini (bundled) |
| `panel` | ANI against curated panel (skani) | panel tier |
| `mash_refseq` | MinHash distance (mash) | mash tier |
| `sourmash` | sourmash GTDB gather | sourmash GTDB db |
| `all` | Run all available methods + arbitrate | as needed |

In `all` mode a failure in one method does not abort the run: the error is
recorded in that method's JSON entry and the remaining methods continue.

### Examples

```bash
# Target gene identification (fastest, no extra database needed)
gside species contigs.fna --mode marker

# ANI identification (precise)
gside species contigs.fna --mode panel

# MinHash distance against RefSeq
gside species contigs.fna --mode mash_refseq

# GTDB gather classification
gside species contigs.fna --mode sourmash

# Multi-method arbitration
gside species contigs.fna --mode all

# Specify database directory
gside species contigs.fna --mode panel --db-dir /path/to/dbs
```

## gside db

Database management subcommand.

### Subcommands

| Subcommand | Description |
|------------|-------------|
| `status` | Show database readiness status |
| `list` | List available tiers |
| `setup` | Install/update databases |

### setup Options

| Option | Default | Description |
|--------|---------|-------------|
| `--tier` | `panel` | mini / panel / mash / sourmash / all |
| `--source` | — | Copy from an existing local database directory instead of downloading |

### Examples

```bash
gside db status
gside db setup --tier all
gside db setup --tier panel --source /path/to/existing/data/db
```

## gside validate

Assembly validation with marker genes, optionally CheckM2 + GTDB-Tk.

### Usage

```bash
gside validate <contigs.fasta> [options]
```

### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `contigs` | positional | — | Assembled contigs FASTA file |
| `--mode` | option | `simple` | Validation depth: `simple` (marker genes) or `standard` (adds CheckM2 + GTDB-Tk) |
| `--output-dir` | option | auto | Directory for `validation.json` |

`standard` mode degrades gracefully: missing CheckM2/GTDB-Tk databases or
binaries produce warnings and partial results, never a crash. Configure via
`CHECKM2DB` / `GTDBTK_DATA_PATH` (or `GTDBDB`) environment variables.

### Examples

```bash
gside validate contigs.fna --mode simple
gside validate contigs.fna --mode standard --output-dir ./taxonomy
```

## gside --version

```bash
gside --version
# gside 0.1.0
```
