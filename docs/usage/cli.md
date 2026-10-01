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
| `--db-dir` | option | auto | Database root directory |

### Identification Modes

| Mode | Description | Database Required |
|------|-------------|-------------------|
| `marker` | Target gene combination rules (blastn) | mini (bundled) |
| `panel` | ANI against curated panel (skani) | panel tier |
| `mash_refseq` | MinHash distance (mash) | mash tier |
| `sourmash` | sourmash GTDB gather | sourmash GTDB db |
| `all` | Run all available methods + arbitrate | as needed |

### Examples

```bash
# Target gene identification (fastest, no extra database needed)
gside species contigs.fna --mode marker

# ANI identification (precise)
gside species contigs.fna --mode panel

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
| `--tier` | `panel` | mini / panel / mash / all |
| `--source` | — | Copy from existing bacmap databases |

### Examples

```bash
gside db status
gside db setup --tier all
gside db setup --tier panel --source ~/bacmap/data/db
```

## gside --version

```bash
gside --version
# gside 0.1.0
```
