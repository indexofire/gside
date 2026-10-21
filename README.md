# gside — Genome Species IDentification Engine

[![CI](https://github.com/indexofire/gside/actions/workflows/test.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![Docs](https://github.com/indexofire/gside/actions/workflows/docs.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![PyPI](https://img.shields.io/pypi/v/gside.svg)](https://pypi.org/project/gside/)

A standalone CLI for multi-layer species identification from assembled contigs and metagenomic profiling from FASTQ reads.

```bash
pip install gside

# Single-genome species identification
gside species contigs.fna --mode all

# Metagenomic profiling
gside meta reads_R1.fq reads_R2.fq --mode all
```

`pip install` ships the Python package and the bundled marker database, but the
aligners (`blastn`, `skani`, `mash`, `sourmash`, `kraken2`, `bracken`, `sylph`)
come from conda/BioConda. Either use the pixi env (`pixi install`), or point
`GSIDE_PIXI_BIN` at a directory that contains them:

```bash
export GSIDE_PIXI_BIN=/path/to/conda/env/bin
gside db status          # shows which tiers are ready
```

## Species Identification (`gside species`)

| Mode | Engine | Method | Database |
|---|---|---|---|
| `marker` (default) | blastn | Target-gene combination rules (38 rules, 80 reference sequences) | Bundled (mini tier) |
| `panel` | skani | ANI against curated reference panel (291 genomes) | ~130 MB |
| `mash_refseq` | mash | RefSeq MinHash distance | ~179 MB |
| `sourmash` | sourmash | GTDB gather | ~3.9 GB |
| `gtdbtk` | GTDB-Tk | Phylogenetic placement | External |
| `all` | All above | Sequential execution + layered arbitration | As needed |

### Layered arbitration

| Layer | Methods | Rationale |
|---|---|---|
| L3 | gtdbtk | Phylogenetic placement — most authoritative |
| L2 | panel, mash_refseq, sourmash | Comparative genomics (ANI / sketch) |
| L1 | marker | Target-gene presence/absence |

Higher layers always outrank lower ones; within a layer, higher confidence wins.

### Quick start

```bash
pixi install && pip install -e .

# Identify species (marker mode works out of the box)
gside species contigs.fna

# Install reference databases
gside db setup --tier panel

# Multi-method identification
gside species contigs.fna --mode all

# Parallel batch processing with streaming output
gside species *.fna --mode all -j 8 -t 4 -o results.tsv
```

## Metagenomic Profiling (`gside meta`)

| Mode | Engine | Input | Best for |
|---|---|---|---|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | Fast screening, contig support, abundance |
| `sylph` | sylph | FASTQ only | ANI-anchored precision, species-level calls |
| `all` | Both | FASTQ (both) / FASTA (kraken2) | Merged profile |

```bash
gside meta reads_R1.fq reads_R2.fq --mode all
gside meta contigs.fna --mode kraken2
```

## Agent-Friendly Design

gside is built for AI agents as much as humans:

```bash
# JSON on stdout (always exit 0; errors in methods.<mode>.error)
gside species genome.fna --mode all

# Structured missing-database errors
# → methods.panel.suggested_command = "gside db setup --tier panel"
# → methods.panel.missing_tier = "panel"
# → methods.panel.estimated_size = "~130MB sketch / ~1GB genomes"

# Auto-install missing databases
gside species genome.fna --mode panel --auto-setup

# Machine-readable database status
gside db status -f json
gside db list -f json
```

### MCP Server

gside includes an MCP (Model Context Protocol) stdio server:

```bash
gside mcp    # or: gside-mcp
```

Exposes `species`, `db_status`, `db_list`, `db_setup` as MCP tools for Claude Desktop or any MCP client.

```json
{
  "mcpServers": {
    "gside": { "command": "gside-mcp" }
  }
}
```

## CLI Reference

```
gside species <contigs.fna> [options]
  -m, --mode         marker|panel|mash_refseq|sourmash|gtdbtk|all
  -D, --db-dir       Database root
  -f, --format       json|tsv|md
  -o, --output       Write to file (tsv/md stream per input)
  -j, --jobs         Parallel input files
  -t, --threads      Threads per aligner
  -a, --auto-setup   Install missing databases automatically
  -n, --no-header    Suppress TSV header row

gside meta <reads/files> [options]
  -m, --mode         kraken2|sylph|all
  -D, --db-dir       Database root
  -f, --format       json|tsv
  -o, --output       Write to file
  -t, --threads      Threads per tool

gside db <subcommand>
  status             Show database readiness (--format table|tsv|csv|json)
  setup              Install databases (--tier mini|panel|mash|sourmash|all)
  list               List tiers with readiness (--format table|tsv|csv|json)

gside validate <contigs.fna>
  --mode simple|standard

gside mcp            MCP stdio server
gside --version      Version info
gside -V -J          Version as JSON
```

## Databases

| Tier | Size | Engine | Status |
|---|---|---|---|
| mini | ~1 MB | blastn (marker) | ✅ Bundled |
| panel | ~130 MB | skani | Downloadable |
| mash | ~179 MB | mash | Downloadable |
| sourmash | ~3.9 GB | sourmash | Downloadable |
| meta | ~5.5 GB | kraken2 | Manual setup |
| meta_sylph | ~14 GB | sylph | Manual setup |

GTDB-Tk manages its own external database (`GTDBTK_DATA_PATH`).

## Installation

```bash
# With pixi (recommended — includes blastn, skani, mash, sourmash)
curl -fsSL https://pixi.sh/install.sh | bash
pixi install
pip install -e .

# From PyPI (Python 3.12, binaries not included)
pip install gside
```

## Documentation

- [Quick Start](https://indexofire.github.io/gside/)
- [CLI Reference](https://indexofire.github.io/gside/usage/cli/)
- [Identification Modes](https://indexofire.github.io/gside/usage/modes/)
- [Metagenomic Profiling](https://indexofire.github.io/gside/usage/meta/)
- [JSON Output Contract](https://indexofire.github.io/gside/reference/json-contract/)
- [中文文档](README.zh.md)

## License

MIT

## Ecosystem

| Tool | Purpose |
|---|---|
| [gside](https://github.com/indexofire/gside) | Species identification + metagenomic profiling |
| [gapit](https://github.com/indexofire/gapit) | Gene screening + phenotype designation (serogroup, toxigenicity) |
| [gmlst](https://github.com/indexofire/gmlst) | MLST / cgMLST typing |
