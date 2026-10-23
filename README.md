# gside — Genome Species IDentification Engine

[![CI](https://github.com/indexofire/gside/actions/workflows/test.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![Docs](https://github.com/indexofire/gside/actions/workflows/docs.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![PyPI](https://img.shields.io/pypi/v/gside.svg)](https://pypi.org/project/gside/)

A standalone CLI for multi-layer species identification from assembled contigs and metagenomic profiling from FASTQ reads.

```bash
pip install gside

# Single-genome species identification
gside scan contigs.fna -m all

# Metagenomic profiling
gside scan --meta reads_R1.fq reads_R2.fq -m all
```

`scan` infers FASTA vs FASTQ from the extension and rejects a mode that cannot
consume that input (e.g. `-m panel` on reads) before any tool runs. `--meta`
selects the community workflow and restricts modes to `kraken2`/`sylph`/`all`.
The old `gside species` and `gside meta` commands remain as hidden aliases.

`pip install` ships the Python package and the bundled marker database, but the
aligners (`blastn`, `skani`, `mash`, `sourmash`, `kraken2`, `bracken`, `sylph`)
come from conda/BioConda. Either use the pixi env (`pixi install`), or point
`GSIDE_PIXI_BIN` at a directory that contains them:

```bash
export GSIDE_PIXI_BIN=/path/to/conda/env/bin
gside db status          # shows which tiers are ready
```

## Species Identification (`gside scan`)

| Mode | Engine | Method | Database |
|---|---|---|---|
| `marker` (default) | blastn | Target-gene combination rules (42 rules, 139 reference sequences) | Bundled (mini tier) |
| `panel` | skani | ANI against curated reference panel (291 genomes) | ~130 MB |
| `mash_refseq` | mash | RefSeq MinHash distance | ~179 MB |
| `sourmash` | sourmash | GTDB gather | ~3.9 GB |
| `gtdbtk` | GTDB-Tk | Phylogenetic placement | External |
| `all` | All above | Sequential execution + layered arbitration | As needed |

`-m` accepts `+` combinations (`marker+sourmash`); `all` runs every genome
tier. Inputs: `marker`, `mash_refseq` and `sourmash` accept FASTA or FASTQ;
`panel`/`gtdbtk` require assemblies.

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

# One genome (marker works out of the box)
gside scan contigs.fna
gside scan contigs.fna -m marker+sourmash
gside scan contigs.fna -m all --auto-setup

# Isolate sequencing reads (R1/R2 are one sample)
gside scan reads.fq -m marker
gside scan sample_R1.fq sample_R2.fq -m marker

# Metagenomic community
gside scan reads_R1.fq reads_R2.fq --meta

# Genome quality (CheckM2 + GTDB-Tk)
gside validate contigs.fna -o taxonomy

# Batch + streaming output
gside scan *.fna -m all -j 8 -t 4 -o results.tsv
```

## Metagenomic Profiling (`gside scan --meta`)

| Mode | Engine | Input | Best for |
|---|---|---|---|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | Fast screening, contig support, abundance |
| `sylph` | sylph | FASTQ only | ANI-anchored precision, species-level calls |
| `all` | Both | FASTQ only (use `-m kraken2` for FASTA) | Merged profile |

```bash
gside scan --meta reads_R1.fq reads_R2.fq            # kraken2 + sylph
gside scan --meta reads_R1.fq reads_R2.fq -m kraken2
gside scan --meta contigs.fna -m kraken2
```

## Genome Quality (`gside validate`)

`validate` answers a different question from `scan`: **is the assembly good?**
It runs CheckM2 (completeness/contamination) + GTDB-Tk (classification) and
writes `validation.json` with tool versions, parameters, per-tool status and a
quality verdict. Output directory is required:

```bash
gside validate contigs.fna -o taxonomy
```

CheckM2 (~1.7GB) and GTDB-Tk (~57GB download, ~140GB RAM for bacteria) are
external databases — see [Databases](docs/installation/databases.md).

## Agent-Friendly Design

gside is built for AI agents as much as humans:

```bash
# JSON on stdout (always exit 0; errors in methods.<mode>.error)
gside scan genome.fna -m all

# Structured missing-database errors
# → methods.panel.suggested_command = "gside db setup --tier panel"
# → methods.panel.missing_tier = "panel"
# → methods.panel.estimated_size = "~130MB sketch / ~1GB genomes"

# Auto-install missing databases
gside scan genome.fna -m panel --auto-setup

# Machine-readable database status
gside db status -f json
gside db list -f json
```

### MCP Server

gside includes an MCP (Model Context Protocol) stdio server:

```bash
gside mcp    # or: gside-mcp
```

Exposes `species`, `meta`, `db_status`, `db_list`, `db_setup` as MCP tools for Claude Desktop or any MCP client.

```json
{
  "mcpServers": {
    "gside": { "command": "gside-mcp" }
  }
}
```

## CLI Reference

```bash
gside scan <files...> [options]
  -m, --mode         marker|panel|mash_refseq|sourmash|gtdbtk|kraken2|sylph|all
                     (+ combinations: marker+sourmash)
  --meta             Community workflow (kraken2|sylph|all)
  -D, --db-dir       Database root
  -f, --format       json|tsv|md|csv
  -o, --output       Write to file (row formats stream)
  -j, --jobs         Parallel samples (paired mates count as one)
  -t, --threads      Threads per tool
  -n, --no-header    Suppress the header row
  -a, --auto-setup   Install missing databases, retry once
  -s, --samples      --meta batch manifest

gside validate <contigs.fna> -o <dir> [options]
  -o, --output-dir   Required; writes validation.json (directory)
  -t, --threads      CheckM2 + GTDB-Tk threads

gside db <subcommand>
  status             Readiness (--format table|tsv|csv|json)
  setup              Install a tier (--tier mini|panel|mash|sourmash|meta|meta_sylph|all)
  list               Tiers + readiness (--format table|tsv|csv|json)

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
| meta | ~5.5 GB | kraken2 | Auto-download |
| meta_sylph | ~14 GB | sylph | Auto-download |

`gside validate` additionally needs CheckM2 (~1.7GB) and GTDB-Tk
(`GTDBTK_DATA_PATH`, ~57GB download) — see the docs.

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
