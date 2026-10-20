# gside — Genome Species IDentification Engine

[![CI](https://github.com/indexofire/gside/actions/workflows/test.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![PyPI](https://img.shields.io/pypi/v/gside.svg)](https://pypi.org/project/gside/)

A standalone CLI for multi-layer species identification from assembled contigs and metagenomic profiling from FASTQ reads.

## Overview

gside wraps standard bioinformatics tools (blastn, skani, mash, sourmash, GTDB-Tk, Kraken2, sylph) behind data-driven rules, arbitrates across methods, and prints a single JSON verdict per run. It is **agent-friendly by design**: structured error fields, auto-setup, MCP server, machine-readable database status.

```
Bioinformatics tools   blastn / skani / mash / sourmash / gtdbtk / kraken2 / sylph
                         ↑
gside engine            rules + backends + arbitration + metagenomic profiling
```

## Species Identification Layers

| Layer | CLI mode | Engine | Method | Notes |
|-------|----------|--------|--------|-------|
| **L3** | `gtdbtk` | GTDB-Tk | Phylogenetic placement | Most authoritative; external DB |
| **L2** | `panel` | skani | ANI against curated reference panel (291 genomes) | Requires panel tier |
| **L2** | `mash_refseq` | mash | MinHash distance (RefSeq sketch) | Requires mash tier |
| **L2** | `sourmash` | sourmash | GTDB gather (LCA classification) | Requires GTDB database |
| **L1** | `marker` | blastn | Target gene combination rules (38 rules, 80 sequences) | Mini database bundled |
| **Arbitration** | `all` | — | Run all methods, layer consensus (L3 > L2 > L1) | Best of all layers |

## Metagenomic Profiling

| Mode | Engine | Input | Best for |
|------|--------|-------|----------|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | Fast screening, contig support, abundance |
| `sylph` | sylph | FASTQ only | ANI-anchored precision, species-level calls |
| `all` | Both | FASTQ (both) / FASTA (kraken2) | Merged profile |

## Quick Start

```bash
# Install
pixi install
pip install -e .

# Species identification (marker mode works out of the box)
gside species contigs.fna

# Multi-method identification with auto-setup
gside species contigs.fna --mode all --auto-setup

# Parallel batch processing
gside species *.fna --mode all -j 8 -t 4 -o results.tsv

# Metagenomic profiling
gside meta reads_R1.fq reads_R2.fq --mode all

# Database management
gside db status
gside db setup --tier all
```

## Agent-Friendly Design

- **JSON on stdout** — species always exits 0; per-method errors in `methods.<mode>.error`
- **Structured missing-database errors** — `suggested_command`, `missing_tier`, `estimated_size`
- **Auto-setup** — `--auto-setup` installs missing tiers, progress on stderr
- **Machine-readable status** — `db status --format json`, `db list --format json`
- **MCP server** — `gside mcp` exposes tools to Claude Desktop / MCP clients
- **Streaming output** — tsv/md rows flushed per input as they complete

## Output

Every mode prints unified JSON to stdout. The top level carries
`analysis_type`, `tool`, `version`, `contigs`, `methods`, and `verdict`.
See the [JSON Output Contract](reference/json-contract.md) for full examples.

## Documentation

- [Quick Start](installation/quick-start.md)
- [Databases](installation/databases.md)
- [Marker Database](installation/marker-database.md)
- [CLI Reference](usage/cli.md)
- [Identification Modes](usage/modes.md)
- [Metagenomic Profiling](usage/meta.md)
- [Layered Design](architecture/layers.md)
- [JSON Output Contract](reference/json-contract.md)
- 中文版可通过页面语言切换器访问

## License

MIT
