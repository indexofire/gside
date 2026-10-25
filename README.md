# gside

[![CI](https://github.com/indexofire/gside/actions/workflows/test.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![Docs](https://github.com/indexofire/gside/actions/workflows/docs.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![PyPI](https://img.shields.io/pypi/v/gside.svg)](https://pypi.org/project/gside/)
[![Python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)

**English** | [简体中文](README.zh.md)

`gside` is a Python library and command-line tool for **multi-layer microbial
species identification** and **metagenomic profiling**. It wraps established
bioinformatics tools behind reproducible, data-driven rules, arbitrates across
methods, and returns a single structured result per sample.

The library exposes each identification engine independently. The `gside` CLI
adds database provisioning, layered consensus, batch processing, and an MCP
server for agent workflows.

```text
contigs / reads
      │
      ├── marker BLAST ──┐
      ├── ANI (skani) ───┤
      ├── mash ──────────┼── layered arbitration ── species verdict
      ├── sourmash ──────┤
      ├── GTDB-Tk ───────┘
      └── kraken2 / sylph ── community profile
```

## Highlights

- **Library first.** Import individual engines through `gside.api`; use the
  result objects directly in pipelines, notebooks, and tests.
- **Five isolate modes.** Marker genes, curated ANI panel, RefSeq mash,
  sourmash/GTDB, and GTDB-Tk, with L3 > L2 > L1 arbitration.
- **Two community modes.** Kraken2 + Bracken and sylph, or both together.
- **Assembly QC.** CheckM2 completeness/contamination plus GTDB-Tk taxonomy.
- **Machine-friendly.** Stable JSON envelopes, structured missing-database
  errors, `--auto-setup`, machine-readable database status, and an MCP server.
- **Optional gseek engines.** `GSIDE_ENGINE=legacy|gseek|shadow` can route
  marker BLAST, read screening, mash, sourmash, and panel/skani through
  [gseek](https://github.com/indexofire/gseek). Legacy stays the default.

## Installation

`gside` requires **Python 3.12**.

```bash
pip install gside
```

The Python package and the bundled marker database are included. The aligners
and reference tools are not Python dependencies; install them with conda,
Bioconda, or pixi:

```text
blastn, makeblastdb, minimap2, skani, mash, sourmash,
kraken2, bracken, sylph, GTDB-Tk, CheckM2
```

For a complete development environment:

```bash
pixi install
pixi run gside -V
```

Optional gseek engine:

```bash
pip install "gside[gseek]"
# or use the dev pixi environment
```

If tools live outside the active environment, point gside at their directory:

```bash
export GSIDE_PIXI_BIN=/path/to/conda/env/bin
gside db status
```

## Python library

The stable import surface is `gside.api`. Each function returns a typed result
object or a list of dictionaries; `MethodEnvelope`-based results expose
`to_dict()` for JSON serialization.

### Marker-gene identification

```python
from gside.api import identify_multigene

result = identify_multigene("contigs.fna", threads=4)
print(result.species, result.confidence)
print(result.to_dict()["detected_markers"])
```

### ANI against the curated panel or RefSeq mash

```python
from gside.api import identify_by_ani

panel = identify_by_ani("contigs.fna", mode="panel", threads=8)
mash = identify_by_ani("contigs.fna", mode="mash_refseq", threads=8)

print(panel.to_dict()["result"]["species"])
```

`mode` is `"panel"` (skani ANI) or `"mash_refseq"` (mash distance). Pass
`db_dir=` to use a non-default database root.

### Sourmash / GTDB-Tk

```python
from gside.api import identify_by_gtdbtk, identify_by_sourmash

sourmash = identify_by_sourmash("contigs.fna").to_dict()
gtdbtk = identify_by_gtdbtk("contigs.fna", threads=8).to_dict()
```

GTDB-Tk requires a configured `GTDBTK_DATA_PATH`. Sourmash uses the GTDB gather
database from the `sourmash` tier.

### Read screening and metagenomic profiling

```python
from pathlib import Path

from gside.api import identify_by_kraken2, identify_by_sylph, screen_reads
from gside.config import SPECIES_DB_DIR

markers = SPECIES_DB_DIR / "D1_marker" / "markers.fasta"
hits = screen_reads(("sample_R1.fq", "sample_R2.fq"), markers, threads=8)

kraken = identify_by_kraken2(("sample_R1.fq", "sample_R2.fq"), threads=8)
sylph = identify_by_sylph(("sample_R1.fq", "sample_R2.fq"), threads=8)
```

### Assembly quality

```python
from gside.api import validate_genome

validation = validate_genome("contigs.fna", "qc/", threads=8)
print(validation.status)
```

`validate_genome` writes `qc/validation.json` and also returns the result.

### Result contract

| Result family | Shape |
| --- | --- |
| ANI, sourmash, GTDB-Tk, kraken2, sylph | `MethodEnvelope` with `.to_dict()` |
| Marker identification | `MultiGeneResult.to_dict()` |
| Read screening | `list[dict]` rows |
| Assembly validation | `ValidateResult.to_dict()` |

`MethodEnvelope.to_dict()` returns:

```json
{
  "analysis_type": "species_identification",
  "method": "panel",
  "database": {"name": "D2_ani", "version": "..."},
  "result": {"species": "...", "confidence": "high"}
}
```

The CLI JSON contract, including error fields and exit behavior, is documented
in [JSON Output Contract](docs/reference/json-contract.md).

## CLI

The CLI uses the same engines and adds database management and layer
arbitration.

```bash
# Version; no separate "version" subcommand
gside -V
gside --version
gside -V -J              # JSON version object

# One genome (marker mode is bundled and works immediately)
gside scan contigs.fna
gside scan contigs.fna -m marker+sourmash
gside scan contigs.fna -m all --auto-setup

# Isolate reads; paired mates are grouped as one sample
gside scan reads.fq -m marker
gside scan sample_R1.fq sample_R2.fq -m marker

# Metagenomic community
gside scan --meta reads_R1.fq reads_R2.fq

# Assembly quality
gside validate contigs.fna -o qc/

# Batch + streaming table output
gside scan *.fna -m panel+mash_refseq -j 8 -t 4 -o results.tsv
```

### Modes and arbitration

| Mode | Engine | Input | Reference |
| --- | --- | --- | --- |
| `marker` (default) | blastn | FASTA / FASTQ | Bundled marker DB |
| `panel` | skani | FASTA | Curated panel tier |
| `mash_refseq` | mash | FASTA / FASTQ | RefSeq mash tier |
| `sourmash` | sourmash | FASTA / FASTQ | GTDB gather tier |
| `gtdbtk` | GTDB-Tk | FASTA | External GTDB-Tk data |
| `kraken2` | Kraken2 + Bracken | FASTA / FASTQ | Kraken2 tier |
| `sylph` | sylph | FASTQ | sylph tier |
| `all` | all applicable | depends | all required tiers |

`-m` accepts `+` combinations such as `marker+sourmash`. `all` in isolate mode
runs every genome layer; `--meta -m all` merges Kraken2 and sylph.

Arbitration always prefers the higher-confidence layer:

```text
L3  GTDB-Tk                 phylogenetic placement
L2  panel / mash / sourmash comparative genomics (ANI / sketches)
L1  marker                  target-gene presence/absence
```

### Database commands

Human-readable lists use bordered Rich tables; machine formats are stable.

```bash
gside db status
gside db list
gside db status -f json
gside db list -f tsv
gside db setup --tier panel
```

| Tier | Contents | Approx. size |
| --- | --- | --- |
| `mini` | Marker rules + sequences (bundled) | ~1 MB |
| `panel` | Curated skani reference panel | ~130 MB sketch / ~1 GB genomes |
| `mash` | RefSeq MinHash sketch | ~179 MB |
| `sourmash` | GTDB gather database | ~3.9 GB |
| `meta` | Kraken2 Standard-8 index | ~5.5 GB |
| `meta_sylph` | sylph GTDB-r226 database | ~14 GB |

`gside validate` additionally needs CheckM2 (~1.7 GB) and GTDB-Tk
(`GTDBTK_DATA_PATH`, ~57 GB download).

## Output and automation

- `gside scan` emits one JSON document per sample by default; `-f tsv`,
  `-f md`, and `-f csv` are available for tables.
- `gside db status/list` default to a human table; use `-f json|tsv|csv` for
  automation.
- `gside validate` writes `validation.json` and exits non-zero on a failed
  quality verdict.
- Structured database errors include `suggested_command`, `missing_tier`, and
  `estimated_size`.
- `gside mcp` (or `gside-mcp`) starts an MCP stdio server exposing species,
  meta, and database operations.

## Documentation

- [Documentation home](https://indexofire.github.io/gside/)
- [Quick start](docs/installation/quick-start.md)
- [Python API](docs/usage/python-api.md)
- [CLI reference](docs/usage/cli.md)
- [Identification modes](docs/usage/modes.md)
- [Metagenomic profiling](docs/usage/meta.md)
- [JSON contract](docs/reference/json-contract.md)
- [中文文档](docs/index.zh.md)

## Development

```bash
pixi install
pixi run check          # ruff check/format
pixi run typecheck      # pyright
pixi run pytest         # unit tests
pixi run docs-build     # strict MkDocs build
```

The test suite covers the CLI, database formatting, layer arbitration, result
contracts, and optional gseek shadow paths. External-tool tests skip when the
required binary is absent.

## Project layout

```text
src/gside/
  api.py                 public Python facade
  cli.py                 Click CLI
  analysis/              identification engines
  engine/                tool backends and optional gseek bridge
  db.py                  database tiers/downloads
  mcp.py                 MCP server
docs/                    bilingual English/Chinese documentation
data/db/D1_marker/       bundled marker database
tests/                   unit and parity tests
```

## License and ecosystem

MIT. Part of the `indexofire` microbial genomics stack:

| Project | Purpose |
| --- | --- |
| [gseek](https://github.com/indexofire/gseek) | Reproducible DNA alignment infrastructure |
| [gside](https://github.com/indexofire/gside) | Species identification + metagenomic profiling |
| [gapit](https://github.com/indexofire/gapit) | Gene screening + phenotype designation |
| [gmlst](https://github.com/indexofire/gmlst) | MLST / cgMLST typing |
