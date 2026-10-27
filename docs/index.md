# gside documentation

`gside` is a Python library and CLI for multi-layer microbial species
identification and metagenomic profiling. This documentation covers the public
Python API, the `gside` command line, reference databases, identification
modes, and the JSON contract.

## Start here

| I want to… | Read |
| --- | --- |
| Install gside and run a first sample | [Quick start](installation/quick-start.md) |
| Call identification engines from Python | [Python API](usage/python-api.md) |
| Learn every CLI flag | [CLI reference](usage/cli.md) |
| Understand marker/ANI/sketch/GTDB modes | [Identification modes](usage/modes.md) |
| Profile a metagenomic community | [Metagenomic profiling](usage/meta.md) |
| Prepare databases | [Databases](installation/databases.md) |
| Use my own marker database | [Custom marker databases](usage/custom-databases.md) |
| Consume JSON output safely | [JSON output contract](reference/json-contract.md) |
| Migrate an engine to gseek | [gseek migration](development/gseek-migration.md) |

## What gside does

```text
assembled contigs ── marker / panel / mash / sourmash / GTDB-Tk ── verdict
FASTQ reads ──────── marker screening / kraken2 / sylph ─────────── profile
assembled contigs ── CheckM2 + GTDB-Tk ─────────────────────────── assembly QC
```

Two interfaces share the same engines:

- **Python library** — `gside.api` exposes individual identification functions
  and typed results for notebooks, pipelines, and tests.
- **CLI** — `gside scan` adds layered arbitration, database provisioning,
  batch execution, table/JSON output, and an MCP server.

## Installation at a glance

```bash
pip install gside          # Python 3.12; bundled marker DB
pixi install               # dev environment with all external tools
```

External tools such as `blastn`, `skani`, `mash`, `sourmash`, `kraken2`,
`sylph`, and GTDB-Tk are installed separately. See
[Quick start](installation/quick-start.md) and [Databases](installation/databases.md).

## Library example

```python
from gside.api import identify_by_ani, identify_multigene

marker = identify_multigene("contigs.fna", threads=4)
panel = identify_by_ani("contigs.fna", mode="panel", threads=8)

print(marker.species, marker.confidence)
print(panel.to_dict()["result"]["species"])
```

The full function/result table is in the [Python API](usage/python-api.md).

## CLI example

```bash
gside -V
gside scan contigs.fna -m all --auto-setup
gside scan --meta reads_R1.fq reads_R2.fq
gside validate contigs.fna -o qc/
gside db status
gside db list
```

`gside db status` and `gside db list` render bordered Rich tables in a terminal.
Use `-f json`, `-f tsv`, or `-f csv` when another program will read them.

## Design principles

- **Rules over ad-hoc thresholds.** Marker combinations, ANI confidence bands,
  and layer arbitration are data-driven and testable.
- **One structured result.** Every method emits a JSON envelope; the CLI adds
  consensus without discarding per-method evidence.
- **Fail clearly.** Missing tools/databases produce actionable errors, including
  a suggested setup command and estimated download size.
- **Reproducible by default.** Database versions, tool versions, and method
  parameters are recorded in results.
- **gseek by default.** Marker BLAST, read screening, mash, sourmash, and
  panel/skani run through gseek; `GSIDE_ENGINE=legacy|shadow` is kept for
  debugging and parity checks.

## Project layout

```text
src/gside/
  api.py                 public Python facade
  cli.py                 Click CLI
  analysis/              identification engines
  engine/                backends + gseek bridge
  db.py                  database tiers and downloads
  mcp.py                 MCP server
docs/                    English/Chinese documentation
data/db/D1_marker/       bundled marker database
tests/                   unit and parity tests
```

## License

MIT. See [README](https://github.com/indexofire/gside) for the project overview
and ecosystem links.
