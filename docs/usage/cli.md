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
| `-m`, `--mode` | option | `marker` | Identification mode (see below) |
| `-D`, `--db-dir` | option | auto | Database root (`$GSIDE_DB_DIR` or `data/db`) |
| `-f`, `--format` | option | `json` | Output format: `json` / `tsv` / `md` |
| `-n`, `--no-header` | option | off | Suppress the TSV header row (tsv only; md/json keep their structure) |
| `-a`, `--auto-setup` | option | off | Automatically install missing database tiers before identification (progress on stderr) |
| `-j`, `--jobs N` | option | `1` | Number of input files to identify concurrently (output order is always input order) |
| `-t`, `--threads N` | option | `4` | Threads per aligner invocation (blastn `-num_threads` / skani `-t` / mash `-p`) |
| `-o`, `--output PATH` | option | stdout | Write the report to PATH instead of stdout (truncates; tsv/md rows stream to the file as each input completes, flushed per row — a mid-batch failure still persists the prefix; json is written once at the end; stdout then carries no data, a notice goes to stderr) |

### Identification Modes

Six modes are available. `marker` is the default.

| Mode | Description | Database Required |
|------|-------------|-------------------|
| `marker` | Target gene combination rules (blastn) | mini (bundled) |
| `panel` | ANI against curated panel (skani) | panel tier |
| `mash_refseq` | MinHash distance (mash) | mash tier |
| `sourmash` | sourmash GTDB gather | sourmash GTDB db |
| `gtdbtk` | GTDB-Tk phylogenetic placement (classify_wf) | external GTDB-Tk install (`GTDBTK_DATA_PATH`) |
| `all` | Run all available methods + arbitrate | as needed |

Progress bar: multi-file runs show an stderr progress bar (spinner, M/N, elapsed) when stderr is a TTY; suppressed under `--quiet`, single-file, or piped stderr.

In `all` mode a failure in one method does not abort the run: the error is
recorded in that method's JSON entry and the remaining methods continue.
Note that `gtdbtk` is slow (minutes per genome); in `all` mode it runs too —
use an explicit `-m` if you want to skip it.

### Structured error fields (agent-friendly)

When a method fails because its database tier is missing, that method's JSON
entry carries structured fields alongside the `error` string — no regex
parsing needed:

```json
{
  "method": "sourmash",
  "error": "sourmash GTDB database not found in ... (run: gside db setup --tier sourmash)",
  "suggested_command": "gside db setup --tier sourmash",
  "missing_tier": "sourmash",
  "estimated_size": "~3.9GB"
}
```

With `--auto-setup`, gside runs `suggested_command` itself (progress on
stderr), retries that single method once, and recomputes the verdict.

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

# GTDB-Tk phylogenetic placement (slow, needs external gtdbtk + GTDBTK_DATA_PATH)
gside species contigs.fna --mode gtdbtk

# Multi-method arbitration
gside species contigs.fna --mode all

# Specify database directory
gside species contigs.fna --mode panel --db-dir /path/to/dbs

# Auto-install missing database tiers, then identify
gside species contigs.fna --mode all --auto-setup

# Write the report to a file (stdout stays empty; notice goes to stderr)
gside species contigs.fna -m all -o report.json
```

## gside db

Database management subcommand.

### Subcommands

| Subcommand | Description |
|------------|-------------|
| `status` | Show database readiness status |
| `list` | List available tiers |
| `setup` | Install/update databases |

`status` and `list` accept `-f`/`--format [table|tsv|json]` (default `table`).
`tsv` and `json` are machine-readable (e.g. `db status --format tsv` emits a
`COMPONENT\tREADY\tTIER\tPATH` table with `ready`/`missing` booleans).

### setup Options

| Option | Default | Description |
|--------|---------|-------------|
| `-t`, `--tier` | `panel` | mini / panel / mash / sourmash / meta / meta_sylph / all |
| `-s`, `--source` | — | Copy from an existing local database directory instead of downloading |
| `-f`, `--format` | `text` | Output format: `text` / `json` (json routes progress to stderr and prints the result dict as a single JSON document on stdout; exit-code contract unchanged) |

`all` installs **panel + mash only** (~310MB). The sourmash GTDB database
(~3.9GB) is never part of `all`; provision it separately:

```bash
gside db setup --tier sourmash
```

The metagenomic tiers `meta` (Kraken2 index, ~5.5GB) and `meta_sylph`
(sylph GTDB database, ~14GB) are **manual setup** — `db setup` prints a
pointer to the upstream downloads instead of fetching them. See
[Metagenomic Profiling](meta.md) for placement instructions.

### Examples

```bash
gside db status
gside db status --format tsv    # machine-readable
gside db setup --tier all
gside db setup --tier sourmash   # not included in --tier all (~3.9GB)
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

## gside meta

Metagenomic species profiling from FASTQ reads or assembled contigs
(Kraken2+Bracken and sylph). See [Metagenomic Profiling](meta.md).

### Usage

```bash
gside meta <reads.fastq> [reads_2.fastq] [options]
```

### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `reads` | positional | — | 1 file = single-end reads or contigs FASTA; 2 files = paired FASTQ |
| `-m`, `--mode` | option | `all` | Profiling method: `kraken2` / `sylph` / `all` (all merges the species lists) |
| `-D`, `--db-dir` | option | auto | Database root (`$GSIDE_DB_DIR` or `data/db`) |
| `-f`, `--format` | option | `json` | Output format: `json` / `tsv` |
| `-o`, `--output PATH` | option | stdout | Write the report to PATH instead of stdout (json written once at the end; stdout then carries no data, a notice goes to stderr) |
| `-t`, `--threads N` | option | `4` | Threads per tool invocation |

### Examples

```bash
gside meta sample_1.fastq sample_2.fastq
gside meta sample_1.fastq sample_2.fastq -m kraken2 -f tsv -o profile.tsv
gside meta contigs.fna          # kraken2 only (sylph requires reads)
gside meta reads.fastq -m sylph
```

## gside --version

```bash
gside --version
# gside, version 0.1.4

gside --version --json     # or -V -J
# {"tool": "gside", "version": "0.1.4"}
```
