# CLI Reference

## gside scan

Single entry point for both workflows: identify the species in one genome
(default) or profile a metagenomic community (`--meta`).

### Usage

```bash
gside scan <files...> [options]
gside scan --meta <reads...> [options]
gside scan --meta --samples samples.tsv [options]
```

### Mode / input matrix

| Mode | `--meta` | FASTA | FASTQ | Notes |
|---|---|---|---|---|
| `marker` | no | ✅ | ✅ | FASTQ: reads mapped with minimap2 (`sr`/`map-ont`); per-gene breadth is the union of read alignments — ~1x genome coverage suffices |
| `panel` | no | ✅ | ❌ | skani ANI |
| `mash_refseq` | no | ✅ | ✅ | reads are sketched with `mash sketch -m 2 -r`; species-level at >=10x, strain calls unreliable |
| `sourmash` | no | ✅ | ✅ | gather on read sketches (~4 min for a 10x 1.6Mb set); containment-ANI based |
| `gtdbtk` | no | ✅ | ❌ | phylogenetic placement |
| `all` | no | ✅ | ❌ | every species tier + arbitration (panel/gtdbtk are assembly-only) |
| `kraken2` | yes | ✅ | ✅ | reads or contigs |
| `sylph` | yes | ❌ | ✅ | reads only |
| `all` | yes | ❌ | ✅ | kraken2+sylph; for contigs use `-m kraken2` |

`-m` accepts a `+`-joined combination (e.g. `marker+sourmash`,
`kraken2+sylph`); duplicates collapse and order is preserved. `all` means
every mode for the current context and cannot be combined with others.

The input type is inferred from the extension (`.fna/.fa/.fasta` vs
`.fq/.fastq`, `.gz` allowed) and validated **before** any tool runs: a
mismatch exits 2 and names the offending mode. `--meta` likewise rejects
isolate modes, and isolate mode rejects `kraken2`/`sylph`. Because `all`
includes sylph, `--meta -m all` needs FASTQ; use `-m kraken2` for contigs.

### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `files` | positional | — | One or more FASTA/FASTQ files (isolate: a batch; `--meta`: one sample of 1-2 files) |
| `-m`, `--mode` | option | `marker` / `all` | Method (defaults to `marker`, or `all` with `--meta`) |
| `--meta` | flag | off | Treat the input as a metagenomic community |
| `-D`, `--db-dir` | option | auto | Database root (`$GSIDE_DB_DIR` or `data/db`) |
| `-f`, `--format` | option | `json` | `json` / `tsv` / `md` / `csv` |
| `-o`, `--output PATH` | option | stdout | Write the report to PATH (row formats stream as inputs complete; json once at the end) |
| `-j`, `--jobs N` | option | `1` | Concurrent isolate inputs, or manifest samples with `--meta --samples` |
| `-t`, `--threads N` | option | `4` | Threads per tool invocation |
| `-n`, `--no-header` | option | off | Suppress the header row (tsv/csv/md) |
| `-a`, `--auto-setup` | option | off | Install a missing database tier, then retry that method once |
| `-s`, `--samples PATH` | option | — | `--meta` batch manifest: `sample_id<TAB>reads [<TAB>reads2]` |

### Examples

```bash
# One genome (isolate or MAG)
gside scan contigs.fna                      # marker
gside scan contigs.fna -m panel -f csv
gside scan contigs.fna -m marker+sourmash
gside scan *.fna -m marker -j 8 -o results.tsv

# Long-read isolate data (marker treats reads as sequences)
gside scan ont_reads.fastq -m marker

# Metagenomic community
gside scan reads_R1.fq reads_R2.fq --meta                 # kraken2 + sylph
gside scan contigs.fna --meta -m kraken2 -f md
gside scan --meta --samples samples.tsv -j 4 -f csv -o profiles.csv
```

## gside species (alias)

`gside species <files...>` is a hidden compatibility alias for `gside scan`
without `--meta`: same options and output, with mode choices limited to
`marker|panel|mash_refseq|sourmash|gtdbtk|all`.

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

Assembly **quality control**: CheckM2 (completeness / contamination) plus
GTDB-Tk classification. It does not run marker rules — use `scan` for
identification.

### Usage

```bash
gside validate <contigs.fasta> -o <output-dir> [options]
```

### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `contigs` | positional | — | Assembled contigs FASTA (must exist) |
| `-o`, `--output-dir` | option | **required** | Directory for `validation.json` (created if missing; a directory, not a file) |
| `-t`, `--threads` | option | `4` | Threads for CheckM2 `--threads` and GTDB-Tk `--cpus` |

Both tools are required for a complete report. Missing databases/binaries are
recorded per tool instead of crashing:

| `status` | Meaning | Exit code |
|----------|---------|-----------|
| `ok` | CheckM2 and GTDB-Tk both ran | 0 |
| `partial` | exactly one tool ran | 1 |
| `error` | neither tool ran | 1 |

A low-quality genome (e.g. `WARNING`) is a **result**, not a command failure:
it still exits 0 when both tools ran. Configure the tools with `CHECKM2DB` and
`GTDBTK_DATA_PATH` (or `GTDBDB`).

### `validation.json`

```json
{
  "schema": "gside.validate/1",
  "genome": {"name": "GCF_000006945.2_ASM694v2_genomic", "path": "/abs/contigs.fna"},
  "parameters": {"threads": 4, "checkm2_db": "/db/checkm2", "gtdbtk_db": "/db/gtdb"},
  "software": {"gside": "0.1.5", "checkm2": "1.0.2", "gtdbtk": "2.4.0"},
  "tools": {"checkm2": "ok", "gtdbtk": "ok"},
  "quality": {"completeness": 99.1, "contamination": 0.5, "verdict": "PASS"},
  "taxonomy": {"classification": "d__Bacteria;p__...;s__Salmonella enterica", "note": ""},
  "status": "ok",
  "errors": [],
  "interpretation": "Genome quality: PASS (completeness=99.1%, contamination=0.5%); GTDB-Tk: ...",
  "created_at": "2026-10-21T14:30:00+00:00"
}
```

### Examples

```bash
gside validate contigs.fna -o ./taxonomy
gside validate contigs.fna -o ./taxonomy -t 8
```

## gside meta (alias)

`gside meta <reads...>` is a hidden compatibility alias for
`gside scan --meta`: same options (including `--samples`), with mode choices
`kraken2|sylph|all`.

## gside --version

```bash
gside --version
# gside, version 0.1.5

gside --version --json     # or -V -J
# {"tool": "gside", "version": "0.1.5"}
```
