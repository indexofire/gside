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

**Paired-end FASTQ is one sample.** Positional files whose names carry a mate
marker (`_R1`/`_R2`, `_1`/`_2`, `.1`/`.2`, `-1`/`-2`, optional `_001`,
case-insensitive, `.gz` allowed) are grouped: both mates go to `marker`
(one minimap2 run), `mash_refseq` and `sourmash` (one sketch each), and the
payload reports the first mate as `contigs`. Unmatched files stay separate
samples, so `gside scan A_R1.fq A_R2.fq B_R1.fq B_R2.fq -m marker` is **two**
samples, not four. FASTA files are never paired. In `--meta` mode exactly two
positional files are always one sample; with more files the same mate grouping
applies.

### Reading the verdict

```json
"verdict": {
  "species": "Campylobacter_jejuni",
  "display_name": "Campylobacter jejuni",
  "strain": "subsp. jejuni NCTC 11168",
  "confidence": "high",
  "basis": ["panel"],
  "supporting": ["marker", "sourmash"],
  "agreement": "3/3"
}
```

| Field | Meaning |
|---|---|
| `species` | Canonical `Genus_species` (the machine contract) |
| `display_name` / `strain` | Human form and strain detail from the winning method |
| `confidence` | high / medium / low (per-method thresholds) |
| `basis` | The single method that produced the call (winners only) |
| `supporting` | Other successful methods that agree on `species` |
| `agreement` | `agreeing/successful` — `2/3` means one tier disagreed or failed |

`methods.<mode>.error` explains a failed method; `species`/`display_name` stay
`Unknown` when nothing was called. Field-by-field details:
[JSON Output Contract](../reference/json-contract.md).

### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `files` | positional | — | One or more FASTA/FASTQ files; FASTQ mates are auto-grouped (`_R1/_R2`, `_1/_2`, `.1/.2`); unmatched files are separate samples (`--meta`: 2 files = one sample) |
| `-m`, `--mode` | option | `marker` / `all` | Method (defaults to `marker`, or `all` with `--meta`) |
| `--meta` | flag | off | Treat the input as a metagenomic community |
| `-D`, `--db-dir` | option | auto | Database root (`$GSIDE_DB_DIR` or `data/db`) |
| `-f`, `--format` | option | `json` | `json` / `tsv` / `md` / `csv` |
| `-o`, `--output PATH` | option | stdout | Write the report to PATH (row formats stream as samples complete; json once at the end) |
| `-j`, `--jobs N` | option | `1` | Concurrent samples (isolate batch, positional `--meta` batch, or `--samples`) |
| `-t`, `--threads N` | option | `4` | Threads per tool invocation |
| `-n`, `--no-header` | option | off | Suppress the header row (tsv/csv/md) |
| `-a`, `--auto-setup` | option | off | Install a missing database tier, then retry that method once |
| `-s`, `--samples PATH` | option | — | `--meta` batch manifest: `sample_id<TAB>reads [<TAB>reads2]` |

With more than one sample on a TTY, a progress bar is drawn on stderr. Streamed
rows are handed to that bar: it is suspended for each write and redrawn
afterwards, so JSON/TSV/CSV/MD output appears immediately and stays intact
instead of being overwritten by the next refresh. When stdout is piped or
redirected, no bar is shown and stdout contains data only.

### Examples

```bash
# One genome (isolate or MAG)
gside scan contigs.fna                      # marker
gside scan contigs.fna -m panel -f csv
gside scan contigs.fna -m marker+sourmash
gside scan *.fna -m marker -j 8 -o results.tsv

# Isolate reads: mates are auto-grouped into ONE sample
gside scan sample_R1.fq sample_R2.fq -m marker
gside scan A_R1.fq A_R2.fq B_R1.fq B_R2.fq -m marker     # two samples

# Long-read isolate data (marker treats reads as sequences)
gside scan ont_reads.fastq -m marker

# Metagenomic community
gside scan reads_R1.fq reads_R2.fq --meta                 # kraken2 + sylph
gside scan contigs.fna --meta -m kraken2 -f md
gside scan --meta --samples samples.tsv -j 4 -f csv -o profiles.csv
```

## gside db

Database management. Running `gside db` without a subcommand prints this
help (the same as `gside db -h`) and exits 0.

| Subcommand | Purpose |
|---|---|
| `status` | Readiness of every component (name, tier, path) |
| `list` | Tiers with readiness, size and description |
| `setup` | Install or update one tier |

`status` and `list` accept `-f/--format table|tsv|csv|json` (default `table`);
`setup` accepts `-f/--format text|json` (default `text`). Human output uses
bordered Rich tables; use `json`/`tsv`/`csv` for scripts.

### status

```bash
gside db status
# ┌─────────────┬───────┬────────────┬──────────────────────────┐
# │ COMPONENT   │ READY │ TIER       │ PATH                     │
# ├─────────────┼───────┼────────────┼──────────────────────────┤
# │ D1_marker   │ ready │ mini       │ data/db/D1_marker        │
# │ D4_sourmash │ missing │ sourmash │ data/db/D4_sourmash      │
# └─────────────┴───────┴────────────┴──────────────────────────┘

gside db status -f json    # {component: {ready, path, tier}}
gside db status -f csv     # COMPONENT,READY,TIER,PATH
```

### list

```bash
gside db list              # bordered tier table
gside db list -f tsv       # TIER	READY	SIZE	DESCRIPTION
```

### setup

| Option | Default | Description |
|--------|---------|-------------|
| `-t`, `--tier` | `panel` | `mini|panel|mash|sourmash|meta|meta_sylph|all` |
| `-s`, `--source` | — | Copy an existing local database directory instead of downloading |
| `-f`, `--format` | `text` | `text` / `json` (json routes progress to stderr) |

- `mini` ships with the package (always ready).
- `all` = **panel + mash**; `sourmash` is separate (~3.9GB).
- `meta` (Kraken2 Standard-8 + Bracken, ~5.5GB) and `meta_sylph`
  (GTDB r226 `.syl2db`, ~14GB) download automatically with checksum
  verification.
- CheckM2 and GTDB-Tk are **not** tiers; see
  [Databases](../installation/databases.md#genome-quality-databases).

```bash
gside db setup --tier panel
gside db setup --tier all
gside db setup --tier sourmash
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

### Reading the verdict

```json
"verdict": {
  "species": "Campylobacter_jejuni",
  "display_name": "Campylobacter jejuni",
  "strain": "subsp. jejuni NCTC 11168",
  "confidence": "high",
  "basis": ["panel"],
  "supporting": ["marker", "sourmash"],
  "agreement": "3/3"
}
```

| Field | Meaning |
|---|---|
| `species` | Canonical `Genus_species` (the machine contract) |
| `display_name` / `strain` | Human form and strain detail from the winning method |
| `confidence` | high / medium / low (per-method thresholds) |
| `basis` | The single method that produced the call (winners only) |
| `supporting` | Other successful methods that agree on `species` |
| `agreement` | `agreeing/successful` — `2/3` means one tier disagreed or failed |

`methods.<mode>.error` explains a failed method; `species`/`display_name` stay
`Unknown` when nothing was called. Field-by-field details:
[JSON Output Contract](../reference/json-contract.md).

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
`GTDBTK_DATA_PATH` (or `GTDBDB`); install commands are in
[Databases](../installation/databases.md#genome-quality-databases).

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

## gside mcp

MCP (Model Context Protocol) stdio server: exposes gside as a tool provider
for agents. JSON-RPC 2.0 over stdin/stdout; stdout carries protocol frames
only (progress and notices go to stderr).

```bash
gside mcp        # or the console script: gside-mcp
```

Tools: `species`, `meta`, `db_status`, `db_list`, `db_setup`. Each mirrors its
CLI twin and returns the same JSON text; failures become `isError` envelopes
instead of crashing the server.

## gside species / gside meta (aliases)

Hidden compatibility aliases: `gside species` = `gside scan` without `--meta`;
`gside meta` = `gside scan --meta`. Mode choices stay context-scoped
(isolate / meta). New scripts should use `scan`.

## Environment variables

| Variable | Used by | Meaning |
|---|---|---|
| `GSIDE_DATA_DIR` | all | Data root (default `<repo>/data`, else bundled package data) |
| `GSIDE_DB_DIR` | all | Database root (default `<data>/db`) |
| `GSIDE_PIXI_BIN` | all | Extra `bin/` for `blastn`/`minimap2`/`skani`/`mash`/… |
| `GSIDE_ENGINE` | `scan` | `legacy` (default, retained), `gseek`, or `shadow` for marker/mash/sourmash/panel paths; `shadow` returns the legacy result |
| `GSIDE_SOURMASH_SIG_SHA256` | `db setup` | Override the pinned sourmash signature digest |
| `CHECKM2DB` | `validate` | Path to CheckM2's `uniref100.KO.1.dmnd` |
| `GTDBTK_DATA_PATH` / `GTDBDB` | `validate`, `scan -m gtdbtk` | GTDB-Tk database directory |

## Exit codes

| Command | Code | Meaning |
|---|---|---|
| `scan` | 0 | Always; per-method failures live in `methods.<mode>.error` |
| `scan` | 2 | Usage error (unknown mode, incompatible input, bad flag) |
| `validate` | 0 | Both CheckM2 and GTDB-Tk ran (quality verdict does not matter) |
| `validate` | 1 | `partial` or `error` (a tool/database is missing) |
| `db setup` | 0 / 1 | 1 only when a result message starts with `ERROR` |
| `db status` / `db list` / `--version` | 0 | Informational |

## Upgrading from 0.1.4

- `gside species` and `gside meta` still work as hidden aliases for
  `gside scan` (without / with `--meta`).
- Species labels are now canonical `Genus_species` (`verdict.species`,
  `methods.<mode>.species`); the display form moved to
  `verdict.display_name`, strain detail to `verdict.strain`. The verdict also
  gained `supporting` and `agreement`.
- Meta rows use the same canonical `name`, with `display_name`/`strain` and a
  fraction (0-1) `abundance`.
- `validate` no longer has `-m/--mode`; `-o/--output-dir` is required and its
  `validation.json` gained schema/genome/parameters/software/tools/status/
  errors/created_at.
- `--meta -m all` requires FASTQ (sylph); use `-m kraken2` for contigs.

## gside --version

There is no separate `gside version` subcommand. Use the root `-V/--version`
flag; `-J` adds a JSON envelope.

```bash
gside -V
# gside, version 0.1.5

gside --version --json     # or -V -J
# {"tool": "gside", "version": "0.1.5"}
```
