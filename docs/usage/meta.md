# Metagenomic Profiling

`gside scan --meta` profiles a metagenomic sample — FASTQ reads or assembled
contigs — reporting the species present, their abundances, and per-method
confidence. It is the community workflow of `gside scan` (add `--meta`); the
isolate workflow identifies a **single** organism and answers "which species
is this?", while `--meta` answers "what is in this sample, and how much of
each?".

Two methods are available:

| Method | Engine | Input | Evidence |
|--------|--------|-------|----------|
| `kraken2` | Kraken2 + Bracken | FASTQ (single/paired) or contigs FASTA | k-mer read classification + abundance re-estimation |
| `sylph` | sylph | FASTQ reads only | containment ANI against GTDB genomes |

In `all` mode (the default) both run and the verdict merges their species
lists: duplicates keep the row from the higher-confidence method (sylph ANI
≥ 95% outranks a kraken2 fraction), and rows are ordered by abundance.

## Usage

```bash
gside scan --meta <reads.fastq> [reads_2.fastq] [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `reads` (positional) | — | 1 file = single-end reads or contigs FASTA; 2 files = paired FASTQ |
| `-m`, `--mode` | `all` | `kraken2` / `sylph` / `all` |
| `-D`, `--db-dir` | auto | Database root (`$GSIDE_DB_DIR` or `data/db`) |
| `-f`, `--format` | `json` | Output format: `json` / `tsv` / `csv` / `md` |
| `-o`, `--output PATH` | stdout | Write the report to PATH instead of stdout |
| `-s`, `--samples PATH` | — | Batch manifest (`sample_id<TAB>reads[<TAB>reads2]`) |
| `-j`, `--jobs N` | `1` | Concurrent samples when using `--samples` |
| `-n`, `--no-header` | off | Suppress the header row (tsv/csv/md) |
| `-a`, `--auto-setup` | off | Install a missing database tier, then retry that method once |
| `-t`, `--threads N` | `4` | Threads per tool invocation |

## Batch samples

For multiple samples, pass a manifest instead of positional reads:

```bash
gside scan --meta --samples samples.tsv -j 4 -f tsv -o profiles.tsv
```

`samples.tsv` is tab-separated, one sample per line (header optional):

```
sample_id	reads_1	reads_2
sampleA	sampleA_R1.fq	sampleA_R2.fq
sampleB	sampleB.fq
```

Output order always follows the manifest. `-j/--jobs` sets how many samples
are profiled concurrently (each one still uses `-t` threads internally).
JSON output becomes an array, each object carrying a `sample` key; `tsv`/`csv`
rows gain a leading `sample` column.

Input type is detected automatically: two files are treated as paired FASTQ,
a single `.fna`/`.fa`/`.fasta` file as assembled contigs, anything else as
single-end FASTQ.

## Confidence

- **kraken2**: per-species `fraction_total_reads` ≥ 0.10 → `high`, ≥ 0.01 →
  `medium`, else `low` (Bracken re-estimated reads).
- **sylph**: adjusted containment ANI ≥ 95 → `high`, ≥ 93 → `medium`, else
  `low` (same thresholds as the panel ANI layer of `gside scan`).

## Output

JSON (default) — a single document:

```json
{
  "analysis_type": "metagenomic_profiling",
  "tool": "gside",
  "version": "0.1.5",
  "input": {"files": ["reads_1.fastq", "reads_2.fastq"], "type": "fastq_paired"},
  "methods": {
    "kraken2": {"method": "kraken2_bracken", "database": {...}, "result": {
      "top_species": [{"name": "...", "taxid": 562, "reads": 850, "fraction": 0.85, "confidence": "high"}],
      "total_reads": 10000, "classified_reads": 9870, "n_species": 4
    }},
    "sylph": {"method": "sylph_ani", "database": {...}, "result": {
      "species": [{"name": "Escherichia coli str. K-12", "ani": 97.53, "abundance": 0.7812, "coverage": 264.0, "confidence": "high"}],
      "n_species": 3
    }}
  },
  "verdict": {
    "top_species": [{
      "name": "Escherichia_coli", "display_name": "Escherichia coli",
      "strain": "str. K-12", "species_raw": "Escherichia coli str. K-12",
      "method": "sylph", "abundance": 0.7812, "confidence": "high", "ani": 97.53
    }],
    "richness": 1,
    "basis": ["kraken2", "sylph"]
  }
}
```

`abundance` is always a **fraction (0-1)** for both engines (sylph's percentage
is normalised at parse time), and `name` is the same canonical `Genus_species`
label used by `gside scan`; strain detail moves to `strain`, and the raw
tool string to `species_raw`. Species-group placeholders (`Campylobacter sp.`)
are dropped from the merged verdict.

TSV (`-f tsv`) — one row per merged species (`ani` is blank for kraken2 rows):

```
method	name	abundance	confidence	ani
sylph	Escherichia_coli	0.7812	high	97.53
kraken2	Salmonella_enterica	0.20	medium	
```

Like the isolate workflow, `gside scan --meta` always exits 0: a failing method records
its error under `methods.<mode>.error` and the remaining methods continue.
FASTA input is rejected by `sylph` (it requires reads) with a structured
error pointing at `-m kraken2`.

## Choosing a method

| Input | Recommended | Why |
|---|---|---|
| FASTQ (single/paired) | `-m all` (default) | kraken2 + sylph complement each other |
| FASTQ, speed | `-m kraken2` | seconds; sylph adds ANI on top |
| FASTQ, ANI precision | `-m sylph` | containment ANI, species-level |
| Contigs FASTA | `-m kraken2` | sylph cannot read FASTA |

`--meta -m all` requires FASTQ (it includes sylph). For a FASTA sample use
`-m kraken2`, otherwise the command exits 2 before running anything.

## Examples

```bash
# Paired-end reads, both methods (default)
gside scan --meta sample_1.fastq sample_2.fastq

# Kraken2 only, TSV to a file
gside scan --meta sample_1.fastq sample_2.fastq -m kraken2 -f tsv -o profile.tsv

# Contigs / assembled FASTA (kraken2 handles FASTA)
gside scan --meta contigs.fna -m kraken2

# Single-end reads with sylph ANI profiling
gside scan --meta reads.fastq -m sylph

# Batch manifest, 4 samples at once
gside scan --meta --samples samples.tsv -j 4 -f csv -o profiles.csv

# Custom database root
gside scan --meta reads.fastq -D /path/to/db
```

## Performance

| Input | kraken2 | sylph |
|---|---|---|
| 50k paired 150bp reads | ~1s | ~1-2s |
| 10x 1.6Mb genome | ~1s | ~2s |

## Troubleshooting

| Symptom | Fix |
|---|---|
| `sylph requires FASTQ reads` | Use `-m kraken2`, or supply reads |
| `Kraken2 database not found ... (run: gside db setup --tier meta)` | `gside db setup --tier meta`, or add `--auto-setup` |
| `Bracken distribution for 150bp reads not found` | Reinstall the tier; the archive must include `database150mers.kmer_distrib` |
| `--jobs requires --samples` | Add `--samples`, or drop `-j` |
| `duplicate sample id` | Manifest IDs must be unique |
| `paired input must be two FASTQ files` | Do not mix FASTQ and FASTA in one sample |

## Databases

Both metagenomic tiers download automatically and are checksum-verified:

```bash
gside db setup --tier meta          # Kraken2 Standard-8 + Bracken (~5.5GB)
gside db setup --tier meta_sylph    # sylph GTDB r226 .syl2db (~14GB)
```

- **D5_kraken2** (`meta`): Ben Langmead's Standard-8 archive, which bundles
  the Bracken 50-300-mer distributions. Readiness requires `hash.k2d`,
  `taxo.k2d`, `opts.k2d` and `database150mers.kmer_distrib` (gside uses 150bp
  Bracken).
- **D6_sylph** (`meta_sylph`): two-stage `gtdb-r226-c200-dbv2.syl2db`,
  MD5-pinned. Both `.syldb` and `.syl2db` are accepted; `.syl2db` is preferred.
- `gside scan --meta --auto-setup` installs the missing tier on first use.
- `--source <dir>` copies an existing local database instead.
- `gside db status` shows `D5_kraken2` / `D6_sylph` readiness.
