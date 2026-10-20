# Metagenomic Profiling

`gside meta` profiles a metagenomic sample — FASTQ reads or assembled
contigs — reporting the species present, their abundances, and per-method
confidence. It complements `gside species`, which identifies a **single**
organism from contigs; `meta` answers "what is in this sample, and how much
of each?".

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
gside meta <reads.fastq> [reads_2.fastq] [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `reads` (positional) | — | 1 file = single-end reads or contigs FASTA; 2 files = paired FASTQ |
| `-m`, `--mode` | `all` | `kraken2` / `sylph` / `all` |
| `-D`, `--db-dir` | auto | Database root (`$GSIDE_DB_DIR` or `data/db`) |
| `-f`, `--format` | `json` | Output format: `json` / `tsv` |
| `-o`, `--output PATH` | stdout | Write the report to PATH instead of stdout |
| `-t`, `--threads N` | `4` | Threads per tool invocation |

Input type is detected automatically: two files are treated as paired FASTQ,
a single `.fna`/`.fa`/`.fasta` file as assembled contigs, anything else as
single-end FASTQ.

## Confidence

- **kraken2**: per-species `fraction_total_reads` ≥ 0.10 → `high`, ≥ 0.01 →
  `medium`, else `low` (Bracken re-estimated reads).
- **sylph**: adjusted containment ANI ≥ 95 → `high`, ≥ 93 → `medium`, else
  `low` (same thresholds as the panel ANI layer of `gside species`).

## Output

JSON (default) — a single document:

```json
{
  "analysis_type": "metagenomic_profiling",
  "tool": "gside",
  "version": "0.1.2",
  "input": {"files": ["reads_1.fastq", "reads_2.fastq"], "type": "fastq_paired"},
  "methods": {
    "kraken2": {"method": "kraken2_bracken", "database": {...}, "result": {
      "top_species": [{"name": "...", "taxid": 562, "reads": 850, "fraction": 0.85, "confidence": "high"}],
      "total_reads": 10000, "classified_reads": 100, "n_species": 4
    }},
    "sylph": {"method": "sylph_ani", "database": {...}, "result": {
      "species": [{"name": "...", "ani": 97.53, "abundance": 78.12, "coverage": 264.0, "confidence": "high"}],
      "n_species": 3
    }}
  },
  "verdict": {
    "top_species": [{"name": "...", "method": "sylph", "abundance": 78.12, "confidence": "high", "ani": 97.53}],
    "richness": 3,
    "basis": ["kraken2", "sylph"]
  }
}
```

TSV (`-f tsv`) — one row per merged species (`ani` is blank for kraken2 rows):

```
method	name	abundance	confidence	ani
sylph	Escherichia coli	78.12	high	97.53
kraken2	Salmonella enterica	0.20	medium	
```

Like `gside species`, `gside meta` always exits 0: a failing method records
its error under `methods.<mode>.error` and the remaining methods continue.
FASTA input is rejected by `sylph` (it requires reads) with a structured
error pointing at `--mode kraken2`.

## Examples

```bash
# Paired-end reads, both methods (default)
gside meta sample_1.fastq sample_2.fastq

# Kraken2 only, TSV to a file
gside meta sample_1.fastq sample_2.fastq -m kraken2 -f tsv -o profile.tsv

# Contigs / assembled FASTA (kraken2 handles FASTA, sylph reports the structured error)
gside meta contigs.fna

# Single-end reads with sylph ANI profiling
gside meta reads.fastq -m sylph

# Custom database root
gside meta reads.fastq -D /path/to/db
```

## Databases

Both metagenomic tiers are **manual setup** for now — `gside db setup`
prints an actionable pointer instead of downloading:

```bash
gside db setup --tier meta          # → ERROR: not auto-downloadable yet — see the Kraken2 index page
gside db setup --tier meta_sylph    # → ERROR: not auto-downloadable yet — see the sylph docs
```

### D5_kraken2 (tier `meta`, ~5.5GB)

Download a prebuilt Kraken2 index (e.g. Standard-8) from
<https://benlangmead.github.io/aws-indexes/k2> and extract it into
`data/db/D5_kraken2/` so that `hash.k2d`, `opts.k2d`, and `taxo.k2d` sit
there. Bracken needs the same database plus its `hash.k2d` / `database100mer.kmer_distrib` files.

### D6_sylph (tier `meta_sylph`, ~14GB)

Download a prebuilt sylph `.syl2db` database (e.g. GTDB) from
<https://sylph-docs.github.io/> and place it under `data/db/D6_sylph/`
(the first `*.syl2db` file found is used).

Check readiness with `gside db status` (components `D5_kraken2` /
`D6_sylph`).
