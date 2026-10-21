# Quick Start

From zero to a species call, a community profile, or a genome-quality
report. The marker tier is bundled, so the first identification takes
seconds; optional tiers add ANI/sketch/read-based methods.

## Workflow summary

| Goal | Command |
|---|---|
| One genome, fastest | `gside scan contigs.fna` |
| One genome, best evidence | `gside scan contigs.fna -m all` |
| Isolate sequencing reads | `gside scan reads.fq -m marker` |
| Community from reads | `gside scan R1.fq R2.fq --meta` |
| Community from contigs | `gside scan contigs.fna --meta -m kraken2` |
| Assembly quality (QC) | `gside validate contigs.fna -o taxonomy` |
| Batch with streaming output | `gside scan *.fna -m all -j 8 -o calls.tsv` |

## 1. Prerequisites

- Linux x86_64, Python ≥ 3.12
- [pixi](https://pixi.sh) — provides `blastn`, `minimap2`, `skani`, `mash`,
  `sourmash`, `kraken2`, `bracken`, `sylph` (the Python package alone cannot
  install these binaries)

## 2. Install

```bash
git clone https://github.com/indexofire/gside.git
cd gside

pixi install          # bioinformatics binaries
pip install -e .      # gside itself (or: pip install gside)
```

If you installed gside from PyPI and keep the binaries elsewhere, point
`GSIDE_PIXI_BIN` at that `bin/` directory.

## 3. Verify

```bash
gside -V              # gside, version 0.1.5
gside db status       # mini (D1_marker) must be ✅
```

`db status` shows every component, its tier and path. Anything ❌ is missing
and can be installed with `gside db setup --tier <tier>`.

## 4. First identification (seconds)

```bash
gside scan contigs.fna
```

- Default mode is `marker` (bundled mini tier, blastn).
- `scan` infers FASTA vs FASTQ from the extension and rejects a mode that
  cannot consume it — before running anything.
- Output is one JSON document on stdout: `methods` per method plus a
  `verdict` (canonical `Genus_species`, confidence, basis, supporting,
  agreement).

Common single-genome runs:

```bash
gside scan contigs.fna -m marker                 # target genes
gside scan contigs.fna -m marker+sourmash        # combine methods
gside scan contigs.fna -m all                    # every genome tier + arbitration
gside scan contigs.fna -m all -j 8 -o run.tsv    # batch + streaming TSV
```

## 5. Reads

```bash
# Pure culture / isolate sequencing reads
gside scan reads.fq -m marker                    # minimap2, ~1x coverage is enough
gside scan reads.fq -m marker+mash_refseq

# Metagenomic community (add --meta)
gside scan reads_R1.fq reads_R2.fq --meta        # kraken2 + sylph
gside scan contigs.fna --meta -m kraken2         # metagenomic contigs
```

`--meta` switches to the community workflow; only `kraken2`, `sylph` and
`all` are accepted there. Without `--meta`, isolate modes run instead.

## 6. Databases

| Tier | Contents | Size | Install |
|------|----------|------|---------|
| `mini` | marker rules + sequences + BLAST db | ~1MB | bundled (always ready) |
| `panel` | curated 291-genome panel (skani) | ~130MB sketch / ~1GB genomes | `gside db setup --tier panel` |
| `mash` | RefSeq MinHash sketch | ~179MB | `gside db setup --tier mash` |
| `sourmash` | GTDB r226 reps + lineages | ~3.9GB | `gside db setup --tier sourmash` |
| `meta` | Kraken2 Standard-8 + Bracken | ~5.5GB | `gside db setup --tier meta` |
| `meta_sylph` | sylph GTDB r226 (`.syl2db`) | ~14GB | `gside db setup --tier meta_sylph` |
| `all` | panel + mash | ~310MB | `gside db setup --tier all` |

Downloads are checksum-verified (panel SHA256, mash MD5, sourmash
SHA256, Kraken2 per-file MD5, sylph MD5). `--auto-setup` installs a missing
tier on first use and retries the method once:

```bash
gside scan contigs.fna -m panel --auto-setup
```

## 7. Genome quality (`gside validate`)

`validate` answers a different question: **is the assembly good, and where
does GTDB-Tk place it?** It runs CheckM2 + GTDB-Tk (no marker rules):

```bash
gside validate contigs.fna -o ./taxonomy          # required output directory
gside validate contigs.fna -o ./taxonomy -t 8
```

Prerequisites (external, not `db setup` tiers):

- CheckM2 binary + database (~1.7GB) — set `CHECKM2DB` to the
  `uniref100.KO.1.dmnd` file.
- GTDB-Tk binary + R232 database (~57GB download, ~98GB disk, ~140GB RAM
  for bacteria) — set `GTDBTK_DATA_PATH`.

The report (`validation.json`) carries genome name/path, tool versions,
parameters, CheckM2 completeness/contamination + verdict,
GTDB-Tk classification, per-tool status, errors and a timestamp.
Exit code is 0 only when both tools ran; a low-quality genome is still
exit 0. See [Databases](databases.md#genome-quality-databases) for the
install commands.

## 8. Batch and formats

`scan` accepts many inputs and four output formats:

```bash
gside scan *.fna -m marker -j 8 -f tsv -o calls.tsv
gside scan *.fna -m marker -f csv
gside scan contigs.fna -m marker -f md
```

- `-j/--jobs` parallelizes isolate inputs; output order is input order.
- `-o` writes the report to a file (row formats stream per input).
- `-n/--no-header` suppresses the header row.
- JSON is a single object for one input, an array for several; `meta`
  batches (below) add a `sample` key.

Metagenomic batches use a manifest:

```bash
gside scan --meta --samples samples.tsv -j 4 -f tsv -o profiles.tsv
```

```
sample_id	reads_1	reads_2
sampleA	sampleA_R1.fq	sampleA_R2.fq
sampleB	sampleB.fq
```

## 9. Troubleshooting

| Symptom | Fix |
|---|---|
| `blastn not found` | Run inside `pixi run -e dev …`, or set `GSIDE_PIXI_BIN` |
| `mode 'panel' does not accept FASTQ input` | Use `--meta -m kraken2`, or an isolate-capable mode (`marker`, `mash_refseq`, `sourmash`) |
| `mode 'kraken2' is not available without --meta` | Add `--meta` |
| `--meta -m all` complains about sylph/FASTA | `all` = kraken2+sylph; use `-m kraken2` for contigs |
| `db status` shows ❌ | `gside db setup --tier <tier>` (or `--auto-setup` on `scan`) |
| `validate` exits 1 with `partial` | One of CheckM2/GTDB-Tk is missing; read `tools`/`errors` in `validation.json` |

## Next Steps

- [CLI Reference](../usage/cli.md) — every command, option, exit code
- [Identification Modes](../usage/modes.md) — per-mode principles and inputs
- [Metagenomic Profiling](../usage/meta.md) — community workflow
- [JSON Output Contract](../reference/json-contract.md) — field-by-field
