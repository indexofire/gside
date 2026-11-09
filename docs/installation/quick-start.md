# Quick Start

From zero to a species call or a genome-quality report. The marker tier is
bundled, so the first identification takes seconds; optional tiers add
ANI/sketch/read-based methods.

## Workflow summary

| Goal | Command |
|---|---|
| One genome, fastest | `gside scan contigs.fna` |
| One genome, best evidence | `gside scan contigs.fna -m all` |
| Isolate sequencing reads | `gside scan reads.fq -m marker` |
| Isolate paired-end reads | `gside scan S_R1.fq S_R2.fq -m marker` (one sample) |
| Assembly quality (QC) | `gside validate contigs.fna -o taxonomy` |
| Batch with streaming output | `gside scan *.fna -m all -j 8 -o calls.tsv` |

## 1. Prerequisites

- Linux x86_64, Python ≥ 3.12
- [pixi](https://pixi.sh) — provides `blastn`, `minimap2`, `skani`, `mash`,
  `sourmash`, GTDB-Tk (the Python package alone cannot
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
gside -V              # gside, version 0.2.4
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

# Paired-end: R1/R2 (and _1/_2, .1/.2) are ONE sample, not two
gside scan S_R1.fq S_R2.fq -m marker
gside scan A_R1.fq A_R2.fq B_R1.fq B_R2.fq -m marker   # two samples

