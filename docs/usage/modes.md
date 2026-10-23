# Identification Modes

## Layered Architecture

gside's identification is layered by evidence strength — higher layers
override lower ones:

```
L3  gtdbtk (GTDB-Tk phylogenetic placement)                        [authority=3]
L2  panel (skani) / mash_refseq (mash) / sourmash (GTDB gather)   [authority=2]
L1  marker target gene combination rules                          [authority=1]
```

## Combining modes

`-m` accepts one mode or a `+`-joined combination; `all` expands to every
mode for the current context. `--meta` switches to the community set
(`kraken2`/`sylph`), where `all` means both.

```bash
gside scan contigs.fna -m marker+sourmash     # just those two, arbitrated
gside scan contigs.fna -m all                 # every genome tier
gside scan reads.fq --meta -m kraken2+sylph   # both community methods
```

| Mode | FASTA | FASTQ | Database |
|---|---|---|---|
| `marker` | ✅ | ✅ | mini (bundled) |
| `panel` | ✅ | ❌ | panel (~130MB sketch) |
| `mash_refseq` | ✅ | ✅ | mash (~179MB) |
| `sourmash` | ✅ | ✅ | sourmash (~3.9GB) |
| `gtdbtk` | ✅ | ❌ | external GTDB-Tk DB |
| `all` (isolate) | ✅ | ❌ | every genome tier |
| `kraken2` (`--meta`) | ✅ | ✅ | meta (~5.5GB) |
| `sylph` (`--meta`) | ❌ | ✅ | meta_sylph (~14GB) |
| `all` (`--meta`) | ❌ | ✅ | kraken2+sylph |

`all` includes `panel`/`gtdbtk` (assembly-only), so a read set needs an
explicit combination such as `-m marker+sourmash` or `--meta`.

## L1: marker (Target Gene Combination)

### Principle

**FASTA**: one BLAST scan of the contigs against the markers database
(139 sequences). **FASTQ**: reads are mapped with minimap2 (`-x sr` for short
reads, `-x map-ont` for long) and each gene's breadth is the union of the
reference intervals covered by its primary alignments, divided by gene
length (the same aggregate-breadth recipe as gapit's reads mode). Detected
genes are then matched against the 42 combination rules
(data-driven, `marker_rules.yaml`).

Paired-end mates are detected from the file names (`_R1`/`_R2`, `_1`/`_2`,
`.1`/`.2`) and mapped in one minimap2 run, so per-gene breadth is the union of
both mates; they count as **one** sample.

Reads only need ~1x genome coverage: 150 bp reads at 1x already reach
60-80% breadth for 400-1000 bp marker genes (10x reaches ~99%). Marker
genes are not assembled first, so a read set is never turned into contigs.

### Thresholds

| Gate | Value | Description |
|------|-------|-------------|
| Global minimum identity | 85% | Hits below this are discarded outright |
| Global minimum coverage | 60% | Per-gene breadth (union of aligned intervals, reads or contigs) |
| Rule minimum identity | 90% default | Per-rule override via `min_identity` |
| High confidence | average identity ≥ 90% | `confidence=high`, otherwise `medium` |

### Rule Mechanism

```yaml
- species: Campylobacter_jejuni
  genes: [mapa, hipo, cadf]     # candidate markers
  min_hits: 2                     # at least 2 must hit
  min_identity: 90                # each gene ≥90% identity
  exclude_genes: [ceue]           # exclusion guard: ceuE present → prefer C_coli
```

### Cross-Reaction Guards

- **exclude_genes**: exclusion rules prevent near-relative misidentification
  (e.g., C.jejuni vs C.coli)
- **tlh near-relative guard**: single-gene <90% hit → suppressed to Unknown
  (V.alginolyticus case study)
- **coverage gate**: global ≥60% (prevents short-HSP high-identity false fires)

### Use Cases

- Rapid initial screening (minutes)
- No large database needed (mini tier bundled)
- Species-rich genera (Salmonella/E.coli/Campylobacter) should confirm with L2

## L2: panel (ANI Curated Panel)

### Principle

skani performs whole-genome ANI comparison of contigs against a curated
reference panel (291 genomes, 30 enabled pathogen species).

### Thresholds

| ANI | Aligned fraction | Call |
|-----|------------------|------|
| ≥95% | ≥0.65 | species=matched reference, confidence=high |
| 93-95% | ≥0.65 | confidence=medium (review recommended) |
| <93% | — | Unknown |

### Use Cases

- Precise species confirmation (subspecies resolution)
- Near-relative discrimination (when marker layer cannot distinguish)
- Requires panel tier database (~130MB)

## L2: mash_refseq (MinHash Distance)

### Inputs: reads or assemblies

Assemblies are compared directly. For FASTQ, gside runs
`mash sketch -m 2 -r` first (both mates of a pair in one sketch) (drop single-copy k-mers, estimate genome size
from k-mer content) and compares that read sketch to the reference sketch.
The Mash paper reports correct **species** calls from unassembled reads
(including noisy Nanopore), but the best-hit **strain** was often wrong;
mash's own docs recommend 10x-100x coverage, and the distance inflates with
coverage/error rate (D ≈ 1 - ANI is calibrated on assemblies). Treat
`mash_refseq` on reads as a species-level screen, not an assembly-grade ANI.
Containment screening of reads is `mash screen`, a different tool.

### Principle

mash computes MinHash distance between the query (assembly, or a read sketch
for FASTQ) and the entire RefSeq sketch (~179MB).

### Thresholds

| Mash identity | Call |
|---------------|------|
| ≥0.97 | confidence=high |
| 0.90-0.97 | confidence=medium |
| <0.90 | Unknown |

### Use Cases

- Broad-spectrum screening (covers all RefSeq prokaryotes)
- Fast (seconds)
- Lower precision than panel ANI (MinHash is an estimate)

## L2: sourmash (GTDB gather)

### Inputs: reads or assemblies

`sourmash sketch` accepts FASTQ, and `gather` is explicitly designed for
metagenomes — "assembled and unassembled" reads both work (sourmash docs;
Irber et al. 2022). On reads this is a **containment profile**: the fraction
of read k-mers in the best reference, not an assembly-vs-assembly ANI.
Coverage gaps and sequencing errors lower sensitivity, and short reads are
not individually classified (gather + read mapping is the recommended
follow-up). Measured: a 10x 150bp read set for a 1.65Mb genome takes
~4 minutes (vs ~40s for the assembly) and returns the same call/ANI as the
assembly (98.2% containment ANI).

### Principle

sourmash gather builds a minimum set cover (containment) over the GTDB
taxonomic framework.

### Thresholds

Confidence comes from the gather containment ANI (`average_containment_ani`),
using the same 95/93 bands as the panel layer:

| Containment ANI | Call |
|-----------------|------|
| ≥95 | confidence=high (unless a mixture is flagged) |
| 93-95 | confidence=medium |
| <93 | Mixed/Unknown |

If the gather CSV has no ANI column (old exports), the fallback bands are
`f_unique_weighted` ≥0.90 → high, ≥0.70 → medium.

When two or more **distinct species** each reach `f_unique_weighted` ≥ 0.10
with ANI ≥ 93, the result carries a `possible_mixture` flag and the
high-confidence branch is withheld. Multiple strains of one species do not
trigger it.

The result also exposes `ani` and a per-match `gather_partition`
(`species`, `f_unique_weighted`, `ani`).

### Use Cases

- Alignment with GTDB taxonomy
- Mixed samples (gather decomposes the query and can flag `possible_mixture`)
- Requires the sourmash GTDB database

## L3: gtdbtk (GTDB-Tk Phylogenetic Placement)

### Principle

GTDB-Tk `classify_wf` places the assembly into the GTDB reference tree
(`classify_wf --genome_dir ... --out_dir ... --database $GTDBTK_DATA_PATH`).
The classification string is parsed from `gtdbtk.bac120.summary.tsv`
(fallback: `gtdbtk.ar122.summary.tsv`); the species is the last `s__` rank.

### Thresholds

| Classification | Call |
|----------------|------|
| contains an `s__` rank | species extracted, confidence=high (authoritative) |
| no `s__` rank | Unknown, confidence=low |

GTDB-Tk is treated as the most authoritative layer: a `gtdbtk` call beats any
L2 ANI/MinHash call regardless of the lower layer's confidence.

### Requirements

- External GTDB-Tk install (not provisioned by `gside db setup`): the
  `gtdbtk` binary on PATH plus the reference database via `GTDBTK_DATA_PATH`
  (or `GTDBDB`) — see https://ecogenomics.github.io/GTDBTk/
- Slow: minutes per genome (MAFFT + pplacer); pass `-t/--threads` through

### Use Cases

- Authoritative species calls aligned to the GTDB tree
- Resolving disagreements between L1/L2 methods
- Publication-grade taxonomic assignment

## all Mode (Multi-Method Arbitration)

Runs all available methods and arbitrates. A single-method failure does not
abort the run: the error goes into that method's JSON entry and arbitration
uses the methods that succeeded. Missing-database errors carry structured
fields (`missing_tier`, `suggested_command`, `estimated_size`); with
`--auto-setup` gside installs the tier itself and retries that method once.

```
Layer 3 (gtdbtk) beats layer 2 (panel/mash_refseq/sourmash) beats layer 1 (marker);
within the same layer, high > medium > low confidence.
```

`gtdbtk` is included in `all` but is slow (minutes per genome) — use an
explicit `-m` to skip it when speed matters.

```bash
gside scan contigs.fna -m all
```

Output includes all method results plus final verdict:

```json
{
  "verdict": {
    "species": "Campylobacter_jejuni",
    "display_name": "Campylobacter jejuni",
    "strain": "subsp. jejuni NCTC 11168",
    "confidence": "high",
    "basis": ["panel"],
    "supporting": ["marker"],
    "agreement": "2/3"
  },
  "methods": {
    "marker": { "species": "Campylobacter_jejuni", "..." : "..." },
    "panel": { "result": { "species": "Campylobacter_jejuni", "species_raw":
      "Campylobacter jejuni subsp. jejuni NCTC 11168", "ani": 100.0, "..." : "..." } },
    "mash_refseq": { "result": { "species": "...", "..." : "..." } },
    "sourmash": { "result": { "species": "...", "..." : "..." } },
    "gtdbtk": { "result": { "species": "...", "classification": "...", "..." : "..." } }
  }
}
```

`verdict.species` is always the canonical `Genus_species` label
(`verdict.display_name` carries the spaced form, `verdict.strain` the strain).
`supporting` lists other methods that agree, `agreement` summarises
`agreeing/successful`; `basis` stays the single winning method.

## Method Selection Guide

| Scenario | Recommended mode | Rationale |
|----------|-----------------|-----------|
| Rapid screening | `marker` | Seconds, zero external dependencies |
| Precise identification | `panel` | ANI gold standard |
| Broad-spectrum scan | `mash_refseq` | Covers all RefSeq |
| GTDB-aligned call | `sourmash` | gather/containment under GTDB |
| Authoritative placement | `gtdbtk` | Phylogenetic placement in the GTDB tree (slow) |
| Comprehensive | `all` | Multi-method complementary + arbitration |
| Near-relative discrimination | `panel` or `all` | marker may cross-react |

## Metagenomic Profiling (`gside scan --meta`)

For mixed-community samples, `gside scan --meta` provides taxonomic profiling using
two complementary classifiers. See [Metagenomic Profiling](meta.md) for details.

| Mode | Engine | Input | Confidence basis |
|------|--------|-------|-----------------|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | Bracken fraction (≥10% high, ≥1% medium) |
| `sylph` | sylph | FASTQ reads only | Adjusted ANI (≥95 high, ≥93 medium) |
| `all` | Both | FASTQ only (use `-m kraken2` for FASTA) | Merged profile |
