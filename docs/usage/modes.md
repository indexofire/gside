# Identification Modes

## Layered Architecture

gside's identification is layered by evidence strength — higher layers
override lower ones:

```
L3  gtdbtk (GTDB-Tk phylogenetic placement)                        [authority=3]
L2  panel (skani) / mash_refseq (mash) / sourmash (GTDB gather)   [authority=2]
L1  marker target gene combination rules                          [authority=1]
```

## L1: marker (Target Gene Combination)

### Principle

Performs a single BLAST scan of contigs against the markers database
(81 sequences), then matches detected genes against 38 combination rules
(data-driven, `marker_rules.yaml`).

### Thresholds

| Gate | Value | Description |
|------|-------|-------------|
| Global minimum identity | 85% | Hits below this are discarded outright |
| Global minimum coverage | 60% | Prevents short-HSP high-identity false fires |
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

### Principle

mash computes MinHash distance between contigs and the entire RefSeq sketch
(~179MB).

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

### Principle

sourmash gather performs LCA classification based on the GTDB taxonomic
framework.

### Thresholds

| f_unique_weighted | Call |
|-------------------|------|
| ≥0.90 | confidence=high (unless a mixture is flagged) |
| 0.70-0.90 | confidence=medium |
| <0.70 | Mixed/Unknown |

When two or more gathered lineages each reach `f_unique_weighted` ≥ 0.10,
the result carries a `possible_mixture` flag and the high-confidence branch
is withheld.

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
within the same layer, high confidence beats lower confidence.
```

`gtdbtk` is included in `all` but is slow (minutes per genome) — use an
explicit `-m` to skip it when speed matters.

```bash
gside species contigs.fna --mode all
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
| GTDB-aligned call | `sourmash` | LCA classification under GTDB |
| Authoritative placement | `gtdbtk` | Phylogenetic placement in the GTDB tree (slow) |
| Comprehensive | `all` | Multi-method complementary + arbitration |
| Near-relative discrimination | `panel` or `all` | marker may cross-react |

## Metagenomic Profiling (`gside meta`)

For mixed-community samples, `gside meta` provides taxonomic profiling using
two complementary classifiers. See [Metagenomic Profiling](meta.md) for details.

| Mode | Engine | Input | Confidence basis |
|------|--------|-------|-----------------|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | Bracken fraction (≥10% high, ≥1% medium) |
| `sylph` | sylph | FASTQ reads only | Adjusted ANI (≥95 high, ≥93 medium) |
| `all` | Both | FASTQ (both) / FASTA (kraken2) | Merged profile |
