# Identification Modes

## Layered Architecture

gside's identification is layered by evidence strength — higher layers
override lower ones:

```
L2  panel (skani) / mash_refseq (mash) / sourmash (GTDB gather)   [authority=2]
L1  marker target gene combination rules                          [authority=1]
```

## L1: marker (Target Gene Combination)

### Principle

Performs a single BLAST scan of contigs against the markers_v2 database
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

## all Mode (Multi-Method Arbitration)

Runs all available methods and arbitrates. A single-method failure does not
abort the run: the error goes into that method's JSON entry and arbitration
uses the methods that succeeded.

```
Layer 2 (panel/mash_refseq/sourmash) beats layer 1 (marker);
within the same layer, high confidence beats lower confidence.
```

```bash
gside species contigs.fna --mode all
```

Output includes all method results plus final verdict:

```json
{
  "verdict": {
    "species": "Campylobacter jejuni subsp. jejuni NCTC 11168",
    "confidence": "high",
    "basis": ["panel"]
  },
  "methods": {
    "marker": { "species": "Campylobacter_jejuni", "..." : "..." },
    "panel": { "result": { "species": "...", "ani": 100.0, "..." : "..." } },
    "mash_refseq": { "result": { "species": "...", "..." : "..." } }
  }
}
```

## Method Selection Guide

| Scenario | Recommended mode | Rationale |
|----------|-----------------|-----------|
| Rapid screening | `marker` | Seconds, zero external dependencies |
| Precise identification | `panel` | ANI gold standard |
| Broad-spectrum scan | `mash_refseq` | Covers all RefSeq |
| GTDB-aligned call | `sourmash` | LCA classification under GTDB |
| Comprehensive | `all` | Multi-method complementary + arbitration |
| Near-relative discrimination | `panel` or `all` | marker may cross-react |
