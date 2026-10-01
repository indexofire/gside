# Identification Modes

## Layered Architecture

gside's identification is layered by evidence strength — higher layers
override lower ones (recorded, not fatal):

```
L3  GTDB-Tk (standard mode, requires heavy environment)     [authority=3]
L2  panel skani / mash_refseq / sourmash                    [authority=2]
L1  marker target gene combination rules                    [authority=1]
```

## L1: marker (Target Gene Combination)

### Principle

Performs a single BLAST scan of contigs against the markers_v2 database
(81 sequences), then matches detected genes against 38 combination rules
(data-driven, `marker_rules.yaml`).

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
- **coverage threshold**: global ≥60% (prevents short-HSP high-identity false fires)

### Use Cases

- Rapid initial screening (minutes)
- No large database needed (mini tier bundled)
- Species-rich genera (Salmonella/E.coli/Campylobacter) should confirm with L2

## L2: panel (ANI Curated Panel)

### Principle

skani performs whole-genome ANI comparison of contigs against a curated
reference panel (291 genomes, 30 enabled pathogen species).

### Thresholds

| ANI | Coverage | Call |
|-----|----------|------|
| ≥95% | ≥70% | species=matched reference, confidence=high |
| 90-95% | ≥70% | confidence=medium (review recommended) |
| <90% | — | Unknown |

### Use Cases

- Precise species confirmation (subspecies resolution)
- Near-relative discrimination (when marker layer cannot distinguish)
- Requires panel tier database (~88MB)

## L2: mash_refseq (MinHash Distance)

### Principle

mash computes MinHash distance between contigs and the entire RefSeq sketch
(~179MB).

### Use Cases

- Broad-spectrum screening (covers all RefSeq prokaryotes)
- Fast (seconds)
- Lower precision than panel ANI (MinHash is an estimate)

## L2: sourmash (GTDB gather)

### Principle

sourmash gather performs LCA classification based on the GTDB taxonomic
framework.

### Use Cases

- Alignment with GTDB taxonomy
- Mixed samples (gather can decompose)
- Requires 3.7GB sourmash GTDB database

## all Mode (Multi-Method Arbitration)

Runs all available methods and arbitrates by layer:

```
ANI layer (panel/mash/sourmash) high-confidence hit > marker layer hit
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
| Comprehensive | `all` | Multi-method complementary + arbitration |
| Near-relative discrimination | `panel` or `all` | marker may cross-react |
