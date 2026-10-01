# gside — Genome Species IDentification Engine

An independent pathogen species identification CLI tool, extracted from the
hermes-bacmap architecture.

## Overview

gside (**G**enome **S**pecies **ID**entification **E**ngine) provides
multi-layer species identification capabilities as a standalone command-line
tool. It can be used independently or integrated into agent platforms
(e.g., hermes-bacmap).

```
Tools        blastn / skani / mash / sourmash / mmseqs2
               ↑
gside        Species ID engine: marker rules + ANI backends + arbitration
               ↑
hermes-bacmap  Hermes bridge: GOM storage / smk orchestration / gbrain memory
```

## Identification Layers

| Layer | CLI mode | Engine | Method | Validated |
|-------|----------|--------|--------|-----------|
| **L1** | `marker` | blastn | Target gene combination rules (38 rules, 81 sequences) | 17/17 |
| **L2** | `panel` | skani | ANI against curated reference panel (291 genomes) | 17/17 |
| **L2** | `mash_refseq` | mash | MinHash distance (RefSeq sketch) | 17/17 |
| **L2** | `sourmash` | sourmash | GTDB gather (LCA classification) | DB required |
| **Arbitration** | `all` | — | Layer consensus (ANI > marker) | ✓ |

## Quick Start

```bash
# Install
pixi install && pip install -e .

# Species identification
gside species contigs.fna --mode marker
gside species contigs.fna --mode panel
gside species contigs.fna --mode all

# Database management
gside db status
gside db setup --tier all
```

## Output

All modes produce a unified JSON format (GOM-compatible contract) containing
`analysis_type`, `method`, `database`, `result`, and `verdict` fields —
directly consumable by downstream systems.

## License

MIT
