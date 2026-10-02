# gside — Genome Species IDentification Engine

A standalone CLI tool for pathogen species identification from assembled contigs.

## Overview

gside (**G**enome **S**pecies **ID**entification **E**ngine) provides
multi-layer species identification as a self-contained command-line tool.
It wraps standard bioinformatics tools (blastn, skani, mash, sourmash) behind
data-driven rules, arbitrates across methods, and prints a single JSON
verdict per run.

```
Bioinformatics tools   blastn / skani / mash / sourmash / mmseqs2
                         ↑
gside engine           marker rules + ANI backends + arbitration
```

## Identification Layers

| Layer | CLI mode | Engine | Method | Notes |
|-------|----------|--------|--------|-------|
| **L1** | `marker` | blastn | Target gene combination rules (38 rules, 81 sequences) | Mini database bundled |
| **L2** | `panel` | skani | ANI against curated reference panel (291 genomes) | Requires panel tier |
| **L2** | `mash_refseq` | mash | MinHash distance (RefSeq sketch) | Requires mash tier |
| **L2** | `sourmash` | sourmash | GTDB gather (LCA classification) | Requires GTDB database |
| **Arbitration** | `all` | — | Run all methods, layer consensus (ANI layer beats marker) | Best of both layers |

## Quick Start

```bash
# Install
pixi install
pip install -e .

# Species identification
gside species contigs.fna --mode marker

# Database management
gside db status
gside db setup --tier all
```

## Output

Every mode prints unified JSON to stdout. The top level carries
`analysis_type`, `tool`, `version`, `contigs`, `methods`, and `verdict`.
Inside `methods`, the marker result is flat (species and confidence sit at
the top level of the method object), while the ANI modes (`panel`,
`mash_refseq`) and `sourmash` nest their findings under a `result` object.
See the [JSON Output Contract](reference/json-contract.md) for full examples.

## License

MIT
