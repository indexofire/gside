# Layered Design

## Architecture Position

gside's position in the overall ecosystem:

```
┌─────────────────────────────────────────────────────────┐
│                  Hermes Agent (LLM)                     │
├─────────────────────────────────────────────────────────┤
│                                                         │
│  hermes-bacmap (bridge layer)                           │
│  42 tools · GOM · smk orchestration · gbrain · skills   │
│                                                         │
├──────────────────────┬──────────────────────────────────┤
│                      │                                  │
│  gside (this tool)   │  gapit · gmlst                   │
│  Species ID engine   │  Gene screening · MLST typing    │
│  marker + ANI + arb  │                                  │
│                      │                                  │
├──────────────────────┴──────────────────────────────────┤
│                                                         │
│  Bioinformatics tools                                   │
│  blastn · skani · mash · sourmash · mmseqs2 · minimap2  │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

## Internal Structure

```
gside/
├── src/gside/
│   ├── cli.py              CLI entry point (species / db subcommands)
│   ├── config.py           Runtime config (data dirs / binary discovery)
│   ├── db.py               Database management (status/setup/list)
│   ├── engine/             Algorithm abstraction layer
│   │   ├── hits.py         Unified Hit dataclass (BLAST tabular + PAF)
│   │   ├── registry.py     Backend registry (name→class, lazy loading)
│   │   ├── read_mapper.py  ReadMapper facade (BWA/Minimap2)
│   │   └── backends/       Swappable backends
│   │       ├── blast.py    blastn/blastp/blastx/tblastn
│   │       ├── skani.py    ANI search
│   │       ├── kmer.py     mash + sourmash backends
│   │       ├── minimap2.py PAF alignment
│   │       ├── kma.py      Exact alignment
│   │       └── mmseqs2.py  Linear-time clustering
│   └── analysis/
│       ├── multigene_identifier.py   L1 target gene combination
│       ├── ani_identifier.py         L2 ANI (panel/skani_gtdb/mash)
│       ├── sourmash_identifier.py    L2 sourmash gather
│       ├── taxonomic_validator.py    L3 GTDB-Tk (code migrated)
│       └── species_canon.py          Species name canonicalization
├── data/
│   ├── reference/species/  marker_rules.yaml + markers_v2.fasta
│   ├── panel_manifest/     panel_accessions.tsv + metadata.tsv
│   └── db/                 Large databases (gitignored)
│       ├── refseq_panel/   skani sketch (291 genomes)
│       └── mash_refseq/    mash sketch (RefSeq)
└── tests/
```

## Design Principles

| Principle | Implementation |
|-----------|---------------|
| **Rules as data** | 38 species rules in YAML; new species without code changes |
| **Swappable backends** | Engine abstraction layer; blastn↔minimap2 by parameter |
| **JSON contract** | All output in unified GOM-compatible format |
| **Graceful degradation** | Missing database/binary returns clear error, never crashes |
| **Independent operation** | No hermes-bacmap dependency; standalone install and use |

## Relationship with hermes-bacmap

gside was extracted from hermes-bacmap's architecture (V0.9 refactor):

| Capability | Original location | Current location |
|-----------|-------------------|-----------------|
| Target gene identification | bacmap/analysis/multigene | gside |
| ANI identification | bacmap/analysis/ani_identifier | gside |
| sourmash identification | bacmap/analysis/sourmash | gside |
| GTDB-Tk | bacmap/analysis/taxonomic_validator | gside |
| Consensus arbitration | bacmap/analysis/species_consensus | bacmap (reads GOM version) |
| GOM storage | bacmap | bacmap |
| smk orchestration | bacmap | bacmap |

bacmap's snakemake rules (species.smk) now **prefer calling gside CLI**, with
fallback to internal modules when gside is not installed.
