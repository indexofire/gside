# Layered Design

## Architecture Position

gside sits between standard bioinformatics tools and your workflows:

```
┌─────────────────────────────────────────────────────────┐
│                                                         │
│  gside (this tool)                                      │
│  Species ID engine: marker rules + ANI backends +       │
│  arbitration, exposed as one CLI with JSON output       │
│                                                         │
├─────────────────────────────────────────────────────────┤
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
│       ├── ani_identifier.py         L2 ANI (panel/mash_refseq)
│       ├── sourmash_identifier.py    L2 sourmash gather
│       ├── taxonomic_validator.py    L3 GTDB-Tk (present in the source
│       │                             tree but not yet wired into the CLI)
│       └── species_canon.py          Species name canonicalization
│                                   (present but not yet wired into the CLI)
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
| **JSON contract** | All output in one unified JSON format on stdout |
| **Graceful degradation** | Missing database/binary returns clear error, never crashes |
| **Independent operation** | Standalone install and use, no external service required |
