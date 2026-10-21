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
│  blastn · skani · mash · sourmash · kraken2 · sylph     │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

## Internal Structure

```
gside/
├── src/gside/
│   ├── cli.py              CLI entry point (species / meta / db / validate / mcp)
│   ├── config.py           Runtime config (data dirs / binary discovery)
│   ├── db.py               Database management (status/setup/list)
│   ├── engine/             Alignment backends used by the analysis layer
│   │   ├── _env.py         Binary discovery + checked subprocess execution
│   │   ├── hits.py         Unified Hit dataclass (BLAST tabular)
│   │   └── backends/
│   │       ├── blast.py    blastn (marker layer)
│   │       ├── skani.py    ANI search (panel layer)
│   │       └── kmer.py     mash + sourmash backends
│   └── analysis/
│       ├── multigene_identifier.py   L1 target gene combination
│       ├── ani_identifier.py         L2 ANI (panel/mash_refseq)
│       ├── sourmash_identifier.py    L2 sourmash gather
│       ├── gtdbtk_identifier.py      L3 GTDB-Tk placement
│       ├── kraken2_identifier.py     meta: kraken2 + bracken
│       ├── sylph_identifier.py       meta: sylph ANI profiling
│       ├── taxonomic_validator.py    CheckM2 + GTDB-Tk QC (`gside validate`)
│       ├── _envelope.py              Shared result-dict envelope
│       └── _thresholds.py            Shared ANI confidence cutoffs
├── data/
│   └── db/                 Per-scheme directories (one folder per method)
│       ├── D1_marker/      marker rules + fasta + BLAST db (tracked in git)
│       ├── D2_ani/         skani sketch + genomes (gitignored) + manifests/ (tracked recipe)
│       ├── D3_mash/        mash sketch (gitignored)
│       ├── D4_sourmash/    GTDB reps k=31 + lineages (via db setup)
│       ├── D5_kraken2/     Kraken2 Standard-8 index (manual / gitignored)
│       └── D6_sylph/       sylph GTDB-r226 .syldb (manual / gitignored)
└── tests/
```

## Design Principles

| Principle | Implementation |
|-----------|---------------|
| **Rules as data** | 38 species rules in YAML; new species without code changes |
| **Backends injected per method** | The analysis layer picks blastn / skani / mash / sourmash directly |
| **JSON contract** | All output in one unified JSON format on stdout |
| **Graceful degradation** | Missing database/binary returns clear error, never crashes |
| **Independent operation** | Standalone install and use, no external service required |
