# Python API

`gside.api` is the supported library entry point. It exposes the individual
identification engines without the CLI's database provisioning or layered
arbitration.

```python
from gside.api import identify_multigene
```

The API is pre-1.0: new optional arguments may be added, but existing
positional arguments and result shapes are kept stable.

## Requirements

Each function needs the same external tools and reference data as its CLI mode:

| Function | External tools | Data |
| --- | --- | --- |
| `identify_multigene` | `blastn`, `makeblastdb` | bundled `D1_marker` |
| `identify_by_ani(mode="panel")` | `skani` | panel tier |
| `identify_by_ani(mode="mash_refseq")` | `mash` | mash tier |
| `identify_by_sourmash` | `sourmash` | sourmash tier |
| `identify_by_gtdbtk` | `gtdbtk` | `GTDBTK_DATA_PATH` |
| `identify_by_kraken2` | `kraken2`, `bracken` | meta tier |
| `identify_by_sylph` | `sylph` | meta_sylph tier |
| `screen_reads` | `minimap2` | marker FASTA |
| `validate_genome` | `checkm2`, `gtdbtk` | external databases |

Set `GSIDE_DB_DIR` or pass `db_dir=` to use a non-default database root. Tools
are resolved from `PATH` and `GSIDE_PIXI_BIN`; see
[Configuration](https://github.com/indexofire/gside).

## Marker genes

```python
from gside.api import identify_multigene

result = identify_multigene("contigs.fna", threads=4)
# custom marker DB:
# result = identify_multigene("contigs.fna", marker_db="my-markers")

print(result.species)
print(result.confidence)          # high / medium / low
print(result.detected_markers)    # list[dict]
print(result.to_dict())           # JSON-ready
```

`identify_multigene` accepts one path or a list of contig files and applies the
bundled marker-combination rules. It returns a `MultiGeneResult`.

## ANI: panel and RefSeq mash

```python
from gside.api import identify_by_ani

panel = identify_by_ani("contigs.fna", mode="panel", threads=8)
mash = identify_by_ani("contigs.fna", mode="mash_refseq", threads=8)

print(panel.to_dict()["result"]["species"])
print(panel.to_dict()["result"]["top_hits"][:3])
```

`mode` must be `"panel"` or `"mash_refseq"`. The return value is an
`AniIdResult`; its `result` contains `species`, `confidence`, and `top_hits`
with ANI/identity and aligned fraction.

## Sourmash and GTDB-Tk

```python
from gside.api import identify_by_gtdbtk, identify_by_sourmash

sourmash = identify_by_sourmash("contigs.fna").to_dict()
gtdbtk = identify_by_gtdbtk("contigs.fna", threads=8).to_dict()
```

`identify_by_sourmash` performs a GTDB gather and returns `SourmashIdResult`.
`identify_by_gtdbtk` returns `GtdbkResult` and requires a GTDB-Tk database
configured through `GTDBTK_DATA_PATH`.

## Reads and metagenomes

### Marker screening

```python
from pathlib import Path

from gside.api import screen_reads
from gside.config import SPECIES_DB_DIR

markers = SPECIES_DB_DIR / "D1_marker" / "markers.fasta"
rows = screen_reads(("sample_R1.fq", "sample_R2.fq"), markers, threads=8)
```

Paired mates are mapped in one invocation. The result is a list of aggregated
per-gene rows, not a `MethodEnvelope`.

### Kraken2/Bracken and sylph

```python
from gside.api import identify_by_kraken2, identify_by_sylph

kraken = identify_by_kraken2(("sample_R1.fq", "sample_R2.fq"), threads=8)
sylph = identify_by_sylph(("sample_R1.fq", "sample_R2.fq"), threads=8)

print(kraken.to_dict()["result"]["top_species"])
```

Kraken2 accepts FASTA or FASTQ; sylph requires FASTQ. Both return
`MethodEnvelope` subclasses.

## Assembly validation

```python
from gside.api import validate_genome

result = validate_genome("contigs.fna", "qc/", threads=8)
print(result.status)       # ok / warning / fail (tool-dependent)
print(result.to_dict())
```

`validate_genome` runs CheckM2 and GTDB-Tk, writes `qc/validation.json`, and
returns `ValidateResult`.

## Result shapes

Most engines return a `MethodEnvelope` subclass:

```python
{
    "analysis_type": "species_identification",
    "method": "panel",
    "database": {"name": "D2_ani", "version": "..."},
    "result": {...},
}
```

- `identify_multigene` returns `MultiGeneResult`; use `.to_dict()`.
- `screen_reads` returns `list[dict]`.
- `validate_genome` returns `ValidateResult`; use `.to_dict()`.

See [JSON output contract](../reference/json-contract.md) for CLI-level fields
and error semantics.

## Errors and exit behavior

Library functions raise Python exceptions (`ValueError` for invalid arguments,
`RuntimeError` for missing/invalid tools or databases). The CLI converts these
to structured JSON errors and exit codes; it does not hide per-method failures.

## Concurrency

Each function processes one sample. Parallelize at the sample level in your
application (threads/processes) and set `threads=` to the per-tool thread count.
Do not share mutable result objects between workers.

## Engine selection

Set `GSIDE_ENGINE` to select the internal implementation for migrated paths:

```bash
GSIDE_ENGINE=legacy   # default
GSIDE_ENGINE=gseek    # use gseek where implemented
GSIDE_ENGINE=shadow   # run both, compare, return legacy
```

Scope and validation details are in [gseek migration](../development/gseek-migration.md).
