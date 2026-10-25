# Custom marker databases

**English** | [简体中文](https://indexofire.github.io/gside/zh/usage/custom-databases/)

`gside` ships one marker database (`D1_marker`) and exposes the same marker
algorithm for user-provided sequence databases. The algorithm is unchanged:
BLAST/minimap2 hits are normalized, gene-level thresholds are applied, and
`marker_rules.yaml` decides the species call. Only the marker sequences and
rules change.

## Built-in tiers vs custom marker databases

| Concept | `mini` tier / `D1_marker` | Custom marker database |
| --- | --- | --- |
| What it is | bundled database | directory you own |
| Installation | ships with gside | no download; point gside at it |
| Managed by | `gside db setup` | you (or your pipeline) |
| Used by | `-m marker` default | `--marker-db PATH` or `$GSIDE_MARKER_DB` |
| Applicable mode | marker | marker only |

`db list` / `db status` describe the built-in tiers and components. They do not
enumerate every custom database; validate a custom database explicitly with
`gside db validate`.

## Directory layout

```text
my-markers/
├── markers.fasta          # required
├── marker_rules.yaml      # required
├── metadata.yaml          # optional
└── markers_blastdb.*       # optional; built on first use
```

### `markers.fasta`

Headers must keep the existing three-part contract:

```text
>markers~~~<gene>~~~<accession> <description>
```

Example:

```text
>markers~~~tlh~~~MY001 thermolabile hemolysin
ATGCGT...
>markers~~~tdh~~~MY002 thermostable direct hemolysin
ATGCGT...
```

Requirements:

- IDs are unique.
- The gene token is non-empty and is matched case-insensitively against rules.
- `db validate` rejects malformed headers and rules referencing unknown genes.

### `marker_rules.yaml`

Same schema as the bundled database:

```yaml
rules:
  - species: Vibrio_parahaemolyticus
    genes: [tlh, tdh]
    min_hits: 2
    min_identity: 90
    min_coverage: 60
    exclude_genes: [some_confounder]
    note: custom scheme note
```

Rules with more matched genes win; `priority_over` can make a species dominate
non-priority candidates. See [Marker rules](https://github.com/indexofire/gside/blob/main/docs/reference/marker-rules.md)
for the full semantics.

### `metadata.yaml` (optional)

```yaml
schema: gside.marker-db/1
name: my_vibrio_panel
description: In-house Vibrio marker panel
version: 2026.1
# Optional relative/absolute BLAST prefix:
# blast_db: markers_blastdb
```

Without `metadata.yaml`, the database name is the directory name and the
version is the first 8 hex characters of a SHA-256 over FASTA + rules.

## Validation

```bash
gside db validate my-markers
gside db validate my-markers -f json
```

`db validate` checks:

- required files;
- FASTA header shape and duplicate IDs;
- YAML syntax;
- `rules` shape;
- unknown gene references;
- empty FASTA/rules.

Exit code is `0` for a valid DB and `1` when errors are reported.

## BLAST database

The legacy marker path uses BLAST. If `markers_blastdb.nhr` / `.nin` is missing,
gside builds it in the database directory on first use:

```bash
makeblastdb -in my-markers/markers.fasta -dbtype nucl -out my-markers/markers_blastdb
```

You can prebuild it (recommended for read-only databases):

```bash
makeblastdb \
  -in my-markers/markers.fasta \
  -dbtype nucl \
  -out my-markers/markers_blastdb
```

If the directory is read-only and the database is missing, gside raises an
actionable error instead of writing outside it.

## CLI usage

```bash
# explicit path
gside scan sample.fna -m marker --marker-db ./my-markers

# environment variable
GSIDE_MARKER_DB=./my-markers gside scan sample.fna -m marker

# a database root containing D1_marker
gside scan sample.fna -m marker --db-dir ./my-db-root
```

`--marker-db` affects only the marker method when used with a combined mode
such as `-m marker+sourmash`; the other methods keep using their tiers.

`--auto-setup` never installs custom databases. The custom DB is user-managed;
`db validate` is the supported check.

## Python usage

```python
from gside.api import identify_multigene
from gside.marker_db import MarkerDatabase

db = MarkerDatabase.from_path("my-markers")
errors = db.validate()
if errors:
    raise SystemExit("\n".join(errors))

result = identify_multigene("sample.fna", marker_db=db)
print(result.database_name, result.database_version)
print(result.to_dict())
```

You can also pass a path directly:

```python
result = identify_multigene("sample.fna", marker_db="my-markers")
```

The result's `database` object reports the custom `name` and `version`, so
downstream JSON remains traceable.

## Limits

- Custom databases are supported for the `marker` mode only.
- They do not participate in `db list` / `db status` or `db setup`.
- The gseek shadow path uses the same custom FASTA; the legacy BLAST path uses
  the custom BLAST index.
- Marker rules cannot create new algorithm behavior; they only select species
  from gene hits.
