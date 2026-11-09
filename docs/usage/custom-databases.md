# Custom databases

**English** | [简体中文](https://indexofire.github.io/gside/zh/usage/custom-databases/)

`gside` modes draw on databases. Besides the bundled/setup tiers
(`D1_marker`, `D2_ani`, `D3_mash`, `D4_sourmash`), every database-capable
mode — `marker`, `panel`, `mash_refseq`, `sourmash` — accepts additional
named databases, and one identification run can select several of them.

## Selecting databases

`scan` takes repeatable `-d/--db [MODE:]NAME` values. A value resolves as:

1. an existing directory path → used directly;
2. a name under `data/db/<family>/<name>/` (see below);
3. the reserved name `default` → the mode's bundled tier.

```bash
# one named marker database (single mode: the MODE: prefix is optional)
gside scan sample.fna -m marker -d food

# two marker databases in one run: rules and hits are pooled and arbitrated
gside scan sample.fna -m marker -d food -d lab2026

# combined modes need qualified names
gside scan sample.fna -m marker+panel -d marker:food -d panel:lab2026

# a directory path also works (validated per family layout)
gside scan sample.fna -m panel -d /opt/dbs/my-panel
```

`--marker-db PATH` remains as a compatibility alias for `-d marker:PATH`.

### Family directories

Named databases live under the mode's family directory inside the database
root (`data/db/` unless `--db-dir` moves it):

```text
data/db/
├── D1_marker/               # marker default tier (bundled)
├── D2_ani/                  # panel default tier (db setup --tier panel)
├── D3_mash/                 # mash_refseq default tier
├── D4_sourmash/             # sourmash default tier
├── marker/<name>/           # custom marker DBs (layout below)
├── panel/<name>/            # panel.sketch/ + metadata.tsv
├── mash_refseq/<name>/      # mash.msh
└── sourmash/<name>/         # gtdb-reps-k31.zip + lineages.csv
```

`gside db list` enumerates them ("named databases" section, and a
`custom_databases` key in `--format json`). Modes without local databases
(`gtdbtk`) reject `-d` with an actionable error.

### Merge semantics

Selecting several databases of one mode runs **all** of them and merges:

- `marker`: rules from every database compete over pooled hits — equivalent
  to one merged database. Control genes and reads defaults merge first-wins
  (selection order).
- `panel` / `mash_refseq` / `sourmash`: each database is searched; hits are
  pooled and re-ranked, and the verdict is derived once from the global best
  hit. Taxa maps merge first-wins.

The method's `database` field reports the joined provenance
(`{"name": "food+lab2026", ...}`) and a `databases` array lists each
contributing database (name, version, path). A single selection keeps the
exact single-database output shape.

## Marker database layout

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

The default gseek marker path aligns against `markers.fasta` directly. The
legacy BLAST path uses `markers_blastdb`; if `markers_blastdb.nhr` / `.nin` is
missing, gside builds it in the database directory on first use:

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

# named / multi-database selection (see "Selecting databases" above)
gside scan sample.fna -m marker -d food -d lab2026
```

`--marker-db` affects only the marker method when used with a combined mode
such as `-m marker+sourmash`; the other methods keep using their tiers. The
same holds for `--db marker:food`.

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

You can also pass a path directly, or a list to run and merge several
databases:

```python
result = identify_multigene("sample.fna", marker_db="my-markers")
merged = identify_multigene("sample.fna", marker_db=["my-markers", "other"])
```

The result's `database` object reports the custom `name` and `version` (joined
with `+` when several databases ran), so downstream JSON remains traceable.

## Limits

- Named databases are supported for `marker`, `panel`, `mash_refseq`, and
  `sourmash`; `gtdbtk` (external GTDB-Tk database) does not accept `-d`.
- Custom databases do not participate in `db setup` / `db status`; `db list`
  shows named databases, `db validate` checks marker layouts.
- The default gseek path and the shadow comparison use the same custom FASTA;
  the legacy BLAST path uses the custom BLAST index.
- Marker rules cannot create new algorithm behavior; they only select species
  from gene hits.
