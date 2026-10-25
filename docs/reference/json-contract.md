# JSON Output Contract

All gside output is unified JSON printed to stdout.

## Single-Method Output (-m marker/panel/mash_refseq/sourmash)

Example with `-m marker` (marker results are flat; no `result` wrapper):

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.2.0",
  "contigs": "/path/to/contigs.fasta",
  "methods": {
    "marker": {
      "species": "Campylobacter_jejuni",
      "confidence": "high",
      "method": "marker",
      "database": {
        "name": "species_markers_v2",
        "version": "a1b2c3d4"
      },
      "detected_markers": [
        {
          "gene": "hipo",
          "identity": 100.0,
          "coverage": 100.0
        }
      ],
      "matched_rule": "hipo+cadf+mapa (3/2 required)",
      "notes": [],
      "interpretation": "Identified as Campylobacter jejuni (confidence: high)."
    }
  },
  "verdict": {
    "species": "Campylobacter_jejuni",
    "display_name": "Campylobacter jejuni",
    "strain": null,
    "confidence": "high",
    "basis": ["marker"],
    "supporting": [],
    "agreement": "1/1"
  }
}
```

The ANI modes (`panel`, `mash_refseq`) and `sourmash` instead nest their
findings under a `result` object, as the `all` example below shows.

## Multi-Method Arbitration Output (-m all)

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.2.0",
  "contigs": "/path/to/contigs.fasta",
  "methods": {
    "marker": {},
    "panel": {
      "method": "panel",
      "database": {
        "name": "D2_ani",
        "version": "e5f6g7h8"
      },
      "result": {
        "species": "Campylobacter_jejuni",
        "species_raw": "Campylobacter jejuni subsp. jejuni NCTC 11168",
        "strain": "subsp. jejuni NCTC 11168",
        "confidence": "high",
        "ani": 100.0,
        "aligned_fraction": 95.2,
        "top_hits": [
          {
            "genome": "GCF_000009085.1_ASM908v1_genomic.fna",
            "species": "Campylobacter jejuni",
            "ani": 100.0,
            "af": 95.2
          }
        ]
      }
    },
    "mash_refseq": {}
  },
  "verdict": {
    "species": "Campylobacter_jejuni",
    "display_name": "Campylobacter jejuni",
    "strain": "subsp. jejuni NCTC 11168",
    "confidence": "high",
    "basis": ["panel"],
    "supporting": ["mash_refseq"],
    "agreement": "2/2"
  }
}
```

## Multiple Inputs (Samples)

`gside scan` accepts several inputs in one run. FASTQ mates (`_R1`/`_R2`,
`_1`/`_2`, `.1`/`.2`) are grouped into **one sample**; unmatched files are
separate samples. The top-level JSON shape depends on the sample count:

- **1 sample**: a single JSON object, as in the examples above
- **2+ samples**: a JSON **array** of those objects, one per sample, in input
  order

For a paired sample, `contigs` is the first mate (both were analyzed).

The `--format tsv` and `--format md` options are alternatives to the default
`json`: they render one row per input file (TSV or Markdown table) instead of
nested objects.

## Field Reference

| Field | Type | Description |
|-------|------|-------------|
| `analysis_type` | str | Always `species_identification` |
| `tool` | str | Always `gside` |
| `version` | str | gside version |
| `contigs` | str | Input path; for paired FASTQ, the first mate |
| `methods` | dict | Results per method (key = mode name) |
| `methods.<mode>.species` | str | Canonical `Genus_species` call for this method |
| `methods.<mode>.species_raw` | str | Original tool label, present only when it differs from `species` |
| `methods.<mode>.strain` | str | Strain / accession detail, when the method resolved one |
| `methods.<mode>.confidence` | str | high / medium / low |
| `methods.<mode>.error` | str | Method execution error (only on failure) |
| `methods.<mode>.suggested_command` | str | Fix-up command when the error is a missing database |
| `methods.<mode>.missing_tier` | str | Tier to provision (`panel`, `mash`, `sourmash`, `meta`, `meta_sylph`) |
| `methods.<mode>.estimated_size` | str | Approximate download size for that tier |
| `methods.<mode>.help_url` | str | Upstream landing page for downloadable tiers (`meta`, `meta_sylph`) |
| `verdict.species` | str | Final arbitrated call (canonical `Genus_species`) |
| `verdict.display_name` | str | Same call with spaces (`Escherichia coli`) |
| `verdict.strain` | str \| null | Strain of the winning method, if any |
| `verdict.confidence` | str | Final confidence |
| `verdict.basis` | list[str] | The single method that produced the verdict |
| `verdict.supporting` | list[str] | Other successful methods that agree on `verdict.species` |
| `verdict.agreement` | str | `agreeing/successful`, e.g. `2/3` (exact label match) |

`agreement` compares canonical labels exactly. Marker rules now emit
species-level labels (genus-level `Shigella` only as the fallback when no
species marker is detected), so marker and the ANI/sketch tiers agree on the
same string. Pathotype/serotype typing (DEC, EIEC, Shigella serotypes) is
gapit's job, not gside's.

## Species name contract

`species` is always `Genus_species` with underscores, matching the labels in
`data/db/D1_marker/marker_rules.yaml` and the downstream bacmap registry:

| Tier | Raw tool label | `species` |
|---|---|---|
| marker | `Streptococcus_agalactiae` | `Streptococcus_agalactiae` (unchanged) |
| panel | `Streptococcus agalactiae 2603V/R` | `Streptococcus_agalactiae` |
| mash_refseq | `Streptococcus_agalactiae_GCF_949788565.1` | `Streptococcus_agalactiae` |
| sourmash | `...;s__Vibrio parahaemolyticus` | `Vibrio_parahaemolyticus` |
| gtdbtk | `Vibrio cholerae` | `Vibrio_cholerae` |

Marker registry labels are already in contract form: `Salmonella_enterica`,
`Vibrio_parahaemolyticus`, `Escherichia_coli`, `Shigella`, `Shigella_flexneri`,
`Shigella_sonnei`, `Shigella_dysenteriae`, `Shigella_boydii`; genus-level
`Shigella` is used only when no species marker is detected. Strain/accession text is never glued
into `species` (a dot would break the downstream matcher); it goes to `strain`.
Species-group placeholders (`Campylobacter sp.`) and `Unknown` are emitted as
`Unknown`.

## ANI-Specific Fields (panel / mash_refseq)

| Field | Type | Description |
|-------|------|-------------|
| `result.ani` | float | Highest ANI value (0-100) |
| `result.aligned_fraction` | float | Alignment coverage fraction |
| `result.top_hits` | list | Top N hits (genome/species/ani/af) |

## sourmash-Specific Fields

| Field | Type | Description |
|-------|------|-------------|
| `result.ani` | float \| null | Containment ANI of the best match (0-100) |
| `result.f_unique_weighted` | float | Fraction of query hashes in the best match |
| `result.gather_partition` | list | Per-match `species` / `f_unique_weighted` / `ani` |
| `result.flags` | list[str] | `possible_mixture` when >=2 distinct species each reach 10% with ANI >= 93 |

Confidence follows the containment ANI (95/93 bands, same as panel). If the
gather CSV has no ANI column, the hash-fraction fallback (0.90/0.70) applies.

## Metagenomic Output (`gside scan --meta`)

```json
{
  "analysis_type": "metagenomic_profiling",
  "tool": "gside",
  "version": "0.2.0",
  "input": {"files": ["reads_R1.fq", "reads_R2.fq"], "type": "fastq_paired"},
  "methods": {
    "kraken2": {"method": "kraken2_bracken", "database": {"name": "D5_kraken2", "path": "..."},
      "result": {"top_species": [{"name": "Escherichia coli", "taxid": 562, "reads": 850,
        "fraction": 0.85, "confidence": "high"}], "total_reads": 10000,
        "classified_reads": 9870, "n_species": 4}},
    "sylph": {"method": "sylph_ani", "database": {"name": "D6_sylph", "path": "..."},
      "result": {"species": [{"name": "Escherichia coli str. K-12", "ani": 99.9,
        "abundance": 0.86, "coverage": 120.0, "confidence": "high"}], "n_species": 1}}
  },
  "verdict": {
    "top_species": [{"name": "Escherichia_coli", "display_name": "Escherichia coli",
      "strain": "str. K-12", "method": "sylph", "abundance": 0.86,
      "confidence": "high", "ani": 99.9, "species_raw": "Escherichia coli str. K-12"}],
    "richness": 1,
    "basis": ["kraken2", "sylph"]
  }
}
```

| Field | Type | Description |
|---|---|---|
| `input.files` / `input.type` | list / str | Inputs and detected type (`fastq_paired`, `fastq_single`, `fasta`) |
| `methods.<mode>.result` | dict | Per-method output (`top_species` for kraken2, `species` for sylph) |
| `verdict.top_species[].name` | str | Canonical `Genus_species` (same contract as `scan`) |
| `verdict.top_species[].display_name` | str | Spaced form |
| `verdict.top_species[].strain` | str \| null | Strain detail from the winning method |
| `verdict.top_species[].method` | str | `kraken2` or `sylph` (higher confidence wins; sylph wins ties) |
| `verdict.top_species[].abundance` | float | Fraction 0-1 (sylph percentages are normalised) |
| `verdict.top_species[].confidence` | str | high / medium / low |
| `verdict.top_species[].ani` | float \| null | sylph ANI when available |
| `verdict.top_species[].species_raw` | str | Original tool label, only when it differs |
| `verdict.richness` | int | Number of merged species |
| `verdict.basis` | list[str] | Methods that contributed |

Batch (`--samples manifest.tsv`) emits a **JSON array**, each object with a
`sample` key; `tsv`/`csv`/`md` rows gain a leading `sample` column.

## Genome Quality Output (`gside validate`)

```json
{
  "schema": "gside.validate/1",
  "genome": {"name": "GCF_000006945.2_ASM694v2_genomic", "path": "/abs/contigs.fna"},
  "parameters": {"threads": 4, "checkm2_db": "/db/checkm2", "gtdbtk_db": "/db/gtdb"},
  "software": {"gside": "0.2.0", "checkm2": "1.0.2", "gtdbtk": "2.4.0"},
  "tools": {"checkm2": "ok", "gtdbtk": "ok"},
  "quality": {"completeness": 99.1, "contamination": 0.5, "verdict": "PASS"},
  "taxonomy": {"classification": "d__Bacteria;p__...;s__Salmonella enterica", "note": ""},
  "status": "ok",
  "errors": [],
  "interpretation": "Genome quality: PASS (...); GTDB-Tk: ...",
  "created_at": "2026-10-21T14:30:00+00:00"
}
```

| Field | Type | Description |
|---|---|---|
| `schema` | str | `gside.validate/1` |
| `genome` | dict | `name` (file stem) and absolute `path` |
| `parameters` | dict | `threads`, `checkm2_db`, `gtdbtk_db` used for this run |
| `software` | dict | `gside`, `checkm2`, `gtdbtk` versions (`""` when unavailable) |
| `tools` | dict | Per-tool status: `ok` / `missing` / `failed` |
| `quality` | dict | CheckM2 `completeness`, `contamination`, `verdict` (`PASS`/`CHECK`/`WARNING`/`UNKNOWN`) |
| `taxonomy` | dict | GTDB-Tk `classification` (full lineage) and `note` |
| `status` | str | `ok` (both tools) / `partial` (one) / `error` (none) |
| `errors` | list[str] | Missing/failed tool messages |
| `created_at` | str | UTC ISO-8601 timestamp |

The same document is written to `<output-dir>/validation.json`. Exit code is
0 only when `status == "ok"`; a low-quality verdict is still exit 0.

## Consuming the output

The JSON is self-describing, so consuming it from Python needs nothing
beyond the standard library:

```python
import json
import subprocess

result = subprocess.run(
    ["gside", "scan", "contigs.fna", "-m", "all"],
    capture_output=True,
    text=True,
)
payload = json.loads(result.stdout)

verdict = payload["verdict"]
print(verdict["species"], verdict["confidence"])

for mode, entry in payload["methods"].items():
    if "error" in entry:
        print(f"{mode}: failed — {entry['error']}")
```

A method that failed carries an `error` key instead of a species call, plus the
structured fix-up fields above when the cause is a missing database, and
in `all` mode the remaining methods still contribute to the verdict.

## Exit Codes

`gside scan` always exits 0, even when one or more methods fail. Judge
success by the JSON content (`methods.<mode>.error` and `verdict`), not by
the exit code.

| Command | Exit code | Behavior |
|---------|-----------|----------|
| `gside scan` | always 0 | Per-method failures go to `methods.<mode>.error`, never the exit code |
| `gside db setup` | 0 / 1 | 1 iff any setup result message contains `ERROR`, else 0 |
| `gside db` (unknown subcommand) | 1 | Prints `unknown subcommand: ...` |
| `gside validate` | 0 / 1 | 0 only when CheckM2 and GTDB-Tk both ran; `partial`/`error` exit 1 (quality verdict does not affect the code) |
