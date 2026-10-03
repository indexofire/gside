# JSON Output Contract

All gside output is unified JSON printed to stdout.

## Single-Method Output (--mode marker/panel/mash_refseq/sourmash)

Example with `--mode marker` (marker results are flat; no `result` wrapper):

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.1.0",
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
      "all_hits": [],
      "notes": [],
      "rule_evidence": "hipo+cadf+mapa (3/2 required)"
    }
  },
  "verdict": {
    "species": "Campylobacter_jejuni",
    "confidence": "high",
    "basis": ["marker"]
  }
}
```

The ANI modes (`panel`, `mash_refseq`) and `sourmash` instead nest their
findings under a `result` object, as the `all` example below shows.

## Multi-Method Arbitration Output (--mode all)

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.1.0",
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
        "species": "Campylobacter jejuni subsp. jejuni NCTC 11168",
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
    "species": "Campylobacter jejuni subsp. jejuni NCTC 11168",
    "confidence": "high",
    "basis": ["panel"]
  }
}
```

## Multiple Input Files

`gside species` accepts several contigs files in one run. The top-level JSON
shape depends on how many:

- **1 file**: a single JSON object, as in the examples above
- **2+ files**: a JSON **array** of those objects, one per input file, in
  input order

The `--format tsv` and `--format md` options are alternatives to the default
`json`: they render one row per input file (TSV or Markdown table) instead of
nested objects.

## Field Reference

| Field | Type | Description |
|-------|------|-------------|
| `analysis_type` | str | Always `species_identification` |
| `tool` | str | Always `gside` |
| `version` | str | gside version |
| `methods` | dict | Results per method (key = mode name) |
| `methods.<mode>.species` | str | Species called by this method |
| `methods.<mode>.confidence` | str | high / medium / low |
| `methods.<mode>.error` | str | Method execution error (only on failure) |
| `verdict.species` | str | Final arbitrated call |
| `verdict.confidence` | str | Final confidence |
| `verdict.basis` | list[str] | Methods that contributed to the verdict |

## ANI-Specific Fields (panel / mash_refseq)

| Field | Type | Description |
|-------|------|-------------|
| `result.ani` | float | Highest ANI value (0-100) |
| `result.aligned_fraction` | float | Alignment coverage fraction |
| `result.top_hits` | list | Top N hits (genome/species/ani/af) |

## Consuming the output

The JSON is self-describing, so consuming it from Python needs nothing
beyond the standard library:

```python
import json
import subprocess

result = subprocess.run(
    ["gside", "species", "contigs.fna", "--mode", "all"],
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

A method that failed carries an `error` key instead of a species call, and
in `all` mode the remaining methods still contribute to the verdict.

## Exit Codes

`gside species` always exits 0, even when one or more methods fail. Judge
success by the JSON content (`methods.<mode>.error` and `verdict`), not by
the exit code.

| Command | Exit code | Behavior |
|---------|-----------|----------|
| `gside species` | always 0 | Per-method failures go to `methods.<mode>.error`, never the exit code |
| `gside db setup` | 0 / 1 | 1 iff any setup result message contains `ERROR`, else 0 |
| `gside db` (unknown subcommand) | 1 | Prints `unknown subcommand: ...` |
| `gside validate` | 0 / non-zero | No per-method error guard; failures raise and exit non-zero |
