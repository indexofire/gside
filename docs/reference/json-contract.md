# JSON Output Contract

All gside output is unified JSON, compatible with GOM (Genome Object Model).

## Single-Method Output (--mode marker/panel/mash_refseq)

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
        "name": "refseq_panel",
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

## GOM Integration

Output format maps directly to a GOM ANALYSIS object's `payload`:

```python
import json, subprocess

result = subprocess.run(
    ["gside", "species", "contigs.fna", "--mode", "all"],
    capture_output=True, text=True,
)
payload = json.loads(result.stdout)

# Store directly as GOM payload
gos.create(GenomeObject(
    object_type=ObjectType.ANALYSIS,
    payload=payload,
    ...
))
```
