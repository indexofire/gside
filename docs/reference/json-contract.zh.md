# JSON 输出契约

所有 gside 输出均为统一 JSON 格式，与 GOM（Genome Object Model）兼容。

## 单方法输出（--mode marker/panel/mash_refseq）

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
      "all_hits": [...],
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

## 多方法仲裁输出（--mode all）

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.1.0",
  "contigs": "/path/to/contigs.fasta",
  "methods": {
    "marker": { ... },
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
    "mash_refseq": { ... }
  },
  "verdict": {
    "species": "Campylobacter jejuni subsp. jejuni NCTC 11168",
    "confidence": "high",
    "basis": ["panel"]
  }
}
```

## 字段说明

| 字段 | 类型 | 说明 |
|---|---|---|
| `analysis_type` | str | 固定为 `species_identification` |
| `tool` | str | 固定为 `gside` |
| `version` | str | gside 版本号 |
| `methods` | dict | 每个运行的方法的结果（key = mode 名） |
| `methods.<mode>.species` | str | 该方法判定的物种（Unknown 表示无法判定） |
| `methods.<mode>.confidence` | str | high / medium / low |
| `methods.<mode>.error` | str | 方法执行错误（仅失败时出现） |
| `verdict.species` | str | 最终仲裁判定 |
| `verdict.confidence` | str | 最终置信度 |
| `verdict.basis` | list[str] | 判定依据的方法列表 |

## ANI 方法特有字段（panel / mash_refseq）

| 字段 | 类型 | 说明 |
|---|---|---|
| `result.ani` | float | 最高 ANI 值（0-100） |
| `result.aligned_fraction` | float | 比对覆盖比例 |
| `result.top_hits` | list | 前 N 个命中（genome/species/ani/af） |

## GOM 兼容

输出格式可直接映射到 GOM ANALYSIS 对象的 `payload`：

```python
import json, subprocess

result = subprocess.run(
    ["gside", "species", "contigs.fna", "--mode", "all"],
    capture_output=True, text=True,
)
payload = json.loads(result.stdout)

# 直接作为 GOM payload 入库
gos.create(GenomeObject(
    object_type=ObjectType.ANALYSIS,
    payload=payload,
    ...
))
```
