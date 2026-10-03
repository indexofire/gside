# JSON 输出契约

所有 gside 输出为统一 JSON（stdout）。

`marker` 方法的结果为扁平结构（无 `result` 包装）；`panel`、`mash_refseq`、
`sourmash` 的结果嵌套在 `result` 字段内。

## 单方法输出（--mode marker）

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
    "mash_refseq": { ... },
    "sourmash": { ... }
  },
  "verdict": {
    "species": "Campylobacter jejuni subsp. jejuni NCTC 11168",
    "confidence": "high",
    "basis": ["panel"]
  }
}
```

## 多输入文件

`gside species` 一次可接受多个 contigs 文件，顶层 JSON 结构随文件数量变化：

- **1 个文件**：单个 JSON 对象，即上文示例的结构
- **多个文件**：JSON **数组**，每个输入文件对应一个对象，按输入顺序排列

`--format tsv` 与 `--format md` 是默认 `json` 之外的两种输出格式：每个输入
文件渲染一行（TSV 或 Markdown 表格），而非嵌套对象。

## 字段说明

| 字段 | 类型 | 说明 |
|---|---|---|
| `analysis_type` | str | 固定为 `species_identification` |
| `tool` | str | 固定为 `gside` |
| `version` | str | gside 版本号 |
| `contigs` | str | 输入 contigs 文件路径 |
| `methods` | dict | 每个运行的方法的结果（key = mode 名） |
| `methods.<mode>.species` | str | 该方法判定的物种（Unknown 表示无法判定） |
| `methods.<mode>.confidence` | str | high / medium / low |
| `methods.<mode>.result` | dict | ANI/sourmash 方法的嵌套结果（marker 无此层） |
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

## 消费输出

用 subprocess 调用 gside，再用标准库 json 解析 stdout：

```python
import json, subprocess

result = subprocess.run(
    ["gside", "species", "contigs.fna", "--mode", "all"],
    capture_output=True, text=True,
)
payload = json.loads(result.stdout)

print(payload["verdict"]["species"])      # 最终判定
print(payload["verdict"]["confidence"])   # 最终置信度
for name, method in payload["methods"].items():
    if "error" in method:
        print(f"{name} failed: {method['error']}")
```

## 退出码

`gside species` 恒以退出码 0 结束，即使一个或多个方法失败。判断鉴定成败
应读取 JSON 内容（`methods.<mode>.error` 与 `verdict`），而不是依赖退出码。

| 命令 | 退出码 | 行为 |
|---|---|---|
| `gside species` | 恒为 0 | 单法失败写入 `methods.<mode>.error`，不影响退出码 |
| `gside db setup` | 0 / 1 | 任一安装结果消息含 `ERROR` 时为 1，否则为 0 |
| `gside db`（未知子命令） | 1 | 输出 `unknown subcommand: ...` |
| `gside validate` | 0 / 非零 | 无单法错误保护，失败直接抛出并以非零退出 |
