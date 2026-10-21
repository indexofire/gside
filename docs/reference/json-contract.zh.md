# JSON 输出契约

所有 gside 输出为统一 JSON（stdout）。

`marker` 方法的结果为扁平结构（无 `result` 包装）；`panel`、`mash_refseq`、
`sourmash` 的结果嵌套在 `result` 字段内。

## 单方法输出（--mode marker）

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.1.4",
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
    "display_name": "Campylobacter jejuni",
    "strain": null,
    "confidence": "high",
    "basis": ["marker"],
    "supporting": [],
    "agreement": "1/1"
  }
}
```

## 多方法仲裁输出（--mode all）

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.1.4",
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
    "mash_refseq": { ... },
    "sourmash": { ... }
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
| `methods.<mode>.species` | str | 该方法判定的规范名 `Genus_species`（Unknown 表示无法判定） |
| `methods.<mode>.species_raw` | str | 工具原始标签，仅当与 `species` 不同时出现 |
| `methods.<mode>.strain` | str | 该方法解析出的菌株/登录号细节 |
| `methods.<mode>.confidence` | str | high / medium / low |
| `methods.<mode>.result` | dict | ANI/sourmash 方法的嵌套结果（marker 无此层） |
| `methods.<mode>.error` | str | 方法执行错误（仅失败时出现） |
| `methods.<mode>.suggested_command` | str | 缺库时的修复命令 |
| `methods.<mode>.missing_tier` | str | 需要安装的档位（`panel` / `mash` / `sourmash` / `meta` / `meta_sylph`） |
| `methods.<mode>.estimated_size` | str | 该档位的近似下载体积 |
| `methods.<mode>.help_url` | str | 仅能手动安装档位（`meta` / `meta_sylph`）的上游页面 |
| `verdict.species` | str | 最终仲裁判定（规范名 `Genus_species`） |
| `verdict.display_name` | str | 同一判定的可读形式（空格分隔，如 `Escherichia coli`） |
| `verdict.strain` | str \| null | 胜出方法给出的菌株信息 |
| `verdict.confidence` | str | 最终置信度 |
| `verdict.basis` | list[str] | 产生该判定的方法（单个） |
| `verdict.supporting` | list[str] | 其他成功且与 `verdict.species` 一致的方法 |
| `verdict.agreement` | str | `一致数/成功数`，如 `2/3`（按标签精确匹配） |

`agreement` 按规范标签精确比对。marker 的致病型标签（`DEC`、`Shigella_EIEC`）
有意与所属二名法区分：同一个 E. coli 基因组可能 marker 报 `DEC`、panel 报
`Escherichia_coli`，此时 `agreement` 为 `2/3`，但两者指向同一物种。**不要**把
DEC 折叠成 `Escherichia_coli`——这个区分正是下游致病型分型的路由依据。

## 物种名契约

`species` 一律为下划线分隔的 `Genus_species`，与
`data/db/D1_marker/marker_rules.yaml` 及下游 bacmap 注册表一致：

| 档位 | 工具原始标签 | `species` |
|---|---|---|
| marker | `Streptococcus_agalactiae` | `Streptococcus_agalactiae`（原样） |
| panel | `Streptococcus agalactiae 2603V/R` | `Streptococcus_agalactiae` |
| mash_refseq | `Streptococcus_agalactiae_GCF_949788565.1` | `Streptococcus_agalactiae` |
| sourmash | `...;s__Vibrio parahaemolyticus` | `Vibrio_parahaemolyticus` |
| gtdbtk | `Vibrio cholerae` | `Vibrio_cholerae` |

带血清群/致病型后缀的 marker 注册标签（`Vibrio_cholerae_O1`、
`Shigella_EIEC`、`DEC`、`Salmonella`）本身即契约形式，原样透传。菌株/登录号
文本绝不拼进 `species`（点号会破坏下游匹配），统一放入 `strain`。
种组占位（`Campylobacter sp.`）与 `Unknown` 一律输出 `Unknown`。

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
