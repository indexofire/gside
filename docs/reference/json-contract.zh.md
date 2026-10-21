# JSON 输出契约

所有 gside 输出为统一 JSON（stdout）。

`marker` 方法的结果为扁平结构（无 `result` 包装）；`panel`、`mash_refseq`、
`sourmash` 的结果嵌套在 `result` 字段内。

## 单方法输出（-m marker）

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.1.5",
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

## 多方法仲裁输出（-m all）

```json
{
  "analysis_type": "species_identification",
  "tool": "gside",
  "version": "0.1.5",
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

## 多输入（样本）

`gside scan` 一次可接受多个 contigs 文件，顶层 JSON 结构随文件数量变化：

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
| `methods.<mode>.help_url` | str | 可下载档位（`meta` / `meta_sylph`）的上游页面 |
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

## sourmash 特有字段

| 字段 | 类型 | 说明 |
|---|---|---|
| `result.ani` | float \| null | 最佳匹配的 containment ANI（0-100） |
| `result.f_unique_weighted` | float | 最佳匹配占查询哈希的比例 |
| `result.gather_partition` | list | 逐匹配的 `species` / `f_unique_weighted` / `ani` |
| `result.flags` | list[str] | 两个及以上不同物种各自 ≥10% 且 ANI ≥93 时为 `possible_mixture` |

置信度按 containment ANI 的 95/93 分界（与 panel 一致）；若 gather CSV 无 ANI
列，则回退到哈希比例 0.90/0.70。

## 宏基因组输出（`gside scan --meta`）

```json
{
  "analysis_type": "metagenomic_profiling",
  "tool": "gside",
  "version": "0.1.5",
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

| 字段 | 类型 | 说明 |
|---|---|---|
| `input.files` / `input.type` | list / str | 输入与识别类型（`fastq_paired`、`fastq_single`、`fasta`） |
| `methods.<mode>.result` | dict | 各方法结果（kraken2 为 `top_species`，sylph 为 `species`） |
| `verdict.top_species[].name` | str | 规范 `Genus_species`（与 `scan` 同一契约） |
| `verdict.top_species[].display_name` | str | 空格形式 |
| `verdict.top_species[].strain` | str \| null | 胜出方法给出的菌株信息 |
| `verdict.top_species[].method` | str | `kraken2` 或 `sylph`（置信度高者优先，同置信 sylph 胜） |
| `verdict.top_species[].abundance` | float | 分数 0-1（sylph 的百分比已归一化） |
| `verdict.top_species[].confidence` | str | high / medium / low |
| `verdict.top_species[].ani` | float \| null | sylph 的 ANI（若有） |
| `verdict.top_species[].species_raw` | str | 原始工具标签，仅在与规范名不同时出现 |
| `verdict.richness` | int | 合并后的物种数 |
| `verdict.basis` | list[str] | 参与合并的方法 |

批量（`--samples manifest.tsv`）输出 **JSON 数组**，每个对象带 `sample` 键；
`tsv`/`csv`/`md` 每行前置一列 `sample`。

## 基因组质量输出（`gside validate`）

```json
{
  "schema": "gside.validate/1",
  "genome": {"name": "GCF_000006945.2_ASM694v2_genomic", "path": "/abs/contigs.fna"},
  "parameters": {"threads": 4, "checkm2_db": "/db/checkm2", "gtdbtk_db": "/db/gtdb"},
  "software": {"gside": "0.1.5", "checkm2": "1.0.2", "gtdbtk": "2.4.0"},
  "tools": {"checkm2": "ok", "gtdbtk": "ok"},
  "quality": {"completeness": 99.1, "contamination": 0.5, "verdict": "PASS"},
  "taxonomy": {"classification": "d__Bacteria;p__...;s__Salmonella enterica", "note": ""},
  "status": "ok",
  "errors": [],
  "interpretation": "Genome quality: PASS (...); GTDB-Tk: ...",
  "created_at": "2026-10-21T14:30:00+00:00"
}
```

| 字段 | 类型 | 说明 |
|---|---|---|
| `schema` | str | `gside.validate/1` |
| `genome` | dict | `name`（文件名主干）与绝对 `path` |
| `parameters` | dict | 本次运行的 `threads`、`checkm2_db`、`gtdbtk_db` |
| `software` | dict | `gside`、`checkm2`、`gtdbtk` 版本（不可用时为 `""`） |
| `tools` | dict | 逐工具状态：`ok` / `missing` / `failed` |
| `quality` | dict | CheckM2 的 `completeness`、`contamination`、`verdict`（`PASS`/`CHECK`/`WARNING`/`UNKNOWN`） |
| `taxonomy` | dict | GTDB-Tk 的 `classification`（完整 lineage）与 `note` |
| `status` | str | `ok`（两工具）/ `partial`（其一）/ `error`（都无） |
| `errors` | list[str] | 缺失/失败的工具信息 |
| `created_at` | str | UTC ISO-8601 时间戳 |

同一份文档会写入 `<output-dir>/validation.json`。仅当 `status == "ok"` 时退出
0；质量结论为低质量仍退出 0。

## 消费输出

用 subprocess 调用 gside，再用标准库 json 解析 stdout：

```python
import json, subprocess

result = subprocess.run(
    ["gside", "scan", "contigs.fna", "-m", "all"],
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

`gside scan` 恒以退出码 0 结束，即使一个或多个方法失败。判断鉴定成败
应读取 JSON 内容（`methods.<mode>.error` 与 `verdict`），而不是依赖退出码。

| 命令 | 退出码 | 行为 |
|---|---|---|
| `gside scan` | 恒为 0 | 单法失败写入 `methods.<mode>.error`，不影响退出码 |
| `gside db setup` | 0 / 1 | 任一安装结果消息含 `ERROR` 时为 1，否则为 0 |
| `gside db`（未知子命令） | 1 | 输出 `unknown subcommand: ...` |
| `gside validate` | 0 / 非零 | 无单法错误保护，失败直接抛出并以非零退出 |
