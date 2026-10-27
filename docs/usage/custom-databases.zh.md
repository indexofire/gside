# 自定义 marker 数据库

[English](https://indexofire.github.io/gside/usage/custom-databases/) | **简体中文**

`gside` 自带一个 marker 数据库（`D1_marker`），同时允许用户用自己的序列库复用
同一套 marker 算法。算法本身不变：BLAST/minimap2 命中被规范化，应用基因阈值，
再由 `marker_rules.yaml` 决定物种判定；变化的只是 marker 序列和规则。

## 内置 tier 与自定义 marker 数据库

| 概念 | `mini` tier / `D1_marker` | 自定义 marker 数据库 |
| --- | --- | --- |
| 是什么 | 随包数据库 | 用户自己的目录 |
| 安装 | 随 gside 分发 | 无需下载，直接指向 |
| 管理 | `gside db setup` | 用户/流水线自己管理 |
| 使用 | `-m marker` 默认 | `--marker-db PATH` 或 `$GSIDE_MARKER_DB` |
| 适用模式 | marker | 仅 marker |

`db list` / `db status` 描述内置 tier 与组件，不会枚举所有自定义数据库；自定义库
请用 `gside db validate` 显式校验。

## 目录结构

```text
my-markers/
├── markers.fasta          # 必需
├── marker_rules.yaml      # 必需
├── metadata.yaml          # 可选
└── markers_blastdb.*       # 可选；首次使用时构建
```

### `markers.fasta`

FASTA header 必须保留三段式契约：

```text
>markers~~~<gene>~~~<accession> <description>
```

示例：

```text
>markers~~~tlh~~~MY001 thermolabile hemolysin
ATGCGT...
>markers~~~tdh~~~MY002 thermostable direct hemolysin
ATGCGT...
```

要求：

- ID 唯一；
- gene token 非空，匹配规则时大小写不敏感；
- `db validate` 会拒绝 malformed header 和引用未知 gene 的规则。

### `marker_rules.yaml`

与内置库 schema 相同：

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

命中 gene 更多的规则优先；`priority_over` 可让某物种压制非 priority 候选。
完整语义见 [Marker 规则](https://github.com/indexofire/gside/blob/main/docs/reference/marker-rules.zh.md)。

### `metadata.yaml`（可选）

```yaml
schema: gside.marker-db/1
name: my_vibrio_panel
description: In-house Vibrio marker panel
version: 2026.1
# 可选，相对或绝对 BLAST prefix：
# blast_db: markers_blastdb
```

没有 `metadata.yaml` 时，数据库名取目录名，版本取 FASTA + rules 的 SHA-256 前
8 位十六进制。

## 校验

```bash
gside db validate my-markers
gside db validate my-markers -f json
```

`db validate` 检查：

- 必需文件；
- FASTA header 形状与重复 ID；
- YAML 语法；
- `rules` 结构；
- 未知 gene 引用；
- FASTA/rules 为空。

有效库退出码 `0`，有错误退出 `1`。

## BLAST 数据库

默认的 gseek marker 路径直接对 `markers.fasta` 做比对。legacy BLAST 路径使用
`markers_blastdb`；若 `markers_blastdb.nhr` / `.nin` 不存在，gside 会在首次使用时
于数据库目录构建：

```bash
makeblastdb -in my-markers/markers.fasta -dbtype nucl -out my-markers/markers_blastdb
```

只读数据库建议预先构建：

```bash
makeblastdb \
  -in my-markers/markers.fasta \
  -dbtype nucl \
  -out my-markers/markers_blastdb
```

如果目录只读且数据库缺失，gside 会给出可操作错误，而不会写到目录之外。

## CLI 用法

```bash
# 显式路径
gside scan sample.fna -m marker --marker-db ./my-markers

# 环境变量
GSIDE_MARKER_DB=./my-markers gside scan sample.fna -m marker

# 数据库根目录，其中包含 D1_marker
gside scan sample.fna -m marker --db-dir ./my-db-root
```

在 `-m marker+sourmash` 这类组合模式中，`--marker-db` 只影响 marker 方法，其他
方法仍使用各自 tier。

`--auto-setup` 不会安装自定义数据库；自定义库由用户管理，支持用 `db validate`
检查。

## Python 用法

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

也可以直接传路径：

```python
result = identify_multigene("sample.fna", marker_db="my-markers")
```

结果的 `database` 对象会报告自定义 `name` 与 `version`，下游 JSON 仍可追溯。

## 限制

- 自定义数据库仅支持 `marker` 模式。
- 不参与 `db list` / `db status` / `db setup`。
- 默认 gseek 路径与 shadow 对比使用同一个自定义 FASTA；legacy BLAST 路径使用自定义索引。
- marker rules 不能改变算法行为，只负责从 gene 命中里选择物种。
