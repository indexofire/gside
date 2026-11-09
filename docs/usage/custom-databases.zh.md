# 自定义数据库

[English](https://indexofire.github.io/gside/usage/custom-databases/) | **简体中文**

`gside` 的各模式依赖数据库。除内置/安装的 tier（`D1_marker`、`D2_ani`、
`D3_mash`、`D4_sourmash`）外,所有支持本地数据库的模式 —— `marker`、
`ani`、`mash`、`sourmash` —— 都可以添加命名数据库,并且一次鉴定
可以选择其中的多个。

## 选择数据库

`scan` 接受可重复的 `-d/--db [MODE:]NAME` 参数,取值按以下顺序解析:

1. 已存在的目录路径 → 直接使用;
2. `data/db/<family>/<name>/` 下的命名库(见下);
3. 保留名 `default` → 该模式的内置 tier。

```bash
# 单模式时 MODE: 前缀可省略
gside scan sample.fna -m marker -d food

# 一次跑两个 marker 库:规则与命中合并后统一仲裁
gside scan sample.fna -m marker -d food -d lab2026

# 组合模式必须带前缀
gside scan sample.fna -m marker+panel -d marker:food -d panel:lab2026

# 目录路径也可以(按 family 布局校验)
gside scan sample.fna -m panel -d /opt/dbs/my-panel
```

`--marker-db PATH` 保留为 `-d marker:PATH` 的兼容别名。

### family 目录

命名数据库放在数据库根目录(`data/db/`,`--db-dir` 可整体移动)下各模式的
family 子目录中:

```text
data/db/
├── D1_marker/               # marker tier:一个子目录一个数据库
│   ├── default/             # 随包默认库(markers.fasta + 规则 + blastdb)
│   └── cpin/                # 随包发行的 CPIN 扩展库
├── D2_ani/                  # panel tier(default/ 存 sketch;genomes/ 为构建输入)
│   └── <name>/              # 你的命名 panel 库
├── D3_mash/                 # mash tier(default/mash.msh;<name>/ 放更多)
└── D4_sourmash/             # sourmash tier(default/;<name>/ 放更多)
```

每个模式独占一个 tier 目录,数据库放在其**内部**——一个子目录一个库,
`default` 加任意命名库。目录内名字唯一,因此 (模式, 库名) 天然不重复。
`-d MODE:` 前缀接受模式名(`panel:lab`)或别名(`ani:lab`),二者等价。

**一个模式多个 bundled 数据库的存放约定。** 其余 bundled 库与用户库
同格式,作为 `default/` 的兄弟子目录发行——`D1_marker/cpin` 即随包示例
(git 跟踪,同一套布局校验)。你的库放在旁边(`D1_marker/my-panel/`),
或 `-d` 直接指向任意目录路径。

`gside db list` / `gside db status` 会枚举它们:`-d` 可选的每个数据库都
以 `MODE / DATABASE` 行出现(JSON 为 `modes`,另附可安装的 `tiers`)。
没有本地数据库的模式
(`gtdbtk`)传入 `-d` 会得到可操作的报错。

### 合并语义

同一模式选中多个数据库时会**全部运行**并合并:

- `marker`:各库规则在合并后的命中上共同竞争 —— 等价于一个合并数据库。
  对照基因与 reads 默认值按选择顺序 first-wins 合并。
- `panel` / `mash_refseq` / `sourmash`:逐库检索,命中合并重排,判定由
  全局最优命中一次导出。taxa 映射 first-wins 合并。

方法的 `database` 字段报告合并后的溯源(`{"name": "food+lab2026", ...}`),
`databases` 数组列出每个参与的库(name/version/path)。单选时输出与单库
形态完全一致。

## marker 数据库目录结构

```text
my-markers/
├── markers.fasta          # 必需
├── marker_rules.yaml      # 必需
├── metadata.yaml          # 可选
└── markers_blastdb.*       # 可选；首次使用时构建
```

### `markers.fasta`

FASTA header 采用 gspec 风格的竖线契约：

```text
>gside|db=<database>|gene=<gene>|acc=<accession>|func=<role>|note=<text>
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

# 命名/多库选择（见上文「选择数据库」）
gside scan sample.fna -m marker -d food -d lab2026
```

在 `-m marker+sourmash` 这类组合模式中，`--marker-db` 只影响 marker 方法，其他
方法仍使用各自 tier。`--db marker:food` 同理。

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

也可以直接传路径，或传列表一次运行并合并多个库：

```python
result = identify_multigene("sample.fna", marker_db="my-markers")
merged = identify_multigene("sample.fna", marker_db=["my-markers", "other"])
```

结果的 `database` 对象会报告自定义 `name` 与 `version`（多库时以 `+` 连接），
下游 JSON 仍可追溯。

## 限制

- 命名数据库支持 `marker`、`ani`、`mash`、`sourmash`；`gtdbtk`
  （外部 GTDB-Tk 库）不接受 `-d`。
- 自定义库不参与 `db setup` / `db status`；`db list` 会列出命名库，
  `db validate` 校验 marker 布局。
- 默认 gseek 路径与 shadow 对比使用同一个自定义 FASTA；legacy BLAST 路径使用自定义索引。
- marker rules 不能改变算法行为，只负责从 gene 命中里选择物种。
