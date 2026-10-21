# gside — 基因组物种鉴定引擎

[![CI](https://github.com/indexofire/gside/actions/workflows/test.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![Docs](https://github.com/indexofire/gside/actions/workflows/docs.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![PyPI](https://img.shields.io/pypi/v/gside.svg)](https://pypi.org/project/gside/)

独立 CLI：从组装 contigs 做多层物种鉴定，从 FASTQ reads 做宏基因组谱分析。

```bash
pip install gside

# 单基因组物种鉴定
gside scan contigs.fna --mode all

# 宏基因组谱分析
gside scan --meta reads_R1.fq reads_R2.fq --mode all
```

`scan` 由扩展名推断 FASTA/FASTQ，并在任何工具启动前拒绝不匹配的模式
（例如对 reads 用 `-m panel`）。`--meta` 切换到群落工作流，模式限定为
`kraken2`/`sylph`/`all`。旧的 `gside species` 与 `gside meta` 保留为隐藏别名。

## 物种鉴定（`gside scan`）

| 模式 | 引擎 | 原理 | 数据库 |
|---|---|---|---|
| `marker`（默认） | blastn | 靶基因组合规则（38 条规则，80 条参考序列） | 随仓库分发 |
| `panel` | skani | ANI 比对精选参考面板（291 基因组） | ~130 MB |
| `mash_refseq` | mash | RefSeq MinHash 距离 | ~179 MB |
| `sourmash` | sourmash | GTDB gather | ~3.9 GB |
| `gtdbtk` | GTDB-Tk | 系统发育位置（最权威） | 外部管理 |
| `all` | 以上全部 | 顺序执行 + 分层仲裁 | 按需 |

### 分层仲裁

| 层 | 方法 | 依据 |
|---|---|---|
| L3 | gtdbtk | 系统发育位置——最权威 |
| L2 | panel, mash_refseq, sourmash | 比较基因组学（ANI / sketch） |
| L1 | marker | 靶基因存在/缺失 |

高层始终优先于低层；同层内高置信度优先。

### 快速开始

```bash
pixi install && pip install -e .

# 物种鉴定（marker 模式开箱即用）
gside scan contigs.fna

# 安装参考数据库
gside db setup --tier panel

# 多方法鉴定
gside scan contigs.fna --mode all

# 并行批量处理 + 流式输出
gside scan *.fna --mode all -j 8 -t 4 -o results.tsv
```

## 宏基因组谱分析（`gside scan --meta`）

| 模式 | 引擎 | 输入 | 适用 |
|---|---|---|---|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | 快速筛查、contig 支持、丰度 |
| `sylph` | sylph | 仅 FASTQ | ANI 精度、种级判定 |
| `all` | 双引擎 | FASTQ（双）/ FASTA（kraken2） | 合并谱 |

```bash
gside scan --meta reads_R1.fq reads_R2.fq --mode all
gside scan --meta contigs.fna --mode kraken2
```

## 智能体友好设计

gside 为 AI 智能体和人类用户同等设计：

```bash
# stdout 纯 JSON（恒 exit 0；错误在 methods.<mode>.error）
gside scan genome.fna --mode all

# 结构化缺库错误
# → methods.panel.suggested_command = "gside db setup --tier panel"
# → methods.panel.missing_tier = "panel"
# → methods.panel.estimated_size = "~130MB sketch / ~1GB genomes"

# 自动安装缺失数据库
gside scan genome.fna --mode panel --auto-setup

# 机器可读数据库状态
gside db status -f json
gside db list -f json
```

### MCP 服务器

gside 内置 MCP（Model Context Protocol）stdio 服务器：

```bash
gside mcp    # 或: gside-mcp
```

向 Claude Desktop 或任何 MCP 客户端暴露 `species`、`db_status`、`db_list`、`db_setup` 四个工具。

```json
{
  "mcpServers": {
    "gside": { "command": "gside-mcp" }
  }
}
```

## CLI 参考

```
gside scan <contigs.fna> [options]
  -m, --mode         marker|panel|mash_refseq|sourmash|gtdbtk|all
  -D, --db-dir       数据库根目录
  -f, --format       json|tsv|md
  -o, --output       写入文件（tsv/md 逐输入流式）
  -j, --jobs         并行输入文件数
  -t, --threads      每次比对线程数
  -a, --auto-setup   自动安装缺失数据库
  -n, --no-header    隐藏 TSV 表头

gside scan --meta <reads/files> [options]
  -m, --mode         kraken2|sylph|all
  -D, --db-dir       数据库根目录
  -f, --format       json|tsv
  -o, --output       写入文件
  -t, --threads      每工具线程数

gside db <subcommand>
  status             数据库就绪状态（--format table|tsv|csv|json）
  setup              安装数据库（--tier mini|panel|mash|sourmash|all）
  list               层级列表+状态（--format table|tsv|csv|json）

gside validate <contigs.fna> -o <目录>
  -o, --output-dir   必填；输出 validation.json
  -t, --threads      CheckM2 + GTDB-Tk 线程数

gside mcp            MCP stdio 服务器
gside --version      版本信息
gside -V -J          版本 JSON
```

## 数据库

| 层级 | 大小 | 引擎 | 状态 |
|---|---|---|---|
| mini | ~1 MB | blastn (marker) | ✅ 随包分发 |
| panel | ~130 MB | skani | 可下载 |
| mash | ~179 MB | mash | 可下载 |
| sourmash | ~3.9 GB | sourmash | 可下载 |
| meta | ~5.5 GB | kraken2 | 手动安装 |
| meta_sylph | ~14 GB | sylph | 手动安装 |

GTDB-Tk 自管外部数据库（`GTDBTK_DATA_PATH`）。

## 安装

```bash
# pixi（推荐——包含 blastn、skani、mash、sourmash）
curl -fsSL https://pixi.sh/install.sh | bash
pixi install
pip install -e .

# PyPI（Python 3.12，不含二进制）
pip install gside
```

## 文档

- [快速开始](https://indexofire.github.io/gside/)
- [CLI 参考](https://indexofire.github.io/gside/usage/cli/)
- [鉴定模式](https://indexofire.github.io/gside/usage/modes/)
- [宏基因组](https://indexofire.github.io/gside/usage/meta/)
- [JSON 输出契约](https://indexofire.github.io/gside/reference/json-contract/)
- [English](README.md)

## License

MIT

## 生态

| 工具 | 用途 |
|---|---|
| [gside](https://github.com/indexofire/gside) | 物种鉴定 + 宏基因组谱分析 |
| [gapit](https://github.com/indexofire/gapit) | 基因筛查 + 表型指定（血清群、产毒） |
| [gmlst](https://github.com/indexofire/gmlst) | MLST / cgMLST 分型 |
