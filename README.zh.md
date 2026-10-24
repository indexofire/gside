# gside — 基因组物种鉴定引擎

[![CI](https://github.com/indexofire/gside/actions/workflows/test.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![Docs](https://github.com/indexofire/gside/actions/workflows/docs.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![PyPI](https://img.shields.io/pypi/v/gside.svg)](https://pypi.org/project/gside/)

独立 CLI：从组装 contigs 做多层物种鉴定，从 FASTQ reads 做宏基因组谱分析。

```bash
pip install gside

# 单基因组物种鉴定
gside scan contigs.fna -m all

# 宏基因组谱分析
gside scan --meta reads_R1.fq reads_R2.fq
```

`scan` 由扩展名推断 FASTA/FASTQ，并在任何工具启动前拒绝不匹配的模式
（例如对 reads 用 `-m panel`）。`--meta` 切换到群落工作流，模式限定为
`kraken2`/`sylph`/`all`。旧的 `gside species` 与 `gside meta` 保留为隐藏别名。

## 物种鉴定（`gside scan`）

| 模式 | 引擎 | 原理 | 数据库 |
|---|---|---|---|
| `marker`（默认） | blastn | 靶基因组合规则（42 条规则，139 条参考序列） | 随仓库分发 |
| `panel` | skani | ANI 比对精选参考面板（291 基因组） | ~130 MB |
| `mash_refseq` | mash | RefSeq MinHash 距离 | ~179 MB |
| `sourmash` | sourmash | GTDB gather | ~3.9 GB |
| `gtdbtk` | GTDB-Tk | 系统发育位置（最权威） | 外部管理 |
| `all` | 以上全部 | 顺序执行 + 分层仲裁 | 按需 |

`-m` 支持 `+` 组合（如 `marker+sourmash`）；`all` 运行全部基因组层。输入方面：
`marker`、`mash_refseq`、`sourmash` 接受 FASTA 或 FASTQ；`panel`/`gtdbtk`
仅接受组装。

可选引擎迁移：`GSIDE_ENGINE=legacy|gseek|shadow` 可把 marker BLAST/FASTQ screening、mash 距离与 sourmash gather 切到 [gseek](https://github.com/indexofire/gseek)，同时保留
legacy 结果。范围与验证见 [gseek 迁移](docs/development/gseek-migration.zh.md)。

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
gside scan contigs.fna -m all

# 并行批量处理 + 流式输出
gside scan *.fna -m all -j 8 -t 4 -o results.tsv
```

## 宏基因组谱分析（`gside scan --meta`）

| 模式 | 引擎 | 输入 | 适用 |
|---|---|---|---|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | 快速筛查、contig 支持、丰度 |
| `sylph` | sylph | 仅 FASTQ | ANI 精度、种级判定 |
| `all` | 双引擎 | 仅 FASTQ（FASTA 请用 `-m kraken2`） | 合并谱 |

```bash
gside scan --meta reads_R1.fq reads_R2.fq
gside scan --meta contigs.fna -m kraken2
```

## 基因组质量（`gside validate`）

`validate` 回答与 `scan` 不同的问题：**这个组装质量如何？** 它运行 CheckM2
（完整度/污染度）+ GTDB-Tk（分类），并把工具版本、参数、逐工具状态与质量结论
写入 `validation.json`。输出目录必填：

```bash
gside validate contigs.fna -o taxonomy
```

CheckM2（约 1.7GB）与 GTDB-Tk（下载约 57GB，细菌运行约 140GB 内存）属于外部
数据库，见[数据库](docs/installation/databases.md)。

## 智能体友好设计

gside 为 AI 智能体和人类用户同等设计：

```bash
# stdout 纯 JSON（`scan` 恒 exit 0；错误在 methods.<mode>.error）
gside scan genome.fna -m all

# 结构化缺库错误
# → methods.panel.suggested_command = "gside db setup --tier panel"
# → methods.panel.missing_tier = "panel"
# → methods.panel.estimated_size = "~130MB sketch / ~1GB genomes"

# 自动安装缺失数据库
gside scan genome.fna -m panel --auto-setup

# 机器可读数据库状态
gside db status -f json
gside db list -f json
```

### MCP 服务器

gside 内置 MCP（Model Context Protocol）stdio 服务器：

```bash
gside mcp    # 或: gside-mcp
```

向 Claude Desktop 或任何 MCP 客户端暴露 `species`、`meta`、`db_status`、`db_list`、`db_setup` 五个工具。

```json
{
  "mcpServers": {
    "gside": { "command": "gside-mcp" }
  }
}
```

## CLI 参考

```bash
gside scan <files...> [选项]
  -m, --mode         marker|panel|mash_refseq|sourmash|gtdbtk|kraken2|sylph|all
                     （支持组合：marker+sourmash）
  --meta             群落工作流（kraken2|sylph|all）
  -D, --db-dir       数据库根目录
  -f, --format       json|tsv|md|csv
  -o, --output       写入文件（行式格式流式写入）
  -j, --jobs         并行样本数（双端 mate 算一个样本）
  -t, --threads      每工具线程数
  -n, --no-header    隐藏表头行
  -a, --auto-setup   自动安装缺失数据库并重试一次
  -s, --samples      --meta 批量清单

gside validate <contigs.fna> -o <目录> [选项]
  -o, --output-dir   必填；输出 validation.json（目录）
  -t, --threads      CheckM2 + GTDB-Tk 线程数

gside db <子命令>
  status             就绪状态（--format table|tsv|csv|json）
  setup              安装 tier（--tier mini|panel|mash|sourmash|meta|meta_sylph|all）
  list               层级 + 状态（--format table|tsv|csv|json）

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
| meta | ~5.5 GB | kraken2 | 自动下载 |
| meta_sylph | ~14 GB | sylph | 自动下载 |

`gside validate` 还需要 CheckM2（约 1.7GB）与 GTDB-Tk
（`GTDBTK_DATA_PATH`，下载约 57GB），详见文档。

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
