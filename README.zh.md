# gside

[![CI](https://github.com/indexofire/gside/actions/workflows/test.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![Docs](https://github.com/indexofire/gside/actions/workflows/docs.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![PyPI](https://img.shields.io/pypi/v/gside.svg)](https://pypi.org/project/gside/)
[![Python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)

[English](README.md) | **简体中文**

`gside` 是一个用于**多层微生物物种鉴定**和**宏基因组谱分析**的 Python 库与命令行工具。它把成熟的生信工具封装在可复现、数据驱动的规则之后，跨方法仲裁，并为每个样本返回一个结构化结果。

库本身按引擎提供独立函数；`gside` CLI 额外提供数据库安装、分层共识、批处理和面向智能体工作流的 MCP 服务器。

```text
contigs / reads
      │
      ├── marker BLAST ──┐
      ├── ANI (skani) ───┤
      ├── mash ──────────┼── 分层仲裁 ── 物种判定
      ├── sourmash ──────┤
      ├── GTDB-Tk ───────┘
      └── kraken2 / sylph ── 群落谱
```

## 特点

- **库优先。** 通过 `gside.api` 直接调用各引擎，结果对象可用于流水线、notebook 和测试。
- **五种分离株模式。** marker 基因、精选 ANI 面板、RefSeq mash、sourmash/GTDB、GTDB-Tk，按 L3 > L2 > L1 仲裁。
- **两种群落模式。** Kraken2 + Bracken、sylph，或两者合并。
- **组装质量评估。** CheckM2 完整度/污染度 + GTDB-Tk 分类。
- **面向自动化。** 稳定 JSON 信封、结构化缺库错误、`--auto-setup`、机器可读数据库状态和 MCP 服务器。
- **可选 gseek 引擎。** `GSIDE_ENGINE=legacy|gseek|shadow` 可把 marker BLAST、read screening、mash、sourmash 和 panel/skani 切到 [gseek](https://github.com/indexofire/gseek)；默认仍是 legacy。

## 安装

`gside` 要求 **Python 3.12**。

```bash
pip install gside
```

Python 包和随包分发的 marker 数据库会一起安装。比对器与参考工具不是 Python 依赖，需要用 conda、Bioconda 或 pixi 安装：

```text
blastn, makeblastdb, minimap2, skani, mash, sourmash,
kraken2, bracken, sylph, GTDB-Tk, CheckM2
```

完整开发环境：

```bash
pixi install
pixi run gside -V
```

可选 gseek 引擎：

```bash
pip install "gside[gseek]"
# 或使用 dev pixi 环境
```

如果工具不在当前环境，可指定其目录：

```bash
export GSIDE_PIXI_BIN=/path/to/conda/env/bin
gside db status
```

## Python 库

稳定导入入口是 `gside.api`。每个函数返回带类型的结果对象或字典列表；基于 `MethodEnvelope` 的结果可用 `to_dict()` 序列化为 JSON。

### Marker 基因鉴定

```python
from gside.api import identify_multigene

result = identify_multigene("contigs.fna", threads=4)
print(result.species, result.confidence)
print(result.to_dict()["detected_markers"])
```

### 自定义 marker 数据库

```python
from gside.api import identify_multigene
from gside.marker_db import MarkerDatabase

db = MarkerDatabase.from_path("my-markers")
errors = db.validate()          # 有效时为 []
result = identify_multigene("contigs.fna", marker_db=db)
```

或使用 CLI：

```bash
gside db validate my-markers
gside scan contigs.fna -m marker --marker-db my-markers
```

自定义 marker 数据库与 `mini` 使用同一算法，只替换 `markers.fasta`、
`marker_rules.yaml` 和可选 `metadata.yaml`。见[自定义 marker 数据库](docs/usage/custom-databases.zh.md)。

### 精选面板 / RefSeq mash 的 ANI

```python
from gside.api import identify_by_ani

panel = identify_by_ani("contigs.fna", mode="panel", threads=8)
mash = identify_by_ani("contigs.fna", mode="mash_refseq", threads=8)

print(panel.to_dict()["result"]["species"])
```

`mode` 可为 `"panel"`（skani ANI）或 `"mash_refseq"`（mash 距离）。用
`db_dir=` 指定非默认数据库根目录。

### Sourmash / GTDB-Tk

```python
from gside.api import identify_by_gtdbtk, identify_by_sourmash

sourmash = identify_by_sourmash("contigs.fna").to_dict()
gtdbtk = identify_by_gtdbtk("contigs.fna", threads=8).to_dict()
```

GTDB-Tk 需要配置 `GTDBTK_DATA_PATH`；sourmash 使用 `sourmash` 层的 GTDB gather 数据库。

### Read screening 与宏基因组谱分析

```python
from pathlib import Path

from gside.api import identify_by_kraken2, identify_by_sylph, screen_reads
from gside.config import SPECIES_DB_DIR

markers = SPECIES_DB_DIR / "D1_marker" / "markers.fasta"
hits = screen_reads(("sample_R1.fq", "sample_R2.fq"), markers, threads=8)

kraken = identify_by_kraken2(("sample_R1.fq", "sample_R2.fq"), threads=8)
sylph = identify_by_sylph(("sample_R1.fq", "sample_R2.fq"), threads=8)
```

### 组装质量

```python
from gside.api import validate_genome

validation = validate_genome("contigs.fna", "qc/", threads=8)
print(validation.status)
```

`validate_genome` 会写入 `qc/validation.json`，同时返回结果对象。

### 结果契约

| 结果类型 | 形状 |
| --- | --- |
| ANI、sourmash、GTDB-Tk、kraken2、sylph | `MethodEnvelope`，有 `.to_dict()` |
| Marker 鉴定 | `MultiGeneResult.to_dict()` |
| Read screening | `list[dict]` 行 |
| 组装质量 | `ValidateResult.to_dict()` |

`MethodEnvelope.to_dict()` 形如：

```json
{
  "analysis_type": "species_identification",
  "method": "panel",
  "database": {"name": "D2_ani", "version": "..."},
  "result": {"species": "...", "confidence": "high"}
}
```

CLI 的 JSON 契约、错误字段和退出行为见
[JSON 输出契约](docs/reference/json-contract.zh.md)。

## CLI

CLI 使用同一批引擎，并增加数据库管理和分层仲裁。

```bash
# 版本；没有单独的 version 子命令
gside -V
gside --version
gside -V -J              # JSON 版本对象

# 单基因组（marker 随包分发，开箱即用）
gside scan contigs.fna
gside scan contigs.fna -m marker+sourmash
gside scan contigs.fna -m all --auto-setup

# 分离株 reads；paired mates 会合并为一个样本
gside scan reads.fq -m marker
gside scan sample_R1.fq sample_R2.fq -m marker

# 宏基因组群落
gside scan --meta reads_R1.fq reads_R2.fq

# 组装质量
gside validate contigs.fna -o qc/

# 批处理 + 流式表格输出
gside scan *.fna -m panel+mash_refseq -j 8 -t 4 -o results.tsv
```

### 模式与仲裁

| 模式 | 引擎 | 输入 | 参考库 |
| --- | --- | --- | --- |
| `marker`（默认） | blastn | FASTA / FASTQ | 随包 marker DB |
| `panel` | skani | FASTA | 精选面板层 |
| `mash_refseq` | mash | FASTA / FASTQ | RefSeq mash 层 |
| `sourmash` | sourmash | FASTA / FASTQ | GTDB gather 层 |
| `gtdbtk` | GTDB-Tk | FASTA | 外部 GTDB-Tk 数据 |
| `kraken2` | Kraken2 + Bracken | FASTA / FASTQ | Kraken2 层 |
| `sylph` | sylph | FASTQ | sylph 层 |
| `all` | 所有适用方法 | 视输入 | 所有必需层 |

`-m` 支持 `+` 组合，例如 `marker+sourmash`。分离株模式下 `all` 运行所有基因组层；`--meta -m all` 合并 Kraken2 与 sylph。

仲裁始终偏向更高置信层：

```text
L3  GTDB-Tk                 系统发育位置
L2  panel / mash / sourmash 比较基因组学（ANI / sketch）
L1  marker                  靶基因存在/缺失
```

### 数据库命令

人类可读列表使用带边框的 Rich 表格；机器格式保持稳定。

```bash
gside db status
gside db list
gside db status -f json
gside db list -f tsv
gside db setup --tier panel
```

| 层 | 内容 | 约大小 |
| --- | --- | --- |
| `mini` | Marker 规则 + 序列（随包） | ~1 MB |
| `panel` | 精选 skani 参考面板 | ~130 MB sketch / ~1 GB genomes |
| `mash` | RefSeq MinHash sketch | ~179 MB |
| `sourmash` | GTDB gather 数据库 | ~3.9 GB |
| `meta` | Kraken2 Standard-8 索引 | ~5.5 GB |
| `meta_sylph` | sylph GTDB-r226 数据库 | ~14 GB |

`gside validate` 还需要 CheckM2（约 1.7 GB）和 GTDB-Tk
（`GTDBTK_DATA_PATH`，下载约 57 GB）。

## 输出与自动化

- `gside scan` 默认每个样本输出一个 JSON 文档；表格可用 `-f tsv`、`-f md`、`-f csv`。
- `gside db status/list` 默认人类可读表格；自动化用 `-f json|tsv|csv`。
- `gside validate` 写 `validation.json`，质量判定失败时以非零退出。
- 结构化缺库错误包含 `suggested_command`、`missing_tier` 和 `estimated_size`。
- `gside mcp`（或 `gside-mcp`）启动 MCP stdio 服务器，暴露物种鉴定、宏基因组和数据库操作。

## 文档

- [文档主页](https://indexofire.github.io/gside/)
- [快速开始](docs/installation/quick-start.zh.md)
- [Python API](docs/usage/python-api.zh.md)
- [CLI 参考](docs/usage/cli.zh.md)
- [鉴定模式](docs/usage/modes.zh.md)
- [宏基因组谱分析](docs/usage/meta.zh.md)
- [JSON 契约](docs/reference/json-contract.zh.md)
- [English docs](docs/index.md)

## 开发

```bash
pixi install
pixi run check          # ruff check/format
pixi run typecheck      # pyright
pixi run pytest         # 单元测试
pixi run docs-build     # 严格 MkDocs 构建
```

测试覆盖 CLI、数据库输出、分层仲裁、结果契约以及可选 gseek shadow 路径。缺少外部工具时，相关集成测试会自动跳过。

## 目录结构

```text
src/gside/
  api.py                 公共 Python facade
  cli.py                 Click CLI
  analysis/              各鉴定引擎
  engine/                工具后端与可选 gseek 桥
  db.py                  数据库层/下载
  mcp.py                 MCP 服务器
docs/                    中英文双语文档
data/db/D1_marker/       随包 marker 数据库
tests/                   单元与 parity 测试
```

## 许可证与生态

MIT。属于 `indexofire` 微生物基因组学工具栈：

| 项目 | 用途 |
| --- | --- |
| [gseek](https://github.com/indexofire/gseek) | 可复现的 DNA 比对基础设施 |
| [gside](https://github.com/indexofire/gside) | 物种鉴定 + 宏基因组谱分析 |
| [gapit](https://github.com/indexofire/gapit) | 基因筛查 + 表型判定 |
| [gmlst](https://github.com/indexofire/gmlst) | MLST / cgMLST 分型 |
