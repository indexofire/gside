# gside：基因组物种鉴定引擎

[![CI](https://github.com/indexofire/gside/actions/workflows/test.yml/badge.svg)](https://github.com/indexofire/gside/actions)
[![PyPI](https://img.shields.io/pypi/v/gside.svg)](https://pypi.org/project/gside/)

独立 CLI：从组装 contigs 做多层物种鉴定，从 FASTQ reads 做宏基因组谱分析。

## 概述

gside 将标准生信工具（blastn、skani、mash、sourmash、GTDB-Tk、Kraken2、sylph）封装在数据驱动规则之后，跨方法仲裁，输出单一 JSON 判定。**智能体优先设计**：结构化错误字段、自动安装、MCP 服务器、机器可读数据库状态。

```
底层生信工具   blastn / skani / mash / sourmash / gtdbtk / kraken2 / sylph
                    ↑
gside            物种鉴定引擎 + 宏基因组谱分析
```

## 物种鉴定层次

| 层次 | CLI mode | 底层引擎 | 原理 | 说明 |
|---|---|---|---|---|
| **L3** | `gtdbtk` | GTDB-Tk | 系统发育位置 | 最权威；外部数据库 |
| **L2** | `panel` | skani | ANI 精选面板（291 基因组） | 需 panel 层 |
| **L2** | `mash_refseq` | mash | MinHash 距离（RefSeq） | 需 mash 层 |
| **L2** | `sourmash` | sourmash | GTDB gather | 需 GTDB 库 |
| **L1** | `marker` | blastn | 靶基因规则（38 条，80 序列） | mini 层随包分发 |
| **仲裁** | `all` | — | 全方法层级共识（L3 > L2 > L1） | 取最优 |

## 宏基因组谱分析

| 模式 | 引擎 | 输入 | 适用 |
|---|---|---|---|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | 快速筛查、contig 支持、丰度 |
| `sylph` | sylph | 仅 FASTQ | ANI 精度、种级判定 |
| `all` | 双引擎 | FASTQ / FASTA | 合并谱 |

## 基因组质量（`gside validate`）

CheckM2（完整度/污染度）+ GTDB-Tk（分类），写入 `validation.json`，含工具版本、
参数与逐工具状态：

```bash
gside validate contigs.fna -o taxonomy
```

两个数据库均为外部依赖（CheckM2 约 1.7GB；GTDB-Tk 下载约 57GB，细菌运行约
140GB 内存）。见[数据库](installation/databases.md)。

## 快速开始

```bash
# 安装
pixi install
pip install -e .

# 物种鉴定（marker 模式开箱即用）
gside scan contigs.fna

# 多方法鉴定 + 自动安装
gside scan contigs.fna -m all --auto-setup

# 并行批量处理
gside scan *.fna -m all -j 8 -t 4 -o results.tsv

# 宏基因组
gside scan --meta reads_R1.fq reads_R2.fq

# 数据库管理
gside db status
gside db setup --tier all
```

## 智能体友好设计

- **stdout 纯 JSON** — species 恒 exit 0；单法错误在 `methods.<mode>.error`
- **结构化缺库错误** — `suggested_command`、`missing_tier`、`estimated_size`
- **自动安装** — `--auto-setup` 安装缺失层级，进度走 stderr
- **机器可读状态** — `db status --format json`、`db list --format json`
- **MCP 服务器** — `gside mcp` 向 Claude Desktop / MCP 客户端暴露工具
- **流式输出** — tsv/md 逐输入即时写盘

## 输出

所有模式向 stdout 输出统一 JSON，顶层字段为 `analysis_type`、`tool`、`version`、
`contigs`、`methods`、`verdict`。详见
[JSON 输出契约](reference/json-contract.zh.md)。

## 文档

- [快速开始](installation/quick-start.zh.md)
- [数据库](installation/databases.zh.md)
- [标记基因库](installation/marker-database.zh.md)
- [CLI 参考](usage/cli.zh.md)
- [鉴定模式](usage/modes.zh.md)
- [宏基因组](usage/meta.zh.md)
- [分层设计](architecture/layers.zh.md)
- [JSON 输出契约](reference/json-contract.zh.md)
- English version available via the language switcher

## License

MIT
