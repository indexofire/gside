# gside 文档

`gside` 是一个用于多层微生物物种鉴定与宏基因组谱分析的 Python 库和 CLI。本文档覆盖公共 Python API、`gside` 命令行、参考数据库、鉴定模式与 JSON 契约。

## 从这里开始

| 我想… | 阅读 |
| --- | --- |
| 安装并跑第一个样本 | [快速开始](installation/quick-start.zh.md) |
| 在 Python 中调用鉴定引擎 | [Python API](usage/python-api.zh.md) |
| 查所有 CLI 参数 | [CLI 参考](usage/cli.zh.md) |
| 理解 marker/ANI/sketch/GTDB 模式 | [鉴定模式](usage/modes.zh.md) |
| 做宏基因组群落分析 | [宏基因组谱分析](usage/meta.zh.md) |
| 准备数据库 | [数据库](installation/databases.zh.md) |
| 安全消费 JSON 输出 | [JSON 输出契约](reference/json-contract.zh.md) |
| 把引擎迁移到 gseek | [gseek 迁移](development/gseek-migration.zh.md) |

## gside 做什么

```text
组装 contigs ── marker / panel / mash / sourmash / GTDB-Tk ── 判定
FASTQ reads ─── marker screening / kraken2 / sylph ─────────── 群落谱
组装 contigs ── CheckM2 + GTDB-Tk ─────────────────────────── 组装质量
```

两个接口共享同一批引擎：

- **Python 库** —— `gside.api` 暴露各鉴定函数与带类型结果，适合 notebook、流水线和测试。
- **CLI** —— `gside scan` 增加分层仲裁、数据库安装、批处理、表格/JSON 输出与 MCP 服务器。

## 安装概览

```bash
pip install gside          # Python 3.12；含 marker DB
pixi install               # 含全部外部工具的 dev 环境
```

`blastn`、`skani`、`mash`、`sourmash`、`kraken2`、`sylph`、GTDB-Tk 等外部工具需另行安装。见[快速开始](installation/quick-start.zh.md)与[数据库](installation/databases.zh.md)。

## 库示例

```python
from gside.api import identify_by_ani, identify_multigene

marker = identify_multigene("contigs.fna", threads=4)
panel = identify_by_ani("contigs.fna", mode="panel", threads=8)

print(marker.species, marker.confidence)
print(panel.to_dict()["result"]["species"])
```

完整函数/结果表见 [Python API](usage/python-api.zh.md)。

## CLI 示例

```bash
gside -V
gside scan contigs.fna -m all --auto-setup
gside scan --meta reads_R1.fq reads_R2.fq
gside validate contigs.fna -o qc/
gside db status
gside db list
```

`gside db status` 和 `gside db list` 在终端中输出带边框的 Rich 表格；供程序读取时使用 `-f json`、`-f tsv` 或 `-f csv`。

## 设计原则

- **规则优先。** marker 组合、ANI 置信区间与分层仲裁都是数据驱动、可测试的。
- **单一结构化结果。** 每个方法输出 JSON 信封；CLI 做共识但不丢弃逐方法证据。
- **清晰失败。** 缺工具/数据库时给出可操作错误，包括建议安装命令和预估大小。
- **默认可复现。** 结果记录数据库版本、工具版本与方法参数。
- **可选加速。** `GSIDE_ENGINE=gseek` 可替换部分子进程路径，同时保留 legacy 与 shadow。

## 目录结构

```text
src/gside/
  api.py                 公共 Python facade
  cli.py                 Click CLI
  analysis/              鉴定引擎
  engine/                后端 + 可选 gseek 桥
  db.py                  数据库层与下载
  mcp.py                 MCP 服务器
docs/                    中英文文档
data/db/D1_marker/       随包 marker 数据库
tests/                   单元与 parity 测试
```

## 许可证

MIT。项目概览与生态链接见 [README](https://github.com/indexofire/gside)。
