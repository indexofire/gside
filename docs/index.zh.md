# gside — 基因组物种鉴定引擎

独立的病原微生物物种鉴定 CLI 工具，从 hermes-bacmap 架构拆分而来。

## 概述

gside（**G**enome **S**pecies **ID**entification **E**ngine）提供多层次的物种鉴定能力，
可作为独立命令行工具使用，也可被智能体平台（如 hermes-bacmap）集成调用。

```
底层工具     blastn / skani / mash / sourmash / mmseqs2
                ↑
gside          物种鉴定引擎：靶基因组合规则 + ANI 多后端 + 层级仲裁
                ↑
hermes-bacmap  Hermes 桥：GOM 入库 / smk 编排 / gbrain 记忆
```

## 鉴定层次

| 层次 | CLI mode | 底层引擎 | 原理 | 验证 |
|---|---|---|---|---|
| **L1** | `marker` | blastn | 靶基因组合规则（38 条，81 序列） | 17/17 |
| **L2** | `panel` | skani | ANI 比对精选参考面板（291 基因组） | 17/17 |
| **L2** | `mash_refseq` | mash | MinHash 距离（RefSeq sketch） | 17/17 |
| **L2** | `sourmash` | sourmash | GTDB gather（LCA 分类） | 需下库 |
| **仲裁** | `all` | — | 层级共识（ANI 层 > marker 层） | ✓ |

## 快速开始

```bash
# 安装
pixi install && pip install -e .

# 物种鉴定
gside species contigs.fna --mode marker
gside species contigs.fna --mode panel
gside species contigs.fna --mode all

# 数据库管理
gside db status
gside db setup --tier all
```

## 输出

所有模式输出统一的 JSON 格式（GOM 兼容契约），包含 `analysis_type`、
`method`、`database`、`result`、`verdict` 字段，可直接被下游系统消费。

## License

MIT
