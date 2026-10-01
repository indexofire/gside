# 分层设计

## 架构定位

gside 在整体生态中的位置：

```
┌─────────────────────────────────────────────────────────┐
│                  Hermes Agent (LLM)                     │
├─────────────────────────────────────────────────────────┤
│                                                         │
│  hermes-bacmap（桥层）                                   │
│  42 tools · GOM · smk 编排 · gbrain · skills            │
│                                                         │
├──────────────────────┬──────────────────────────────────┤
│                      │                                  │
│  gside（本工具）      │  gapit · gmlst                   │
│  物种鉴定引擎         │  基因筛查 · MLST 分型             │
│  marker + ANI + 仲裁  │                                  │
│                      │                                  │
├──────────────────────┴──────────────────────────────────┤
│                                                         │
│  底层生信工具                                            │
│  blastn · skani · mash · sourmash · mmseqs2 · minimap2  │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

## 内部结构

```
gside/
├── src/gside/
│   ├── cli.py              命令行入口（species / db 子命令）
│   ├── config.py           运行时配置（数据目录/二进制发现）
│   ├── db.py               数据库管理（status/setup/list）
│   ├── engine/             算法抽象层
│   │   ├── hits.py         统一 Hit 数据类（BLAST tabular + PAF）
│   │   ├── registry.py     后端注册表（名→类，懒加载）
│   │   ├── read_mapper.py  ReadMapper 门面（BWA/Minimap2）
│   │   └── backends/       可换后端
│   │       ├── blast.py    blastn/blastp/blastx/tblastn
│   │       ├── skani.py    ANI 搜索
│   │       ├── mash.py     MinHash 距离
│   │       ├── kmer.py     mash + sourmash 后端
│   │       ├── minimap2.py PAF 比对
│   │       ├── kma.py      精确比对
│   │       └── mmseqs2.py  线性时间聚类
│   └── analysis/
│       ├── multigene_identifier.py   L1 靶基因组合
│       ├── ani_identifier.py         L2 ANI（panel/skani_gtdb/mash）
│       ├── sourmash_identifier.py    L2 sourmash gather
│       ├── taxonomic_validator.py    L3 GTDB-Tk（代码已迁）
│       └── species_canon.py          物种名规范化
├── data/
│   ├── reference/species/  marker_rules.yaml + markers_v2.fasta
│   └── db/                 大型数据库（gitignored）
│       ├── refseq_panel/   skani sketch（291 基因组）
│       └── mash_refseq/    mash sketch（RefSeq）
└── tests/
```

## 设计原则

| 原则 | 实现 |
|---|---|
| **规则即数据** | 38 条物种规则存 YAML，新增物种不改代码 |
| **后端可换** | engine 抽象层，blastn↔minimap2 一参数切换 |
| **JSON 契约** | 所有输出统一 GOM 兼容格式，下游零解析成本 |
| **优雅降级** | 缺数据库/二进制时返回明确错误，不崩溃 |
| **独立运行** | 不依赖 hermes-bacmap，可独立安装使用 |

## 与 hermes-bacmap 的关系

gside 从 hermes-bacmap 架构拆分而来（V0.9 架构重构）：

| 能力 | 原位置 | 现位置 |
|---|---|---|
| 靶基因鉴定 | bacmap/analysis/multigene | gside |
| ANI 鉴定 | bacmap/analysis/ani_identifier | gside |
| sourmash 鉴定 | bacmap/analysis/sourmash | gside |
| GTDB-Tk | bacmap/analysis/taxonomic_validator | gside |
| 共识仲裁 | bacmap/analysis/species_consensus | bacmap 保留（读 GOM 版） |
| GOM 入库 | bacmap | bacmap |
| smk 编排 | bacmap | bacmap |

bacmap 的 smk 规则（species.smk）已改为**优先调 gside CLI**，无 gside 时 fallback 到内部模块。
