# 分层设计

## 架构定位

gside 是独立的物种鉴定引擎，直接调用底层生信工具：

```
底层生信工具   blastn · skani · mash · sourmash · mmseqs2 · minimap2
                   ↑
gside            物种鉴定引擎：marker 组合规则 + ANI 多后端 + 层级仲裁
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
│       ├── ani_identifier.py         L2 ANI（panel/mash）
│       ├── sourmash_identifier.py    L2 sourmash gather
│       ├── taxonomic_validator.py    L3 GTDB-Tk（已通过 gside validate 接入）
│       └── species_canon.py          物种名规范化（未接入 CLI）
├── data/
│   └── db/                 按方案分目录（一法一文件夹）
│       ├── D1_marker/      标记规则 + fasta + BLAST 库（git 跟踪）
│       ├── D2_ani/         skani sketch + 基因组（gitignored）+ manifests/（跟踪的配方）
│       ├── D3_mash/        mash sketch（gitignored）
│       └── D4_sourmash/    GTDB reps k=31 + lineages（db setup 安装）
└── tests/
```

注：`taxonomic_validator.py` 已通过 `gside validate` 接入（simple 模式开箱即用，
standard 模式需 CheckM2/GTDB-Tk 环境）；`species_canon.py` 仍未接入 CLI。
当前 `species` 命令提供 marker / panel / mash_refseq / sourmash / all 五种模式。

## 设计原则

| 原则 | 实现 |
|---|---|
| **规则即数据** | 38 条物种规则存 YAML，新增物种不改代码 |
| **后端可换** | engine 抽象层，blastn↔minimap2 一参数切换 |
| **JSON 契约** | 所有输出统一 JSON 格式，下游零解析成本 |
| **优雅降级** | 缺数据库/二进制时返回明确错误，不崩溃 |
| **独立运行** | 无平台依赖，可独立安装使用 |
