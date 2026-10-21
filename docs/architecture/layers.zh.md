# 分层设计

## 架构定位

gside 是独立的物种鉴定引擎，直接调用底层生信工具：

```
底层生信工具   blastn · skani · mash · sourmash · kraken2 · sylph
                   ↑
gside            物种鉴定引擎：marker 组合规则 + ANI 多后端 + 层级仲裁
```

## 内部结构

```
gside/
├── src/gside/
│   ├── cli.py              命令行入口（species / meta / db / validate / mcp）
│   ├── config.py           运行时配置（数据目录/二进制发现）
│   ├── db.py               数据库管理（status/setup/list）
│   ├── engine/             分析层实际调用的比对后端
│   │   ├── _env.py         二进制发现 + 受检子进程执行
│   │   ├── hits.py         统一 Hit 数据类（BLAST tabular）
│   │   └── backends/
│   │       ├── blast.py    blastn（marker 层）
│   │       ├── skani.py    ANI 搜索（panel 层）
│   │       └── kmer.py     mash + sourmash 后端
│   └── analysis/
│       ├── multigene_identifier.py   L1 靶基因组合
│       ├── ani_identifier.py         L2 ANI（panel/mash_refseq）
│       ├── sourmash_identifier.py    L2 sourmash gather
│       ├── gtdbtk_identifier.py      L3 GTDB-Tk 定位
│       ├── kraken2_identifier.py     meta：kraken2 + bracken
│       ├── sylph_identifier.py       meta：sylph ANI 定量
│       ├── taxonomic_validator.py    CheckM2 + GTDB-Tk 质控（`gside validate`）
│       ├── _envelope.py              共享结果字典外壳
│       └── _thresholds.py            共享 ANI 置信度阈值
├── data/
│   └── db/                 按方案分目录（一法一文件夹）
│       ├── D1_marker/      标记规则 + fasta + BLAST 库（git 跟踪）
│       ├── D2_ani/         skani sketch + 基因组（gitignored）+ manifests/（跟踪的配方）
│       ├── D3_mash/        mash sketch（gitignored）
│       ├── D4_sourmash/    GTDB reps k=31 + lineages（db setup 安装）
│       ├── D5_kraken2/     Kraken2 Standard-8 索引（手动/gitignored）
│       └── D6_sylph/       sylph GTDB-r226 .syldb（手动/gitignored）
└── tests/
```

注：`taxonomic_validator.py` 由 `gside validate` 驱动，做 CheckM2 + GTDB-Tk
组装质控（需 `CHECKM2DB` 与 `GTDBTK_DATA_PATH` 环境）。当前 `species` 命令提供
marker / panel / mash_refseq / sourmash / gtdbtk / all 六种模式。

## 设计原则

| 原则 | 实现 |
|---|---|
| **规则即数据** | 38 条物种规则存 YAML，新增物种不改代码 |
| **后端按方法注入** | 分析层直接选用 blastn / skani / mash / sourmash |
| **JSON 契约** | 所有输出统一 JSON 格式，下游零解析成本 |
| **优雅降级** | 缺数据库/二进制时返回明确错误，不崩溃 |
| **独立运行** | 无平台依赖，可独立安装使用 |
