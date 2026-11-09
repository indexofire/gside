# 快速安装

从零到完成一次物种鉴定或基因组质量报告。marker 层随包分发，
首次鉴定只需数秒；可选层再叠加 ANI/sketch/读段方法。

## 工作流速查

| 目标 | 命令 |
|---|---|
| 单基因组，最快 | `gside scan contigs.fna` |
| 单基因组，最强证据 | `gside scan contigs.fna -m all` |
| 单菌株测序读段 | `gside scan reads.fq -m marker` |
| 单菌株双端读段 | `gside scan S_R1.fq S_R2.fq -m marker`（一个样本） |
| 组装质量（QC） | `gside validate contigs.fna -o taxonomy` |
| 批量 + 流式输出 | `gside scan *.fna -m all -j 8 -o calls.tsv` |

## 1. 前置条件

- Linux x86_64，Python ≥ 3.12
- [pixi](https://pixi.sh) — 提供 `blastn`、`minimap2`、`skani`、`mash`、
  `sourmash`、GTDB-Tk（只装 Python 包无法获得这些二进制）

## 2. 安装

```bash
git clone https://github.com/indexofire/gside.git
cd gside

pixi install          # 生信二进制
pip install -e .      # gside 本体（或直接 pip install gside）
```

若从 PyPI 安装且二进制在别处，设置 `GSIDE_PIXI_BIN` 指向那个 `bin/` 目录。

## 3. 验证

```bash
gside -V              # gside, version 0.2.4
gside db status       # mini（D1_marker）必须是 ✅
```

`db status` 会列出每个组件、所属 tier 与路径；❌ 表示缺失，可用
`gside db setup --tier <tier>` 安装。

## 4. 第一次鉴定（秒级）

```bash
gside scan contigs.fna
```

- 默认模式是 `marker`（随包 mini 层，blastn）
- `scan` 会按扩展名推断 FASTA/FASTQ，并在**运行前**拒绝不匹配的模式
- 输出是 stdout 上的一份 JSON：每个方法的 `methods` 加上 `verdict`
  （规范 `Genus_species`、置信度、basis、supporting、agreement）

常见单基因组用法：

```bash
gside scan contigs.fna -m marker                 # 靶基因
gside scan contigs.fna -m marker+sourmash        # 组合方法
gside scan contigs.fna -m all                    # 全部基因组层 + 仲裁
gside scan contigs.fna -m all -j 8 -o run.tsv    # 批量 + 流式 TSV
```

## 5. 读段

```bash
# 纯培养 / 单菌株测序读段
gside scan reads.fq -m marker                    # minimap2，约 1× 覆盖即可
gside scan reads.fq -m marker+mash_refseq

# 双端：R1/R2（以及 _1/_2、.1/.2）是**一个**样本，不是两个
gside scan S_R1.fq S_R2.fq -m marker
gside scan A_R1.fq A_R2.fq B_R1.fq B_R2.fq -m marker   # 两个样本

