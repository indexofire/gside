# 快速安装

从零到完成一次物种鉴定、群落分析或基因组质量报告。marker 层随包分发，
首次鉴定只需数秒；可选层再叠加 ANI/sketch/读段方法。

## 工作流速查

| 目标 | 命令 |
|---|---|
| 单基因组，最快 | `gside scan contigs.fna` |
| 单基因组，最强证据 | `gside scan contigs.fna -m all` |
| 单菌株测序读段 | `gside scan reads.fq -m marker` |
| 单菌株双端读段 | `gside scan S_R1.fq S_R2.fq -m marker`（一个样本） |
| 读段宏基因组 | `gside scan R1.fq R2.fq --meta` |
| Contigs 宏基因组 | `gside scan contigs.fna --meta -m kraken2` |
| 组装质量（QC） | `gside validate contigs.fna -o taxonomy` |
| 批量 + 流式输出 | `gside scan *.fna -m all -j 8 -o calls.tsv` |

## 1. 前置条件

- Linux x86_64，Python ≥ 3.12
- [pixi](https://pixi.sh) — 提供 `blastn`、`minimap2`、`skani`、`mash`、
  `sourmash`、`kraken2`、`bracken`、`sylph`（只装 Python 包无法获得这些二进制）

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
gside -V              # gside, version 0.2.0
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

# 宏基因组群落（加 --meta）
gside scan reads_R1.fq reads_R2.fq --meta        # kraken2 + sylph
gside scan contigs.fna --meta -m kraken2         # 宏基因组 contigs
```

`--meta` 切换到群落工作流，只接受 `kraken2`、`sylph`、`all`；不加 `--meta`
则运行单菌株模式。

## 6. 数据库

| Tier | 内容 | 大小 | 安装 |
|---|---|---|---|
| `mini` | 标记规则 + 序列 + BLAST 库 | ~1MB | 随包分发（始终就绪） |
| `panel` | 精选 291 基因组面板（skani） | ~130MB sketch / ~1GB 基因组 | `gside db setup --tier panel` |
| `mash` | RefSeq MinHash sketch | ~179MB | `gside db setup --tier mash` |
| `sourmash` | GTDB r226 reps + lineages | ~3.9GB | `gside db setup --tier sourmash` |
| `meta` | Kraken2 Standard-8 + Bracken | ~5.5GB | `gside db setup --tier meta` |
| `meta_sylph` | sylph GTDB r226（`.syl2db`） | ~14GB | `gside db setup --tier meta_sylph` |
| `all` | panel + mash | ~310MB | `gside db setup --tier all` |

所有下载都做校验（panel SHA256、mash MD5、sourmash SHA256、Kraken2 逐文件
MD5、sylph MD5）。`--auto-setup` 会在首次使用时安装缺失层并重试该方法一次：

```bash
gside scan contigs.fna -m panel --auto-setup
```

## 7. 基因组质量（`gside validate`）

`validate` 回答另一个问题：**组装质量如何、GTDB-Tk 把它放在哪里？**
它运行 CheckM2 + GTDB-Tk（不跑 marker 规则）：

```bash
gside validate contigs.fna -o ./taxonomy          # 输出目录必填
gside validate contigs.fna -o ./taxonomy -t 8
```

前置条件（外部工具，不属于 `db setup` 的 tier）：

- CheckM2 二进制 + 数据库（约 1.7GB）——`CHECKM2DB` 指向
  `uniref100.KO.1.dmnd` 文件
- GTDB-Tk 二进制 + R232 数据库（下载约 57GB、磁盘约 98GB、细菌运行约
  140GB 内存）——设置 `GTDBTK_DATA_PATH`

报告（`validation.json`）包含基因组名/路径、工具版本、参数、CheckM2
完整度/污染度与结论、GTDB-Tk 分类、逐工具状态、错误与时间戳。
只有两工具都跑通才退出 0；低质量本身仍退出 0。安装命令见
[数据库](databases.md)。

## 8. 批量与格式

`scan` 支持多输入与四种输出格式：

```bash
gside scan *.fna -m marker -j 8 -f tsv -o calls.tsv
gside scan *.fna -m marker -f csv
gside scan contigs.fna -m marker -f md
```

- `-j/--jobs` 并行处理样本（双端 mate 算一个样本）；输出顺序等于样本顺序
- `-o` 写入文件（行式格式随样本完成流式写入）
- `-n/--no-header` 隐藏表头
- 单样本 JSON 为对象，多样本为数组；`meta` 批量（见下）会额外带 `sample` 键

宏基因组批量支持位置参数分组或清单：

```bash
gside scan --meta A_R1.fq A_R2.fq B_R1.fq B_R2.fq -j 2 -f tsv   # 两个样本
gside scan --meta --samples samples.tsv -j 4 -f tsv -o profiles.tsv
```

```
sample_id	reads_1	reads_2
sampleA	sampleA_R1.fq	sampleA_R2.fq
sampleB	sampleB.fq
```

## 9. 故障排查

| 现象 | 处理 |
|---|---|
| `blastn not found` | 在 `pixi run -e dev …` 中运行，或设置 `GSIDE_PIXI_BIN` |
| `mode 'panel' does not accept FASTQ input` | 改用 `--meta -m kraken2`，或读段可用模式（`marker`、`mash_refseq`、`sourmash`） |
| `mode 'kraken2' is not available without --meta` | 加 `--meta` |
| `--meta -m all` 报 sylph/FASTA | `all` = kraken2+sylph；contigs 请用 `-m kraken2` |
| `db status` 显示 ❌ | `gside db setup --tier <tier>`（或 `scan` 加 `--auto-setup`） |
| `validate` 退出 1 且 `partial` | CheckM2/GTDB-Tk 缺一；查看 `validation.json` 的 `tools`/`errors` |

## 下一步

- [CLI 参考](../usage/cli.md) — 全部命令、参数、退出码
- [鉴定模式](../usage/modes.md) — 各模式原理与输入
- [宏基因组分析](../usage/meta.md) — 群落工作流
- [JSON 输出契约](../reference/json-contract.md) — 逐字段说明
