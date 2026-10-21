# 宏分类学分析（metagenomic profiling）

`gside meta` 对宏基因组样本——FASTQ 测序数据或组装 contigs——进行物种
组成分析，报告样本中存在的物种、丰度以及各方法的置信度。它与
`gside species` 互补：后者从 contigs 鉴定**单一**物种，而 `meta`
回答"样本里有什么、各占多少"。

可用两种方法：

| 方法 | 引擎 | 输入 | 证据 |
|---|---|---|---|
| `kraken2` | Kraken2 + Bracken | FASTQ（单端/双端）或 contigs FASTA | k-mer 读段分类 + Bracken 丰度重估 |
| `sylph` | sylph | 仅 FASTQ 读段 | 相对 GTDB 基因组的 containment ANI |

`all` 模式（默认）下两者都运行，verdict 合并两份物种列表：重名物种保留
置信度更高方法的一行（sylph ANI ≥ 95% 优先于 kraken2 丰度比例），并按
丰度排序。

## 用法

```bash
gside meta <reads.fastq> [reads_2.fastq] [选项]
```

| 选项 | 默认 | 说明 |
|---|---|---|
| `reads`（位置参数） | 无 | 1 个文件 = 单端读段或 contigs FASTA；2 个文件 = 双端 FASTQ |
| `-m`，`--mode` | `all` | `kraken2` / `sylph` / `all` |
| `-D`，`--db-dir` | 自动 | 数据库根目录（`$GSIDE_DB_DIR` 或 `data/db`） |
| `-f`，`--format` | `json` | 输出格式：`json` / `tsv` / `csv` |
| `-o`，`--output PATH` | stdout | 报告写入 PATH 而非 stdout |
| `-s`，`--samples PATH` | — | 批量清单（`sample_id<TAB>reads[<TAB>reads2]`） |
| `-j`，`--jobs N` | `1` | 使用 `--samples` 时的并发样本数 |
| `-n`，`--no-header` | 关 | 隐藏 TSV/CSV 表头 |
| `-t`，`--threads N` | `4` | 每次工具调用的线程数 |

## 批量样本

多样本时改用清单文件而不是位置参数：

```bash
gside meta --samples samples.tsv -j 4 -f tsv -o profiles.tsv
```

`samples.tsv` 为制表符分隔，每行一个样本（表头可省略）：

```
sample_id	reads_1	reads_2
sampleA	sampleA_R1.fq	sampleA_R2.fq
sampleB	sampleB.fq
```

输出顺序始终与清单一致。`-j/--jobs` 控制并发处理的样本数（每个样本内部仍用
`-t` 指定的线程数）。JSON 输出为数组，每个对象带 `sample` 键；`tsv`/`csv`
每行前置一列 `sample`。

输入类型自动检测：两个文件视为双端 FASTQ；单个 `.fna`/`.fa`/`.fasta`
文件视为组装 contigs；其余视为单端 FASTQ。

## 置信度

- **kraken2**：每物种 `fraction_total_reads` ≥ 0.10 → `high`，≥ 0.01 →
  `medium`，否则 `low`（Bracken 重估读段数）。
- **sylph**：校正 containment ANI ≥ 95 → `high`，≥ 93 → `medium`，否则
  `low`（与 `gside species` 的 panel ANI 层同阈值）。

## 输出

JSON（默认）——单文档：

```json
{
  "analysis_type": "metagenomic_profiling",
  "tool": "gside",
  "version": "0.1.2",
  "input": {"files": ["reads_1.fastq", "reads_2.fastq"], "type": "fastq_paired"},
  "methods": {
    "kraken2": {"method": "kraken2_bracken", "database": {}, "result": {
      "top_species": [{"name": "...", "taxid": 562, "reads": 850, "fraction": 0.85, "confidence": "high"}],
      "total_reads": 10000, "classified_reads": 100, "n_species": 4
    }},
    "sylph": {"method": "sylph_ani", "database": {}, "result": {
      "species": [{"name": "Escherichia coli str. K-12", "ani": 97.53, "abundance": 0.7812, "coverage": 264.0, "confidence": "high"}],
      "n_species": 3
    }}
  },
  "verdict": {
    "top_species": [{
      "name": "Escherichia_coli", "display_name": "Escherichia coli",
      "strain": "str. K-12", "species_raw": "Escherichia coli str. K-12",
      "method": "sylph", "abundance": 0.7812, "confidence": "high", "ani": 97.53
    }],
    "richness": 3,
    "basis": ["kraken2", "sylph"]
  }
}
```

TSV（`-f tsv`）——每个合并后物种一行（kraken2 行的 `ani` 为空）：

```
method	name	abundance	confidence	ani
sylph	Escherichia_coli	0.7812	high	97.53
kraken2	Salmonella_enterica	0.20	medium	
```

与 `gside species` 一致，`gside meta` 始终以 0 退出：某方法失败时错误
记录在 `methods.<mode>.error`，其余方法继续。FASTA 输入会被 `sylph`
拒绝（其仅接受读段），并给出指向 `--mode kraken2` 的结构化错误。

## 示例

```bash
# 双端读段，两种方法（默认）
gside meta sample_1.fastq sample_2.fastq

# 仅 kraken2，TSV 写入文件
gside meta sample_1.fastq sample_2.fastq -m kraken2 -f tsv -o profile.tsv

# contigs / 组装 FASTA（kraken2 支持 FASTA，sylph 报结构化错误）
gside meta contigs.fna

# 单端读段 sylph ANI 分析
gside meta reads.fastq -m sylph

# 自定义数据库根目录
gside meta reads.fastq -D /path/to/db
```

## 数据库

两个宏基因组 tier 目前均需**手动安装**——`gside db setup` 会打印指引
而非直接下载：

```bash
gside db setup --tier meta          # → ERROR: not auto-downloadable yet（附 Kraken2 索引页链接）
gside db setup --tier meta_sylph    # → ERROR: not auto-downloadable yet（附 sylph 文档链接）
```

### D5_kraken2（tier `meta`，约 5.5GB）

从 <https://benlangmead.github.io/aws-indexes/k2> 下载预建 Kraken2 索引
（如 Standard-8），解压到 `data/db/D5_kraken2/`，使 `hash.k2d`、
`opts.k2d`、`taxo.k2d` 位于该目录下。Bracken 还需要同一数据库及
`hash.k2d` / `database100mer.kmer_distrib` 文件。

### D6_sylph（tier `meta_sylph`，约 14GB）

从 <https://sylph-docs.github.io/> 下载预建 sylph `.syl2db` 数据库
（如 GTDB），放入 `data/db/D6_sylph/`（使用找到的第一个 `*.syl2db`
文件）。

用 `gside db status` 检查就绪状态（组件 `D5_kraken2` / `D6_sylph`）。
