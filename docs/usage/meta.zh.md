# 宏分类学分析（metagenomic profiling）

`gside scan --meta` 对宏基因组样本——FASTQ 测序数据或组装 contigs——进行物种
组成分析，报告样本中存在的物种、丰度以及各方法的置信度。它对应
`gside scan` 的群落工作流（加 `--meta`）；单菌株工作流回答"这是什么物种"，
`--meta` 回答"样本里有什么、各占多少"。

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
gside scan --meta <reads.fastq> [reads_2.fastq] [选项]
```

| 选项 | 默认 | 说明 |
|---|---|---|
| `reads`（位置参数） | 无 | 1 个文件 = 单端读段或 contigs FASTA；2 个文件 = 双端 FASTQ |
| `-m`，`--mode` | `all` | `kraken2` / `sylph` / `all` |
| `-D`，`--db-dir` | 自动 | 数据库根目录（`$GSIDE_DB_DIR` 或 `data/db`） |
| `-f`，`--format` | `json` | 输出格式：`json` / `tsv` / `csv` / `md` |
| `-o`，`--output PATH` | stdout | 报告写入 PATH 而非 stdout |
| `-s`，`--samples PATH` | — | 批量清单（`sample_id<TAB>reads[<TAB>reads2]`） |
| `-j`，`--jobs N` | `1` | 使用 `--samples` 时的并发样本数 |
| `-n`，`--no-header` | 关 | 隐藏表头行（tsv/csv/md） |
| `-a`，`--auto-setup` | 关 | 自动安装缺失档位并重试该方法一次 |
| `-t`，`--threads N` | `4` | 每次工具调用的线程数 |

## 批量样本

多样本时改用清单文件而不是位置参数：

```bash
gside scan --meta --samples samples.tsv -j 4 -f tsv -o profiles.tsv
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
  `low`（与 `gside scan` 的 panel ANI 层同阈值）。

## 输出

JSON（默认）——单文档：

```json
{
  "analysis_type": "metagenomic_profiling",
  "tool": "gside",
  "version": "0.1.5",
  "input": {"files": ["reads_1.fastq", "reads_2.fastq"], "type": "fastq_paired"},
  "methods": {
    "kraken2": {"method": "kraken2_bracken", "database": {}, "result": {
      "top_species": [{"name": "...", "taxid": 562, "reads": 850, "fraction": 0.85, "confidence": "high"}],
      "total_reads": 10000, "classified_reads": 9870, "n_species": 4
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
    "richness": 1,
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

与 `gside scan` 一致，`gside scan --meta` 始终以 0 退出：某方法失败时错误
记录在 `methods.<mode>.error`，其余方法继续。FASTA 输入会被 `sylph`
拒绝（其仅接受读段），并给出指向 `-m kraken2` 的结构化错误。

## 方法选择

| 输入 | 推荐 | 理由 |
|---|---|---|
| FASTQ（单端/双端） | `-m all`（默认） | kraken2 + sylph 互补 |
| FASTQ，追求速度 | `-m kraken2` | 秒级；sylph 额外给 ANI |
| FASTQ，追求 ANI 精度 | `-m sylph` | containment ANI，种级 |
| contigs FASTA | `-m kraken2` | sylph 不读 FASTA |

`--meta -m all` 需要 FASTQ（含 sylph）。FASTA 样本请用 `-m kraken2`，否则
命令在运行前就以退出码 2 结束。

## 示例

```bash
# 双端读段，两种方法（默认）
gside scan --meta sample_1.fastq sample_2.fastq

# 仅 kraken2，TSV 写入文件
gside scan --meta sample_1.fastq sample_2.fastq -m kraken2 -f tsv -o profile.tsv

# 组装 contigs（kraken2 支持 FASTA）
gside scan --meta contigs.fna -m kraken2

# 单端读段，sylph ANI
gside scan --meta reads.fastq -m sylph

# 批量清单，一次 4 个样本
gside scan --meta --samples samples.tsv -j 4 -f csv -o profiles.csv

# 自定义数据库根目录
gside scan --meta reads.fastq -D /path/to/db
```

## 性能

| 输入 | kraken2 | sylph |
|---|---|---|
| 5 万条双端 150bp 读段 | 约 1 秒 | 约 1-2 秒 |
| 10× 1.6Mb 基因组 | 约 1 秒 | 约 2 秒 |

## 故障排查

| 现象 | 处理 |
|---|---|
| `sylph requires FASTQ reads` | 改用 `-m kraken2`，或提供读段 |
| `Kraken2 database not found ... (run: gside db setup --tier meta)` | 执行 `gside db setup --tier meta`，或加 `--auto-setup` |
| `Bracken distribution for 150bp reads not found` | 重装该 tier；归档必须含 `database150mers.kmer_distrib` |
| `--jobs requires --samples` | 加 `--samples`，或去掉 `-j` |
| `duplicate sample id` | 清单中的样本 ID 必须唯一 |
| `paired input must be two FASTQ files` | 同一个样本不要混用 FASTQ 与 FASTA |

## 数据库

两个宏基因组 tier 都会自动下载并校验：

```bash
gside db setup --tier meta          # Kraken2 Standard-8 + Bracken（约 5.5GB）
gside db setup --tier meta_sylph    # sylph GTDB r226 .syl2db（约 14GB）
```

- **D5_kraken2**（`meta`）：Ben Langmead 的 Standard-8 包，已内含 Bracken
  50–300mer 分布。就绪条件：`hash.k2d`、`taxo.k2d`、`opts.k2d` 与
  `database150mers.kmer_distrib`（gside 使用 150bp Bracken）。
- **D6_sylph**（`meta_sylph`）：两阶段 `gtdb-r226-c200-dbv2.syl2db`，MD5 固定。
  `.syldb` 与 `.syl2db` 都接受，优先 `.syl2db`。
- `gside scan --meta --auto-setup` 首次使用时安装缺失 tier
- `--source <dir>` 改为复制本地已有数据库
- `gside db status` 查看 `D5_kraken2` / `D6_sylph` 就绪状态
