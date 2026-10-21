# CLI 参考

## gside species

物种鉴定主命令。

### 用法

```bash
gside species <contigs.fasta> [选项]
```

### 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `contigs` | 位置参数 | 无 | 组装后的 contigs FASTA 文件 |
| `-m`，`--mode` | 选项 | `marker` | 鉴定模式（见下） |
| `-D`，`--db-dir` | 选项 | 自动 | 数据库根目录 |
| `-f`，`--format` | 选项 | `json` | 输出格式：`json` / `tsv` / `md` |
| `-n`，`--no-header` | 选项 | 关 | 隐藏 TSV 表头行（仅 tsv；md/json 保持结构） |
| `-a`，`--auto-setup` | 选项 | 关 | 鉴定前自动安装缺失的数据库 tier（进度走 stderr） |
| `-j`，`--jobs N` | 选项 | `1` | 并发鉴定的输入文件数（输出顺序恒为输入顺序） |
| `-t`，`--threads N` | 选项 | `4` | 每次比对调用的线程数（blastn `-num_threads` / skani `-t` / mash `-p`） |
| `-o`，`--output PATH` | 选项 | stdout | 报告写入 PATH 而非 stdout（截断写入；tsv/md 逐输入流式落盘、逐行 flush——批次中途失败已写前缀仍保留；json 为单文档收尾一次写入；stdout 不再输出数据，提示信息走 stderr） |

### 鉴定模式

共 6 种模式，默认 `marker`：

| Mode | 说明 | 数据库需求 |
|---|---|---|
| `marker` | 靶基因组合规则（blastn） | mini（内置） |
| `panel` | ANI 比对精选面板（skani） | panel tier |
| `mash_refseq` | MinHash 距离（mash） | mash tier |
| `sourmash` | sourmash GTDB gather | sourmash GTDB 库 |
| `gtdbtk` | GTDB-Tk 系统发育 placement（classify_wf） | 外部安装的 GTDB-Tk（`GTDBTK_DATA_PATH`） |
| `all` | 运行全部 5 种方法 + 仲裁 | 按需 |

注意 `gtdbtk` 较慢（每基因组分钟级）；`all` 模式也会运行它——如需跳过请显式指定 `-m`。

### 结构化错误字段（面向 agent）

当某方法因数据库 tier 缺失而失败时，该方法的 JSON 条目在 `error` 字符串
之外还携带结构化字段，无需正则解析：

```json
{
  "method": "sourmash",
  "error": "sourmash GTDB database not found in ... (run: gside db setup --tier sourmash)",
  "suggested_command": "gside db setup --tier sourmash",
  "missing_tier": "sourmash",
  "estimated_size": "~3.9GB"
}
```

加上 `--auto-setup` 后，gside 会自行执行 `suggested_command`（进度走
stderr），对该方法重试一次，并重算 verdict。

### 示例

```bash
# 靶基因鉴定（最快，无需额外数据库）
gside species contigs.fna --mode marker

# ANI 鉴定（精确）
gside species contigs.fna --mode panel

# sourmash GTDB gather
gside species contigs.fna --mode sourmash

# GTDB-Tk 系统发育 placement（慢，需外部 gtdbtk + GTDBTK_DATA_PATH）
gside species contigs.fna --mode gtdbtk

# 多法仲裁（marker + panel + mash_refseq + sourmash + gtdbtk）
gside species contigs.fna --mode all

# 指定数据库目录
gside species contigs.fna --mode panel --db-dir /path/to/dbs

# 自动安装缺失数据库 tier 后鉴定
gside species contigs.fna --mode all --auto-setup

# 报告写入文件（stdout 保持为空，提示信息走 stderr）
gside species contigs.fna -m all -o report.json
```

## gside db

数据库管理子命令。

### 子命令

| 子命令 | 说明 |
|---|---|
| `status` | 显示数据库就绪状态 |
| `list` | 列出可用 tier |
| `setup` | 安装/更新数据库 |

`status` 与 `list` 支持 `-f`/`--format [table|tsv|json]`（默认 `table`）。
`tsv` 与 `json` 为机器可读格式（如 `db status --format tsv` 输出
`COMPONENT\tREADY\tTIER\tPATH` 表，READY 列为 `ready`/`missing`）。

### setup 选项

| 选项 | 默认 | 说明 |
|---|---|---|
| `-t`，`--tier` | `panel` | mini / panel / mash / sourmash / meta / meta_sylph / all |
| `-s`，`--source` | 无 | 从本地已有数据目录复制 |
| `-f`，`--format` | `text` | 输出格式：`text` / `json`（json 模式下进度信息走 stderr，结果字典以单个 JSON 文档输出到 stdout；退出码契约不变） |

`all` **仅安装 panel + mash**（约 310MB）。sourmash GTDB 库（约 3.9GB）不含
在 `all` 内，需单独安装：

```bash
gside db setup --tier sourmash
```

宏基因组 tier `meta`（Kraken2 索引，约 5.5GB）与 `meta_sylph`（sylph GTDB
库，约 14GB）需**手动安装**——`db setup` 只打印上游下载指引而不下载。
放置方法见[宏分类学分析](meta.zh.md)。

### 示例

```bash
gside db status
gside db status --format tsv    # 机器可读
gside db setup --tier all
gside db setup --tier sourmash   # 不含在 --tier all 内（~3.9GB）
gside db setup --tier panel --source /path/to/existing/db
```

## gside validate

标记基因组装验证，可选叠加 CheckM2 + GTDB-Tk。

### 用法

```bash
gside validate <contigs.fasta> [选项]
```

### 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `contigs` | 位置参数 | — | 组装后的 contigs FASTA 文件 |
| `--mode` | 选项 | `simple` | 验证深度：`simple`（标记基因）或 `standard`（叠加 CheckM2 + GTDB-Tk） |
| `--output-dir` | 选项 | 自动 | `validation.json` 输出目录 |

`standard` 模式优雅降级：缺 CheckM2/GTDB-Tk 数据库或二进制时只告警并返回
部分结果，不崩溃。通过 `CHECKM2DB` / `GTDBTK_DATA_PATH`（或 `GTDBDB`）
环境变量配置。

### 示例

```bash
gside validate contigs.fna --mode simple
gside validate contigs.fna --mode standard --output-dir ./taxonomy
```

## gside meta

宏基因组物种分析（Kraken2+Bracken 与 sylph），输入 FASTQ 读段或组装
contigs。详见[宏分类学分析](meta.zh.md)。

### 用法

```bash
gside meta <reads.fastq> [reads_2.fastq] [选项]
```

### 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `reads` | 位置参数 | — | 1 个文件 = 单端读段或 contigs FASTA；2 个文件 = 双端 FASTQ |
| `-m`，`--mode` | 选项 | `all` | 分析方法：`kraken2` / `sylph` / `all`（all 合并两份物种列表） |
| `-D`，`--db-dir` | 选项 | 自动 | 数据库根目录（`$GSIDE_DB_DIR` 或 `data/db`） |
| `-f`，`--format` | 选项 | `json` | 输出格式：`json` / `tsv` / `csv` |
| `-o`，`--output PATH` | 选项 | stdout | 报告写入 PATH 而非 stdout（json 为单文档收尾一次写入；stdout 不再输出数据，提示信息走 stderr） |
| `-s`，`--samples PATH` | 选项 | — | 批量清单 `sample_id<TAB>reads[<TAB>reads2]`（表头可省略） |
| `-j`，`--jobs N` | 选项 | `1` | 并发处理的清单样本数（输出顺序保持不变） |
| `-n`，`--no-header` | 选项 | 关 | 隐藏 TSV/CSV 表头行（json 保持结构） |
| `-a`，`--auto-setup` | 选项 | 关 | 自动安装缺失档位并重试该方法一次 |
| `-t`，`--threads N` | 选项 | `4` | 每次工具调用的线程数 |

### 示例

```bash
gside meta sample_1.fastq sample_2.fastq
gside meta sample_1.fastq sample_2.fastq -m kraken2 -f tsv -o profile.tsv
gside meta contigs.fna          # 仅 kraken2（sylph 只接受读段）
gside meta reads.fastq -m sylph
gside meta --samples samples.tsv -j 4 -f csv -o profiles.csv
```

## gside --version

```bash
gside --version
# gside, version 0.1.4

gside --version --json     # 或 -V -J
# {"tool": "gside", "version": "0.1.4"}
```
