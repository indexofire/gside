# CLI 参考

## gside scan

两种工作流的统一入口：默认识别单基因组物种，加 `--meta` 分析宏基因组群落。

### 用法

```bash
gside scan <files...> [选项]
gside scan --meta <reads...> [选项]
gside scan --meta --samples samples.tsv [选项]
```

### 模式与输入矩阵

| 模式 | `--meta` | FASTA | FASTQ | 说明 |
|---|---|---|---|---|
| `marker` | 否 | ✅ | ✅ | FASTQ：minimap2（`sr`/`map-ont`）比对，基因 breadth 取所有读段比对区间的并集——约 1× 基因组覆盖即可 |
| `panel` | 否 | ✅ | ❌ | skani ANI |
| `mash_refseq` | 否 | ✅ | ✅ | 读段先 `mash sketch -m 2 -r`；≥10× 可做种级判定，菌株级不可靠 |
| `sourmash` | 否 | ✅ | ✅ | 对读段签名做 gather（1.6Mb 的 10× 数据约 4 分钟）；基于 containment ANI |
| `gtdbtk` | 否 | ✅ | ❌ | 系统发育 placement |
| `all` | 否 | ✅ | ❌ | 全部物种层 + 仲裁（panel/gtdbtk 仅接受组装） |
| `kraken2` | 是 | ✅ | ✅ | 读段或 contigs |
| `sylph` | 是 | ❌ | ✅ | 仅读段 |
| `all` | 是 | ❌ | ✅ | kraken2+sylph；contigs 请用 `-m kraken2` |

`-m` 支持 `+` 组合（如 `marker+sourmash`、`kraken2+sylph`）；重复项自动去重，
顺序保留。`all` 表示当前上下文的全部模式，不能与其他模式组合。

输入类型由扩展名推断（`.fna/.fa/.fasta` 与 `.fq/.fastq`，支持 `.gz`），并在
**任何工具启动前**校验：不匹配时退出码 2 并指出具体模式。`--meta` 同样拒绝
单菌株模式，单菌株模式也拒绝 `kraken2`/`sylph`。由于 `all` 含 sylph，
`--meta -m all` 需要 FASTQ；contigs 请用 `-m kraken2`。

### 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `files` | 位置参数 | — | 一个或多个 FASTA/FASTQ（单菌株：批量；`--meta`：单个样本 1–2 个文件） |
| `-m`，`--mode` | 选项 | `marker` / `all` | 方法（默认 `marker`，`--meta` 时默认 `all`） |
| `--meta` | 开关 | 关 | 输入视为宏基因组群落 |
| `-D`，`--db-dir` | 选项 | 自动 | 数据库根目录 |
| `-f`，`--format` | 选项 | `json` | `json` / `tsv` / `md` / `csv` |
| `-o`，`--output PATH` | 选项 | stdout | 报告写入 PATH（行式格式随输入完成流式写入；json 最后一次性写入） |
| `-j`，`--jobs N` | 选项 | `1` | 单菌株并发输入数，或 `--meta --samples` 的并发样本数 |
| `-t`，`--threads N` | 选项 | `4` | 每次工具调用的线程数 |
| `-n`，`--no-header` | 选项 | 关 | 隐藏表头行（tsv/csv/md） |
| `-a`，`--auto-setup` | 选项 | 关 | 自动安装缺失档位并重试该方法一次 |
| `-s`，`--samples PATH` | 选项 | — | `--meta` 批量清单：`sample_id<TAB>reads [<TAB>reads2]` |

### 示例

```bash
# 单基因组（isolate 或 MAG）
gside scan contigs.fna                      # marker
gside scan contigs.fna -m panel -f csv
gside scan contigs.fna -m marker+sourmash
gside scan *.fna -m marker -j 8 -o results.tsv

# 长读长单菌株（marker 把读段当序列比对）
gside scan ont_reads.fastq -m marker

# 宏基因组群落
gside scan reads_R1.fq reads_R2.fq --meta                 # kraken2 + sylph
gside scan contigs.fna --meta -m kraken2 -f md
gside scan --meta --samples samples.tsv -j 4 -f csv -o profiles.csv
```

## gside species（别名）

`gside species <files...>` 是 `gside scan`（不带 `--meta`）的隐藏兼容别名：
选项与输出完全一致，模式限定为 `marker|panel|mash_refseq|sourmash|gtdbtk|all`。

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

## gside meta（别名）

`gside meta <reads...>` 是 `gside scan --meta` 的隐藏兼容别名：选项（含
`--samples`）与输出一致，模式为 `kraken2|sylph|all`。

## gside --version

```bash
gside --version
# gside, version 0.1.5

gside --version --json     # 或 -V -J
# {"tool": "gside", "version": "0.1.5"}
```
