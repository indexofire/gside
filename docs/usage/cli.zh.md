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

**双端 FASTQ 算一个样本。** 位置参数中带配对标记的文件（`_R1`/`_R2`、
`_1`/`_2`、`.1`/`.2`、`-1`/`-2`，可带 `_001`，大小写不敏感，支持 `.gz`）会
自动合并：两条 mate 一起送入 `marker`（一次 minimap2）、`mash_refseq` 与
`sourmash`（各一次 sketch），输出中的 `contigs` 取第一条 mate。未配对的文件
各自成样本，因此 `gside scan A_R1.fq A_R2.fq B_R1.fq B_R2.fq -m marker`
是 **2** 个样本而不是 4 个。FASTA 永不配对。`--meta` 下两个位置参数文件始终
视为一个样本；文件更多时同样按 mate 名分组。

### 解读 verdict

```json
"verdict": {
  "species": "Campylobacter_jejuni",
  "display_name": "Campylobacter jejuni",
  "strain": "subsp. jejuni NCTC 11168",
  "confidence": "high",
  "basis": ["panel"],
  "supporting": ["marker", "sourmash"],
  "agreement": "3/3"
}
```

| 字段 | 含义 |
|---|---|
| `species` | 规范 `Genus_species`（机器契约） |
| `display_name` / `strain` | 可读名与胜出方法给出的菌株信息 |
| `confidence` | high / medium / low（按方法阈值） |
| `basis` | 产生判定的单个方法（仅胜出者） |
| `supporting` | 其他成功且与 `species` 一致的方法 |
| `agreement` | `一致数/成功数`——`2/3` 表示有一层不一致或失败 |

`methods.<mode>.error` 说明失败的方法；没有任何判定时 `species` 为
`Unknown`。逐字段说明见 [JSON 输出契约](../reference/json-contract.md)。

### 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `files` | 位置参数 | — | 一个或多个 FASTA/FASTQ；FASTQ 的 mate 自动合并（`_R1/_R2`、`_1/_2`、`.1/.2`）；未配对文件各自成样本（`--meta` 下 2 个文件 = 一个样本） |
| `-m`，`--mode` | 选项 | `marker` / `all` | 方法（默认 `marker`，`--meta` 时默认 `all`） |
| `--meta` | 开关 | 关 | 输入视为宏基因组群落 |
| `-D`，`--db-dir` | 选项 | 自动 | 数据库根目录 |
| `-f`，`--format` | 选项 | `json` | `json` / `tsv` / `md` / `csv` |
| `-o`，`--output PATH` | 选项 | stdout | 报告写入 PATH（行式格式随样本完成流式写入；json 最后一次性写入） |
| `-j`，`--jobs N` | 选项 | `1` | 并发样本数（单菌株批量、位置参数 `--meta` 批量或 `--samples`） |
| `-t`，`--threads N` | 选项 | `4` | 每次工具调用的线程数 |
| `-n`，`--no-header` | 选项 | 关 | 隐藏表头行（tsv/csv/md） |
| `-a`，`--auto-setup` | 选项 | 关 | 自动安装缺失档位并重试该方法一次 |
| `-s`，`--samples PATH` | 选项 | — | `--meta` 批量清单：`sample_id<TAB>reads [<TAB>reads2]` |

多样本且 stderr 为 TTY 时会显示进度条。流式输出行会交给进度条协调：每次写入前
暂停进度条、写完再恢复，因此 JSON/TSV/CSV/MD 结果会实时出现且不会被下一次刷新
覆盖。stdout 被管道或重定向时不显示进度条，stdout 只有数据。

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

## gside db

数据库管理。`gside db` 不带子命令时显示本帮助（等价于 `gside db -h`）并以 0 退出。

| 子命令 | 用途 |
|---|---|
| `status` | 每个组件的就绪状态（名称、tier、路径） |
| `list` | tier 列表（就绪、大小、说明） |
| `setup` | 安装/更新单个 tier |
| `validate` | 校验自定义 marker 数据库 |

`status` 与 `list` 接受 `-f/--format table|tsv|csv|json`（默认 `table`）；
`setup` 接受 `-f/--format text|json`（默认 `text`）。人类可读输出使用带边框的
Rich 表格；脚本使用 `json`/`tsv`/`csv`。

### status

```bash
gside db status
# ┌─────────────┬─────────┬────────────┬──────────────────────────┐
# │ COMPONENT   │ READY   │ TIER       │ PATH                     │
# ├─────────────┼─────────┼────────────┼──────────────────────────┤
# │ D1_marker   │ ready   │ mini       │ data/db/D1_marker        │
# │ D4_sourmash │ missing │ sourmash   │ data/db/D4_sourmash      │
# └─────────────┴─────────┴────────────┴──────────────────────────┘

gside db status -f json    # {component: {ready, path, tier}}
gside db status -f csv     # COMPONENT,READY,TIER,PATH
```

### list

```bash
gside db list              # 带边框的 tier 表格
gside db list -f tsv       # TIER	READY	SIZE	DESCRIPTION
```

### setup

| 选项 | 默认 | 说明 |
|---|---|---|
| `-t`，`--tier` | `panel` | `mini|panel|mash|sourmash|meta|meta_sylph|all` |
| `-s`，`--source` | — | 从本地已有数据库目录复制，而不是下载 |
| `-f`，`--format` | `text` | `text` / `json`（json 把进度写到 stderr） |

- `mini` 随包分发（始终就绪）
- `all` = **panel + mash**；`sourmash` 需单独安装（约 3.9GB）
- `meta`（Kraken2 Standard-8 + Bracken，约 5.5GB）与 `meta_sylph`
  （GTDB r226 `.syl2db`，约 14GB）会自动下载并校验
- CheckM2 与 GTDB-Tk **不是** tier；见
  [数据库](../installation/databases.md)

```bash
gside db setup --tier panel
gside db setup --tier all
gside db setup --tier sourmash
gside db setup --tier panel --source /path/to/existing/data/db
```

## gside validate

组装**质量评估**：CheckM2（完整度/污染度）+ GTDB-Tk 分类。不再运行 marker
规则——物种鉴定请用 `scan`。

### 用法

```bash
gside validate <contigs.fasta> -o <输出目录> [选项]
```

### 解读 verdict

```json
"verdict": {
  "species": "Campylobacter_jejuni",
  "display_name": "Campylobacter jejuni",
  "strain": "subsp. jejuni NCTC 11168",
  "confidence": "high",
  "basis": ["panel"],
  "supporting": ["marker", "sourmash"],
  "agreement": "3/3"
}
```

| 字段 | 含义 |
|---|---|
| `species` | 规范 `Genus_species`（机器契约） |
| `display_name` / `strain` | 可读名与胜出方法给出的菌株信息 |
| `confidence` | high / medium / low（按方法阈值） |
| `basis` | 产生判定的单个方法（仅胜出者） |
| `supporting` | 其他成功且与 `species` 一致的方法 |
| `agreement` | `一致数/成功数`——`2/3` 表示有一层不一致或失败 |

`methods.<mode>.error` 说明失败的方法；没有任何判定时 `species` 为
`Unknown`。逐字段说明见 [JSON 输出契约](../reference/json-contract.md)。

### 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `contigs` | 位置参数 | — | 组装 contigs FASTA（必须存在） |
| `-o`，`--output-dir` | 选项 | **必填** | `validation.json` 的输出目录（不存在会自动创建；是目录，不是文件） |
| `-t`，`--threads` | 选项 | `4` | 传给 CheckM2 `--threads` 与 GTDB-Tk `--cpus` 的线程数 |

两项工具都成功才算完整报告；缺失数据库/二进制会按工具记录，不会崩溃：

| `status` | 含义 | 退出码 |
|---|---|---|
| `ok` | CheckM2 与 GTDB-Tk 都成功 | 0 |
| `partial` | 只有一个成功 | 1 |
| `error` | 两个都未成功 | 1 |

低质量（如 `WARNING`）是**结果**而非命令失败：只要两工具跑通仍退出 0。
通过 `CHECKM2DB` 与 `GTDBTK_DATA_PATH`（或 `GTDBDB`）配置；安装命令见
[数据库](../installation/databases.md)。

### `validation.json`

```json
{
  "schema": "gside.validate/1",
  "genome": {"name": "GCF_000006945.2_ASM694v2_genomic", "path": "/abs/contigs.fna"},
  "parameters": {"threads": 4, "checkm2_db": "/db/checkm2", "gtdbtk_db": "/db/gtdb"},
  "software": {"gside": "0.2.0", "checkm2": "1.0.2", "gtdbtk": "2.4.0"},
  "tools": {"checkm2": "ok", "gtdbtk": "ok"},
  "quality": {"completeness": 99.1, "contamination": 0.5, "verdict": "PASS"},
  "taxonomy": {"classification": "d__Bacteria;p__...;s__Salmonella enterica", "note": ""},
  "status": "ok",
  "errors": [],
  "interpretation": "Genome quality: PASS (completeness=99.1%, contamination=0.5%); GTDB-Tk: ...",
  "created_at": "2026-10-21T14:30:00+00:00"
}
```

### 示例

```bash
gside validate contigs.fna -o ./taxonomy
gside validate contigs.fna -o ./taxonomy -t 8
```

## gside mcp

MCP（Model Context Protocol）stdio 服务器：把 gside 作为工具提供给 agent。
stdin/stdout 上走 JSON-RPC 2.0；stdout 只输出协议帧（进度与提示走 stderr）。

```bash
gside mcp        # 或控制台脚本：gside-mcp
```

工具：`species`、`meta`、`db_status`、`db_list`、`db_setup`。每个都镜像对应
CLI 命令并返回同样的 JSON 文本；失败返回 `isError` 包络，不会让服务器崩溃。

## gside species / gside meta（别名）

隐藏兼容别名：`gside species` = 不带 `--meta` 的 `gside scan`；
`gside meta` = `gside scan --meta`。模式选择仍按上下文（单菌株/宏数据）限定。
新脚本请使用 `scan`。

## 环境变量

| 变量 | 使用者 | 含义 |
|---|---|---|
| `GSIDE_DATA_DIR` | 全部 | 数据根目录（默认 `<repo>/data`，否则用包内数据） |
| `GSIDE_DB_DIR` | 全部 | 数据库根目录（默认 `<data>/db`） |
| `GSIDE_PIXI_BIN` | 全部 | 额外 `bin/`（`blastn`/`minimap2`/`skani`/`mash`…） |
| `GSIDE_ENGINE` | `scan` | marker/mash/sourmash/panel 路径引擎：`legacy`（默认，保留）/`gseek`/`shadow`；`shadow` 返回 legacy 结果 |
| `GSIDE_SOURMASH_SIG_SHA256` | `db setup` | 覆盖 sourmash 签名摘要 |
| `CHECKM2DB` | `validate` | CheckM2 的 `uniref100.KO.1.dmnd` 路径 |
| `GTDBTK_DATA_PATH` / `GTDBDB` | `validate`、`scan -m gtdbtk` | GTDB-Tk 数据库目录 |

## 退出码

| 命令 | 退出码 | 含义 |
|---|---|---|
| `scan` | 0 | 始终为 0；单方法失败写入 `methods.<mode>.error` |
| `scan` | 2 | 用法错误（未知模式、输入不匹配、参数错误） |
| `validate` | 0 | CheckM2 与 GTDB-Tk 都运行（质量结论不影响退出码） |
| `validate` | 1 | `partial` 或 `error`（工具/数据库缺失） |
| `db setup` | 0 / 1 | 仅当结果消息以 `ERROR` 开头时为 1 |
| `db status` / `db list` / `--version` | 0 | 信息性输出 |

## 从 0.1.4 升级

- `gside species` / `gside meta` 仍作为 `gside scan`（不带/带 `--meta`）的
  隐藏别名可用。
- 物种标签改为规范 `Genus_species`（`verdict.species`、
  `methods.<mode>.species`）；展示名移到 `verdict.display_name`，菌株移到
  `verdict.strain`。verdict 新增 `supporting` 与 `agreement`。
- 宏基因组行使用同一规范 `name`，并带 `display_name`/`strain` 与分数
  （0-1）形式的 `abundance`。
- `validate` 不再有 `-m/--mode`；`-o/--output-dir` 为必填，`validation.json`
  新增 schema/genome/parameters/software/tools/status/errors/created_at。
- `--meta -m all` 需要 FASTQ（含 sylph）；contigs 请用 `-m kraken2`。

## gside --version

没有单独的 `gside version` 子命令。使用根命令的 `-V/--version`；`-J` 输出
JSON。

```bash
gside -V
# gside, version 0.2.0

gside --version --json     # 或 -V -J
# {"tool": "gside", "version": "0.2.0"}
```
