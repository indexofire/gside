# 数据库管理

gside 的参考数据库通过 `gside db` 子命令管理，支持分层数据库安装。

## 数据库分层

| Tier | 包含 | 大小 | 用途 |
|---|---|---|---|
| **mini** | marker_rules.yaml + markers.fasta + 预建 BLAST 库 | ~1MB | 靶基因组合鉴定（L1），随仓库分发 |
| **panel** | D2_ani（291 基因组 skani sketch） | ~130MB | ANI 鉴定（L2） |
| **mash** | D3_mash（RefSeq MinHash sketch） | ~179MB | 距离鉴定（L2） |
| **sourmash** | D4_sourmash（GTDB reps k=31 + lineages） | ~3.9GB | gather 分类（sourmash） |
| **meta** | D5_kraken2（Kraken2 Standard-8 索引） | ~5.5GB | 宏基因组分析（`gside scan --meta`） |
| **meta_sylph** | D6_sylph（sylph GTDB 数据库） | ~14GB | 基于 ANI 的宏基因组分析（`gside scan --meta`） |
| **all** | panel + mash | ~310MB | 全层次能力（不含 sourmash） |

### 目录命名

`data/db/` 下目录按鉴定层次命名：

| 旧名 | 新名 | 对应方法 |
|---|---|---|
| `refseq_panel` | `D2_ani` | panel（skani ANI） |
| `mash_refseq` | `D3_mash` | mash_refseq（MinHash） |
| `sourmash_gtdb` | `D4_sourmash` | sourmash（GTDB gather） |

已有旧目录按如下迁移：

```bash
mv data/db/refseq_panel data/db/D2_ani
mv data/db/mash_refseq data/db/D3_mash
mkdir -p data/db/D4_sourmash
```

注意：GitHub Release 预建包内层 `panel.sketch/` 目录名不变，下载解压流程不受影响。

## 命令

### 查看状态

```bash
gside db status
```

```
gside database status
──────────────────────────────────────────
  ✅ D1_marker          tier=mini   path=data/db/D1_marker
  ✅ D2_ani        tier=panel  path=data/db/D2_ani
  ✅ D3_mash         tier=mash   path=data/db/D3_mash
  ✅ D4_sourmash     tier=sourmash path=data/db/D4_sourmash
```

### 安装/更新

```bash
gside db setup --tier panel          # 仅面板（ANI），默认
gside db setup --tier mash           # 仅 MinHash
gside db setup --tier sourmash       # GTDB gather 库（~3.9GB）
gside db setup --tier meta           # Kraken2 + Bracken（~5.5GB）
gside db setup --tier meta_sylph     # sylph GTDB 库（~14GB）
gside db setup --tier all            # panel + mash（不含 sourmash）
```

不带 `--tier` 时默认安装 `panel`。

`--tier all` 仅安装 panel + mash。sourmash GTDB 库（约 3.9GB）不含在 `all`
内，需另行执行 `gside db setup --tier sourmash` 单独安装。

宏基因组 tier 现在支持自动下载：

```bash
gside db setup --tier meta         # Kraken2 Standard-8 + Bracken（约 5.5GB）
gside db setup --tier meta_sylph   # sylph GTDB r226，.syl2db（约 14GB）
```

- **meta** 下载 Ben Langmead 的 Standard-8 包（Kraken2 索引 + Bracken
  50–300mer 分布）；解压后按官方 `.md5` 清单逐文件校验；就绪条件包含
  `database150mers.kmer_distrib`（gside 使用 150bp 读长）。
- **meta_sylph** 从 HTTPS 镜像下载两阶段 `gtdb-r226-c200-dbv2.syl2db`，
  并用 GCS 公布的 MD5 校验。

两者都是一次性大体积下载。首次使用可加 `gside scan --meta --auto-setup` 自动安装
缺失档位（进度走 stderr）并重试该方法一次。`--source` 仍可从本地已有
数据库目录复制。

### 校验策略

`gside db setup` 对所有能固定的产物做摘要校验：

| 档位 | 摘要 | 说明 |
|---|---|---|
| panel | SHA256（`PANEL_SHA256`） | GitHub Release 资产 |
| mash | MD5（`MASH_MD5`） | 上游 Zenodo 仅提供 MD5 |
| sourmash | lineages SHA256 + 签名 SHA256（`SOURMASH_SIG_SHA256`） | 签名值来自已验证下载 |

sourmash 签名 zip（约 3.9GB）上游没有摘要。gside 记录了一次本地人工校验后的值
（体积精确、zip 完整、143,385 条 k=31 DNA 签名、CC0）：

```
SOURMASH_SIG_SHA256 = 6f3cc74f0bdcf6f84466ba91f87d6f2bf30081d5814bdf006cafd1e9cbbeada9
```

若使用不同镜像，可在安装时覆盖：

```bash
export GSIDE_SOURMASH_SIG_SHA256=<你的 zip 的 sha256>
gside db setup --tier sourmash
```

### 从本地数据库目录复制

若磁盘上已有现成数据库（例如另一台机器或旧 checkout），可跳过下载直接复制：

```bash
gside db setup --tier all --source /path/to/existing/data/db
```

这两个 tier 同样支持 `--source`，可直接从本地已有数据库目录复制。

### 查看可用 tier

```bash
gside db list
```

## 基因组质量数据库

`gside validate` 需要两个**不属于** `db setup` tier 的外部数据库：CheckM2
（完整度/污染度）与 GTDB-Tk（分类）。它们体积大、与程序版本强耦合，且
GTDB-Tk 内存需求高，因此 gside 刻意交给官方安装器处理。

| 数据库 | 下载 | 解压后 | 运行需求 | 环境变量 |
|---|---|---|---|---|
| CheckM2（v1.1.0） | 约 1.7GB | 约 3GB | 较低（DIAMOND） | `CHECKM2DB` → `.../uniref100.KO.1.dmnd` |
| GTDB-Tk R232 | 约 57GB | 约 98GB | 细菌约 140GB 内存 | `GTDBTK_DATA_PATH` → 数据库目录 |

### CheckM2

```bash
pixi add checkm2          # 或：conda install -c bioconda checkm2
checkm2 database --download --path /path/to/D7_checkm2
export CHECKM2DB=$(find /path/to/D7_checkm2 -name uniref100.KO.1.dmnd | head -1)
```

手动方式（无需 CheckM2 二进制）：

```bash
mkdir -p /path/to/D7_checkm2
wget -c https://zenodo.org/records/14897628/files/checkm2_database.tar.gz
md5sum checkm2_database.tar.gz    # 07c10655620843b517d0df0c160d911f
tar xzf checkm2_database.tar.gz -C /path/to/D7_checkm2
export CHECKM2DB=$(find /path/to/D7_checkm2 -name uniref100.KO.1.dmnd | head -1)
```

数据库必须与安装的 CheckM2 版本兼容；优先使用
`checkm2 database --download`，它会选择匹配的版本。

### GTDB-Tk

```bash
pixi add gtdbtk
download-db.sh /path/to/D8_gtdbtk    # 官方脚本；目标目录必须为空
export GTDBTK_DATA_PATH=/path/to/D8_gtdbtk
```

手动方式：

```bash
mkdir -p /path/to/D8_gtdbtk
wget -c https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
md5sum gtdbtk_data.tar.gz            # R232: 25a59e0352b1fd150c589f56559767d4
tar xvzf gtdbtk_data.tar.gz -C /path/to/D8_gtdbtk --strip 1
export GTDBTK_DATA_PATH=/path/to/D8_gtdbtk
```

| 版本 | 所需 GTDB-Tk | MD5 |
|---|---|---|
| R232 | ≥ 2.7.0 | `25a59e0352b1fd150c589f56559767d4` |
| R226 | 2.4.1 – 2.6.1 | `24b476ea5a4ef30519d461e56cc4a27f` |

硬件：细菌约需 140GB 内存（古菌约 100GB），R232 磁盘约 98GB。内存不足时
`validate` 会把 GTDB-Tk 记为 `missing`，仍写出 CheckM2 结果，并报告
`status: partial`（退出码 1）。

### gside 的报告

`gside db status` 只覆盖 `db setup` 的 tier。`gside validate` 读取上述两个
环境变量，并在 `validation.json` 中记录逐工具状态（`tools.checkm2`、
`tools.gtdbtk`、`errors`）。只有两工具都跑通才退出 0。

## 数据库路径

默认路径可通过环境变量覆盖：

| 变量 | 默认值 | 说明 |
|---|---|---|
| `GSIDE_DATA_DIR` | `<repo>/data` | 标记规则等小型参考数据 |
| `GSIDE_DB_DIR` | `<data>/db` | 大型数据库（panel/mash） |
| `GSIDE_PIXI_BIN` | 自动发现 | 生信二进制目录 |

## 自定义标记规则

编辑 `data/db/D1_marker/marker_rules.yaml`：

```yaml
- species: Mycobacterium_tuberculosis
  genes:
  - mpb64
  - is6110
  min_hits: 2
  min_identity: 90
  exclude_genes: [some_cross_reactive_gene]
```

然后重建 BLAST 库：

```bash
cd data/db/D1_marker
makeblastdb -in markers.fasta -dbtype nucl -out markers_blastdb
```

无需修改代码：规则即数据。
