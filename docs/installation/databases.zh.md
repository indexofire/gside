# 数据库管理

gside 的参考数据库通过 `gside db` 子命令管理，支持分层数据库安装。

## 数据库分层

| Tier | 包含 | 大小 | 用途 |
|---|---|---|---|
| **mini** | marker_rules.yaml + markers.fasta + 预建 BLAST 库 | ~1MB | 靶基因组合鉴定（L1），随仓库分发 |
| **panel** | D2_ani（291 基因组 skani sketch） | ~130MB | ANI 鉴定（L2） |
| **mash** | D3_mash（RefSeq MinHash sketch） | ~179MB | 距离鉴定（L2） |
| **sourmash** | D4_sourmash（GTDB reps k=31 + lineages） | ~3.9GB | gather 分类（sourmash） |
| **all** | panel + mash | ~310MB | 全层次能力（不含 sourmash） |

### 目录命名

`data/db/` 下目录按鉴定层次命名：

| 旧名 | 新名 | 对应方法 |
|---|---|---|
| `refseq_panel` | `D2_ani` | panel（skani ANI） |
| `mash_refseq` | `D3_mash` | mash_refseq（MinHash） |
| `sourmash_gtdb` | `D4_sourmash` | sourmash（GTDB gather，需手动下载） |

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
  ✅ markers          tier=mini   path=data/db/D1_marker
  ✅ D2_ani        tier=panel  path=data/db/D2_ani
  ✅ D3_mash         tier=mash   path=data/db/D3_mash
  ✅ D4_sourmash     tier=sourmash path=data/db/D4_sourmash
```

### 安装/更新

```bash
gside db setup --tier panel          # 仅面板（ANI），默认
gside db setup --tier mash           # 仅 MinHash
gside db setup --tier sourmash       # GTDB gather 库（~3.9GB）
gside db setup --tier all            # panel + mash（不含 sourmash）
```

不带 `--tier` 时默认安装 `panel`。

`--tier all` 仅安装 panel + mash。sourmash GTDB 库（约 3.9GB）不含在 `all`
内，需另行执行 `gside db setup --tier sourmash` 单独安装。

### 查看可用 tier

```bash
gside db list
```

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
