# 数据库管理

gside 的参考数据库通过 `gside db` 子命令管理，支持分层数据库安装。

## 数据库分层

| Tier | 包含 | 大小 | 用途 |
|---|---|---|---|
| **mini** | marker_rules.yaml + markers_v2.fasta + 预建 BLAST 库 | ~1MB | 靶基因组合鉴定（L1），随仓库分发 |
| **panel** | refseq_panel（291 基因组 skani sketch） | ~130MB | ANI 鉴定（L2） |
| **mash** | mash_refseq（RefSeq MinHash sketch） | ~179MB | 距离鉴定（L2） |
| **all** | panel + mash | ~310MB | 全层次能力 |

## 命令

### 查看状态

```bash
gside db status
```

```
gside database status
──────────────────────────────────────────
  ✅ markers_v2       tier=mini   path=data/reference/species
  ✅ refseq_panel     tier=panel  path=data/db/refseq_panel
  ✅ mash_refseq      tier=mash   path=data/db/mash_refseq
```

### 安装/更新

```bash
gside db setup --tier panel          # 仅面板（ANI），默认
gside db setup --tier mash           # 仅 MinHash
gside db setup --tier all            # 全部
```

不带 `--tier` 时默认安装 `panel`。

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

编辑 `data/reference/species/marker_rules.yaml`：

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
cd data/reference/species
makeblastdb -in markers_v2.fasta -dbtype nucl -out markers_v2_blastdb
```

无需修改代码：规则即数据。
