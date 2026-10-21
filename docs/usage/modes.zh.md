# 物种鉴定模式

## 分层架构

gside 的鉴定按证据强度分层，高层优先于低层（记录不致命）：

```
L3  gtdbtk（GTDB-Tk 系统发育 placement）      [layer=3]
L2  panel skani / mash_refseq / sourmash    [layer=2]
L1  marker 靶基因组合规则                    [layer=1]
```

## L1: marker（靶基因组合）

### 原理

对 contigs 执行一次 BLAST（vs markers 库，81 条序列），将命中基因与
38 条组合规则匹配（数据驱动，`marker_rules.yaml`）。

### 判定阈值

| 门槛 | 值 | 说明 |
|---|---|---|
| 全局最低一致性 | 85% | 低于此值的命中不计入 |
| 全局最低覆盖率 | 60% | 防短 HSP 高一致性误报 |
| 高置信分界 | 平均一致性 90% | 命中基因平均一致性 ≥90% 时 confidence=high，低于则为 medium |

### 规则机制

```yaml
- species: Campylobacter_jejuni
  genes: [mapa, hipo, cadf]     # 候选标记
  min_hits: 2                     # 至少命中 2 个
  min_identity: 90                # 每个基因 ≥90% 一致性
  exclude_genes: [ceue]           # 排除守卫：ceuE 存在 → 倾向 C_coli
```

### 交叉反应防护

- **exclude_genes**：排除规则防止近缘种误判（如 C.jejuni vs C.coli）
- **tlh 近缘守卫**：单基因 <90% 命中 → 判定抑制为 Unknown（V.alginolyticus 案例）
- **coverage 门槛**：全局 ≥60%（防短 HSP 高一致性误报）

### 适用场景

- 快速初筛（分钟级）
- 无需大数据库（mini tier 内置）
- 近缘种丰富的属（沙门菌/大肠/弯曲菌等）需配合 L2 确认

## L2: panel（ANI 精选面板）

### 原理

skani 对 contigs 与精选参考面板（291 基因组，30 启用物种）做全基因组 ANI 比对。

### 判定阈值

| ANI | aligned_fraction | 判定 |
|---|---|---|
| ≥95% | ≥0.65 | species=命中参考基因组的物种，confidence=high |
| 93~95% | ≥0.65 | confidence=medium（建议复核） |
| <93% | 无 | Unknown |

### 适用场景

- 精确物种确认（亚种级分辨率）
- 近缘种区分（marker 层无法判别时）
- 需 panel tier 数据库（~130MB）

## L2: mash_refseq（MinHash 距离）

### 原理

mash 对 contigs 与 RefSeq 全库 sketch（mash.msh，~179MB）做 MinHash 距离估算。

### 判定阈值

| mash 相似度（1 - 距离） | 判定 |
|---|---|
| ≥0.97 | 物种判定，confidence=high |
| 0.90~0.97 | confidence=medium（建议复核） |
| <0.90 | Unknown |

### 适用场景

- 广谱筛查（覆盖 RefSeq 全部原核生物）
- 快速（秒级）
- 精度低于 panel ANI（MinHash 是估算）

## L2: sourmash（GTDB gather）

### 原理

sourmash gather 对 contigs 做 LCA 分类，基于 GTDB 分类框架。

### 判定阈值

置信度取自 gather 的 containment ANI（`average_containment_ani`），与 panel
层共用 95/93 分界：

| Containment ANI | 判定 |
|---|---|
| ≥95 | confidence=high（除非标记为混合） |
| 93~95 | confidence=medium |
| <93 | Unknown |

若 gather CSV 没有 ANI 列（旧版导出），回退到 `f_unique_weighted`：
≥0.90 → high，≥0.70 → medium。

当两个及以上**不同物种**各自 `f_unique_weighted` ≥ 0.10 且 ANI ≥ 93 时，
结果携带 `possible_mixture` 标记并降为 medium；同一物种的多个菌株不会触发。

结果同时给出 `ani` 以及逐匹配的 `gather_partition`
（`species`、`f_unique_weighted`、`ani`）。

### 适用场景

- 需要与 GTDB 分类体系对齐时
- 混合样本（gather 可分解）
- 需下载 sourmash GTDB 库

## L3: gtdbtk（GTDB-Tk 系统发育 placement）

### 原理

GTDB-Tk `classify_wf` 将组装放入 GTDB 参考树
（`classify_wf --genome_dir ... --out_dir ... --database $GTDBTK_DATA_PATH`）。
分类串从 `gtdbtk.bac120.summary.tsv` 解析（回退 `gtdbtk.ar122.summary.tsv`），
物种名取最后一个 `s__` 位阶。

### 判定阈值

| 分类串 | 判定 |
|---|---|
| 含 `s__` 位阶 | 提取物种名，confidence=high（权威） |
| 无 `s__` 位阶 | Unknown，confidence=low |

GTDB-Tk 视为最高权威层：`gtdbtk` 的判定优先于任何 L2 ANI/MinHash 判定，
与低层置信度无关。

### 依赖

- 外部安装的 GTDB-Tk（`gside db setup` 不负责安装）：PATH 中的 `gtdbtk`
  二进制 + 通过 `GTDBTK_DATA_PATH`（或 `GTDBDB`）指定的参考库——见
  https://ecogenomics.github.io/GTDBTk/
- 较慢：每基因组分钟级（MAFFT + pplacer）；可用 `-t/--threads` 提线程

### 适用场景

- 需要 GTDB 树上的权威物种判定
- 裁决 L1/L2 方法之间的分歧
- 论文级分类学指定

## all 模式（多法仲裁）

运行全部 5 种方法（marker + panel + mash_refseq + sourmash + gtdbtk），
按层级仲裁取最优：

```
L3（gtdbtk）优先于 L2（panel/mash_refseq/sourmash）优先于 L1（marker）：
先比较层级，再比较置信度。
```

单个方法失败不会终止整体流程，错误信息写入该方法结果
（`methods.<mode>.error`）。缺库错误额外携带结构化字段
（`missing_tier`、`suggested_command`、`estimated_size`）；加上
`--auto-setup` 后 gside 会自行安装该 tier 并对该方法重试一次。

`gtdbtk` 包含在 `all` 内但较慢（每基因组分钟级）——追求速度时请显式
指定 `-m` 跳过。

```bash
gside species contigs.fna --mode all
```

输出包含所有方法结果 + 最终 verdict：

```json
{
  "verdict": {
    "species": "Campylobacter_jejuni",
    "display_name": "Campylobacter jejuni",
    "strain": "subsp. jejuni NCTC 11168",
    "confidence": "high",
    "basis": ["panel"],
    "supporting": ["marker"],
    "agreement": "2/3"
  },
  "methods": {
    "marker": { "species": "Campylobacter_jejuni", ... },
    "panel": { "result": { "species": "Campylobacter_jejuni", "species_raw":
      "Campylobacter jejuni subsp. jejuni NCTC 11168", "ani": 100.0, ... } },
    "mash_refseq": { "result": { "species": "...", ... } },
    "sourmash": { "result": { "species": "...", ... } },
    "gtdbtk": { "result": { "species": "...", "classification": "...", ... } }
  }
}
```

`verdict.species` 一律为规范 `Genus_species` 标签（`verdict.display_name`
为空格可读形式，`verdict.strain` 为菌株）；`supporting` 列出其他一致的方法，
`agreement` 为 `一致数/成功数`，`basis` 仍是单个胜出方法。

## 方法选择建议

| 场景 | 推荐 mode | 理由 |
|---|---|---|
| 快速初筛 | `marker` | 最快，零外部依赖 |
| 精确鉴定 | `panel` | ANI 金标准 |
| 广谱筛查 | `mash_refseq` | 覆盖全 RefSeq |
| GTDB 体系对齐 | `sourmash` | LCA 分类，可提示混合样本 |
| 权威 placement | `gtdbtk` | GTDB 树系统发育 placement（慢） |
| 综合判定 | `all` | 多法互补 + 仲裁 |
| 近缘种区分 | `panel` 或 `all` | marker 可能交叉反应 |

## 宏基因组谱分析（`gside meta`）

混合群落样本可使用 `gside meta` 做分类学谱分析。详见
[宏基因组谱分析](meta.zh.md)。

| 模式 | 引擎 | 输入 | 置信度依据 |
|------|--------|-------|-----------------|
| `kraken2` | Kraken2 + Bracken | FASTQ / FASTA | Bracken 丰度（≥10% high，≥1% medium） |
| `sylph` | sylph | 仅 FASTQ | Adjusted ANI（≥95 high，≥93 medium） |
| `all` | 双引擎 | FASTQ / FASTA | 合并谱 |
