# 标记规则格式

`marker_rules.yaml` 定义靶基因组合鉴定规则。新增物种只需添加规则条目，
无需修改代码。

## 文件位置

```
data/db/D1_marker/marker_rules.yaml
```

## 规则结构

```yaml
rules:
  - species: <物种名>              # 必填：注册表键名（下划线连接）
    genes:                          # 必填：候选标记基因列表
      - <gene1>
      - <gene2>
    min_hits: <int>                 # 可选：最少命中数（默认 1）
    min_identity: <float>           # 可选：最低一致性 %（默认 90）
    min_coverage: <float>           # 可选：最低覆盖率 %（默认 60）
    exclude_genes:                  # 可选：排除守卫
      - <cross_reactive_gene>
    note: <str>                     # 可选：注释（文档性）
    co_detect_note: <str>           # 可选：共检出说明（文档性）
```

## 字段说明

### species

物种名，使用注册表键名（空格转下划线）。如 `Campylobacter_jejuni`、
`Vibrio_parahaemolyticus`。

### genes

候选标记基因名列表（小写）。这些基因必须存在于 `markers.fasta` 中。

### min_hits

最少命中数。设为 2 时，`genes` 中至少 2 个基因需达到阈值才算该物种。

多基因组合可显著降低交叉反应（如 C.jejuni 的 3 选 2 设计）。

### min_identity

每个命中基因的最低一致性百分比（blastn pident）。默认 90。

降低至 85 适用于短序列或保守基因（但会增加假阳性风险）。

### exclude_genes

**排除守卫**：如果这些基因在样本中检测到，则该物种的判定被抑制。

典型用法，弯曲菌种间区分：

```yaml
- species: Campylobacter_coli
  genes: [ceue]
  exclude_genes: [mapa, hipo, cadf]    # jejuni 标记存在 → 倾向 jejuni
- species: Campylobacter_jejuni
  genes: [mapa, hipo, cadf]
  min_hits: 2
  # 无排除：jejuni 基因组天然含 ceuE 同源体
```

### min_coverage

基因计入本规则所需的最低 query 覆盖率（默认 60，与全局门槛一致）。
判定时逐条规则检查该基因的所有原始命中，因此"一致性稍低但覆盖率更高"
的命中仍可满足规则。

仅当规则确实接受部分基因匹配时才调低（例如 50）——两条
*Campylobacter* 规则就是这样配置的。

### note 与 co_detect_note

文档性字段。记录规则设计意图与近缘种共检出提示，仅供规则维护者参考，
不参与判定逻辑。

## 全局阈值

引擎层设有全局门槛（在 `multigene_identifier.py` 中定义）：

| 参数 | 值 | 说明 |
|---|---|---|
| `_MIN_IDENTITY` | 85.0 | 全局最低一致性 |
| `_MIN_COVERAGE` | 60.0 | 全局最低覆盖率（防短 HSP 误报） |
| `_HIGH_CONF` | 90.0 | 高置信度分界 |

规则级 `min_identity` / `min_coverage` 作用于该基因的每条原始命中；
未填写的规则使用 90 / 60。全局 `_MIN_IDENTITY`（85）与 `_MIN_COVERAGE`（60）
仍用于 `detected_markers` 上报与 `exclude_genes` 排除守卫。

`priority_over` 引擎支持（命中的 priority 物种会压制非 priority 物种）。四条志贺菌
种级规则（`Shigella_flexneri`、`Shigella_sonnei`、`Shigella_dysenteriae`、
`Shigella_boydii`）使用 `priority_over: [Shigella, Escherichia_coli]`：
种级决定因子（如 `rfc`、`wbgZ`、`sb1_wzx`）会升级属级判定，而单纯的
`uidA+ipaH` 仍停留在 `Shigella` 属级；所有志贺菌规则均为 lacY 阴性
（`exclude_genes: [lacy]`），EIEC 不会触发志贺菌判定。

## Reads 模式深度门控（FASTQ 输入） {: #reads-mode-depth-gate }

FASTQ 输入时，标记引擎在一致性/覆盖规则之上叠加 gspec 风格的丰度门控。
每个基因上报 `mean_depth`（比对目标区间长度之和 / 基因长度 —— 不是读段数）、
`reads` 和 `mapq_mean`；低于 `min_mapq`（默认 20）的比对在聚合前丢弃。

标记基因需同时满足：

```
mean_depth >= min_depth  且  mean_depth / control_depth >= min_relative_depth
```

默认值（在 `metadata.yaml` 中，调用时可用 `min_depth` / `min_relative_depth` /
`min_mapq` 覆盖）：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `min_depth` | 5.0 | 绝对平均深度下限 |
| `min_relative_depth` | 0.10 | 相对对照深度的比例 |
| `min_mapq` | 20 | 聚合前的 minimap2 MAPQ 下限 |

`control_depth` 是该物种配置的单拷贝对照基因（`metadata.yaml: control_genes`，
如 *Salmonella* 的 `rpoD`、*Campylobacter* 的两个 `lpxA` 等位基因）的
mean depth 中位数。每个基因取最深的等位基因，避免少量污染读段压低对照。
对照基因在 `markers.fasta` 中标记 `role=control`，绝不会出现在
`detected_markers` 中。没有专用对照的物种回退到已检出标记的样本中位深度
（`control_source: "fallback"`）；完全无对照时关闭相对门控，仅保留绝对下限。
组装（BLAST）输入没有深度信息，保持旧行为。

## 新增物种示例

```yaml
# 在 rules 列表中添加：
- species: Mycobacterium_tuberculosis
  genes:
  - mpb64
  - is6110
  - rpoB_mtbc
  min_hits: 2
  min_identity: 90
  exclude_genes: [rpob_maa]    # 排除 M.africanum
  note: MTBC complex identification
```

然后：

1. 将对应基因序列追加到 `markers.fasta`，FASTA 头格式：
   `>markers~~~mpb64~~~ACC 描述 [物种] role=primary`
2. 重建 BLAST 库：
   `makeblastdb -in markers.fasta -dbtype nucl -out markers_blastdb`
3. 运行 `gside scan test.fna -m marker` 验证
