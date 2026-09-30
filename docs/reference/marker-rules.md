# 标记规则格式

`marker_rules.yaml` 定义靶基因组合鉴定规则——新增物种只需添加规则条目，无需修改代码。

## 文件位置

```
data/reference/species/marker_rules.yaml
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
    min_coverage: <float>           # 可选：最低覆盖率 %（默认 50）
    exclude_genes:                  # 可选：排除守卫
      - <cross_reactive_gene>
    priority_over:                  # 可选：优先级声明
      - <other_species>
    note: <str>                     # 可选：注释
```

## 字段说明

### species

物种名，使用注册表键名（空格→下划线）。如 `Campylobacter_jejuni`、`V_parahaemolyticus`。

### genes

候选标记基因名列表（小写）。这些基因必须存在于 `markers_v2.fasta` 中。

### min_hits

最少命中数。设为 2 时，`genes` 中至少 2 个基因需达到阈值才算该物种。

多基因组合可显著降低交叉反应（如 C.jejuni 的 3 选 2 设计）。

### min_identity

每个命中基因的最低一致性百分比（blastn pident）。默认 90。

降低至 85 适用于短序列或保守基因（但会增加假阳性风险）。

### exclude_genes

**排除守卫**：如果这些基因在样本中检测到，则该物种的判定被抑制。

典型用法——弯曲菌种间区分：

```yaml
- species: Campylobacter_coli
  genes: [ceue]
  exclude_genes: [mapa, hipo, cadf]    # jejuni 标记存在 → 倾向 jejuni
- species: Campylobacter_jejuni
  genes: [mapa, hipo, cadf]
  min_hits: 2
  # 无排除：jejuni 基因组天然含 ceuE 同源体
```

### priority_over

声明本物种优先于其他物种（如 DEC 优先于 Shigella_EIEC）。

## 全局阈值

引擎层设有全局门槛（在 `multigene_identifier.py` 中定义）：

| 参数 | 值 | 说明 |
|---|---|---|
| `_MIN_IDENTITY` | 85.0 | 全局最低一致性 |
| `_MIN_COVERAGE` | 60.0 | 全局最低覆盖率（防短 HSP 误报） |
| `_HIGH_CONF` | 90.0 | 高置信度分界 |

规则中的 `min_identity` 可以覆盖全局值（但不会低于 85）。

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
1. 将对应基因序列追加到 `markers_v2.fasta`（`>markers_v2~~~mpb64~~~ACC 描述 [物种] role=primary`）
2. 重建 BLAST 库：`makeblastdb -in markers_v2.fasta -dbtype nucl -out markers_v2_blastdb`
3. 运行 `gside species test.fna --mode marker` 验证
