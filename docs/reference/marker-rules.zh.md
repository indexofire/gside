# 标记规则格式

`marker_rules.yaml` 定义靶基因组合鉴定规则。新增物种只需添加规则条目，
无需修改代码。

## 文件位置

```
data/db/L1_marker/marker_rules.yaml
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
    min_coverage: <float>           # 可选：最低覆盖率 %（当前被引擎忽略，生效全局 60）
    exclude_genes:                  # 可选：排除守卫
      - <cross_reactive_gene>
    note: <str>                     # 可选：注释（文档性）
    co_detect_note: <str>           # 可选：共检出说明（文档性）
```

## 字段说明

### species

物种名，使用注册表键名（空格转下划线）。如 `Campylobacter_jejuni`、
`V_parahaemolyticus`。

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

规则级最低覆盖率百分比。**当前引擎忽略该字段**：覆盖率门槛实际由全局
`_MIN_COVERAGE`（60）生效，在规则中设置 min_coverage 不会改变判定。

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

规则中的 `min_identity` 可以覆盖全局值（但不会低于 85）。
规则级 `min_coverage` 当前被忽略，实际生效全局 60。

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
3. 运行 `gside species test.fna --mode marker` 验证
