# 物种鉴定模式

## 分层架构

gside 的鉴定按证据强度分层，高层覆盖低层（记录不致命）：

```
L3  GTDB-Tk（标准模式，需重型环境）          [层权威=3]
L2  panel skani / mash_refseq / sourmash    [层权威=2]
L1  marker 靶基因组合规则                    [层权威=1]
```

## L1: marker（靶基因组合）

### 原理

对 contigs 执行一次 BLAST（vs markers_v2 库，81 条序列），将命中基因与
38 条组合规则匹配（数据驱动，`marker_rules.yaml`）。

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

| ANI | 覆盖率 | 判定 |
|---|---|---|
| ≥95% | ≥70% | species=命中参考基因组的物种，confidence=high |
| 90-95% | ≥70% | confidence=medium（建议复核） |
| <90% | — | Unknown |

### 适用场景

- 精确物种确认（属级分辨率）
- 近缘种区分（marker 层无法判别时）
- panel tier 数据库（2.7GB）

## L2: mash_refseq（MinHash 距离）

### 原理

mash 对 contigs 与 RefSeq 全库 sketch（159MB）做 MinHash 距离估算。

### 适用场景

- 广谱筛查（覆盖 RefSeq 全部原核生物）
- 快速（秒级）
- 精度低于 panel ANI（MinHash 是估算）

## L2: sourmash（GTDB gather）

### 原理

sourmash gather 对 contigs 做 LCA 分类，基于 GTDB 分类框架。

### 适用场景

- 需要与 GTDB 分类体系对齐时
- 混合样本（gather 可分解）
- 需 3.7GB sourmash GTDB 库

## all 模式（多法仲裁）

运行所有可用方法，按层级仲裁取最优：

```
ANI 层 (panel/mash/sourmash) high-confidence 命中 > marker 层命中
```

```bash
gside species contigs.fna --mode all
```

输出包含所有方法结果 + 最终 verdict：

```json
{
  "verdict": {
    "species": "Campylobacter jejuni subsp. jejuni NCTC 11168",
    "confidence": "high",
    "basis": ["panel"]
  },
  "methods": {
    "marker": { "species": "Campylobacter_jejuni", ... },
    "panel": { "result": { "species": "...", "ani": 100.0, ... } },
    "mash_refseq": { "result": { "species": "...", ... } }
  }
}
```

## 方法选择建议

| 场景 | 推荐 mode | 理由 |
|---|---|---|
| 快速初筛 | `marker` | 秒级，零外部依赖 |
| 精确鉴定 | `panel` | ANI 金标准 |
| 广谱筛查 | `mash_refseq` | 覆盖全 RefSeq |
| 综合判定 | `all` | 多法互补 + 仲裁 |
| 近缘种区分 | `panel` 或 `all` | marker 可能交叉反应 |
