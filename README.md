# gside（Genome Species IDentification Engine）

独立自包含的基因组物种鉴定 CLI：blastn 靶基因组合规则打底，skani / mash / sourmash
三路比对补充，`all` 模式按层仲裁，最终输出单一 JSON 判定。

## 鉴定模式

| mode | 底层引擎 | 原理 |
|---|---|---|
| `marker`（默认） | blastn | 靶基因组合规则（38 条规则，81 条参考序列） |
| `panel` | skani | ANI 比对精选参考面板（291 基因组） |
| `mash_refseq` | mash | RefSeq MinHash 距离 |
| `sourmash` | sourmash | GTDB gather |
| `all` | 上述全部 | 顺序执行 4 种方法并按层仲裁 |

## 快速开始

```bash
pixi install && pip install -e .

# 物种鉴定（--mode 可选 marker/panel/mash_refseq/sourmash/all，默认 marker）
gside species contigs.fna --mode marker
gside species contigs.fna --mode panel
gside species contigs.fna --mode mash_refseq
gside species contigs.fna --mode sourmash
gside species contigs.fna --mode all

# 数据库管理（panel/mash/sourmash 需先下载对应库）
gside db status
gside db setup --tier panel

gside --version
```

## 判定阈值

- marker：命中需同时满足 identity >= 85% 且 coverage >= 60%（全局阈值，逐 hit 过滤）；
  规则字段 `min_hits`（默认 1）、`min_identity`（默认 90）；命中基因平均
  identity >= 90% 判 high，否则 medium。
- panel：ANI >= 95% 且 aligned_fraction >= 0.65 判 high；93% <= ANI < 95% 且
  aligned_fraction >= 0.65 判 medium。
- mash_refseq：identity >= 0.97 判 high；0.90 <= identity < 0.97 判 medium，更低为 low。

## 仲裁（all 模式）

按 marker、panel、mash_refseq、sourmash 顺序执行，随后仲裁：

- panel / mash_refseq / sourmash 为第 2 层，marker 为第 1 层，层高者优先
- 层优先于置信度，同层内 high 优先
- 单法失败不终止其余方法，错误写入该方法在 JSON 中的条目

## 数据目录

| 变量 | 默认 | 内容 |
|---|---|---|
| `GSIDE_DATA_DIR` | `<repo>/data` | 参考数据根目录（规则、fasta、预建库） |
| `GSIDE_DB_DIR` | `<data>/db` | 大型数据库（panel / mash sketch 等） |
| `GSIDE_PIXI_BIN` | 自动发现 | 生信二进制目录（blastn / skani / mash） |

数据分层：

- mini：随仓库分发，含 marker_rules.yaml（38 条规则）、markers.fasta（81 序列）
  与预建 BLAST 库，开箱即用
- panel：skani 精选面板 sketch，约 130MB
- mash：RefSeq MinHash sketch，约 179MB

## 自定义规则

编辑 `data/db/L1_marker/marker_rules.yaml`，字段为 `species` / `genes` /
`min_hits` / `min_identity` / `exclude_genes`。修改 fasta 后重建 BLAST 库即生效，
无需改代码：

```bash
makeblastdb -in markers.fasta -dbtype nucl -out markers_blastdb
```

## 输出

统一 JSON 输出（stdout）。顶层字段：`analysis_type` / `tool` / `version` /
`contigs` / `methods` / `verdict`。`methods` 内 panel、mash_refseq、sourmash 为
`result` 嵌套结构，marker 为扁平结构（无 `result` key）。`gside species` 恒以
退出码 0 结束，单法错误写入 JSON，不影响其余方法。

## 验证状态

- marker 法 17/17 金标准基因组（含 10 株近缘阴性对照）
- 3/3 新增物种真实参考株（C. jejuni / C. coli / S. pyogenes）
- panel 法 17/17（291 基因组精选面板）

## License

MIT
