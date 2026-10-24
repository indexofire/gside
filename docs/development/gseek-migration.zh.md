# gseek 引擎迁移

`gside` 可以通过显式开关，把 marker BLAST 与 FASTQ marker screening 两条路径
切到 [`gseek`](https://github.com/indexofire/gseek)。默认仍走原有子进程实现，
gseek 是可选依赖。

## 开关

```bash
GSIDE_ENGINE=legacy   # 默认：旧后端，不 import gseek
GSIDE_ENGINE=gseek    # 使用 gseek 的 BLAST/minimap2 适配器
GSIDE_ENGINE=shadow   # 两条路径都跑，记录差异，返回 legacy 结果
```

`GSIDE_ENGINE` 每次调用时读取，因此同一进程可在测试中切换模式；非法值会在
任何外部工具启动前报错。

本地安装 gseek：

```bash
pip install -e ../gseek
# 或者不安装，直接注入源码：
PYTHONPATH=../gseek/src pixi run -e dev pytest tests/integration/test_gseek_parity.py -q
```


> mash/sourmash 路径需要 gseek 0.2.0（`external_reference`），已发布 PyPI：
> `pip install gseek==0.2.0`。optional extra 已固定为对应版本；sibling checkout
> 仍可用于本地开发。

## 已迁移路径

| 路径 | legacy | gseek |
| --- | --- | --- |
| Marker BLAST | `BlastBackend` + `markers_blastdb` | `gseek.align(backend="blast", target=markers.fasta)` |
| FASTQ marker screening | minimap2 PAF + 本地聚合 | `gseek.align(backend="minimap2", ...)` + 同样的聚合 |
| mash_refseq（组装/reads） | `MashBackend` 对 `mash.msh` | `gseek.distance(backend="mash", external_reference=True)` |
| sourmash gather | sketch/gather CSV 解析 | `gseek.distance(backend="sourmash", external_reference=True)`，每个 query 文件一次 |

gseek 路径显式保留 gside 语义：

- coverage 仍为 `aligned_length / length`，不是 gseek 的 span coverage；
- screening identity 仍为 `nmatch / alen`，通过 `record.matches` 计算；
- BLAST 使用 `task="megablast"`、`word_size=11`、保留 dust，匹配历史
  `blastn` 默认行为；
- 工具发现仍使用 `GSIDE_PIXI_BIN`/pixi/PATH 的顺序。

## 尚未迁移

- `panel` 仍使用 skani database 目录（`panel.sketch`）。gseek 的
  external-reference 契约接受预构建文件，不接受数据库目录。
- sourmash 的 **rocksdb 目录**（`gtdb-k31.rocksdb`，基准脚本用于加速）不能作为
  gseek external reference。请使用 GTDB zip/signature 文件，或对该覆盖路径使用
  `GSIDE_ENGINE=legacy`/`shadow`。gseek 模式现在会显式报错，而不是静默返回空
  partition。

## 验证

- `tests/unit/test_gseek_bridge_units.py`：记录映射、聚合、模式分派、fake gseek
  参数接线。
- `tests/integration/test_gseek_parity.py`：真实 marker DB/基因组对比 legacy，
  以及 FASTQ screening parity；gseek 不可导入时 skip。

迁移分支必须通过 gside 常规门禁（`ruff`、`pyright`、`pytest tests/unit --cov`）。
parity 测试是额外集成检查，在 gseek 成为 CI 可选依赖前不进入默认 CI。


## 本地验证（2026-10-24）

对 tracked N16961 基因组和由 `inva` marker 生成的模拟 FASTQ 运行
`gside scan -m marker -f json`：

| 输入 | legacy | gseek | shadow |
| --- | --- | --- | --- |
| `GCF_000006745.1` 组装 | `Vibrio_cholerae` / high / `ompw` 94.767 / 100.1 | 完全一致 | 完全一致 |
| 覆盖 `inva` 的 21 条 reads | `Salmonella_enterica` / high / `inva` 100.0 / 98.3 | 完全一致 | 完全一致 |
| K12 组装，`-m mash_refseq` | `Escherichia_coli` / high / identity 1.0 | 完全一致 | 完全一致 |

`shadow` 没有产生差异 warning。单测 737 passed；集成 22 passed / 1 skipped；
gseek parity 5 passed；覆盖率 97.19%。



## 真实数据抽查（`~/temps/test`，2026-10-24）

每个物种取一个 RefSeq 基因组：*Acinetobacter baumannii*、*Campylobacter
coli*、*Escherichia coli*、*Salmonella enterica*、*Shigella flexneri*。分别用
`gside scan <genome> -m marker+mash_refseq -f json` 跑 legacy/gseek/shadow。

- 五个基因组的三种模式 JSON **逐字节一致**。
- 例如 `A. baumannii` marker `rpob-aba` 99.068 / `oxa51` 97.39、mash
  0.978914；`S. flexneri` marker `ipah` 100.0、mash 0.993922。
- 从 *A. baumannii* 生成的 260,000 条 150bp reads（约 10x）三种模式结果一致：
  marker `Acinetobacter_baumannii` high（`rpob-aba` 94.84 / `oxa51` 93.56），
  mash 0.9789136 high。
- 用三个 RefSeq 基因组构造的 mini sourmash signature DB，在 *A. baumannii*
  与 *S. enterica* 上三种模式的 `species`、`confidence`、
  `f_unique_weighted`、`ani`、`gather_partition` 完全一致。
- benchmark 的 **rocksdb 目录**（`gtdb-k31.rocksdb`）不能作为 gseek external
  reference：legacy 仍返回完整 partition；gseek 显式报错；shadow 返回 legacy
  并记录该错误。
