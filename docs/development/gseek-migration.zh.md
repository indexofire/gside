# gseek 引擎迁移

`gside` 可以通过显式开关，把 marker BLAST/FASTQ screening、mash distance、
sourmash gather 与 panel/skani 四条路径切到
[`gseek`](https://github.com/indexofire/gseek)。默认仍走原有子进程实现，
gseek 是可选依赖。

## 开关

```bash
GSIDE_ENGINE=legacy   # 默认：旧后端，不 import gseek
GSIDE_ENGINE=gseek    # 所有已迁移路径走 gseek
GSIDE_ENGINE=shadow   # 两条路径都跑，记录差异，返回 legacy 结果
```

`GSIDE_ENGINE` 每次调用时读取，因此同一进程可在测试中切换模式；非法值会在
任何外部工具启动前报错。选择、shadow 比较与回退由 `gseek.compat.Engine`
提供；legacy 模式仍然完全不 import gseek。

## legacy 保留策略

这是有意的双路径设计，不是迁移遗留：

- `legacy` 是默认值，没有 gseek 的安装行为与之前完全一致；
- 每个已迁移路径都保留 legacy 实现，仍可正常使用；
- `shadow` 同时跑两套实现，记录粗粒度差异，并返回 legacy 结果；
- legacy 是 parity 测试的基准，也是 gseek 不可用时的回退；
- 没有 gseek 适配器的模式（`gtdbtk`、`kraken2`、`sylph`）只能走 legacy。

本地安装 gseek：

```bash
pixi install -e dev            # dev feature 包含 gseek>=0.9.1
# 或：pip install -e ".[dev]"
# 或者不安装 release，直接注入源码：
PYTHONPATH=../gseek/src pixi run -e dev pytest tests/integration/test_gseek_parity.py -q
```


> mash/sourmash 路径需要 gseek 0.9.1（`external_reference`），已发布 PyPI：
> `pip install gseek==0.9.1`。optional extra 已固定为对应版本。开发/测试环境
> 通过 `dev` extra（`pip install -e ".[dev]"`）和 pixi `dev` feature 安装；
> sibling checkout 仅用于修改 gseek 本身。

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

所有由 gseek 后端覆盖的 gside 模式均已迁移（marker、mash、sourmash、panel）。
没有 gseek 适配器的模式（`gtdbtk`、`kraken2`、`sylph`）保持原有实现。

## 验证

- `tests/unit/test_gseek_bridge_units.py`：记录映射、聚合、模式分派、fake gseek
  参数接线。
- `tests/integration/test_gseek_parity.py`：真实 marker DB/基因组对比 legacy，
  以及 FASTQ screening parity；gseek 不可导入时 skip。

`main` 上执行同样的门禁（`ruff`、`pyright`、`pytest tests/unit --cov` 及集成）。
安装 gseek 时，parity 测试会被默认测试套件收集；需要外部二进制的用例在缺失时 skip。


## 本地验证（2026-10-24）

对 tracked N16961 基因组和由 `inva` marker 生成的模拟 FASTQ 运行
`gside scan -m marker -f json`：

| 输入 | legacy | gseek | shadow |
| --- | --- | --- | --- |
| `GCF_000006745.1` 组装 | `Vibrio_cholerae` / high / `ompw` 94.767 / 100.1 | 完全一致 | 完全一致 |
| 覆盖 `inva` 的 21 条 reads | `Salmonella_enterica` / high / `inva` 100.0 / 98.3 | 完全一致 | 完全一致 |
| K12 组装，`-m mash_refseq` | `Escherichia_coli` / high / identity 1.0 | 完全一致 | 完全一致 |

`shadow` 没有产生差异 warning。单测 748 passed；集成 28 passed / 1 skipped（安装 gseek 时；
其中 6 个为 parity，未装 gseek 时 22 个）；覆盖率 96.88%。



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
- panel（skani）在 *A. baumannii* 与 *C. coli* 上：legacy/gseek/shadow JSON
  逐字节一致（A. baumannii 98.06 ANI / 85.79 AF；C. coli 99.72 / 98.09）。
- 用三个 RefSeq 基因组构造的 mini sourmash signature DB，在 *A. baumannii*
  与 *S. enterica* 上三种模式的 `species`、`confidence`、
  `f_unique_weighted`、`ani`、`gather_partition` 完全一致。
- 完整 gside GTDB zip（3.7 GB）对 *A. baumannii*：legacy 与 gseek JSON 逐字节
  一致（`species=Acinetobacter_baumannii`，high，`f_unique_weighted=
  0.5978601252609603`，`ani=98.3`，相同的两成员 partition）。gseek 显式传
  `threshold_bp=50000` 对齐 legacy sourmash CLI 默认值；最初的 50 bp 默认会多出
  低权重 partition。
- benchmark 的 **rocksdb 目录**（`gtdb-k31.rocksdb`）已可作为 gseek 0.9.1 的
  external reference，并按确定性 tree hash 记录。legacy/gseek/shadow 均返回相同
  的完整 partition（`Acinetobacter_baumannii`，high，0.5978601252609603，
  ANI 98.3），shadow 无差异 warning。


## 全量数据基准（`~/temps/test`）

实际数据是 **42 个 label、4,178 个基因组**（包含 `Shigella` 属级 label），
每个基因组都跑 legacy 与 gseek 的四种模式，8 个 worker 并行：

| 模式 | 运行数 | 完全一致 | 错误 | legacy 中位数 | gseek 中位数 | 性能比（中位数） |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| marker | 4,178 | 4,178 | 0 | 0.4381s | 0.4480s | 0.978x |
| panel | 4,178 | 4,178 | 0 | 0.0660s | 0.0748s | 0.882x |
| mash_refseq | 4,178 | 4,178 | 0 | 0.5183s | 0.7013s | 0.739x |
| sourmash | 4,178 | 4,178 | 0 | 1.6415s | 1.6452s | 0.998x |

总耗时：legacy 11,309s，gseek 12,183s（并行基准下 gseek 约慢 1.08 倍）。
0.6.0 复用 ToolResolver probe、按进程缓存 reference digest、按需跳过
manifest/result hash，并使用 mash 直接 FASTA dist。该基准是 8 路并行下的
wall-clock，不是单核隔离对比。

原始行与生成报告：`~/temps/test/results/gseek-full-benchmark-v4/`。
