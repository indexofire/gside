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


> mash/sourmash 路径依赖 gseek 的 `external_reference` 支持，该能力不在 gseek
> 0.1.0 中。请使用 sibling checkout（`PYTHONPATH=../gseek/src`）或包含该能力的
> gseek 构建；optional extra 已固定为 `gseek>=0.2`。

## 已迁移路径

| 路径 | legacy | gseek |
| --- | --- | --- |
| Marker BLAST | `BlastBackend` + `markers_blastdb` | `gseek.align(backend="blast", target=markers.fasta)` |
| FASTQ marker screening | minimap2 PAF + 本地聚合 | `gseek.align(backend="minimap2", ...)` + 同样的聚合 |
| mash_refseq（组装） | `MashBackend` 对 `mash.msh` | `gseek.distance(backend="mash", external_reference=True)` |
| sourmash gather（单 FASTA） | sketch/gather CSV 解析 | `gseek.distance(backend="sourmash", external_reference=True)` |

gseek 路径显式保留 gside 语义：

- coverage 仍为 `aligned_length / length`，不是 gseek 的 span coverage；
- screening identity 仍为 `nmatch / alen`，通过 `record.matches` 计算；
- BLAST 使用 `task="megablast"`、`word_size=11`、保留 dust，匹配历史
  `blastn` 默认行为；
- 工具发现仍使用 `GSIDE_PIXI_BIN`/pixi/PATH 的顺序。

## 尚未迁移

- `panel` 仍使用 skani database 目录（`panel.sketch`）。gseek 的
  external-reference 契约只接受单个文件，不接受数据库目录。
- `mash_refseq` 的 FASTQ reads 仍走 legacy `mash sketch -m 2 -r`；只有组装
  FASTA query 走 gseek。
- `sourmash` 多 query 文件仍走 legacy，因为 gseek sourmash 适配器每次只接受
  一个 query 文件。

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

`shadow` 没有产生差异 warning。单测 732 passed；集成 22 passed / 1 skipped；
gseek parity 4 passed；覆盖率 97.07%。

