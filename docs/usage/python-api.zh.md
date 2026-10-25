# Python API

`gside.api` 是受支持的库入口。它按引擎提供独立鉴定函数，不包含 CLI 的数据库安装与分层仲裁。

```python
from gside.api import identify_multigene
```

API 处于 1.0 之前：可能新增可选参数，但现有位置参数与结果形状保持稳定。

## 依赖

每个函数需要与对应 CLI 模式相同的外部工具和参考数据：

| 函数 | 外部工具 | 数据 |
| --- | --- | --- |
| `identify_multigene` | `blastn`, `makeblastdb` | 随包 `D1_marker` |
| `identify_by_ani(mode="panel")` | `skani` | panel 层 |
| `identify_by_ani(mode="mash_refseq")` | `mash` | mash 层 |
| `identify_by_sourmash` | `sourmash` | sourmash 层 |
| `identify_by_gtdbtk` | `gtdbtk` | `GTDBTK_DATA_PATH` |
| `identify_by_kraken2` | `kraken2`, `bracken` | meta 层 |
| `identify_by_sylph` | `sylph` | meta_sylph 层 |
| `screen_reads` | `minimap2` | marker FASTA |
| `validate_genome` | `checkm2`, `gtdbtk` | 外部数据库 |

用 `GSIDE_DB_DIR` 或 `db_dir=` 指定非默认数据库根目录。工具从 `PATH` 和
`GSIDE_PIXI_BIN` 解析。

## Marker 基因

```python
from gside.api import identify_multigene

result = identify_multigene("contigs.fna", threads=4)
# 自定义 marker DB：
# result = identify_multigene("contigs.fna", marker_db="my-markers")

print(result.species)
print(result.confidence)          # high / medium / low
print(result.detected_markers)    # list[dict]
print(result.to_dict())           # 可直接 JSON 序列化
```

`identify_multigene` 接受单个路径或 contig 文件列表，应用随包 marker 组合规则，返回 `MultiGeneResult`。

## ANI：panel 与 RefSeq mash

```python
from gside.api import identify_by_ani

panel = identify_by_ani("contigs.fna", mode="panel", threads=8)
mash = identify_by_ani("contigs.fna", mode="mash_refseq", threads=8)

print(panel.to_dict()["result"]["species"])
print(panel.to_dict()["result"]["top_hits"][:3])
```

`mode` 必须是 `"panel"` 或 `"mash_refseq"`。返回 `AniIdResult`；`result` 含
`species`、`confidence` 与 `top_hits`（ANI/identity、aligned fraction）。

## Sourmash 与 GTDB-Tk

```python
from gside.api import identify_by_gtdbtk, identify_by_sourmash

sourmash = identify_by_sourmash("contigs.fna").to_dict()
gtdbtk = identify_by_gtdbtk("contigs.fna", threads=8).to_dict()
```

`identify_by_sourmash` 执行 GTDB gather，返回 `SourmashIdResult`。
`identify_by_gtdbtk` 返回 `GtdbkResult`，需要 `GTDBTK_DATA_PATH`。

## Reads 与宏基因组

### Marker screening

```python
from pathlib import Path

from gside.api import screen_reads
from gside.config import SPECIES_DB_DIR

markers = SPECIES_DB_DIR / "D1_marker" / "markers.fasta"
rows = screen_reads(("sample_R1.fq", "sample_R2.fq"), markers, threads=8)
```

Paired mates 在一次调用中映射。结果是按基因聚合的行列表，不是 `MethodEnvelope`。

### Kraken2/Bracken 与 sylph

```python
from gside.api import identify_by_kraken2, identify_by_sylph

kraken = identify_by_kraken2(("sample_R1.fq", "sample_R2.fq"), threads=8)
sylph = identify_by_sylph(("sample_R1.fq", "sample_R2.fq"), threads=8)

print(kraken.to_dict()["result"]["top_species"])
```

Kraken2 接受 FASTA 或 FASTQ；sylph 要求 FASTQ。二者返回 `MethodEnvelope` 子类。

## 组装质量

```python
from gside.api import validate_genome

result = validate_genome("contigs.fna", "qc/", threads=8)
print(result.status)
print(result.to_dict())
```

`validate_genome` 运行 CheckM2 与 GTDB-Tk，写 `qc/validation.json`，返回 `ValidateResult`。

## 结果形状

大多数引擎返回 `MethodEnvelope` 子类：

```python
{
    "analysis_type": "species_identification",
    "method": "panel",
    "database": {"name": "D2_ani", "version": "..."},
    "result": {...},
}
```

- `identify_multigene` 返回 `MultiGeneResult`，用 `.to_dict()`。
- `screen_reads` 返回 `list[dict]`。
- `validate_genome` 返回 `ValidateResult`，用 `.to_dict()`。

CLI 层字段与错误语义见 [JSON 输出契约](../reference/json-contract.zh.md)。

## 错误与退出行为

库函数抛出 Python 异常（参数错误为 `ValueError`，工具/数据库缺失为 `RuntimeError`）。CLI 会转成结构化 JSON 错误和退出码，不隐藏逐方法失败。

## 并发

每个函数处理一个样本。请在应用层按样本并行，并用 `threads=` 控制单工具线程数。不要在工作线程间共享可变结果对象。

## 引擎选择

已迁移路径可用 `GSIDE_ENGINE` 切换内部实现：

```bash
GSIDE_ENGINE=legacy   # 默认
GSIDE_ENGINE=gseek    # 在已实现路径使用 gseek
GSIDE_ENGINE=shadow   # 同时运行、比较并返回 legacy
```

范围与验证见 [gseek 迁移](../development/gseek-migration.zh.md)。
