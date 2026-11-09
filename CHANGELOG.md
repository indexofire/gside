# Changelog

本项目所有显著变更记录于此。
格式基于 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，
版本遵循 [Semantic Versioning](https://semver.org/lang/zh-CN/)。

## [Unreleased]

## [0.2.3] - 2026-11-09

### Added

- **一个模式可对应多个数据库，鉴定时可自行选择**（`-d/--db [MODE:]NAME`，可重复）。

  - `marker` / `panel` / `mash_refseq` / `sourmash` 支持命名数据库，放在
    `data/db/<family>/<name>/`；取值可为家族命名库、已存在目录路径或保留名
    `default`；`--marker-db` 保留为 `-d marker:PATH` 兼容别名。
  - 多选语义：全部运行并合并仲裁 —— marker 为"规则并集 + 命中汇池"的合并库
    视图；panel/mash/sourmash 逐库检索后命中合并重排，判定由全局最优一次
    导出；taxa/对照基因按选择顺序 first-wins。
  - 溯源：`database` 报告合并名（如 `food+lab2026`），多库运行时方法载荷
    额外携带 `databases` 数组（name/version/path）；单选输出与单库形态完全
    一致。
  - `gside db list` 新增 named databases 小节（json 为 `custom_databases` 键）；
    `gtdbtk` 拒绝 `-d`。

### Removed

- **宏基因组（群落）分析整体迁移至 [gmeta](https://github.com/indexofire/gmeta)**。gside 回归单菌物种鉴定的定位。

  - 移除 `gside scan --meta`（kraken2/bracken 与 sylph 群落谱分析）、`--samples` 批量清单、`gside meta` 别名与 MCP `meta` 工具。
  - 移除 `analysis/kraken2_identifier.py`、`analysis/sylph_identifier.py` 及 `gside.api` 中的 `identify_by_kraken2` / `identify_by_sylph`。
  - 移除数据库层 `meta`（D5_kraken2，Kraken2 Standard-8）与 `meta_sylph`（D6_sylph）及其自动下载器；gseek 桥接的 sylph 后端一并移除。
  - pixi 依赖移除 kraken2/bracken/sylph。

  迁移对照：`gside scan --meta …` → `gmeta run -m kraken2|sylph …`（含 `--samples` 批量与 `-J/--jobs` 并行，均已在 gmeta 实现并扩展至 centrifuge/kaiju/metaphlan/kmcp）；`kreport`/bracken/sylph 结果导入 → `gmeta report` / `gmeta net build`（国家致病菌识别网上报格式）。

### Changed

- 覆盖率基线刷新：775 tests / 99.8%（配置 D2_ani+D3_mash 后），除 ani 97%
  外全部模块 100%。
