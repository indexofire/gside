# 未完成事项(2026-11-03 会话记录,同日收尾更新)

来源:gseek 0.10.4 升级 + marker 深度门控 WIP 收尾。

## 已完成(本次会话)

- [x] gseek 0.10.4 升级(依赖 + lock + 文档)
- [x] marker 深度门控 + 志贺菌种级规则:测试对齐新契约并新增 ~40 个用例
- [x] 覆盖率基线重测:D2_ani + D3_mash 配置后 **864 tests / 99.0%**,
      AGENTS.md 基线已更新(multigene / gseek_bridge / marker_db 均 100%)
- [x] gseek 0.10.4 parity 真库验证:marker/mash/sourmash/panel/screening
      全部通过(29 个集成测试 0 跳过);顺手修了 **db.py 预构建 panel 下载
      不发布 taxa map** 的 bug(`_publish_taxa_map`)
- [x] 预存覆盖欠账清零:marker_db 77%→100%,multigene 6 行→0,
      bridge 16 行→0
- [x] 下游确认:bacmap `canonical_from_verdict("Shigella_flexneri")` →
      `"Shigella"`(实测),门控行为不变;gside 新 JSON 字段经 `.get()`
      读取,无严格校验,兼容
- [x] 中文锚点:mkdocs 启用 attr_list,`#reads-mode-depth-gate` 精确锚点

## 仍开放

- [ ] **marker 数据库内容性欠账**(需人工审核参考基因组,不可自动完成):
      14 个 sample_median 回退物种仍无专用对照基因(Vibrio_cholerae、
      V._vulnificus、E._faecium、Cronobacter_sakazakii、Mycoplasma_pneumoniae、
      H._pylori、S._suis、B._pseudomallei、L._interrogans、T._pallidum、
      Yersinia×2、A._hydrophila、C._botulinum)— 各需一条经核实参考基因组的
      rpoB/gyrB 序列(流程见 CONTROL_GENES.md「How to add a new control gene」)
- [ ] 志贺菌 4 个种规则用合成序列验证过规则引擎;真实志贺菌参考基因组
      (S._flexneri 2457T、S._sonnei 53G、S._dysenteriae Sd197、S._boydii CDC
      3083-94)的端到端判定尚未跑(可选;基因组放入 D2_ani/genomes 后
      `pytest tests/integration -k golden` 即可)
- [ ] 预存低覆盖行(与旧基线持平,非回归):ani 92%(main() argparse 分支)、
      sylph 95%、sourmash 96%、kraken2/db/cli 99%
- [ ] 提交后推送 + CI 观察一次(CI 无 gitignore 数据库,预计 ~97%,floor 95)

## 本机最终验证状态(2026-11-03)

ruff check/format ✅;pyright 0 错误 ✅;`mkdocs build --strict` ✅;
**864 tests 全部通过(集成 0 跳过)** ✅;总覆盖 **99.0%**;
D1_marker(bundled)+ D2_ani(panel)+ D3_mash(mash)已配置。
