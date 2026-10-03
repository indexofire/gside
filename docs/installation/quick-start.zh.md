# 快速安装

从零搭建 gside 到跑通物种鉴定。预计耗时 5~15 分钟（含数据库）。

## 前置条件

- Linux x86_64
- [pixi](https://pixi.sh)（管理生信 CLI 工具）
- Python ≥ 3.12

## 安装步骤

### 1. 克隆仓库

```bash
git clone https://github.com/indexofire/gside.git
cd gside
```

### 2. 生信工具环境

```bash
pixi install
```

自动拉取 blast / skani / mash / sourmash / mmseqs2 等。

### 3. Python 包

```bash
pip install -e .
```

### 4. 验证

```bash
gside --version
gside db status     # markers 应为 ✅
```

## 数据库

gside 采用分层数据库设计：

| Tier | 内容 | 大小 | 状态 |
|---|---|---|---|
| mini | 标记规则 + 序列 + BLAST 库 | ~1MB | ✅ 随仓库分发 |
| panel | 精选参考面板（skani sketch，291 基因组） | ~130MB | 需 `gside db setup` |
| mash | RefSeq MinHash sketch | ~179MB | 需 `gside db setup` |
| all | panel + mash | ~310MB | 一次性全装 |

详见[数据库管理](databases.md)。

## 下一步

- [CLI 参考](../usage/cli.md) — 完整命令行参数说明
- [物种鉴定模式](../usage/modes.md) — 各层次原理与适用场景
