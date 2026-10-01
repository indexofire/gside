# CLI 参考

## gside species

物种鉴定主命令。

### 用法

```bash
gside species <contigs.fasta> [选项]
```

### 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `contigs` | 位置参数 | — | 组装后的 contigs FASTA 文件 |
| `--mode` | 选项 | `marker` | 鉴定模式（见下） |
| `--db-dir` | 选项 | 自动 | 数据库根目录 |

### 鉴定模式

| Mode | 说明 | 数据库需求 |
|---|---|---|
| `marker` | 靶基因组合规则（blastn） | mini（内置） |
| `panel` | ANI 比对精选面板（skani） | panel tier |
| `mash_refseq` | MinHash 距离（mash） | mash tier |
| `sourmash` | sourmash GTDB gather | sourmash GTDB 库 |
| `all` | 运行全部可用方法 + 仲裁 | 按需 |

### 示例

```bash
# 靶基因鉴定（最快，无需额外数据库）
gside species contigs.fna --mode marker

# ANI 鉴定（精确）
gside species contigs.fna --mode panel

# 多法仲裁
gside species contigs.fna --mode all

# 指定数据库目录
gside species contigs.fna --mode panel --db-dir /path/to/dbs
```

## gside db

数据库管理子命令。

### 子命令

| 子命令 | 说明 |
|---|---|
| `status` | 显示数据库就绪状态 |
| `list` | 列出可用 tier |
| `setup` | 安装/更新数据库 |

### setup 选项

| 选项 | 默认 | 说明 |
|---|---|---|
| `--tier` | `panel` | mini / panel / mash / all |
| `--source` | — | 从现有 bacmap 数据库复制 |

### 示例

```bash
gside db status
gside db setup --tier all
gside db setup --tier panel --source ~/bacmap/data/db
```

## gside --version

```bash
gside --version
# gside 0.1.0
```
