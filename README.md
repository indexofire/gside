# gside — Genome Species IDentification Engine

独立的多层物种鉴定 CLI 工具（hermes-bacmap 架构拆分产物）。

## 架构定位

```
底层工具     blastn / skani / mash / sourmash / mmseqs2
                ↑
gside          物种鉴定引擎：靶基因组合规则 + ANI 多后端 + 层级仲裁
                ↑
hermes-bacmap  Hermes 桥：GOM 入库 / smk 编排 / gbrain 记忆
```

## 快速开始

```bash
pixi install && pip install -e .

gside species contigs.fna --mode marker          # 靶基因组合规则（38 条）
gside species contigs.fna --mode panel           # ANI（精选面板 skani）
gside species contigs.fna --mode mash_refseq     # MinHash（RefSeq sketch）
gside species contigs.fna --mode all             # 多法 + 仲裁（ANI 层 > marker 层）
```

输出：单 JSON 判定（`analysis_type/method/database/result/verdict`，GOM 兼容契约）。

## 数据目录

| 变量 | 默认 | 内容 |
|---|---|---|
| `GSIDE_DATA_DIR` | `<repo>/data` | marker_rules.yaml + markers_v2.fasta（81 序列）+ BLAST 库 |
| `GSIDE_DB_DIR` | `<data>/db` | panel.sketch / mash.msh（大库，与 bacmap 共享） |
| `GSIDE_PIXI_BIN` | 自动发现 | 生信二进制目录（blastn/skani/mash） |

自定义规则：编辑 `data/reference/species/marker_rules.yaml`（min_hits/min_identity/
exclude_genes 数据化），重建 BLAST 库后即生效——无需改代码。

## 验证状态

- 17/17 金标准基因组（marker 法，含 10 株近缘阴性对照）
- 3/3 新增物种真实参考株（C.jejuni / C.coli / S.pyogenes）
- panel 法 17/17（共享 bacmap 291 基因组面板）

## License

MIT
