# 标记规则格式

`marker_rules.yaml` 定义靶基因组合鉴定规则。新增物种只需添加规则条目，
无需修改代码。

## 文件位置

```
data/db/D1_marker/marker_rules.yaml
```

## 规则结构

```yaml
rules:
  - species: <物种名>              # 必填：注册表键名（下划线连接）
    genes:                          # 必填：候选标记基因列表
      - <gene1>
      - <gene2>
    min_hits: <int>                 # 可选：最少命中数（默认 1）
    min_identity: <float>           # 可选：最低一致性 %（默认 90）
    min_coverage: <float>           # 可选：最低覆盖率 %（默认 60）
    exclude_genes:                  # 可选：排除守卫
      - <cross_reactive_gene>
    note: <str>                     # 可选：注释（文档性）
    co_detect_note: <str>           # 可选：共检出说明（文档性）
```

## 字段说明

### species

物种名，使用注册表键名（空格转下划线）。如 `Campylobacter_jejuni`、
`Vibrio_parahaemolyticus`。

### genes

候选标记基因名列表（小写）。这些基因必须存在于 `markers.fasta` 中。

### min_hits

最少命中数。设为 2 时，`genes` 中至少 2 个基因需达到阈值才算该物种。

多基因组合可显著降低交叉反应（如 C.jejuni 的 3 选 2 设计）。

### min_identity

每个命中基因的最低一致性百分比（blastn pident）。默认 90。

降低至 85 适用于短序列或保守基因（但会增加假阳性风险）。

### exclude_genes

**排除守卫**：如果这些基因在样本中检测到，则该物种的判定被抑制。

典型用法，弯曲菌种间区分：

```yaml
- species: Campylobacter_coli
  genes: [ceue]
  exclude_genes: [mapa, hipo, cadf]    # jejuni 标记存在 → 倾向 jejuni
- species: Campylobacter_jejuni
  genes: [mapa, hipo, cadf]
  min_hits: 2
  # 无排除：jejuni 基因组天然含 ceuE 同源体
```

### min_coverage

基因计入本规则所需的最低 query 覆盖率（默认 60，与全局门槛一致）。
判定时逐条规则检查该基因的所有原始命中，因此"一致性稍低但覆盖率更高"
的命中仍可满足规则。

仅当规则确实接受部分基因匹配时才调低（例如 50）——两条
*Campylobacter* 规则就是这样配置的。

### note 与 co_detect_note

文档性字段。记录规则设计意图与近缘种共检出提示，仅供规则维护者参考，
不参与判定逻辑。

## 全局阈值

引擎层设有全局门槛（在 `multigene_identifier.py` 中定义）：

| 参数 | 值 | 说明 |
|---|---|---|
| `_MIN_IDENTITY` | 85.0 | 全局最低一致性 |
| `_MIN_COVERAGE` | 60.0 | 全局最低覆盖率（防短 HSP 误报） |
| `_HIGH_CONF` | 90.0 | 高置信度分界 |

规则级 `min_identity` / `min_coverage` 作用于该基因的每条原始命中；
未填写的规则使用 90 / 60。全局 `_MIN_IDENTITY`（85）与 `_MIN_COVERAGE`（60）
仍用于 `detected_markers` 上报与 `exclude_genes` 排除守卫。

`priority_over` 引擎支持（命中的 priority 物种会压制非 priority 物种）。四条志贺菌
种级规则（`Shigella_flexneri`、`Shigella_sonnei`、`Shigella_dysenteriae`、
`Shigella_boydii`）使用 `priority_over: [Shigella, Escherichia_coli]` 与
`require: [uida, ipah]`：种级决定因子（如 `rfc`、`wbgZ`、`sb1_wzx`）会升级属级
判定——且 O-抗原决定因子单独在场、缺少 uidA+ipaH 属锚时永不触发种级判定——
而单纯的 `uidA+ipaH` 仍停留在 `Shigella` 属级；所有志贺菌规则均为 lacY 阴性
（`exclude_genes: [lacy]`），EIEC 不会触发志贺菌判定。

## 跨物种共检出与仲裁 {: #cross-species-arbitration }

纯培养基因组也可能同时命中多个物种的标记基因。常见原因：保守管家基因的
共同祖先信号、水平转移（毒力质粒、致病岛、前噬菌体）、多拷贝旁系同源基因，
以及混样/污染组装。标记命中只是证据，不是结论；引擎按固定顺序仲裁：

1. **命中门槛**——低于全局下限（85% 一致性、60% 覆盖率）的命中直接丢弃。
   随后每条规则再对自身基因的每条原始命中应用规则级 `min_identity`
   （默认 90）与 `min_coverage`（默认 60）。
2. **规则法定数**——`min_hits` 规定至少几个基因通过；`require` 中的基因
   必须全部在其中（硬性 AND）。达不到法定数的规则不进入候选。
3. **排除守卫**——当 `exclude_genes` 中的交叉反应基因在全局下限被检出时，
   即使规则自身标记通过，该规则也被抑制（如 hipo/cadf 抑制大肠弯曲菌；
   所有志贺菌种级规则均为 lacY 阴性，EIEC 仍归大肠杆菌）。
4. **优先级**——候选规则声明 `priority_over` 时，无论得分高低都压制所列
   物种（鼻疽 BimA(Bm) 压制类鼻疽；志贺菌种级压制志贺菌属/大肠杆菌；
   邦戈尔沙门菌压制肠道沙门菌）。优先级是解决共享标记竞争的标准手段。
5. **得分**——剩余候选中，命中基因数最多者胜出；平局优先高置信度
   （平均一致性 ≥90%），仍平局则按 YAML 中的规则顺序决定。

输出保留证据：`detected_markers` 列出所有显著基因，`matched_rule` 给出
获胜基因组合，`notes` / `co_detect_note` 解释混合画像。例如只检出炭疽
芽孢杆菌质粒基因、缺染色体标记时，报告蜡样芽孢杆菌群并附 "Anthrax-like"
提示，而不是强行给出种级判定。

marker、ANI 与 GTDB-Tk 结果不一致时（`all` 模式），仲裁按层级优先：
GTDB-Tk（L3）> ANI/mash/sourmash（L2）> marker（L1）；同层高置信度
优先于中置信度，`verdict.supporting` / `verdict.agreement` 记录跨方法
互证。标记命中永不覆盖全基因组 ANI 判定，ANI 判定也不隐藏标记画像——
两者不一致时应同时查看。

解读建议：

- 单个脱靶标记、未达规则法定数：属于背景，不是混合样本。
- 多个外源物种毒力标记：提示水平转移或致病型杂合，需 L2 确认。
- 跨属广泛高一致性命中：提示污染或混样组装，应复核输入。

## CPIN 覆盖块

规则文件中有一段标注的块，覆盖 `data/CPIN.txt` 物种清单。条目分三类：

- **经典靶基因**：按 `locus_tag` 从经核验的 RefSeq 参考基因组提取
  （IS6110、TUL4/FopA、侵袭素+VirF、BimA(Bm)、ompA、hugA、BP26/Omp25 等）。
  规则 note 记录参考株与实测最近邻交叉同一性；阈值设在交叉之上。
- **条码基因**（`Citrobacter`/`Escherichia_albertii` 用 `gyrB`，
  `Vibrio_fluvialis`/`Burkholderia_gladioli` 用 `rpoB+gyrB`，
  `Bartonella` 用 `gltA+rpoB`）：经典靶点无法区分近邻时使用 ——
  肠杆菌科 `rpoB` 交叉高达 97%，因此改由 `gyrB`（交叉 ≤89.5%）承担。
- **真菌 ITS 条码**（`its-*`，来自类型株/临床株记录的 ITS1 或
  ITS1-5.8S-ITS2）：属级判定，92% 阈值 —— Onygenales 各属间交叉最高 90.5%。

每条规则都经参考基因组端到端验证 + 近缘对抗基因组验证
（鼠疫杆菌 vs 假结核耶尔森氏菌、鼻疽伯克霍尔德菌 vs 类鼻疽伯克霍尔德菌），
见 `tests/integration/test_cpin_markers.py`。

BimA(Bm) 规则的阈值为 98% 一致性：100 株临床类鼻疽中有 3 株携带 Bm 样
等位基因（96.2–96.4%），而真正鼻疽菌为 99.6–100%。旧阈值 95% 会让这 3 株
被 `priority_over` 压制、误报为 `Burkholderia_mallei`。两者为同一基因组种
（ANI ~99.5%）；详见下文类鼻疽伯克霍尔德菌条目。

## Reads 模式深度门控（FASTQ 输入） {: #reads-mode-depth-gate }

FASTQ 输入时，标记引擎在一致性/覆盖规则之上叠加 gspec 风格的丰度门控。
每个基因上报 `mean_depth`（比对目标区间长度之和 / 基因长度 —— 不是读段数）、
`reads` 和 `mapq_mean`；低于 `min_mapq`（默认 20）的比对在聚合前丢弃。

标记基因需同时满足：

```
mean_depth >= min_depth  且  mean_depth / control_depth >= min_relative_depth
```

默认值（在 `metadata.yaml` 中，调用时可用 `min_depth` / `min_relative_depth` /
`min_mapq` 覆盖）：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `min_depth` | 5.0 | 绝对平均深度下限 |
| `min_relative_depth` | 0.10 | 相对对照深度的比例 |
| `min_mapq` | 20 | 聚合前的 minimap2 MAPQ 下限 |

`control_depth` 是该物种配置的单拷贝对照基因（`metadata.yaml: control_genes`，
如 *Salmonella* 的 `rpoD`、*Campylobacter* 的两个 `lpxA` 等位基因）的
mean depth 中位数。每个基因取最深的等位基因，避免少量污染读段压低对照。
对照基因在 `markers.fasta` 中标记 `role=control`，绝不会出现在
`detected_markers` 中。没有专用对照的物种回退到已检出标记的样本中位深度
（`control_source: "fallback"`）；完全无对照时关闭相对门控，仅保留绝对下限。
组装（BLAST）输入没有深度信息，保持旧行为。

## 沙门菌种级规则

内置库用三层规则区分沙门菌的两个种（实证基础：26 株沙门菌面板——
*S. enterica* 全部 6 个亚种 + 13 株 *S. bongori*——及 5 株外群；背景文献
Fookes 2011, PLoS Pathog e1002191）：

| 规则 | 基因 | 阈值 | 含义 |
|---|---|---|---|
| `Salmonella_enterica` | `inva`（必命中）+（`ssav` \| `ssec`），排除 `asst`/`zadh` | 3 选 2，≥ 90 | 属锚 AND SPI-2 标记（T3SS-2 装置 `ssav`/转位孔 `ssec`，11 株 *S. enterica* 面板[覆盖 6 亚种]核心、所有邦戈尔与外群缺失）；属级阈值保证 divergent enterica 谱系不漏检 |
| `Salmonella_bongori` | `inva`（必命中）+（`asst` \| `zadh`） | 3 选 2，≥ 90 | 属锚 AND 邦戈尔特异标记（`asst` = 芳基硫酸磺基转移酶、`zadh` = 锌依赖乙醇脱氢酶） |
| `Salmonella` | `inva` | 一致性 ≥ 90 | 属级回退——inva 阳性但种级规则未触发（邦戈尔 draft 缺双标记、SPI-2 阴性 enterica、标记落在 85–90 检出带） |

两个种级规则均 `require: inva`（属信号必须在命中集合内），
`Salmonella_bongori` 另有 `priority_over: [Salmonella, Salmonella_enterica]`：
特异标记单独在场永不触发种级判定，携带标记的邦戈尔基因组必然判对。
已知权衡（有意）：缺失 SPI-2 的 *S. enterica*（罕见，历史 subsp. arizonae
菌株）回退属级标签。血清型/亚种分型仍归 gspec。

## 百日咳鲍特菌：抗 IS 碎片化设计

百日咳鲍特菌携带数百拷贝 IS481，draft 组装高度碎片化（700+ contigs 常见；
20 株面板中 30% 的 RefSeq draft 将 tcfA 断在两个 contig 上）。规则设计吸收
了这一现实：

| 基因 | 角色 |
|---|---|
| `cupin-bp`（681bp，必命中） | 染色体种锚——20 株多样化菌株（含 ptxP3 时代）100% 保守，近缘种最高仅 71%（ptxS1 与副百日咳/霍氏/支气管败血鲍特菌交叉 96–99%、IS481 与霍氏共享，均不能作锚） |
| `tcfa5` / `tcfa3`（各 ~972bp） | 1944bp 的 tcfA 拆半——IS 断点将基因切碎时，每个半段可完整落在单个 contig 上 |
| `ptxs1`（810bp） | 百日咳毒素 S1——毒力位点报告 |

百日咳的毒力在染色体上（无炭疽那样的质粒问题）：Bvg⁻ 相变异株与解毒疫苗
衍生株均可检出——历史上的漏检是 IS 碎片化伪影，而非无毒株漏检。

## 鼠疫耶尔森菌：染色体锚 + 质粒毒力报告

鼠疫菌是假结核耶尔森菌的极新近进化枝（染色体同一性 >97%），属内染色体鉴
别最难。50 株真实面板（覆盖全部三个生物型）× 5 近缘基因组的反向筛选仅
发现 4 个染色体特异基因——该菌高度克隆（所有候选 100% 保守）：

| 基因 | 角色 |
|---|---|
| `spi-yp`（987bp，必命中） | 染色体糖磷酸异构酶种锚——50 株 100% 保守，近缘零同源 |
| `btub-yp`（1878bp） | TonB 依赖 BtuB——第二染色体标记（有文献先例） |
| `caf1`（pMT1） | F1 荚膜抗原——毒力质粒报告 |
| `pla`（pPCP1） | 纤溶酶原激活因子——毒力报告；13% 真实世界菌株丢失小 pPCP1 质粒 |

原 `min_hits: 2` 的 `[caf1, pla]` 设计使所有 pPCP1 阴性株报 Unknown。
引入染色体锚后，质粒丢失不再影响种级判定。

## 小肠结肠炎耶尔森菌：种锚 + 致病性报告

原 `ail`/`yada` 条目是与任何基因组都不匹配的坏序列（规则从未触发——100 株
真实世界基因组中 76 株报 Unknown）。重构后的规则将物种身份与致病性解耦：

| 基因 | 角色 |
|---|---|
| `glk-ye`（1164bp，必命中） | 甘油酸激酶染色体锚——50 株真实面板（覆盖 subsp. enterocolitica/palearctica 与生物型 1A）100% 在场且 ≥96%，对假结核/弗雷德里克森/中间/克里斯滕森耶尔森菌、鼠疫菌与大肠埃希菌零同源 |
| `yddg-ye`（888bp） | 第二染色体标记（冗余随行） |
| `ail`（537bp） | 真实 AIL CDS——仅致病性生物型（1B/2–5）存在；非致病 1A 生物型缺失 |
| `yada`（1269bp） | 真实 pYV 质粒 YadA——毒力报告；质粒易丢失 |

非致病 1A 生物型与 pYV 丢失株因此仍报 `Yersinia_enterocolitica`（仅锚命中），
并附 note 提示无致病性标记。

## 蜡样芽孢杆菌群：种 vs complex

蜡样芽孢杆菌群（炭疽/蜡样/苏云金/蕈状/韦氏芽孢杆菌等）染色体同一性 >99%，
只有炭疽支系可通过 λBa03 前噬菌体在染色体层面区分。内置规则如实反映这一
生物学事实（实证基础：404 基因组面板——150 株炭疽 + 246 株非炭疽群成员 +
8 株外群；ring trial 文献：Wielinga 2013、Bentahir 2018）：

| 规则 | 基因 | 阈值 | 含义 |
|---|---|---|---|
| `Bacillus_anthracis` | `adk`（必命中）+ `ba5345`（必命中），`paga`/`capb` 随行 | 4 选 2，≥ 90 | 群锚 AND λBa03 染色体标记；质粒基因只报告毒力状态、绝不作为种级门槛——脱毒/疫苗株（Sterne）仍是炭疽 |
| `Bacillus_cereus_group` | `adk` | ≥ 90 | complex 级回退：蜡样/苏云金等群成员，以及 λBa03 阴性的炭疽 draft（该前噬菌体短读组装易丢失——RefSeq draft 中占 27%） |

`paga`/`capb` 阳性而 `ba5345` 缺席时（B. cereus biovar anthracis / 交叉菌株
画像，如 Bcbva CI、BcFL2013、03BB87），结果附带 "Anthrax-like" 提示并建议
送参考实验室确认——质粒对单独在场永不触发 `Bacillus_anthracis`。群级标签的
`verdict.scope` 为 `complex`。

## 物种标记参考（含病原生物学）

每条规则的标记集合反映了底层病原生物学。下表解释每个基因测量什么、为什么
选择它，帮助用户深入理解输出中的 `detected_markers` 和 `notes`。

### 沙门菌（肠道沙门菌 vs 邦戈尔沙门菌）

| 基因 | 角色 | 生物学 |
|---|---|---|
| `inva` | 属锚（必命中） | 侵袭蛋白 A；全属保守（invA⁻ 突变体为非侵袭实验室构建）。一致性带区分种：≥93% = 肠道沙门菌，90–93% = 邦戈尔沙门菌（26 基因组面板验证）。 |
| `ssav` | 肠道沙门菌特异 | SPI-2 三型分泌系统装置蛋白。SPI-2 是胞内存活岛，所有邦戈尔沙门菌和外群均缺失。 |
| `ssec` | 肠道沙门菌特异 | SPI-2 转位孔蛋白。同岛独立基因，提供冗余。 |
| `asst` | 邦戈尔沙门菌特异 | 芳基硫酸磺基转移酶，染色体。肠道沙门菌和外群零全长同源。 |
| `zadh` | 邦戈尔沙门菌特异 | 锌依赖乙醇脱氢酶。与 asst 同种特异性。 |

**解读**：`Salmonella_enterica` = 属信号 + SPI-2 种特异标记。
`Salmonella_bongori` = 属信号 + 邦戈尔代谢标记。裸 `Salmonella` = invA
阳性但种级规则未触发（部分组装/分歧谱系）。血清型分型（Typhi、
Enteritidis 等）属于 **gspec**，不属于 gside。

### 蜡样芽孢杆菌群

| 基因 | 角色 | 生物学 |
|---|---|---|
| `adk` | 群锚（炭疽必命中） | 腺苷酸激酶，全群保守（>99%）。区分蜡样芽孢杆菌群与其他芽孢杆菌和外群。 |
| `ba5345` | 炭疽特异（必命中） | λBa03 前噬菌体染色体标记；唯一经 ring trial 验证的炭疽染色体区分标记（Wielinga 2013，134 基因组）。 |
| `paga` | 毒力报告（pXO1 质粒） | 保护性抗原，炭疽毒素组分。报告产毒性；质粒清除/疫苗株（Sterne）无此基因仍为炭疽。 |
| `capb` | 毒力报告（pXO2 质粒） | 荚膜合成。完全毒力需双质粒。 |

**解读**：`Bacillus_anthracis` = 染色体身份确认。`Bacillus_cereus_group`
= 蜡样群成员但无炭疽染色体——若 paga/capb 在场，输出 "Anthrax-like" 提示
（Bcbva/交叉菌株，或 λBa03 阴性 draft）。**scope** 为 `complex`。

### 百日咳鲍特菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `cupin-bp` | 种锚（必命中） | Cupin 域染色体蛋白；多样化菌株 100% 保守。近缘种最高 71%。 |
| `ptxs1` | 毒力报告 | 百日咳毒素 S1 亚基（染色体，Bvg 调控）。与副百日咳/霍氏鲍特菌交叉 96–99%——仅报告，不作种级门槛。 |
| `tcfa5`/`tcfa3` | 毒力报告（拆半） | 气管定植因子 1944bp 拆为两半。IS481（数百拷贝）碎片化组装；每半可完整落在单 contig。 |

**解读**：所有百日咳鲍特菌（含 Bvg⁻ 相变异和解毒疫苗衍生株）携带染色体
锚。ptxs1/tcfa 报告毒素位点存在。IS481 已弃用（与霍氏鲍特菌共享）。

### 鼠疫耶尔森菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `spi-yp` | 种锚（必命中） | 染色体糖磷酸异构酶。鼠疫菌是假结核菌的极近进化枝（>97% 染色体同一性）——仅 4 个染色体基因可区分。 |
| `btub-yp` | 第二染色体标记 | 维生素 B12 受体（有文献先例）。 |
| `caf1` | 毒力报告（pMT1 质粒） | F1 荚膜抗原。13% 真实世界菌株丢失 pMT1。 |
| `pla` | 毒力报告（pPCP1 质粒） | 纤溶酶原激活因子。9.6kb 小质粒易在传代中丢失。 |

**解读**：种级身份来自染色体，产毒性来自质粒状态。pPCP1 阴性株（常见
实验室衍生株）仍报鼠疫耶尔森菌。

### 小肠结肠炎耶尔森菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `glk-ye` | 种锚（必命中） | 甘油酸激酶，染色体。50 株真实面板（生物型 1A–4）100% 保守。近缘种零同源。 |
| `yddg-ye` | 第二染色体标记 | 芳香族氨基酸通透酶。 |
| `ail` | 致病性报告 | 附着侵袭位点，染色体。仅致病性生物型（1B、2–5）存在；非致病 1A 缺失。 |
| `yada` | 毒力报告（pYV 质粒） | 耶尔森菌黏附素 A。70kb 质粒易丢失。 |

**解读**：非致病 1A 生物型（食品调查中常见）仍报小肠结肠炎耶尔森菌，
附注"未检出致病性标记"。致病性生物型额外显示 ail/yada。

### 脑膜炎奈瑟氏菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `sodc` | 种锚（必命中） | 周质 Cu/Zn 超氧化物歧化酶。脑膜炎/淋病奈瑟氏菌鉴别的经典标记（Corless 2001）。所有脑膜炎菌携带（含无荚膜携带株）。 |
| `ctra` | 荚膜报告 | 荚膜转运蛋白。不可分群（NG）和无荚膜株缺失；种级判定不受影响。 |

### 淋病奈瑟氏菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `pora` | 等位基因家族 1（porB1a） | 主要膜孔蛋白，原命名 "porA"。覆盖 porB1a/PIA 等位基因家族。 |
| `porb1b` | 等位基因家族 2（porB1b） | 第二抗原分歧膜孔蛋白等位基因家族（~78% 家族间同一性）。任一家族命中即可。 |

### 肺炎链球菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `lyta5`/`lyta3` | 种标记（拆半） | 自溶素（N-乙酰胞壁酰-L-丙氨酸酰胺酶），肺炎链球菌鉴定金标准。957bp 拆半抗碎片化（高重组率导致 6/100 draft 断裂）。 |

### 流感嗜血杆菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `hpd` | 种标记 | 血红蛋白结合蛋白 D。流感嗜血杆菌特异；副流感/溶血/皮特曼嗜血杆菌无全长同源。 |

### 副流感嗜血杆菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `ribf-hpi` | 种锚（必命中） | GTP 环化水解酶 II（核黄素合成）。确认株 98.6%；近缘最高 83.4%。临床实验室常将副流感误鉴定为流感嗜血杆菌。 |
| `sera-hpi` | 第二标记 | 磷酸甘油酸脱氢酶。 |

### 肺炎支原体

| 基因 | 角色 | 生物学 |
|---|---|---|
| `tuf-mp` | 种锚（必命中） | 延长因子 Tu，1185bp。30 株多样化菌株 100% 保守；生殖/人型支原体零同源。 |
| `gyrb-mp` | 第二染色体标记 | DNA 拓扑异构酶 B。 |
| `p1a`/`p1c` | 细胞黏附报告（4 段中 2 段） | P1 黏附素（4884bp）拆为四段。全长需 810kb 基因组的 0.6% 在单 contig 上；分段抗碎片化。毒力基因，非种锚。 |

### 白喉棒状杆菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `rela-cd` | 种锚（必命中） | GTP 焦磷酸激酶（严紧反应）。多样化菌株 98.3–100%；近缘最高 70.5%。 |
| `mena-cd` | 第二染色体标记 | 邻琥珀酰苯甲酸合酶（维生素 K2 合成）。 |
| `tox` | 产毒报告 | 白喉毒素，β-棒状杆菌噬菌体携带。37% 真实世界菌株非产毒（无噬菌体）。非产毒白喉可致侵袭性疾病（1990 年代以来报告）。 |

### 肉毒梭菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `bonta` | 神经毒素 A 型 | 3891bp，ATCC 3502（A1 型参考株）。BoNT 血清型间仅 30–40% 同一性——每型需独立标记。 |
| `bontb` | 神经毒素 B 型 | 覆盖 B 型菌株（部分亚型）。 |
| `bonte` | 神经毒素 E 型 | 覆盖 E 型菌株（食源性，鱼类相关）。 |
| `bontc` | 神经毒素 C 型 | 噬菌体 c-st 来源。覆盖动物相关 C 型。 |
| `bontd` | 神经毒素 D 型 | 覆盖 D 型及 C/D 嵌合体。 |
| `bontf` | 神经毒素 F 型 | F634（F5 亚型）。覆盖罕见但临床重要的 F 型。 |

**解读**：肉毒梭菌是**毒素表型**而非基因组种——4 个基因组群（I–IV）仅由
BoNT 产毒统一。种级身份基于毒素；染色体背景（肉毒 vs 生胞梭菌，ANI
93.4%）需 L2 ANI 方法（`-m marker+mash`）。非产毒梭菌正确报 Unknown。

### 霍乱弧菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `ompw` | 种标记（一致性 85%） | 外膜蛋白 W。非 O1/非 O139 环境株携带分歧等位基因（86–90%，标记源为 O1 菌株）。最近缘非霍乱弧菌：拟态弧菌 80.8%。 |

**解读**：阈值降至 85% 以容纳非 O1 环境株多样性。O1/O139 大流行谱系
~99–100% 命中。

### 类鼻疽伯克霍尔德菌

| 基因 | 角色 | 生物学 |
|---|---|---|
| `glgb-bp` | 种锚（required） | GlgB 1,4-α-葡聚糖分支酶，2217bp 染色体基因。100 株临床类鼻疽基因组一致性 99.6–100%、覆盖率 100%；最近交叉：泰国伯克霍尔德菌 73.3%、洋葱伯克霍尔德菌 82.5%。 |
| `bima-bp` | 毒力报告 | BimA(Bp) 肌动蛋白成核因子，介导肌动蛋白运动与细胞间传播。95/100 面板基因组携带（部分菌株为 Bm 样或分歧等位基因）；永不参与种级判定。 |
| `sctn-bp` | 毒力报告 | T3SS ATPase SctN（TTS1 簇）。99/100 面板基因组 ≥95% 命中（1 株 TTS 簇缺失）；永不参与种级判定。 |

**解读**：类鼻疽伯克霍尔德菌（类鼻疽）与鼻疽伯克霍尔德菌（鼻疽）属于
同一基因组种（ANI ~99.5%）——后者是丢失约 1.5Mb、宿主受限的克隆。
标记层无法区分两者。内置 `cpin` 库另有一条 BimA(Bm) 等位基因规则，
以 `priority_over` 识别鼻疽；阈值定为 98%，因为 100 株临床类鼻疽中
有 3 株携带 Bm 样等位基因（96.2–96.4%）。GTDB rs226 将整个分支命名为
`s__Burkholderia mallei`，而 NCBI 与临床微生物学沿用 *B. pseudomallei*，
本规则报告后者。marker 与 ANI 结果不一致（pseudomallei vs mallei）反映的是
命名体系差异，不是误判——用 L2 ANI 确认并同时报告两者。

## 新增物种示例

```yaml
# 在 rules 列表中添加：
- species: Mycobacterium_tuberculosis
  genes:
  - mpb64
  - is6110
  - rpoB_mtbc
  min_hits: 2
  min_identity: 90
  exclude_genes: [rpob_maa]    # 排除 M.africanum
  note: MTBC complex identification
```

然后：

1. 将对应基因序列追加到 `markers.fasta`，FASTA 头格式：
   `>gside|db=<database>|gene=<gene>|acc=<accession>|func=<role>|note=<text>`；
   `func` 取 `primary`/`confirm`/`typing`/`virulence`/`control`
2. 重建 BLAST 库：
   `makeblastdb -in markers.fasta -dbtype nucl -out markers_blastdb`
3. 运行 `gside scan test.fna -m marker` 验证
