# 选择器自适应收益文献核查：算子选择的外部先例与"确有收益"的判定框架（2026-09-26）

日期：2026-09-26。状态：**快速文献核查（非系统综述），为选择器设计的"自适应确有收益"前置证据要求提供外部方法先例**。

## 1. 目的与范围

仓库纪律（AGENTS.md"ML 与评价"、[传统基线研究](2026-09-24_cycle03_traditional_baselines.md)、[运行器/选择器设计 v0.1](2026-09-26_eval_runner_and_selector_design_v0.1.md) §B/§C）要求：进入选择器训练前，必须先有固定流程/有限搜索的参照与"自适应确有收益"的证据。本核查把该要求外投到文献，回答两个问题：

1. GPR 处理文献中，"自适应"（参数或算子层面）相对固定基线**被定量证明有收益**的先例以什么证据形态出现？
2. 更一般的算子选择（algorithm selection）文献中，"选择器比最佳单一算法更好"的**可判定条件**是什么？

本次为**快速核查**：全文通读 2 篇，摘要/片段级 4 项，仅登记引用 3 条；阅读范围逐条记录于[本轮台账](2026-09-26_selector_adaptivity_literature_ledger.json)。

## 2. 来源与要点

### S1（全文通读）Xue et al. 2019

*Noise Suppression for GPR Data Based on SVD of Window-Length-Optimized Hankel Matrix*，Sensors 19，[PMC6749384 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6749384/)。阅读范围：摘要至 §3.2 主体（方法、合成例 1/2）；§3.3 之后实测例与结论仅经页面提取片段。

- **参数级自适应的 GPR 先例**：Hankel 矩阵窗长（SVD 去噪的唯一参数）用"奇异值四阶中心矩四次方根（FRFCM）最大化"准则**逐数据自适应选定**；奇异值分界点用差分谱阈值自动确定。
- **证据形态（与本仓库纪律同构）**：自适应方法与两个固定基线（局部能量比准则的常规 SVD、小波阈值）在 gprMax 合成数据（混凝土背景 εr=6/σ=0.01、PEC 圆柱、900 MHz Ricker，加高斯白噪至 SNR=−5 dB）上对比，报 SNR（7.55 vs 4.23/7.08 dB）**并同时报成本**（耗时 4.17 vs 1.9/2.31 s、内存 99 vs 71/39 MB）。
- 关键界限：这是**单一算子内部的参数自适应**，不是"在多个候选算子之间做选择"；收益证明依赖已知真值的合成数据（SNR 可算），实测数据上只有定性/对比证据。

### S2（全文通读）Messelis & De Causmaecker（EJOR 预印本）

*An automatic algorithm selection approach for the multi-mode resource-constrained project scheduling problem*，[KU Leuven LIRIAS 全文](https://lirias.kuleuven.be/retrieve/f285366b-15a6-4928-a9c4-b5a515886cf2/)。阅读范围：摘要至 §3.2 Step 4（问题定义、有用性框架、经验硬度模型构建）；Step 5 之后与附录仅经提取片段。

- **"选择器是否有用"的可判定前置条件**（§3.1）：提出竞争力比 c（两算法各自占优的实例子集之较小占比 ×2，可分解为等势度 e × 覆盖度 r）与潜在收益 i（完美选择器相对最佳单一算法策略的相对性能改进上界）。**只有在实例分布上 c 与 i 都足够大时，选择器才值得建**。
- **反面证据同样重要**：同一对算法在 PSPlib 基准上 c=0.2194、i=0.41% ⇒ 选择器无意义（最佳单一算法在 1−c/2≈89% 实例上已是最优，构成任何可接受选择器的**自然下界**）；在 MMlib 基准上 c=0.7957、i=2.52% ⇒ 选择器显著优于任一单算法。**收益是实例分布的条件属性，不是方法属性**。
- 方法链：实例特征 → 经验硬度模型（逐算法性能预测）→ 选预测最优者（Rice 1976 框架的工程化，Leyton-Brown 六步流程）；训练/验证划分在实例层面（3140/1000）。
- 提及 SATzilla（Xu et al. 2008）为同框架在 SAT 求解器上的成功案例（仅登记）。

### S3（摘要/片段级）GPR 与超声的自适应杂波抑制

- Chen & Fu 2017（*Adaptive Ground Clutter Reduction in GPR Data Based on PCA*，[Semantic Scholar 摘要](https://www.semanticscholar.org/paper/e65423d7b58066a2152ac11df7948fe61f02f5c3)）：SVD 奇异值自适应选择 + 自动数据评价——与 S1 同族的参数级自适应。
- Hou et al. 2025（*An Adaptive SVD-Based Approach to Clutter Suppression for Slow-Moving Targets*，Remote Sensing 17(15):2697，[MDPI 页](https://www.mdpi.com/2072-4292/17/15/2697)，本次未获全文）：A-SVD 解决"SVD 阈值选择问题"——同为参数级。
- arXiv 2405.11105（超声小血管成像自适应杂波滤波，[PDF 片段](https://arxiv.org/pdf/2405.11105v1)）：明确指出**现有"自适应"方法多数仍需预定义初值，且最优值随应用场景逐案经验调整**——对"自适应"宣称的重要保留。
- APSIPA 2022（Cao et al.，[PDF 片段](https://www.apsipa.org/proceedings/2022/APSIPA%202022/ThPM2-7/1570834693.pdf)）：基于 Stein 无偏风险估计（SURE）的自适应阈值去噪——统计风险驱动的参数自适应先例。

### S4（仅登记引用，未读）

- Rice 1976，*The algorithm selection problem*，Advances in Computers 15:65–118（经多个片段与 S2 交叉登记；问题形式化 P/F/A/Y 四空间）。
- Smith-Miles 2009，*Cross-disciplinary perspectives on meta-learning for algorithm selection*，ACM Computing Surveys 41(1):6:1–6:25（经多个片段登记；meta-learning 综述）。
- Xu et al. 2008，SATzilla（JAIR 32:565–606，经 S2 正文提及登记）。

## 3. 对选择器设计的含义（设计建议级，不修改任何既有文档）

1. **"自适应确有收益"可以操作化为可计算的前置检查。** S2 的竞争力比 c 与潜在收益 i 提供了直接模板：在我们的开发组实例（场景族 × 档位）上，先计算 27 项候选中"最佳固定流程在多大比例实例上已最优"（下界 1−c/2）与"完美选择的收益上界 i"；**若 c/i 显示选择空间贫瘠（如某候选在绝大多数实例最优或差异微小），则按 S2 的证据标准应判定"不值得建选择器"**——这正是 AGENTS.md 要求的前置证据的量化形态。建议在将来的选择器评价协议中显式引入这两个量（或其对多算法的推广）。
2. **GPR 内部的自适应先例集中在参数级，证据形态与本仓库纪律一致。** S1/S3 表明该领域"自适应有收益"的证明 = 已知真值合成数据 + 固定基线对照 + 定量指标 + 成本并报；这与运行器设计 §B 的五类比较对象、§C 的"耗时/内存分别报"完全同构，无需改动设计，只需执行时保持该形态。
3. **"自适应"宣称需要保留审查。** 超声片段（S3）显示许多自适应方法仍隐含逐案预定义值；我们的评价协议应要求报告"自适应准则在开发组内是否稳定、是否对每个场景族都要重调"（与损伤阶梯/种子分组对齐）。
4. **算子级选择的 GPR 先例仍缺。** 检索范围内未找到"在 GPR 背景抑制算子集合上做 meta-learning 选择并证明优于最佳固定流程"的直接先例；现有证据是参数级（S1/S3）与跨领域算子级（S2，调度/SAT）。这意味着我们的选择器若建成，其"算子级选择在 GPR 背景抑制上有收益"将是一个**需要自行产出证据的命题**，S2 的 c/i 框架是产出该证据的合适工具。此缺口已同步登记（本段仅作记录，不修改证据材料原文）。

## 4. 限制

- 快速核查非系统综述；2 篇全文通读（S2 为机构库预印本，未经同行评审版核对）、4 项片段级、3 条仅登记引用；检索词覆盖英文主要组合。
- S1 的收益数值（SNR 7.55 dB 等）来自其合成设置，不可迁移为本仓库任何阈值或预期收益；S2 的 c/i 数值属于调度问题，仅框架可迁移。
- 本文不构成对运行器/选择器设计 v0.1 或传统基线研究的修改；选择器仍为纯设计，不训练、不生成标签；`training_eligible=false` 不变。
- 与仓库既有文献工作的关系：本笔记与[采集几何核查](2026-09-26_multitrace_acquisition_literature.md)、[G4 阈值校准核查](2026-09-26_g4_threshold_calibration_literature.md)平行，同属[持续研究任务书](continuous_research_brief.md)下的补充核查。

**起草标注**：本文由 **kimi** 直接调研撰写（2026-09-26，未经 codebuddy 委派），来源链接与阅读范围见[本轮台账](2026-09-26_selector_adaptivity_literature_ledger.json)；未执行仿真、未训练、未修改任何既有研究文件。
