# 多道采集几何文献核查：变偏移距道集 vs 常偏移 B-scan（2026-09-26）

日期：2026-09-26。状态：**快速文献核查（非系统综述），为待敲定决策提供外部证据**。

## 1. 目的与范围

为用户敲定[多道场景族批次设计提案 v0.1](2026-09-26_multitrace_batch_proposal.md)的核心取舍——方案 A（单发多收、变偏移距，1 次运行/母模型）vs 方案 B（常偏移多发多收，33 次运行/母模型）——核查公开文献中的相关实践。问题具体化为两点：

1. 变偏移距道集（CSG，单固定源 + 接收阵列）是否是评价处理算法的合法数据形态？
2. 背景抑制类算子（沿道均值、SVD 低秩等，即本仓库 27 项候选的文献来源）的文献建立在哪种采集几何上？

本次为**快速核查**：正文通读 2 篇，摘要/片段级 4 项，未做系统综述；阅读范围逐条记录于[本轮台账](2026-09-26_multitrace_acquisition_literature_ledger.json)。

## 2. 来源与要点

### S1（正文通读）Roncoroni, Koyan, Forte, Tronicke & Pipan 2025

*A realistic 2D multi-offset, multi-frequency synthetic GPR data set as a benchmark for testing new algorithms*，Scientific Data 11，[PMC11802766](https://pmc.ncbi.nlm.nih.gov/articles/PMC11802766/)（DOI: 10.1038/s41597-024-04300-1）。阅读范围：正文 Background & Summary / Methods / Data Records / Technical Validation 全文；补充材料未读；60 条参考文献未逐条追踪。

- 该文是**同行评审的合成基准数据集**，明确目的即"evaluate and test processing, analysis and inversion techniques"，其采集几何为**多偏移距 CSG**：每剖面 161 个固定接收点 + 161 个移动源点（0.1 m 步进），gprMax v3.1.7，50/100/200 MHz。即"固定接收阵列 + 变源"（与我们方案 A 的"固定源 + 变接收阵列"互为互易几何）。
- 该文明确区分 **CO（common-offset）剖面**与 **MO（multi-offset）道集**，并指出 MO 因采集耗时"not a standard technique"，现今 GPR 数据仍以 CO 为主。
- 其基准面向的处理/反演类别是**叠前深度偏移、AVO、FWI**——不是 B-scan 背景抑制。
- 数据与代码公开（figshare 10.6084/m9.figshare.26508976.v1；GitHub `Giacomo-Roncoroni/MO-GPR_data`）。**本次不下载、不引入**，仅登记为将来可能需要外部基准对照时的候选来源。

### S2（正文通读）Picchi & Brell-Çokcan 2022

*A modified common midpoint approach for GPR radars*，Construction Robotics 6:319–328，[Springer 全文 PDF](https://link.springer.com/content/pdf/10.1007/s41693-022-00086-z.pdf)（DOI: 10.1007/s41693-022-00086-z）。阅读范围：§1–§4（方法、推导、仿真设置与对比）；§5 之后的结论部分只经页面提取片段，未逐句核对。

- 用 gprMax 仿真直接对比**常规 CMP（变偏移距）**与**固定偏移距设备**两种几何：CMP 道集上目标反射到时随偏移距增大而增大（双曲/moveout），固定偏移距设备在目标正上方到时最小——两种几何的事件形态**结构性不同**。
- 该文为变偏移距与常偏移几何到时关系给出了显式方程（其式 1/2：$vt_n=2\sqrt{(x-x_0)^2+(vt_0/2)^2}$ 及含固定偏移 a 的修正形式），与我们提案 §3.1/§2.4 的自算到时方法是同一类几何论证。
- 与本提案直接相关的证据点：**变偏移距几何下平层/目标事件天然带 moveout**（不是逐道对齐的平事件），这正是提案 §2.4"平层=秩一在本几何下不严格成立"论断的外部印证。

### S3（摘要级）背景抑制文献的采集几何前提

- Rashed 2014（*A novel background removal algorithm for GPR data*，JAG 106:154，[ADS 摘要](https://ui.adsabs.harvard.edu/abs/2014JAG...106..154R/abstract)）：BMS（背景矩阵扣除）针对"horizontal background noise"——**水平相干性假设**。
- *Critical Analysis of Background Subtraction Techniques on Real GPR Data*（[CORE 全文 PDF](https://core.ac.uk/download/pdf/333722781.pdf)，本次只读检索片段与公式区）：SVD 背景扣除把 B-scan 矩阵分解为目标分量与杂波分量，其成立前提同样是 B-scan 中背景事件**沿道水平相干**。
- arXiv 2305.18775 与 DR-NTU 论文（检索片段）：常规预处理链"time-zero → background removal using SVD"同样作用于常规 B-scan。
- **结论（摘要级证据）**：均值扣除/SVD/BMS 这一族背景抑制方法的文献形态全部建立在**常偏移 B-scan**上——常偏移下直耦波与平层界面逐道对齐，"背景 = 低秩/水平分量"才近似成立。

### S4（官方文档，片段）gprMax 的 B-scan 生成方式

[gprMax 官方文档 examples_simple_2D](https://docs.gprmax.com/en/latest/examples_simple_2D.html)（检索片段）：物理天线在场景中移动时**每个位置需要一次新的几何构建（新仿真）**，再用 `outputfiles_merge.py` 合并成 B-scan。这与提案 §2.4-B 的成本结构（常偏移 = O(n_traces) 次运行/母模型）一致。

### S5（仅登记引用，未读）

Forte & Pipan 2017，*Review of multi offset GPR applications: data acquisition, processing and analysis*，Signal Processing 132C:210–220（S1 之参考文献 4，DOI: 10.1016/j.sigpro.2016.04.011）。多偏移距 GPR 综述；本次未获取正文，仅登记存在性，不引用其具体结论。

## 3. 对待敲定决策的含义（设计建议，不是结论升级）

1. **方案 A 的数据形态合法。** CSG（单固定源 + 接收阵列）是有同行评审合成基准先例的合法采集几何（S1）；用它解除 G2 的"形态层"阻塞（让 `[sample, trace]` 输入存在、27 项候选通过 `min(shape)≥2` 门禁）不属于伪造数据形态。
2. **几何标注是必要的、且目前充分。** 27 项候选的文献来源（S3）全部建立在常偏移 B-scan 的"背景水平相干"前提上；在变偏移距道集上评价它们时，直耦波与平层界面带 moveout（S2），提案 §2.4 的固有限制声明（结论须带变偏移距几何标注、不得直接外推为实测常偏移 B-scan 性能）与文献证据一致。本次核查**未发现需要推翻"A 先行、B 待定"决策的证据**。
3. **方案 B 的触发条件维持不变。** 仅当 A 批次分析表明变偏移距可能改变方向性结论时，才按提案 §2.4-B 另立小规模常偏移配对子集（建议 C3-BG 与 C3-A0 两母模型）；其"每位置一次新仿真"的成本结构有官方文档依据（S4）。
4. **新登记的缺口**：在检索范围内**未找到**专门研究"在变偏移距道集上评价 B-scan 背景抑制算子（mean/SVD/BMS）有效性"的文献——该具体问题无直接先例。这意味着本批多道结果若用于算子评价，其"变偏移距"标注不是形式主义，而是真实的证据边界；若后续需要强结论，方案 B 或 CO 派生（如从 CSG 重排 CMP/CO 道集）将不可避免。此缺口已同步登记于证据材料 G2 的解除条件解释中（本段仅作记录，不修改证据材料原文）。

## 4. 限制

- 快速核查非系统综述；两篇文章正文通读，其余为摘要/片段级；检索词覆盖英文主要组合，未覆盖德文/法文文献（HAL 法语论文仅命中片段未读）。
- S1/S2 的阅读基于出版方页面全文，未核对补充材料与全部参考文献。
- 本文不构成对多道提案的修改；若用户敲定多道批次，本笔记作为决策参考附件一并生效。
- 与仓库既有文献工作的关系：本笔记是[跨领域 ML 证据](2026-09-23_cross_domain_ml_evidence.md)与[持续研究任务书](continuous_research_brief.md)下的补充核查，不替代其中任何已登记来源。

**起草标注**：本文由 **kimi** 直接调研撰写（2026-09-26，未经 codebuddy 委派），来源链接与阅读范围见[本轮台账](2026-09-26_multitrace_acquisition_literature_ledger.json)；未执行仿真、未修改任何既有研究文件。
