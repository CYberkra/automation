# G4 阈值校准方法文献核查：NDT/POD 框架与 gprMax 参考校准先例（2026-09-26）

日期：2026-09-26。状态：**快速文献核查（非系统综述），为 G4（物理阈值与容差为 null）的解除提供外部方法先例**。

## 1. 目的与范围

G4 的解除条件由[评价契约 v0.2 §9]规定为四步：先用独立网格/边界/源采样敏感性**校准参考** → **损伤阶梯** → **任务容差** → 敏感性检查；阈值为 null 时可行性一律 `undetermined`（见[证据材料 G4 行](2026-09-26_operator_evaluation_sim_evidence_v0.1.md)）。本核查为其中两件方法论问题寻找外部先例：

1. "任务容差/验收阈值"如何以方法论方式设定（即阈值的**形式与推导程序**，而不是任何具体数值）？
2. FDTD/GPR 仿真"参考"的校准在文献中如何做（独立参考、源模型不确定性的处理）？

本次为**快速核查**：全文通读 2 项，摘要/片段级 4 项，仅登记引用 4 条；阅读范围逐条记录于[本轮台账](2026-09-26_g4_threshold_calibration_literature_ledger.json)。

## 2. 来源与要点

### S1（全文通读）Warren, Giannopoulos & Giannakis 2016

*gprMax: Open source software to simulate electromagnetic wave propagation for Ground Penetrating Radar*，Computer Physics Communications 209:163–170，[UWL 机构库接受稿全文](https://repository.uwl.ac.uk/id/eprint/5367/1/gprmax_new%282%29.pdf)（DOI: 10.1016/j.cpc.2016.08.020）。阅读范围：正文摘要至结论全文及参考文献表；无补充材料。

- FDTD 精度结构性约束：计算域离散必须**相对最高频率**确定（引言），时域方法单次仿真覆盖宽频带但代价是全空间离散——这与本仓库"网格依据写入规格书、收敛链校核"的做法同构。
- **定量幅度研究的门槛**：文中明确"许多仿真用理论赫兹偶极子代表真实天线仅在关心远场行为或到时信息时可接受；要研究 GPR 的定量幅度信息，必须开发并使用真实天线的精细 3D FDTD 模型"（§3.1）。这直接支持本仓库"幅度类指标在 G1/G6 未解除前只给方向"的纪律。
- 该文自身的验证路径是**引用专门的独立验证工作**：Lambot et al. 2004（全波形 GPR 建模与实测反演验证，其参考文献 23）、Tran et al. 2014（近场 GPR 建模经全波形反演验证，参考文献 13）、Warren & Giannopoulos 天线实验验证系列（参考文献 33/34/48）。即"gprMax 输出作为参考"的合法性来自**逐模型的独立验证链**，而非软件本身。
- 附带确认（与既有文献核查一致）：B-scan 的每条 A-scan 是一次独立仿真（§2.1 MPI 任务农场、§4.1）。

### S2（全文通读）nde-ed.org POD 教程（Iowa State CNDE）

*Introduction to POD / Outcome of a POD Study*，[NDE Engineering Resource Center](https://www.nde-ed.org/NDEEngineering/POD/introPOD.xhtml)。阅读范围：该教程页全文（Outcome / Two Types / Uncertainty / â vs a 示例 / POD(a) / 置信区间 / 假设检验各节）。

- **任务容差的成熟形式化**：NDT 可靠性用 POD(a) = P(信号 â > 检测门限 â_th | 缺陷尺寸 a) 表述；â vs a 数据拟合线性模型 â = β0 + β1·a + ε（ε 高斯），POD 曲线由噪声分布对门限积分得到。
- **能力表述模板 a90/95**：检出概率 90% 且 95% 置信对应的缺陷尺寸；置信区间"builds conservatism"，用于安全攸关决策——即**阈值/能力声明必须带不确定性界限，且界限取保守方向**。
- **第二产出是误报率（false call rate）**：完整评价须同时报告检出与误报，对应本仓库负控（NC/OFF）窗纪律。
- **假设回检是强制步骤**：残差随机性、Q-Q 图、方差齐性、线性性逐项检查，假设不成立则结论部分失效——对应 v0.2 §9 的"敏感性检查"第四步。
- 数据两类：â vs a（连续信号）与 hit/miss（二值），前者可经门限转后者。

### S3（片段级）MIL-HDBK-1823A 与 ASTM E2862

- [ASTM E2862-23 标准页（检索片段）](https://kpt-bj.com/astm-e2862-23-standard-practice-for-probability-of-detection-analysis-for-hit-miss-data-1957228530.html)：hit/miss 数据 POD 分析的标准实践，源自 MIL-HDBK-1823A；a90/95 定义为"90% 检出概率、95% 置信对应的缺陷尺寸"。
- MIL-HDBK-1823A（2009）本身仅登记引用（经多个片段交叉确认其存在与地位），未获正文。

### S4（片段级）MAPOD：仿真辅助 POD 的三项先例

- [ndt.net #32190（检索片段）](https://www.ndt.net/search/docs.php?id=32190)：胶接接头超声检测的 MAPOD，用 CIVA 仿真预测 0.01–5 mm 缺陷响应，建立 POD——**仿真-only/仿真为主的能力评估在 NDT 是已确立实践**。
- [ECNDT 2023 报告（检索片段）](https://www.ndt.net/article/ecndt2023/presentation/ECNDT2023_PRESENTATION_272.pdf)：工业超声检测的纯仿真 POD 流程，明确以 Berens 1989 与 MIL-HDBK-1823A 为方法基础。
- [APCNDT 2013（检索片段）](https://www.ndt.net/article/apcndt2013/papers/224.pdf)：贝叶斯框架把仿真（FE）数据作为先验、实验数据作更新——**仿真导出的 POD 需要实验锚定或明确的不确定性传递**，不能静默外推。

### S5（仅登记引用，未读）

- Berens 1989，*NDE reliability data analysis*，ASM Handbook Vol. 17:689–701（经 ECNDT 2023 片段之参考文献 1 登记）。
- Lambot et al. 2004，*Modeling of ground-penetrating radar for accurate characterization of subsurface electric properties*，IEEE TGRS 42(11):2555–2568（经 S1 参考文献 23 登记）。
- Tran et al. 2014，*Validation of near-field ground-penetrating radar modeling using full-wave inversion for soil moisture estimation*，IEEE TGRS 52(9):5483–5497（经 S1 参考文献 13 登记）。
- MIL-HDBK-1823A（2009），*Nondestructive Evaluation System Reliability Assessment*（经 S2/S3/S4 多处交叉登记）。

## 3. 对 G4 的含义（设计建议级，不解除任何缺口）

1. **v0.2 §9 的"损伤阶梯 → 任务容差"顺序有成熟外部同构。** POD 框架中"缺陷尺寸阶梯 a × 信号响应模型 â(a) × 检测门限 â_th"与本仓库"损伤阶梯 × 评价指标 × 物理阈值"一一对应；a90/95 提供了"能区分不可接受损伤的阈值与理由"的表述模板（能力值 + 保守置信界限），可直接借用作 G4 阈值文档的结构参照。
2. **仿真-only 阈值合法但须标注域。** MAPOD 先例（S4）表明仿真可承担 POD/阈值推导的主体，但文献中的稳健做法要求实验锚定或显式不确定性传递（贝叶斯更新）。这与本仓库"阈值带仿真域标注、G6/G7 未解除前不外推实测"的纪律一致；无需修改。
3. **"先校准参考"一步有 gprMax 直接先例。** S1 表明 gprMax 输出的参考资格来自逐模型独立验证链（解析/实测反演），且定量幅度断言要求真实源模型——支持 §9 第一步以独立网格/边界/源采样敏感性校准参考（E6 量级参照），并再次确认幅度类指标在 G1/G6 解除前只报方向。
4. **误报率是与检出并列的必报项**（S2）——支持负控能量上限作为独立约束而非可折抵项。
5. **新登记的缺口**：检索范围内**未找到**"针对 GPR B-scan 背景抑制算子评价的仿真 POD 阈值校准"的直接先例——POD/MAPOD 文献集中在超声、涡流等 NDT 模态，迁移到本问题需要显式域映射（缺陷尺寸 a → 损伤阶梯级别；信号 â → 所选评价指标；检测门限 â_th → 最小有效收益/损伤上限）。该映射本身应作为 G4 解除文档的一部分写明，不得隐含。此缺口已同步登记（本段仅作记录，不修改证据材料原文）。

## 4. 限制

- 快速核查非系统综述；2 项全文通读（其中 S2 为教程页非同行评审论文，但内容是 MIL-HDBK-1823A 方法的标准教学转述）、4 项片段级、4 条仅登记引用；检索词覆盖英文主要组合。
- S1 阅读基于机构库接受稿，未核对出版版排版差异；S3/S4 均未获正文，其要点仅为片段级转述。
- 本文不构成对证据材料 v0.1、运行器设计 v0.1 或 v0.2 契约的任何修改；G4 全部阈值保持 null，`physical_acceptance_threshold=null`、`training_eligible=false` 不变。
- 与仓库既有文献工作的关系：本笔记是[持续研究任务书](continuous_research_brief.md)下的补充核查，与[采集几何文献核查](2026-09-26_multitrace_acquisition_literature.md)平行，不替代任何已登记来源。

**起草标注**：本文由 **kimi** 直接调研撰写（2026-09-26，未经 codebuddy 委派），来源链接与阅读范围见[本轮台账](2026-09-26_g4_threshold_calibration_literature_ledger.json)；未执行仿真、未修改任何既有研究文件。
