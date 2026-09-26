# 2D vs 3D GPR 建模差异文献核查：差异的结构性来源与"以 3D 校 2D"先例（2026-09-26）

日期：2026-09-26。状态：**快速文献核查（非系统综述），为 P10 决策与待敲定的 A0 3D 校核契约提供外部参照**。

## 1. 目的与范围

本仓库已有 dep3d_gold_v1 的自产证据（D20m 同族场景 2D/3D 谱相关 0.935/0.958，远场后 0.967/0.979；P10 决策），并有待用户敲定的 [A0 3D 校核契约提案](2026-09-26_a0_3d_validation_contract_proposal.md)（G1：C3 族锚点尚无 3D 运行）。本核查回答两个外部问题：

1. 文献中 2D 与 3D 建模差异的**结构性来源**是什么（哪些差异是先验可预期的）？
2. "用 3D 结果评判 2D 有效域"这一方法论在 GPR 文献中是否有先例？

本次为**快速核查**：全文通读 2 项，摘要/片段级 5 项，仅登记引用 1 条；阅读范围逐条记录于[本轮台账](2026-09-26_2d3d_modeling_literature_ledger.json)。

## 2. 来源与要点

### S1（全文通读）Cheng, Schultz & Toksoz 1995

*A Comparison of Scattering from 2-D and 3-D Rough Interface*，MIT Earth Resources Laboratory 年报论文，[MIT DSpace 全文](http://dspace.mit.edu/bitstream/handle/1721.1/75257/1995.10%20Cheng_Schultz_Toksoz.pdf?sequence=1)。阅读范围：正文全文（摘要至结论及参考文献）。**域标注：声-弹界面地震散射，非 GPR/电磁**，仅结构结论可迁移。

- 用同一 3D 时域有限差分程序比较"几何继承"的 2D（沿 Y 不变）与 3D 粗糙界面：F-K 谱显示 **2D 散射前向/后向能量近似相等，3D 则前向占优、后向减少**；"2D 仿真会高估背向散射量"（摘要原文）。
- 主反射首周期 rms 幅度 2D 与 3D 相近且都小于平界面；计入更多周期（更多散射能量）后 **3D 因面外散射贡献，近法向入射处幅度大于 2D**——差异随计入的散射程数增大。
- 3D 面外散射在 F-K 谱上直接可见；作者自陈这是"第一步"，统计意义需多实现检验——单一实现对比的局限声明与本仓库"单例不升级统计结论"纪律一致。

### S2（全文通读）Giannopoulos 1997 博士论文

*The Investigation of Transmission-Line Matrix and Finite-Difference Time-Domain Methods for the Forward Problem of Ground Probing Radar*，University of York，[White Rose eTheses 全文](https://etheses.whiterose.ac.uk/id/eprint/2443/1/DX203108.pdf)。阅读范围：摘要、目录、图表清单、第 1–2 章正文；第 4–7 章仅经目录/图表清单与检索片段（§6.1 论断"2D 模型的不足在 3D 效应重要时是显然的，尤其对小型局部目标"来自检索片段，非本次正文通读范围）。

- 这是 gprMax 的奠基论文（gprMax 2D/3D 代码源头）。其方法分工即本核查问题的原型：**第 5 章 2D 模型用于"长（二维）目标"**（管线、向斜、断层、暗渠，并与 Fountains Abbey 实测暗渠数据对比）；**第 7 章 3D 模型用于"小型局部目标"**，并专门研究**离线距效应**（§7.4.1：响应随测线相对目标中心的横向偏移变化，图 7.14–7.17）——这是 2D 模型原理上无法表示的 3D 现象。
- **2D 模型的验证链**：图 4.24–4.26 以"线源在半空间上的闭式解"为独立参考验证 2D TLM/FDTD——2D 仿真的合法参考是**同维度的闭式解**，不是 3D 结果；这支持本仓库 2D 收敛链（YFINE/ZFINE/ZFINE2）作为 2D 内部校核、3D 校核另立的分层验证结构。
- 3D 章同样以闭式解（偶极子、PEC 球背向散射）验证——"逐模型独立参考"的实践与 [G4 核查](2026-09-26_g4_threshold_calibration_literature.md) S1 结论互洽。

### S3（摘要/片段级）

- **Belli, Rappaport, Zhan & Wadia-Fascetti 2009**，*Effectiveness of 2-D and 2.5-D FDTD Ground-Penetrating Radar Modeling for Bridge-Deck Deterioration Evaluated by 3-D FDTD*，IEEE TGRS 47(11):3656–3663（[ResearchGate 页](https://www.researchgate.net/publication/224440268_Effectiveness_of_2-D_and_25-D_FDTD_Ground-Penetrating_Radar_Modeling_for_Bridge-Deck_Deterioration_Evaluated_by_3-D_FDTD)，全文未获，经多个引用片段交叉确认书目）。**这是"以 3D FDTD 为裁判评估 2D/2.5D 有效域"在 GPR 中的直接先例**——方法论与 dep3d_gold/A0 校核契约完全同构；片段还提示其结论之一为"T/R 天线特性需要在仿真中考虑"（与 G4 核查 S1 的定量幅度门槛一致）。
- **FZ Jülich 学位论文**（Energie & Umwelt 643，[JuSER 记录](https://juser.fz-juelich.de/record/1032276/files/Energie_Umwelt_643.pdf)，获取失败）：片段给出理论骨架——**3D 球面扩散 A∝1/r，2D 柱面扩散 A∝1/√r**，辐射方向图不同 ⇒ 2D/3D 的**幅度衰减律先验不同**，幅度类量不可跨维直接相减。
- **Bradford 2010**（*2D GPR AVO Response to a 3D NAPL Accumulation*，[Boise State ScholarWorks](https://scholarworks.boisestate.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1044&context=cgiss_facpubs)，获取失败）：片段指出偏差可由"面外极化效应或 2D 模型未包含的非均匀性"引起。
- **Sena et al. 2008**（JGE 5(2):129，[Oxford Academic](https://academic.oup.com/jge/article/5/2/129/5127498)）：片段指出 2D FDTD 建模中源型选择问题（无限小……）。
- **Toronto 学位论文**（[Scholaris](https://utoronto.scholaris.ca/bitstreams/e198de4f-5579-4620-85e2-77a7dae1982f/download)）：片段——"2D 中的点源即 3D 中的线源，与半波偶极子不同"，故 2D 结果与 3D 仿真可显著不同。

### S4（仅登记引用，未读）

Busch et al. 2012（3D→2D 转换方案的近似误差，经 [GeoScienceWorld Geophysics 83(6):H43 片段](https://pubs.geoscienceworld.org/seg/geophysics/article/83/6/H43/566055/Radius-estimation-of-subsurface-cylindrical)登记；原文未获）。

## 3. 对 P10 与 A0 校核契约的含义（设计建议级，不修改任何既有文档）

1. **2D/3D 差异的结构性来源在文献中高度一致**：几何扩散律不同（S3-Jülich：1/√r vs 1/r）、源/辐射方向图不同（S2、S3-Toronto/Sena）、面外散射与离线距效应（S1、S2 §7.4.1、S3-Bradford）。由此**先验预期：到时/相位/谱形类量跨维可比性好于幅度类量**——这与 dep3d_gold 的自产证据方向一致（谱相关 0.935/0.958、远场后更优），也与证据材料"幅度类指标只给方向"的现行纪律一致；本次核查**未发现需要修改 P10 决策或证据材料表述的证据**。
2. **"以 3D 校 2D"有 GPR 同行评审先例**（S3-Belli 2009）：A0 校核契约的方法论结构（3D 运行为 2D 批量规格书的有效域提供裁判）不是自创，其外部合法性成立。注意其片段提示天线特性对结论的影响——A0 契约若未来要量化幅度差，源模型不确定性须单列（已在 G4 核查中登记）。
3. **2D 内部校核与 3D 校核是两层，不可互替**（S2）：2D 仿真的数值合法性由同维度独立参考（收敛链/闭式解）保证；2D→3D 的物理差异只能由 3D 运行量化。这正是本仓库 YFINE/ZFINE/ZFINE2（2D 层）与 A0 契约（3D 层）的分工，文献结构支持维持该分工。
4. **单实现对比的统计局限**（S1 作者自陈）与地震域结论的迁移界限（S1 为声-弹散射）都应写入 A0 校核报告的局限段：A0 仅 2 次运行，结论为方向性/锚点级，不得升级为统计结论。
5. **新登记的缺口**：检索范围内未找到针对本仓库场景型（厚覆盖层、深部目标、SFCW 设备链）的公开 2D/3D 差异定量基准——dep3d_gold 类自产证据仍是唯一来源，A0 校核正是补此缺口；其必要性不因文献而降低。

## 4. 限制

- 快速核查非系统综述；2 项全文通读（S1 为地震声学域、S2 第 4–7 章未经正文通读）、5 项片段级（其中 2 项获取失败仅有检索片段）、1 条仅登记引用；检索词覆盖英文主要组合。
- S1 的定量结论属声-弹散射域，迁移到电磁 GPR 只取结构方向，不取数值；S3-Jülich 的扩散律公式为片段转述，未核对上下文。
- 本文不构成对 P10 决策、dep3d_gold 报告、A0 契约提案或证据材料的修改；不解除 G1；未执行任何仿真。
- 与仓库既有文献工作的关系：本笔记与[采集几何核查](2026-09-26_multitrace_acquisition_literature.md)、[G4 阈值校准核查](2026-09-26_g4_threshold_calibration_literature.md)、[选择器自适应收益核查](2026-09-26_selector_adaptivity_literature.md)平行，同属[持续研究任务书](continuous_research_brief.md)下的补充核查。

**起草标注**：本文由 **kimi** 直接调研撰写（2026-09-26，未经 codebuddy 委派），来源链接与阅读范围见[本轮台账](2026-09-26_2d3d_modeling_literature_ledger.json)；未执行仿真、未修改任何既有研究文件。
