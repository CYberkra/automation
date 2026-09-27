# 2D 批量仿真规格书 v1.1（草案）：固化 A0/ZFINE2/FINE2 收敛与 3D 校核证据

> **codebuddy(glm-5.3-flash) 起草，kimi 验收前为草稿。**
> 草稿来源：CodeBuddy CLI（glm-5.3-flash）按任务书 `docs/research/task_codebuddy_2d_spec_v11.md` 起草；验收状态：**未验收**（kimi 逐条核验前不构成任何口径；验收不过退回重派）。
> 本文只新建此一个文件；v1.0 规格书（`docs/research/2026-09-26_batch_2d_spec_v1.md`）是已执行批次的冻结引用，**永不改写**。
> 本文不是执行授权；全局 gate 当前无任何待执行契约。文末附"待 kimi 裁决清单"。

---

## 1. 版本关系

- v1.1 在 v1.0（`docs/research/2026-09-26_batch_2d_spec_v1.md`）之上**增补**，不 retroactive 改写（【依据】任务书 §纪律）。
- 已执行批次（`batch2d_v1` BASE 29 例、`batch2d_v1_fine2` 8 例、`batch2d_v1_mt` 22 例、`batch2d_v1_co` 22 例）一律继续引用 v1.0 及各自已消耗并关闭的 gate 契约；v1.1 只取代**未来批次**的口径（【依据】`START_HERE.md` 当前状态段与执行契约状态段：上述 gate 均已 `consumed` 并关闭）。
- v1.0 遗留的"待 dep3d 证据最终确认 P10"事项，已由 [dep3d_gold 结果](2026-09-26_dep3d_gold_results.md)、[3D 校核子集决策](2026-09-26_3d_validation_subset_decision.md)与 [A0 3D 校核结果](2026-09-27_a0_3d_results.md)闭环（决策落盘 + C3 锚点 A0 已执行完成）；v1.1 §7 固化其证据含义。v1.0 原文不作回改。

## 2. 网格依据：BASE 25 mm 档的证据链

粗档 BASE（dy=25 mm / dz=25 mm）作为批量主力的依据由三条独立证据链构成（【依据】任务书 §结构 2）：

**（a）ZFINE2 收敛链（垂向细化，全频带复数，TGT−BG）**（【依据】`docs/research/2026-09-25_deep_zfine2_results.md` "收敛链"表）：

|相邻比较|全带相对 L2 变化|最大相位变化|
|---|---:|---:|
|原2.5cm均匀→FINE（均匀1.25cm）|66.6647%|64.8431°|
|FINE→ZFINE（dy12.5/dz6.25mm）|17.2452%|16.2741°|
|ZFINE→ZFINE2（dy12.5/dz3.125mm）|4.1216%|3.8847°|

最近两步比值 17.2452/4.1216≈4.18、16.2741/3.8847≈4.19，为垂向二阶行为的直接指示。**限制原文照录**："此为实测相邻差的比值趋势，未做 Richardson 外推，相邻差也不等于相对精确解的误差"；"这不声明绝对收敛，也不做任意精度无限加密"；`reference_state` 保持 `numerically_unresolved`（同文"结论与资格"节）。**收敛未认证、Richardson 未做。**

**（b）FINE2 方向一致性（BASE↔FINE2 档间，四族锚点配对）**（【依据】`docs/research/2026-09-26_batch2d_v1_fine2_results.md` "P2 方向一致性"节）：四族差分谱形状相关 **0.99813（C8）–0.99955（C1）**、逐频点符号一致率 **88.0%–93.8%**、包络峰到时 FINE2 一致早于 BASE **15.7–19.1 ns**（系统性档间偏移，各族同向）。粗档三项机制结论（带内能量比随覆盖厚度/深度递减、包络峰到时递增、电性对比增强方向差分上升）获细档方向一致支持，无机制按 v1.0 §6-P2 降级。**限制原文照录**："该支持保留为**机制/算子评价证据**，不是物理充分性认证：细档复核只判方向一致，不用于生成绝对幅度/相位结论（规格书 §6-P3），亦不得把方向一致性外推为粗档精度达标（§7.8 禁止事项）"。

**（c）A0 3D 形状同形**（【依据】`docs/research/2026-09-27_a0_3d_results.md` §3/§4）：D10m 锚点 2D/3D 归一化谱形状相关 BASE **0.9686**、FINE2 **0.9729**，与 dep3d D20m 类（0.935/0.958）同类量级且略高——支持"谱形状类指标上 2D 与 3D 的配对差分高度同形"（详见 §7）。

**综合限制（v1.0 §6-P1/P2 原文沿用）**：粗档**物理充分性未认证**；相邻差不等同相对精确解误差；"不声称绝对精度、不宣称网格收敛、不设相位或 L2 合格阈值"；粗档结论须细档子集方向一致方可保留（BASE 机制结论已获满足），结果按网格档分层报告、不跨档合并（§8）。

## 3. PML / 时窗（照录）

以下为 v1.0 §2.1/§2.2 与冻结契约 `configs/research/batch2d_v1/cases.json`（BASE）/ `configs/research/batch2d_v1_fine2/cases.json`（FINE2）的实际值（【依据】cases.json `grid`/`geometry` 字段、`2026-09-26_batch2d_v1_results.md` "h5 元数据"节、`2026-09-26_batch2d_v1_fine2_results.md` "执行记录"节）：

| 项 | BASE（粗档） | FINE2（细档） |
|---|---|---|
| 模式/不变轴 | `#domain_mode: TM`，X = `inf` | 同左 |
| 域 | `inf 32 50` m | 同左 |
| dy / dz | 25 mm / 25 mm | 12.5 mm / 3.125 mm |
| 网格 | 1×1280×2000（2,560,000 单元） | 1×2560×16000（40,960,000 单元） |
| PML 格数 | `[0,40,40,0,40,40]` | `[0,80,320,0,80,320]` |
| PML 物理厚度 | 1.0 m（各向恒定，【依据】v1.0 §2.2） | 同左 |
| PML 公式/CFS | HORIPML；`constant forward 0 0` / `constant forward 1 1` / `quartic forward 0 0.21235349838321013`（【依据】v1.0 §2.1 表，源 `DEP_BG_ZFINE2.in:8,12`） | 同左 |
| dt | 5.896635841874211e-11 s | 1.0112647039820337e-11 s |
| 时窗 | 1200 ns（1.2e-06 s），全批统一不缩短（v1.0 §2.3/§6-P8） | 同左 |
| 时间步 | 20352（实测 `iterations`，ceil 约定） | 118665（实测；契约字段 118663 为 floor 计划值，`iterations_minus_contract_time_steps=2`，【依据】fine2 结果 "iterations 约定"节原文） |
| 源/接收 | `#waveform: impulse 1 1 impulse`；`#hertzian_dipole: x inf 15.35 45 impulse`；`#rx: inf 16.65 45 measurement Ex`（【依据】v1.0 §2.1 表） | 同左 |
| 尾窗诊断 | 固定 200 / 400 ns 两档（BASE `tail_taper_fraction` 200ns=0.16666257186379047、400ns=0.33333742813620953；FINE2 200ns=0.16666175082586127、400ns=0.33332560844063913） | 同左 |
| 分析频点网格 | `linspace(20e6,170e6,501)`（501 频点） | 同左 |

时窗限制沿用 v1.0 §2.3：统一 1200 ns 不因浅部档位缩短；"浅部成组 + 更短窗"如后续提出须另立批次并单独分组。

## 4. 材料表（照录）

按 `configs/research/batch2d_v1/cases.json`（`materials` 字段）与 `configs/research/project_context_v1.json` 照录：

| 材料 | εr | σ (S/m) | 说明 |
|---|---:|---:|---|
| 基岩砂岩 | 9.0 | 0.001 | 研究假设，非场地反演值 |
| 覆盖层粉质粘土 | 16.0 | 0.01 | 研究假设，非场地反演值；场地覆盖层厚度在 project_context 中仍为 `null`（`silty_clay_cover_thickness_m`） |
| 目标（锚点"湿"） | 20.0 | 0.02 | 沿用既有 wet 定义 |
| 目标（弱档） | 12.0 | 0.005 | 弱反射/低对比 |
| 目标（强档） | 28.0 | 0.05 | 高对比上端 |
| 目标（电导单因子） | 20.0 | 0.005 / 0.02 / 0.05 | 固定 εr20，介电对比与吸收分开记账 |
| 零对比负控 NC | 9.0 | 0.001 | 与围岩完全相同 |

场景背景（【依据】`configs/research/project_context_v1.json` `user_confirmed`/`derived_frequency_grid`）：SFCW 20–170 MHz、步长 0.3 MHz，均匀含端点推导 f[k] = 20e6 + 3e5·k（k=0..500）共 501 频点；约 20 m 浅层非显性滑坡（`target_depth_range_m [0,20]`，"a research target, not demonstrated radar penetration"）；航高约 15 m、收发间距 1.3 m（横跨航线）。

**Peplinski 频段限制警告（照录）**："官方 Peplinski 土壤模型范围为 0.3–1.3 GHz，不直接作为本项目标称低频材料真值"（【依据】AGENTS.md "gprMax 训练数据"节）。本批材料为常数 εr/σ，不含弛豫谱；"常数 εr/σ 不含弛豫谱，不能据本批量结果反推场地材料谱"（v1.0 §2.4）；"σ 已含导电损耗时不得重复叠加同一 DC 电导项"（v1.0 §2.4、AGENTS.md）。若未来批次引入 Peplinski 类混合物模型，上述 0.3–1.3 GHz 适用范围警告随引用照录。

## 5. 单例成本（实测，非 ETA 承诺）

以下全部为监督器实测/分析复算值（【依据】`2026-09-26_batch2d_v1_results.md`、`2026-09-26_batch2d_v1_fine2_results.md`、`2026-09-27_a0_3d_results.md`、`2026-09-26_dep3d_gold_results.md`）。**实测非 ETA 承诺**；估算与实测的差异是 v1.0 §3.1 已知限制所述情形。

| 档 | 单例墙钟（实测） | 批量合计 | Job 峰值提交内存（实测） |
|---|---|---|---|
| BASE（25 mm，2D） | 首例 `B2D-C3m-BG` 532.156 s；其余 28 例 29.859–59.218 s | 29 例全批 1400.62 s ≈ **23.4 min** | 2,914,930,688–2,918,010,880 B（约 2.72 GiB）< 4 GiB 上限 |
| FINE2（12.5/3.125 mm，2D） | 1318.437–1323.031 s ≈ **约 1320 s/例** | 8 例串行约 10569 s | 15,246,221,312–15,248,601,088 B（约 14.2 GiB）< 20 GiB 上限 |
| 3D 5CM（各向同性） | A0 3D：BG 909.0 s / TGT 908.9 s；dep3d `DEP3D_5CM`：916.078 / 916.156 s → **约 909–916 s/例** | —（配对 2 例） | A0 7.95 GiB；dep3d 约 7.98 GiB |
| ZFINE2（12.5/3.125 mm，2D，历史参照） | BG 1438.6 s / TGT 1334.8 s | — | 15.28 GB |

- 精度/后端：各批均 CUDA double、无 CPU 回退、每例 1 次 attempt、无重试（各结果文档"执行记录"节）。
- **执行机器标注待裁决**：任务书要求标注"RTX 3060 Laptop, CUDA double"；但仓库证据显示 dim23/ZFINE2 时期（2026-09-25）执行机为 RTX 3060 Laptop（约 6 GiB 显存，【依据】`START_HERE.md` 2026-09-25 段、`docs/research/2026-09-25_2d3d_paired_plan.md:27`），而 [dep3d_gold 设计 §7.1]（`configs/research/dep3d_gold_v1/budget.json:9`）口径为"dim23 3D 反算，RTX 3060 Laptop 实测；**本机 RTX 4090 Laptop** 按 1–2 倍估计"。各批次精确机器归属待 kimi 裁决后定稿（见文末裁决清单 #1）。
- Job 提交内存不是 GPU 显存（v1.0 §3.1 限制原文）；主数组不含 PML、更新系数、CUDA 上下文与主机建模临时数组。
- 未来批次预算仍按 v1.0 §3.1 公式重算落盘，实测值只作排序与上限校准参照，不当 ETA。

## 6. 场景族与种子分组（照录）

**四族定义**（【依据】`configs/research/batch2d_v1/groups.json` `groups` 数组）：以覆盖层厚度为族键，族 = `C{h}`：

| group_id | family | cover_thickness_m | mother_model_hash（SHA-256） |
|---|---|---:|---|
| `B2D-C1m` | C1 | 1.0 | `c68b349aa31b991591c96cbf9737b0aa39be5ba949812adcabb7f47f054fca7d` |
| `B2D-C3m` | C3 | 3.0 | `92f9f5f917fd9136e1f012bc6870c7369bfbf1b0cb9b27a045151489a8fdcc4f` |
| `B2D-C5m` | C5 | 5.0 | `7d4732856d32fde982ddd6432b4c272f29de34051716f8d5e476cbd76e59f531` |
| `B2D-C8m` | C8 | 8.0 | `0e3c8d367c78cfc15c0a98e2ff28e901301ccc0437f22f6d5376d8209e10d161` |

**种子规则**（【依据】groups.json `convention.seed` 原文）："int.from_bytes(sha256('<batch_id>|<group_id>|<case_id>|<variant_tag>')[:8],'big') per spec 4.1"；同一 case_id 重跑必得同一 seed。`group_id` 为族级母模型分组键："family-level mother-model grouping key B2D-C{h}m (spec 1.0 family key + 4.2: same-family derivations share a group); family BG, targets and negative controls of one family share this group_id"。

**开发/测试划分**：开发组 **{C1,C3}** / 测试组 **{C5,C8}**，以 group_id 为单位（【依据】`configs/research/s4_baseline_group_v0.1.json` 冻结契约（SHA-256 `8676f37d…`），`START_HERE.md` S4 行原文："开发组 {B2D-C1m,B2D-C3m}/测试组 {B2D-C5m,B2D-C8m}（group_id 单位，G4 前不参与选择）"）。注意：该划分冻结于 S4 契约，`batch2d_v1/groups.json` 本身不含划分字段（见文末裁决清单 #2）。

## 7. 2D/3D 差异与适用范围（v1.1 新增核心节）

两个 3D 校核锚点：dep3d_gold（D20m 类平层，【依据】`2026-09-26_dep3d_gold_results.md`）与 a0_3d_v1（C3 族锚点 D10m，【依据】`2026-09-27_a0_3d_results.md`）。

### 7.1 已校核锚点的幅度谱同形证据（0.94–0.97 量级）

- D10m 锚点（A0）：归一化谱形状相关 BASE **0.9686** / FINE2 **0.9729**（200 ns 尾窗；400 ns 稳健性同值）。
- D20m 类（dep3d 同格距跨维）：**0.935**（B2D5CM vs DEP3D_5CM）/ **0.958**（B2DANISO vs DEP3D_ANISO）。
- 读法：两锚点支持所分析组合的**幅度谱同形**，不证明实波形NRMSE、负控能量比N_b或算子排序可跨维迁移。A0结果§6已作该项勘误；y=-x可保持幅度谱相关1却产生NRMSE2。2D可在其自身仿真域内开发，跨维评价适用性仍待同指标验证。
- 适用范围限制：本结论仅覆盖已校核的两个锚点场景（D10m/D20m 类）；"C5/C8 锚点（P-A/P-B）与电性档（P-C）的 3D 校核仍未做；本结论不得外推至其他族/其他深度档"（A0 结果 §6 原文）。

### 7.2 到时：方向两锚点相反，不可迁移

- A0（D10m）：包络峰到时 3D 晚于 2D **+2.46 ns**（BASE）/ **+18.13 ns**（FINE2）；FINE2 行与已知 BASE↔FINE2 档间偏移（细档早 15.7 ns）定量自洽（18.13 ≈ 2.46 + 15.67，A0 结果 §4 原文）。
- dep3d（D20m 类）：包络峰到时差 **−71.02 / −53.72 ns**（3D 早）。
- 读法（A0 结果 §4 原文照录）："方向与 dep3d D20m 类（3D 早）相反——不同深度/几何下到时关系不恒定，**不得跨场景迁移**"。
- 涉**到时绝对值**的结论须携带 2D/3D 差异警告（到时差 ns 级~数十 ns 量级且方向不恒定，A0 结果 §6）。

### 7.3 波形逐样本：不可跨维直接比

- 时域差分最大归一化互相关：A0 **0.2907**（BASE，滞后 +25.42 ns）/ **0.0481**（FINE2，+14.35 ns）；dep3d 同格距 **0.738 / 0.732**；3D 纯 dz 细化对 **0.18**（最优时移 21.86 ns）——全范围约 0.05–0.74，且 dep3d 的 dz 细化显示"dz 细化对时域差分波形的影响远大于对谱形状的影响"（dep3d 结果 §"3D 纯 dz 细化对照"，原文"不宣布收敛通过"）。
- 读法（A0 结果 §4 原文）："时域波形级互相关弱（0.29/0.05）：维度/扩散效应在波形层面显著，跨维度可比的稳定量是**归一化谱形状**，不是逐样本波形"。

### 7.4 绝对幅值：永不跨维比较

- 硬限制原文（A0 结果 §5 第 1 条）："不跨维度比较绝对幅值（`time_diff_peak_V_m` 仅带对内激励尺度）"。
- dep3d 补充证据：低频段 2D 与 3D 差分幅值趋势偏离最大（选频点归一化幅度趋势差最大 17.20 / 17.43 dB @ 20 MHz，170 MHz 处收窄至 −0.19 / −0.00 dB）；`cross_dimension_absolute_amplitude_compared: false`。

### 7.5 远场变换：仅探索性

- `D_2D_equiv(f) ∝ D_3D(f) / √(i·2πf)`（常数省略，只比形状）：A0 变换后 **0.9860 / 0.9877**（变换前 0.9686/0.9729）；dep3d 变换后 0.967 / 0.979（变换前 0.935/0.958）。
- 限制原文照录（dep3d 结果）："近场区不适用、常数未定、只比形状；仅作探索性证据单列，不作为结论依据（`far_field_transform_is_conclusion_basis: false`）"。

### 7.6 混淆因素未分离与 3D 校核缺口

- A0 四类混淆因素（维度效应/侧向域宽 32 m vs 10 m/目标 y 位置/格距未对齐）**均未分离**（A0 结果 §4 原文："本分析不声称排除其中任何一个"）；dep3d 三类混淆（维度/侧向域宽/目标 y 位置）"两对差值之间不得再相减当作可加分解"。
- **仍未覆盖的 3D 校核缺口**：C5 锚点（P-B）、C8 锚点（P-A）、电性/电导档（P-C）的 3D 校核未做（【依据】`2026-09-26_3d_validation_subset_decision.md` 依赖链、A0 结果 §6）。

## 8. 不变限制（照录，适用于未来批次）

以下为各已执行批次 `results.json` `hard_limits` 的不变字段，未来批次继续沿用（【依据】`2026-09-26_batch2d_v1_results.md`、`2026-09-26_batch2d_v1_fine2_results.md`、`2026-09-27_a0_3d_results.md` "硬限制"节）：

1. `reference_state = "numerically_unresolved"`——网格收敛未认证（`grid_convergence_certified: false`），20 m 全频带物理真值未认证；
2. 无物理/训练标签：`physical_label_eligible: false`、`training_eligible: false`、`training_labels_generated: false`、`clean_truth_generated: false`；
3. 阈值 null：`physical_acceptance_threshold: null`、`direction_consistency_threshold: null`（"none declared: the numbers are reported as directions/agreement rates only"）；
4. 结论限仿真域：`conclusions_scope` 原文——"mechanism and operator-evidence ... not field performance, not 3D, not absolute accuracy, no maximum-depth claim"（batch2d_v1 版）；"limited to the dep3d_gold_v1 scenario family and the grids analysed; not extrapolated to field data"（dep3d 版）；A0 版限定于"a0_3d_v1 场景族与所分析网格"；
5. `in_band_energy_ratio_is` 原文："diagnostic energy ratio on the 501-point analysis grid; not detectability, not SNR, not a physical threshold"；
6. `fdtd_solver_invoked: false`（分析脚本仅做后处理复算）；
7. 不跨族差分（`cross_family_differencing: false`）、不跨档差分/合并（`cross_tier_differencing: false`、`grid_tiers_merged: false`）、负控不并入目标行（`negative_controls_merged_into_targets: false`）；
8. 相位只做线性拟合形式量不跨维迁移：`phase_not_transferable_across_dimensions: true`，"the equivalent time offset is a formal quantity, not a physical delay"。

## 9. 批量跑启动条件清单（未来新批次）

未来任何新批次启动前，以下条件必须全部满足（【依据】v1.0 §5、`configs/research/gprmax_v4_execution_gate.json` 机制、`START_HERE.md` 执行契约状态段、`2026-09-26_speedup_probe_results.md` 结论——见 `START_HERE.md` 引用）：

1. **冻结契约**：新批次须重新冻结具体契约（`run_id` 清单、输入 SHA-256、资源上限），不重跑已消耗 attempt（AGENTS.md；v1.0 §5.2）；契约字段沿用 `packet_id`、`input_path/input_sha256`、`runtime_identity_path/sha256`、`launcher_sha256`、`supervisor_sha256`、`attempt_record`、`run_directory`、`cuda_identity_path/sha256`、`run_id`，批量形态用 `approved_execution_contracts`(list) + `batch_id`；`max_fdtd_runs`、`retries = 0`、`wall_minutes/job_commit_GiB/output_GiB` 三项预算 > 0。
2. **gate 激活留痕**：`approved_to_simulate` 从 `false` 起步；激活（含 2026-09-26 18:48 用户授权"仿真启动可不经用户逐次批准"的范围）须留痕——该授权免除的是"启动前逐次问用户"，**契约纪律不变**（每 attempt 仅一次、目录不存在语义、墙钟/内存上限、gate 留痕与关闭照旧）；执行后回填 `execution_outcome`（status / attempt_consumed / exit_code / wall_s / peak_job_commit_bytes）并关闭 gate。gate/契约/授权 JSON 的冻结与激活不属委派起草范围（多 agent 协作纪律）。
3. **attempt 单次**：每例单次 attempt、失败即停、不自动重试、不重跑已消耗 attempt；运行目录"不存在"语义照旧。
4. **监督器**：全程监督器实测（`supervision.json`：墙钟、exit、Job 峰值提交内存、进程计数）；执行方式按加速实测结论——未来批次默认**单进程多例**启动（吞吐 ~2.8×）；FP32 仅可作独立新档、不得跨精度配对；单卡多进程并发已否决。
5. **静态检查先行**：逐例 `.in` 生成后先做静态核验（域/格距/单元数/dt/时窗/PML/材料/源逐项自算，含 v1.0 §1.5 网格整除与 PML 余量单列项、`target_fully_in_sandstone` 断言），再写 budget（v1.0 §3.1 公式重算落盘）。
6. **分组与种子**：`groups.json`（group_id=族级母模型键）与 seed 规则随契约冻结；同母模型全部派生同组，任何跨组划分只能以 group_id 为单位；开发/测试划分 {C1,C3}/{C5,C8} 以 S4 冻结契约为准，测试组在 G4 解除前不参与选择。
7. **网格档分层**：结果标注所属网格档（BASE / FINE2 / 未来新档），跨档结果不得直接相减、合并报告或把方向一致性外推为精度达标（v1.0 §7.8）。
8. **分析验收**：分析脚本 r1/r2 两次连跑 `results.json` 逐字节一致；硬限制字段（§8）逐条照录；负控（族 BG/NC/OFF）单列不并入目标行。
9. **归档与登记**：证据归档 `artifacts/research_checks/`（SHA-256 清单入 `results.json` `inputs`）；registry / START_HERE / decision_log 更新后提交推送；vctip 等操作层处置（如触发）逐事件留痕。

---

## 待 kimi 裁决清单

1. **执行机器标注**：任务书要求 §5 标注"RTX 3060 Laptop, CUDA double"；但仓库证据显示 RTX 3060 Laptop 为 dim23/ZFINE2 时期执行机（2026-09-25），dep3d_gold 预算口径（`configs/research/dep3d_gold_v1/budget.json:9`）称"本机 RTX 4090 Laptop"。§5 现按"两台机器证据并列 + 待裁决"书写，定稿时如何标注（按批次逐一归属/统一口径）请裁决。
2. **开发/测试划分出处**：任务书要求"按 groups.json 照录…开发/测试划分（{C1,C3}/{C5,C8}）"，但 `configs/research/batch2d_v1/groups.json` 不含划分字段；划分实际冻结于 `configs/research/s4_baseline_group_v0.1.json`（START_HERE 记载）。§6 现按后者引用并显式说明，请确认。
3. **数值精度口径**：任务书摘要值为舍入值（远场 0.986/0.988 vs 源文件 0.9860/0.9877；波形互相关"0.05–0.74"为 A0 0.0481–dep3d 0.738 的跨锚点范围）。草稿一律采用源文件全精度值并注明范围口径，请确认。
4. **§9 是否引用 2026-09-26 18:48 用户授权**：草稿已将其写入启动条件（免逐次批准、契约纪律不变），表述是否准确请裁决。
5. **范围边界**：任务书未列 MT/CO/事件表/G4 等后续单元；草稿仅在必要处引用（§6 划分、§9 监督器/加速结论），未将其纳入 v1.1 口径，请确认是否需增补。
