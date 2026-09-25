# A0 三维校核契约提案（C3 族锚点 D10m，5 cm 各向同性）：待用户敲定

日期：2026-09-26。状态：**契约提案（初版草稿），不是执行授权**。

本文响应[三维校核子集决策](2026-09-26_3d_validation_subset_decision.md)的决策 D3（首选 3D 校核对象 = C3 族锚点 A0 配对）与 D6（3D 校核批次单独冻结契约、须经用户敲定），为对象 `B2D-C3m-D10m-W4m-T0.5m-E20-S0.02` 的三维对应物起草**待敲定的契约提案**，供用户逐项确认后再按[批量 2D 规格书](2026-09-26_batch_2d_spec_v1.md) §5 冻结为正式契约 JSON。

标注约定（沿用决策文档）：

- 【依据】仓库既有文档/契约/结果中的既有事实或既定参数；
- 【自算】由已引实测值或既有公式做的算术，不是新的测量；
- 【待敲定】本提案提出、须由用户确认后方可写入契约的取值；
- 【解读】方向性/机制性说明，不构成结论升级。

---

## 0. 状态声明（先读）

1. **本文不是执行授权。** 本文只起草契约内容，不启动、不排队、不预生成任何 FDTD 运行。
2. **全局 gate 保持关闭态不动。** `configs/research/gprmax_v4_execution_gate.json` 当前 `approved_to_simulate: false`、`approved_run_ids: []`、`batch_id` 为 `batch2d_v1_fine2_consumed`；本文**不读取修改、不回填、不新增** `approved_execution_contracts`。
3. **本文不生成输入文件、不生成契约 JSON、不修改任何既有文件。** 输入 `.in`、`budget.json`、`cases.json`、`groups.json`、attempt 记录一律在用户敲定后按 §8 流程另行生成并哈希登记。
4. **attempt 未被消费。** 本文所述 2 次运行各自的 1 次 attempt 在契约冻结前处于未消费状态；已消耗 attempt（含 dep3d_gold_v1 的 8 次、batch2d_v1 粗档 29 例、fine2 段 8 例）一律不复用。
5. **用户敲定是冻结的前置条件**；未经敲定，本文任何取值都不得写入 gate 或被执行。

---

## 1. 校核对象

### 1.1 对象定义

| 项 | 取值 |
|---|---|
| 母模型（2D 侧） | `B2D-C3m-D10m-W4m-T0.5m-E20-S0.02`（全批锚点 A0）【依据】规格书 §1.0、决策文档 D3 |
| 族 | C3（覆盖层 3 m，覆盖层底 z = 27、地表 z = 30）【依据】规格书 §1.1、fine2 `cases.json` `cover_bottom_z_m: 27.0` |
| 目标（2D 侧） | y 14–18 m（宽 4 m）、z 19.75–20.25 m（厚 0.5 m、中心深 10 m）、εr20 / σ0.02 S/m【依据】`configs/research/batch2d_v1_fine2/cases.json` `target_box_m: [14.0, 19.75, 18.0, 20.25]` |
| 3D 对象 | 上述母模型的三维对应物，BG + TGT 各 1 例，**共 2 次运行**，5 cm 各向同性档【依据】决策文档 D2/D3 |
| 运行数 / attempt | 2 次 FDTD 运行，各 1 次 attempt，串行，`retries = 0`【依据】决策文档 D3、dep3d 设计书 §3 |

### 1.2 与 dep3d_gold_v1 5 cm 档的对齐项（逐项照录）

对齐基准：`configs/research/dep3d_gold_v1/DEP3D_5CM_BG.in` / `DEP3D_5CM_TGT.in`（[设计书](2026-09-26_dep3d_gold_design.md) §6.1/§6.2；实测已归档，[结果](2026-09-26_dep3d_gold_results.md)）。

| 项 | dep3d 5 cm 档取值 | 本提案取值 | 性质 |
|---|---|---|---|
| 维度 | 3D（无 `#domain_mode`） | 3D（无 `#domain_mode`） | 对齐【依据】dep3d §6.1、§9.1 末条 |
| 域 | 8 × 10 × 50 m | 8 × 10 × 50 m | 对齐【依据】dep3d 决策 2 |
| 格距 | 0.05 / 0.05 / 0.05 m（160×200×1000 = 3 200 万单元） | 同 | 对齐【依据】dep3d §7.1 |
| 时窗 | 1200 ns | 1200 ns | 对齐【依据】规格书 §2.3（全批统一 1200 ns，不缩短） |
| 线程 | `#omp_threads: 8` | 同 | 对齐【依据】dep3d §6.1 |
| PML | `HORIPML`，`pml_cells: 20 20 20 20 20 20`（有限方向各 1 m），`pml_cfs` 与 DEP 2D 链逐字相同 | 同 | 对齐【依据】dep3d §6.1、§8.6 |
| 地表 / 覆盖层 | z = 30；cover `16 0.01 1 0`，box z 27–30 | 同 | 对齐【依据】C3 族覆盖层 3 m 与 DEP 族一致 |
| 基岩 | rock `9 0.001 1 0`，box z 0–30 | 同 | 对齐【依据】dep3d §6.1 |
| 激励 | `#waveform: impulse 1 1 impulse` | 同 | 对齐【依据】dep3d §6.1 |
| 源 | `#hertzian_dipole: x 4 4.35 45 impulse` | 同 | 对齐【依据】dep3d §4（y 向重定位规则） |
| 接收 | `#rx: 4 5.65 45 measurement Ex` | 同 | 对齐【依据】同上；分量 Ex 与 YZ 链一致（规格书 §2.1） |
| 目标电性 | `20 0.02 1 0` | 同 | 对齐【依据】A0 锚点电性 E20-S0.02 |
| 目标 x 范围 | x 3–5 m（2 m 有限、居中；Tx/Rx 在 x = 4） | 同 | 对齐【依据】dep3d §4 |
| 后处理频率表 | 501 频点 20–170 MHz（0.3 MHz 步长），官方 SFCW 后处理 | 同 | 对齐【依据】dep3d §2、规格书 §2.1 |
| 尾窗 | 主 200 ns / 稳健 400 ns | 同 | 对齐【依据】DEP 链固定约定 |

### 1.3 与 dep3d 5 cm 档的差异点（逐项列理由）

| # | 差异项 | dep3d 5 cm 档 | 本提案 | 理由 |
|---|---|---|---|---|
| 1 | **目标 z 范围（对应深度档）** | z 9.75–10.25（中心深 **20 m**） | z **19.75–20.25**（中心深 **10 m**） | 本提案对象是 A0 = **D10m** 锚点；dep3d 的 3D 算例几何是 C3 族的 **D20m** 类（决策文档 §1.2 末条）。这正是 D3 要回答的问题——"深度档改变是否显著改变维度差"，因此该差异是**被研究的量**，不是误差，不得通过改几何消除【依据】决策文档 §3.1、dep3d §2（"维度差异是被测量"的同款处理） |
| 2 | 目标 y 范围的来源 | y 3–7（由 2D y 14–18 按 `y_3D = y_2D − 11` 重定位） | 同值 y 3–7 | 2D 侧目标 y 14–18 与 dep3d 的 DEP 链完全相同，故套用同一中心锚定平移规则得同一 3D 值【依据】dep3d 决策 1、§4；本提案 2D 侧目标 y 见 §1.1 |
| 3 | 场景族归属与命名 | DEP 场景族（`DEP3D_5CM_*`） | C3 族锚点 A0（`run_id` 沿用 `B2D-C3m-*` 母模型 ID + 3D 档后缀） | 保证与 batch2d_v1 的母模型/分组键一致（规格书 §4.2：同一母模型的二维结果与其三维校核子集必须同组）【依据】规格书 §4.2、§6-P10 |
| 4 | 实波形诊断窗 | 440–600 ns（围绕 D20m 自算双程到时 520 ns） | **待敲定**：建议按同一 ±80 ns 规则围绕 D10m 自算双程到时取 240–400 ns | 目标深度改变后原诊断窗不再对准目标到时；须按同一规则重算【自算】 |
| 5 | 单元数上限 / Job 上限的"批次专属"属性 | dep3d 批次专属冻结值（3D 单元上限 6 400 万、Job 20 GiB、墙钟 40 min/例） | 同量级沿用，**但须重新敲定** | 决策文档 §4 明确：dep3d 的批次专属上限不对后续批次自动生效；本提案不预设其继承关系【依据】决策文档 §4、dep3d §7.3 |

**自算双程到时（D10m，供诊断窗与到时方向对照）**：空气段 30 m/c = 100.1 ns + 覆盖层双程 6 m/(c/4) = 80.0 ns + 砂岩双程 14 m/(c/3) = 140.1 ns ≈ **320.2 ns**（目标中心）【自算，c = 299 792 458 m/s，公式与 dep3d §9.3 的 D20m 自算同法；dep3d 同法得 520 ns】。该值仅作诊断窗与方向对照的参考，不作物理时延标定结论。

### 1.4 自算网格对齐与余量（冻结前须由静态检查脚本复核，本文不代替）

以下按 0.05 m 各向同性格距自算，坐标均整除【自算，方法同 dep3d §9.1】：

| 坐标 | 值 (m) | ÷0.05 格索引 | 结论 |
|---|---:|---:|---|
| 域 x / y / z | 8 / 10 / 50 | 160 / 200 / 1000 | 整除，与 dep3d 5 cm 例同网格 |
| 地表 z | 30 | 600 | 整除 |
| 覆盖层底 / 顶 z | 27 / 30 | 540 / 600 | 整除 |
| 天线 z（Tx/Rx 同高） | 45 | 900 | 整除 |
| 目标 z | 19.75 / 20.25 | 395 / 405 | 整除（与 dep3d 的 195/205 同为 10 格厚） |
| 目标 y（重定位后） | 3 / 7 | 60 / 140 | 整除，与 dep3d 一致 |
| 天线 y | 4.35 / 5.65 | 87 / 113 | 整除，与 dep3d 一致 |
| 目标 x / Tx-Rx x | 3 / 5、4 | 60 / 100、80 | 整除，与 dep3d 一致 |
| PML 厚度（有限方向） | 1 | 20 | 整除，`2×厚度 < 该轴格数` 成立 |

余量（自算）：底 PML 内边 z = 1 m，目标下缘 z = 19.75 m → 余量 **18.75 m**；侧向 PML 内边 x = 1 / 7 m（域 x 0–8）→ 目标 x 侧边余量 **2 m**，y 侧边（y = 1 / 9，域 y 0–10）→ 目标 y 侧边余量 **2 m**；顶 PML 内边 z = 49 m，天线 z = 45 m → 余量 4 m，与 dep3d §5 一致【自算，口径同 dep3d §5】。**x 向 2 m 余量沿用 dep3d 的已知横向截断代价，不在本提案中扩域**【依据】dep3d §9.2 决策 2 与其理由 (iii)。

---

## 2. 与已完成 2D 侧的比较设计

### 2.1 比较结构

- **各维内部先做配对差分**：3D 侧 `TGT − BG`（同域同格同源）；2D 侧沿用 batch2d_v1 的"同族同格 BG"差分规则（规格书 §6-P7）。
- **跨维度只比形状**，沿用 dep3d §9.3 的指标集与实现约定【依据】dep3d §9.3、决策文档 §3.1】：

| # | 比较量 | 说明 |
|---|---|---|
| 1 | **归一化谱形状相关** | 各侧内部按自身最大值归一化后求相关（本次比较的主量） |
| 2 | **包络峰到时方向** | 只取方向（3D 相对 2D 更早/更晚），跨档/跨维只报差值(ns)，不取偏移量作误差 |
| 3 | 谱峰频率差 (MHz) | 报数值，不设通过判读 |
| 4 | 时域差分波形最大归一化互相关 + 最优时移 (ns) | 沿用 dep3d 指标定义 |
| 5 | 选频点（20/50/80/110/140/170 MHz）归一化幅度趋势差 (dB) | 沿用 dep3d 指标定义 |
| 6 | 相位差对频率线性拟合：等效时移（形式量）/ 残差 std (°) | 等效时移**不是**物理时延，不跨维迁移【依据】dep3d §9.3 |
| 7 | 远场 3D→2D 变换后的谱形相关 | **探索性证据单列**，`far_field_transform_is_conclusion_basis: false`；本提案列为可选项，是否计算由用户敲定 |

### 2.2 明确不做（逐条写入契约与结果文件）

1. **不做跨维绝对幅值比较**（源归一化与几何扩散不同）：结果文件保留 `cross_dimension_absolute_amplitude_compared: false`【依据】dep3d §9.3、规格书 §7.4】。
2. **不跨格相减**：3D 侧为 5 cm 各向同性；2D 粗档为 dy25/dz25 mm、细档为 dy12.5/dz3.125 mm，三者网格互不相同。故本次比较只在"各自内部配对差分的形状"层面进行，**不得把 3D 与 2D 的差当作可加的网格/维度分解**【依据】决策文档 §3.1、§6；规格书 §7.8】。
3. **不跨档合并、不跨档相减**：2D 粗档（BASE）与细档（FINE2）结果分层报告；本次 3D 结果与两档 2D 结果分别成对照，两对照之间不得再相减【依据】规格书 §6-P1/P2、§7.8】。
4. **不跨族差分**：只与 C3 族 A0 的 2D 结果比较【依据】规格书 §6-P7】。
5. **不做收敛认证、不设通过阈值**：`grid_convergence_certified: false`、`direction_consistency_threshold: null`、`physical_acceptance_threshold: null`【依据】dep3d §9.2 与 fine2 结果 hard_limits】。
6. **不生成任何标签**：`training_labels_generated: false`、`clean_truth_generated: false`、`physical_label_eligible: false`、`training_eligible: false`。

### 2.3 既有 2D 侧参照值（比较基准，只作量级参照）

| 侧 | 带内能量比 DIFF/BG (dB) | 包络峰到时 (ns) | 谱峰 (MHz) | 来源 |
|---|---:|---:|---:|---|
| 2D 粗档 BASE（A0） | 2.8635252336237045e-07 | 331.8626651806806 | 169.7 | 【依据】决策文档 §3.1（batch2d 粗档结果） |
| 2D 细档 FINE2（A0） | 3.092143334046171e-07（−65.09740382883551） | 316.1921349940624 | 170.0 | 【依据】[fine2 结果](2026-09-26_batch2d_v1_fine2_results.md)细档锚点表 |
| 2D 粗→细方向一致性（C3） | 差分谱形状相关 0.9994178181733713；逐频点符号一致率 92.22%；包络峰到时 FINE2−BASE = −15.67053018661818 ns | | | 【依据】fine2 结果 P2 方向一致性表 |

dep3d 在 **D20m 类**实测的跨维量级（本提案的可互引参照，**不是本提案的预期值**）：归一化谱形状相关 0.935（B2D5CM vs DEP3D_5CM）/ 0.958（B2DANISO vs DEP3D_ANISO），远场变换后 0.967 / 0.979；包络峰到时差 −71.02 / −53.72 ns；谱峰频差 −5.70 / −0.30 MHz；时域最大归一化互相关 0.738（时移 0.94 ns）/ 0.732（0.37 ns）；选频点归一化幅度趋势差最大 17.20 / 17.43 dB（均 @20 MHz）【依据】dep3d 结果"2D/3D 同格距跨维比较"】。本提案预期回答的问题（照录决策文档 §3.1）：**锚点深度 D10m 处的 2D→3D 形状差方向与量级，是否与上述 D20m 处的量级处于可互引的同类量级**——即"深度档改变是否显著改变维度差"。本文不预设该问题的答案。

---

## 3. 资源预算（只引 dep3d 实测，不是 ETA 承诺）

**口径声明**：以下全部引用 dep3d_gold_v1 的实测值；dep3d 实测仅覆盖 8×10×50 m 域、5 cm 与 dz 2.5 cm 两档、1200 ns 时窗、本机 RTX 4090 Laptop 的 CUDA double 配置；其他格距/域/机器无任何实测成本依据，**不得外推**【依据】决策文档 §4 口径声明】。

| 量 | 取值 | 来源 |
|---|---|---|
| 单次 3D 运行墙钟（5 cm 档实测） | 916.078 s（BG）/ 916.156 s（TGT）；全批包络 910.6–920.2 s/run | 【依据】dep3d 结果执行记录表（各 run `supervision.json`） |
| 单次 3D 运行 Job 峰值提交内存（5 cm 档实测） | 8,567,218,176 B / 8,565,354,496 B（约 **7.98 GiB**） | 【依据】同上 |
| 求解耗时自报（5 cm 档） | 3 min 21.5987 s / 3 min 21.6987 s；gprMax 进程总时长 3 min 35.8739 s / 3 min 36.0577 s | 【依据】同上 |
| 本提案 2 次运行墙钟包络 | 2 × 约 916 s ≈ **31 min** | 【自算，引决策文档 §4 同款算术】 |
| 本提案 Job 峰值包络 | 约 **7.98 GiB/例**，低于 dep3d 冻结的 20 GiB（21,474,836,480 B）上限 | 【依据】dep3d 结果：全批实测 2,322,751,488–14,979,588,096 B，未触及上限 |
| 3D 单元数 | 3 200 万（160×200×1000） | 【依据】dep3d §7.1 |

**不得下调预算的理由（照录）**：dep3d 8 例墙钟远大于求解耗时，差值（派生量）为 2D 约 897.6–900.0 s、3D 5 cm 约 700 s、3D ANISO 约 338 s，根因在同环境鉴定为 `vctip.exe`（MSVC 编译器遥测）孤儿滞留；batch2d_v1 粗档启用 vctip 处置后 2D 单例墙钟降至 29.9–59.2 s，**但该处置对 3D 运行墙钟的影响没有任何实测记录**，因此 3D 墙钟估算一律沿用 dep3d 的 910.6–920.2 s/run 实测包络，不得据此下调【依据】决策文档 §1.1、dep3d 结果"墙钟构成"与根因补记】。

**上限取值（待敲定，建议沿用 dep3d §7.3 本机值）**：

| 项 | 建议上限 | 出处 |
|---|---|---|
| 墙钟（每例） | 40 min | 【依据】dep3d §7.3（用户指定，同 `deep_zfine2_v1` gate）；实测 916 s 未触及 |
| Job 提交内存（每例） | 20 GiB | 【依据】dep3d 决策 4、§7.3 |
| 预检空闲 RAM | ≥ 24 GiB（以执行当时实测为准） | 【依据】dep3d 决策 8；低于门槛则停止，不降级运行 |
| 预检空闲 VRAM | ≥ 8 GiB | 【依据】dep3d §7.3 |
| 输出（每例） | 1 GiB | 【依据】dep3d §7.3；DEP 2D 完成例实测输出 2 103 634 B |
| retries | 0 | 【依据】dep3d §7.3 |
| 每例 FDTD 运行数 | 1（本批 `max_fdtd_runs` = 2） | 【依据】dep3d §7.3、决策文档 D3 |
| 线程 / 后端 / 精度 / 设备 | 8 / CUDA / double / device 0，**无 CPU 回退** | 【依据】dep3d §7.3、gate `execution_policy` |
| 3D 单元数上限 | 建议沿用 6 400 万量级（本批 3 200 万在其内），**须重新敲定** | 【依据】决策文档 §4：dep3d 批次专属上限不自动继承 |

> **本提案所有耗时/内存为量级参考，不是 ETA 承诺**；实际值以执行时监督器记录为准【依据】决策文档 §7 末条、dep3d §10.6】。预算 `budget.json` 须在冻结时按实测值重新落盘，不沿用设计期速率估算口径充当预测【依据】决策文档 §4 末条】。

---

## 4. attempt 规则与停止条件（待敲定）

1. **2 次运行各 1 次 attempt**，`retries = 0`，`max_fdtd_runs = 2`，已消耗 attempt 不复用、旧会话不重启【依据】决策文档 §7、dep3d 结果末节】。
2. **失败即中止整批**：任一例未 `completed` 即终止，后续 attempt 保持未消费；不跳过、不重试、不重排【依据】规格书 §3.3、dep3d §3】。
3. **串行，不并发**；**BG 失败则不进入 TGT**（沿用 dep3d §3"背景失败则不进入目标"）【依据】dep3d §3】。
4. **无 CPU 回退**：`required_solver_backend: cuda`、`allow_cpu_solver_fallback: false`；执行前重新核验 CUDA double kernel【依据】gate `execution_policy`、dep3d §8.12】。
5. 超出冻结清单的算例属"批量外扩"，须重新敲定【依据】规格书 §6-P4】。

---

## 5. 契约字段（格式参考模板，本提案不生成文件）

沿用 `configs/research/gprmax_v4_execution_gate.json` 的 `approved_execution_contracts[]` 字段与 `configs/research/batch2d_v1_fine2/cases.json` 的批量字段【依据】规格书 §5.2】：

| 字段 | 本提案取值 | 状态 |
|---|---|---|
| `packet_id` | `A0-3D-V1`（建议） | 【待敲定】 |
| `batch_id` | `a0_3d_v1`（建议；与 `batch2d_v1` 分列，因决策 D1 明确 batch2d_v1 不含 3D） | 【待敲定】 |
| `approved_run_ids` | 2 个（BG、TGT），命名沿用 `B2D-C3m-BG-*` / `B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-*` + 3D 档后缀 | 【待敲定】 |
| `input_path` / `input_sha256` | 冻结时生成 `.in` 后逐例计算登记 | 【待敲定，冻结阶段填】 |
| `runtime_identity_path` / `_sha256` | `artifacts/research_checks/2026-09-24_gpu_setup/runtime_identity.json` / `14abf893d12a6db11e264f324ee83ec66194b9983efbcbc46cab5b79d74c9f38` | 【依据】现有 gate；**冻结前重新核验** |
| `cuda_identity_path` / `_sha256` | `artifacts/research_checks/2026-09-24_gpu_setup/cuda_identity.json` / `65c7c240db5146e1bf37f10a171f572f8ac599933d229a2272c44af0b780e9df` | 【依据】现有 gate；**冻结前重新核验** |
| `launcher_sha256` / `supervisor_sha256` | `ae823072743c972c68b2150796d191035ffe7e65feca604726c239e42a6904e8` / `0c111cdd1f8dc1799c85707b9a4de8040d80b0b9f6a0663c9d87dd4275fec9cc` | 【依据】现有 gate；**冻结前重新核验** |
| `attempt_record` / `run_directory` | `artifacts/research_checks/2026-09-26_<RUN>_attempt.json`（建议）/ `artifacts/simulations/2026-09-26_<RUN>` | 【待敲定】 |
| `continuation` | "Serial in frozen order; any failure aborts batch, later attempts untouched. No retries." | 【依据】现有 gate 同款表述 |
| `mother_model_id` / `group_id` | `B2D-C3m-D10m-W4m-T0.5m-E20-S0.02`、`B2D-C3m-BG` / `B2D-C3m`（2D 侧与 3D 校核同组） | 【依据】规格书 §4.2 |
| `variant_tag` | 建议 `3d5cm` | 【待敲定】 |
| `seed` | 按规格书 §4.1：`sha256(f"{batch_id}|{group_id}|{case_id}|{variant}")` 前 8 字节转整数，冻结时计算 | 【待敲定，冻结阶段填】 |
| `grid_tier` | 建议 `3D5CM` | 【待敲定】 |
| `n_cases` / `n_exceptions` | 2 / 0（冻结时核验） | 【待敲定】 |

---

## 6. 验收与归档要求

1. **逐例验收**：`exit_code = 0`、监督器 `supervision.json` 记 `reason = completed`、实测墙钟与 Job 峰值提交内存落 `execution_outcome`（`wall_s`、`peak_job_commit_bytes`）；两项上限（40 min / 20 GiB）未触及【依据】dep3d 结果首段、fine2 结果执行记录】。
2. **attempt 留痕**：每例 1 份 `attempt.json`，由契约 `attempt_record` 字段逐例指向；attempt 已消耗即不可复用【依据】fine2 结果归档节】。
3. **输出完整性校验**（分析前 assert，沿用 dep3d §9.3 / fine2 的做法）：h5 的 `nx_ny_nz`、`dx_dy_dz`、Rx 位置、dtype `float64`、有限值、BG 与 TGT 源样本一致；接收分量 Ex。
4. **哈希登记**：两份 `.in`、生成的 h5、契约 `cases.json` 与 `groups.json`、分析 `results.json`/`arrays.npz` 的 SHA-256 全部登记；历史输出不覆盖，不用新代码哈希改写旧记录【依据】AGENTS.md 版本管理约束、fine2 结果归档节】。
5. **分析 r1/r2 字节一致**：同一输入连跑两次，`results.json` 逐字节一致并登记哈希（dep3d 与 fine2 均已按此执行）【依据】规格书 §4.1.3、dep3d 结果、fine2 结果】；分析脚本为 CPU 后处理复算，字段 `solver_invoked`/`fdtd_solver_invoked` 为 false。
6. **归档目录**：建议 `artifacts/research_checks/2026-09-26_<RUN>/`（含 h5、输入、stdout/stderr、`supervision.json`、attempt 记录、显存采样、record.json 哈希清单）与 `artifacts/research_checks/2026-09-26_a0_3d_analysis_r1|_r2/`；运行工作目录副本在 `artifacts/simulations/`（忽略目录）【依据】dep3d §12、fine2 结果归档节】。
7. **结果文档**：执行后新建结果文档（不覆盖 dep3d / batch2d 既有文档），逐条照录 `hard_limits` 全字段与 §2.2 的"明确不做"清单；未执行/失败算例写明原因【依据】规格书 §5.4】。
8. **关闭 gate 并回填**：执行后按规格书 §5.3 第 7 步回填 `execution_outcome` 并关闭 gate；凭据不入库【依据】AGENTS.md】。
9. **收尾**：更新 registry / `START_HERE.md` / `decision_log.md`，检查暂存范围后提交并推送 `main`（凭据不入库）【依据】AGENTS.md 版本管理、规格书 §5.3 第 9–10 步】。

---

## 7. 冻结流程（沿用规格书 §5，顺序不可跳过）

1. 本文（设计/提案文档）；
2. 逐例生成 `.in` 与静态检查（域/格距/单元数/dt/时窗/PML/材料/源逐项自算，含 §1.4 整除与余量项；沿用 dep3d `static_check.json` 做法）；
3. 写 `budget.json`（按实测值重算并落盘）；
4. 写 gate（`approved_to_simulate: false` + 批量契约 + 资源上限）；
5. **提交用户敲定**（`requires_agreement_before_execution`）；
6. 用户同意后改 `true` 并执行（单次 attempt、串行、失败即停）；
7. 落 `execution_outcome` 并关闭 gate；
8. 归档（§6）；
9. 更新 registry / `START_HERE.md` / `decision_log.md`；
10. 提交并推送。

【依据】规格书 §5.3；gate `requires_agreement_before_execution`：冻结清单外的 FDTD 运行、自动重试或扩大批次均需事先达成一致。

---

## 8. 限制与明确不声称

1. **`reference_state` 保持 `numerically_unresolved`**：本提案及任何后续 3D 校核结果不升级任何参考资格；是否升级由资格评估另行决定【依据】决策文档 §7】。
2. **不生成**物理阈值、训练标签、clean 真值、可探测性结论或最大探测深度结论；`physical_acceptance_threshold: null`、`training_labels_generated: false`、`clean_truth_generated: false`、`physical_label_eligible: false`、`training_eligible: false`【依据】决策文档 §7、fine2 结果 hard_limits】。
3. **3D 网格收敛未认证**：5 cm 在 170 MHz 仅达官方起步准则（覆盖层 8.8 格/波长、砂岩 11.8 格/波长）【依据】dep3d §10.1、决策文档 §1.2】；本提案只给"网格之间/维度之间的差"，不给"距离真值的误差"。
4. **三类混淆因素未排除（照录 dep3d 设计书 §9.3"剩余混杂"）**：即便格距对齐，跨维度比较仍同时含 (i) **维度差异**（被测量：2D 目标沿 x 无限 / 线源 / 柱面扩散 vs 3D 目标沿 x 有限 2 m / 点源 / 球面扩散）；(ii) **侧向域宽差异**（2D 32 m vs 3D 10 m）；(iii) **目标 y 位置差异**（2D y 14–18 vs 3D y 3–7，域内相对位置相同，但离侧向 PML 的绝对距离不同 → 横向截断的贡献不同）。**本提案不声称任何一项已被排除**；不同比较对的差数值之间不得再相减当作可加的分解。
5. **第四类因素（本提案新增，同样未排除）**：本提案的 3D 侧为 5 cm 各向同性，2D 侧为 25 mm（粗档）或 dy12.5/dz3.125 mm（细档），**格距并不对齐**；dep3d 用以剥离格距变量的"同格距 2D 对照"未列入本提案（见 §9）。因此本提案的跨维差同时混入网格差异，只能在"各自内部配对差分的形状"层面报方向与幅度。
6. **相位不可互迁**；跨维绝对幅值不比较【依据】dep3d §9.3】。
7. 材料参数（覆盖层 εr16/σ0.01、砂岩 εr9/σ0.001、目标 εr20/σ0.02）为研究假设，不是场地实测值【依据】dep3d §10.5】。
8. 结论限定于本场景族与所分析网格，不外推到实测性能、其他几何/材料/倾斜界面；在本提案实际执行并落盘前，batch2d_v1 任何结论不得表述为三维校核通过或实测可用，二维幅度—深度趋势不得直接迁移到三维【依据】决策文档 §7、规格书 §6-P10】。
9. 本文所有预算为量级参考，不是 ETA 承诺；实际值以执行时监督器记录为准。

---

## 9. 未列入本提案的项（照录决策文档决定，不在本文扩大）

| 项 | 决定 | 依据 |
|---|---|---|
| C3-D20m 类新增 3D | 不新增（dep3d 已实测该类场景） | 决策文档 D4 |
| NC / OFF 负控的 3D 校核 | 不列入 | 决策文档 D5 |
| dz 2.5 cm 各向异性档（+2 次运行） | **不列入本文**（决策文档 §4 列为预算备选，需另行敲定） | 决策文档 §4、D2 |
| 同格距 2D 对照（5 cm 2D，+2 次运行） | **不列入本文**；其缺失的后果已写入 §8.5 | 决策文档 D3（本批 2 次运行） |
| C8 / C5 族锚点、电性与电导档的 3D 校核 | 待定（P-A / P-B / P-C），不在本文预设 | 决策文档 §2 §5 |

---

## 10. 复核与溯源

- 预算与执行口径唯一来源：`artifacts/research_checks/2026-09-26_DEP3D_5CM_BG|TGT/` 各目录 `supervision.json`；比较数值来源 `artifacts/research_checks/2026-09-26_dep3d_gold_analysis_r1/results.json`（r1/r2 字节一致）【依据】决策文档 §8】。
- 2D 侧参照值来源：`artifacts/research_checks/2026-09-26_batch2d_analysis_r1/results.json`（粗档）与 `2026-09-26_batch2d_fine2_analysis_r1/results.json`（细档，r1/r2 一致，SHA-256 前 16 位 `dd2ba3328a6f24b1`）；契约 `configs/research/batch2d_v1_fine2/cases.json` SHA-256 `1efe6e29f44950f1eb52ecc781b1c6caf378903e00b302d9480189d20f7a49c3`、`groups.json` SHA-256 `2a39f0f24d492a782af2a3728b8146c915739621c31fb74412fd778ae1cbb5fc`【依据】fine2 结果归档节】。
- 本文为纯提案文档：不生成输入文件、不生成契约 JSON、不改动既有文件（含 gate），不执行仿真。

---

**起草标注**：本文由 **codebuddy hy4-preview** 起草（2026-09-26），状态为**待用户敲定的契约提案**。**kimi 已验收**：逐项抽查数值与来源文件一致（fine2 `cases.json`/`groups.json` SHA-256、gate 四项哈希、粗档/细档 A0 参照值、C3 方向一致性数值、dep3d 墙钟 916.078/916.156 s 与 Job 峰值 8,567,218,176/8,565,354,496 B），并修正一处笔误（§1.4 x 向 PML 内边 1/9 → 1/7，余量 2 m 结论不变）。验收后仍须经用户敲定方可按 §7 冻结为正式契约 JSON 并执行。
