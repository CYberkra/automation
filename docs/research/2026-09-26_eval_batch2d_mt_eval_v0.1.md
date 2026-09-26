# batch2d_v1_mt 多道评价运行器 v0.1 首跑：27 候选 × 61 事件行全量执行，可行性维持 undetermined

日期：2026-09-26。入口：用户 15:14 批准"实现评价侧变偏移距道集语义并开跑 27 项候选"。本文是 S2（运行器实现）+ S3（27 候选首跑）的产物报告；评价与标签口径唯一来源为[评价与标签协议 v0.2](2026-09-24_evaluation_and_labels_v0.2.md)，运行器结构按[评价运行器与选择器设计 v0.1](2026-09-26_eval_runner_and_selector_design_v0.1.md) §A。

## 变偏移距道集语义（本轮实现的核心）

- **输入阵列**：每个评价窗取 `[n_window_samples, 33]`（样本轴 × 道轴），33 道为单发多收道集（Tx 固定 y15.35，偏移距 2.70→1.30→5.30 m），**道轴不是空间测线轴，本数据不是常偏移 B-scan**。此语义以 `gather_semantics.version = variable_offset_csg_not_common_offset_bscan_v1` 写入每条记录与 manifest，随数据传递给一切下游。
- **窗映射**：事件窗为冻结样本轴区间（锚点几何自算双程到时 ±80 ns 内收取整），**原样、不加任何道间平移地**施加到全部 33 道。不做逐道对齐——窗在候选运行前已冻结，逐道平移等于重新定义窗（设计 §A.3 门禁 1  spirit、v0.2 §2）。
- **语义后果（记录于每条记录）**：直耦波端到端道间时差约 147 样本（提案 §2.4 自算），事件/背景在道轴上不再平层对齐，"平层=秩一"前提不成立；`mean`/`svd` 的道轴操作看到的是道集而非拉平 B-scan。实测奇异谱平滑衰减（见下），与此前文献核查的标注一致。

## 执行规模与门禁核验

- **门禁 1（事件表冻结）**：整文件 SHA-256 核验等于冻结值 `b0ad1003…a607f56c`，`status=frozen`、43 条。
- **门禁 2（输入形态）**：33 道 ≥ 2，`single_trace_gate_blocked=false`；前次 27/27 被 `min(shape)>=2` 挡下的状态正式结束，**27/27 候选全部可运行**。
- **门禁 3**：BG 序增益前真实中间量必存（756 个 .npy）；GB 序只存 `G⁻¹Y` 审计量（324 个 .npy，标注 `ginv_audit_not_true_intermediate`）；负控行独立成组、不并入目标行。
- **规模**：43 事件 × MT 案例展开为 61 事件行（4 覆盖层界面事件按 `applies_to=all_cases_of_family` 施加到各族全部 22 个 MT 案例；目标 14 行、NC 8 行、BG-NEG 10 行、OFF 1 行）× 27 候选 = **1647 条记录：1458 ran / 189 unavailable，零数值失败**。189 条不可用全部是 7 个粗档有单道无多道的母模型（C1m-D2m/D5m、C3m-OFF、C5m-D20m、C8m-D20m×2、C8m-E28-S0.05），原因枚举 `no_multitrace_case_in_batch2d_v1_mt`。
- **可计算指标（当前证据状态下唯一）**：负控残差 `N_b`（`energy_metrics`，绝对能量 + 输入基准 + 非零基准能量比）；增益风险诊断（曲线、最大增益、截幅样本数、溢出、增益前后能量比）；SVD 奇异谱诊断；资源（墙钟/tracemalloc/dtype/窗维度）。**`D/A/H/ρ` 不可计算**——`reference_state=numerically_unresolved` + `paired_contrast`，逐条记录照录门禁原因 `reference_state_or_scope:numerically_unresolved:paired_contrast:contrast`；**`R_c` 不计算**（无合法纯干扰窗）。
- **可行性**：`configuration_labels` 全 null 界限（G4 未解除）⇒ 全部 gate `undetermined`、状态 `partial_only`；不排名、不选优、无 `unique_admissible`。

## 诊断观察（机制方向，不是质量分）

以下数值全部为**负控窗/绝对残差诊断**，不是可探测性、不是 SNR、不是物理阈值，不构成任何候选的优劣结论：

1. **恒等 sanity**：`B0_G1_BG`（identity）在全部 54 个可评价窗上 `energy_ratio=1.0` 逐位成立，管道无隐性缩放。
2. **负控窗残差能量比中位数**（18 个 nc_zero+bg_absent 窗）：identity 1.0000；mean λ=0.25/0.5/1.0 → 0.6042 / 0.3215 / 0.0953；svd k=1/2/3 → 0.0869 / 0.0320 / 0.0163。背景算子在道集上去除直耦/背景能量的幅度量级在此；因窗内含强直耦且道间不对齐，解释须携带变偏移距语义。
3. **奇异谱平滑衰减**（k=1 截断全部 `status=ok`，无 `svd_cutoff_gap_unresolved`）：例 C3 锚点目标窗 σ/σ1 前 8 = [1.0, 0.2332, 0.1274, 0.0867, 0.0671, 0.0480, 0.0297, 0.0235]；C1 覆盖层界面窗 [1.0, 0.2686, 0.1953, 0.1606, 0.1418, 0.1150, 0.0969, 0.0851]。无尖锐截断间隙，低秩背景模型在本道集上是近似而非精确——这与"平层=秩一"不严格成立的预判一致。
4. **增益安全**：1080 条增益诊断全部 `clipped_samples=0`、`overflow=false`；BG 序 `pre_gain` 已按契约持久化。
5. **NC 窗内零差分**：全部 nc_zero 事件行逐窗核验 `NC−BG max_abs = 0.0`（逐道逐样本），多道负控在本评价口径下成立。

## 确定性与归档

- r1/r2 两遍连跑：剥离每记录 `resource`（墙钟/内存为实测资源量，天然不可复现）后 `records.json` **逐字节一致**；`run_manifest.json` 除 `run_tag`/`total_wall_s` 外一致。增益曲线以 .npy 持久化（记录内带 SHA-256），两遍写出的数组字节一致。
- 记录与 manifest：`artifacts/research_checks/2026-09-26_eval_batch2d_mt_r1/`、`_r2/`（r1 为主口径）。
- 大数组（2160 个 .npy = 1080 个增益前中间量/审计量 + 1080 条增益曲线，约 1.1 GB）：`artifacts/simulations/2026-09-26_eval_batch2d_mt_arrays/`（git 忽略目录；记录内带路径 + SHA-256，可逐数组校验；曲线不内联记录，避免单文件超 git 托管单文件上限）。
- 复现入口：`scripts/run_eval_batch2d_mt.py --output-dir <新目录> --run-tag <tag>`；只读输入（归档 h5、冻结事件表、算子契约），不调用求解器，`fdtd_solver_invoked=false`。

## 边界（照录并新增）

- `hard_limits` 逐字段照录（粗档 16 项 + 多道 4 项），见 `run_manifest.json`；**不产生绝对精度声明、clean 真值、物理阈值、物理/训练标签**；不跨族/跨档差分；负控不并入目标行；参考状态保持 `numerically_unresolved`；网格收敛未认证。
- 运行器自身限制：`r_c_computed=false`、`isolated_event_granted=false`、`cross_family_metrics=false`、`gain_reference_absent=true`、`waveform_metrics_computed=false`。
- N_b 能量比是**绝对残差诊断量**；分母非零才报比值，零基准不加 epsilon；它不构成"背景抑制效果验收"。
- 本跑只覆盖 BASE 档多道数据；FINE2 多道（G5）、3D（G1，A0 待批）未涉及。
- 7 个粗档母模型无多道数据（OFF 等），对应 189 条记录不可用；补多道数据须另冻结批次并经用户批准。

## 下一步

1. **S4 基线参照**：identity 已有（B0 行即 I）；F（冻结固定流程）/L（`local_mean_control` 窗宽 3/5/7）/R（输入奇异谱规则）/O（离线最优上界）现在数据形态全部合法，可进入开发组冻结流程；同等动作集合与调参机会，不削弱基线。
2. **G4 仍是可行性与排序的总阻塞**：阈值校准（v0.2 §9 四步）前，一切候选维持 `undetermined`。
3. **A0 3D 校核**：提案待用户批准，与本文独立。
4. 选择器训练继续禁止（`training_eligible=false`）。
