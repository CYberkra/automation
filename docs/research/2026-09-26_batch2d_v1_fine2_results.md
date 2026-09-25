# batch2d_v1 细档（FINE2）批量结果：8/8 完成，四族锚点方向一致性复核与逐位零差分复现对

按[规格书](2026-09-26_batch_2d_spec_v1.md) §6-P2/P3 执行细档复核段：8 次 FDTD 运行（细档 FINE2，dy=12.5 mm / dz=3.125 mm，ZFINE2 同格）全部串行完成，每例 1 次 attempt、无重试、无 CPU 回退，`reason=completed`、exit 0。每例 40 min 墙钟硬上限（2400 s）与 20 GiB Job 提交内存硬上限（21,474,836,480 B）均未触及（规格书 §6-P5）；全批串行总墙钟按 gate `execution_outcome` 记录约 10569 s（执行窗口 2026-09-25T20:09Z 至 2026-09-25T23:06Z，首例 attempt 记录 `2026-09-25T20:09:40Z`）。冻结契约 `configs/research/batch2d_v1_fine2/cases.json` 共 8 例（`n_cases_contract=8`，4 族各 1 个配对例 = 族锚点目标 + 该族同格 BG），`removed_cases` 为空、`n_exceptions=0`，`unavailable_cases` 为空，8/8 全部 `analyzed`。

以下区分来源：墙钟、exit、Job 峰值提交内存为监督器实测（各 run 目录 `supervision.json`，经分析脚本并入结果文件）；差分、谱、到时与 P2 方向指标全部来自分析脚本对已归档 HDF5 输出的实测复算（`analysis_r1/results.json`），不是求解器声明。r1 与 r2 两次连跑的 `results.json` 已核验逐字节一致（SHA-256 前 16 位 `dd2ba3328a6f24b1`）。

## 执行记录

- 8/8 `exit_code=0`、`reason=completed`，attempt 各 1 次（`attempts_consumed=8`），无重试、无 CPU 回退；单 attempt 已消耗，不可复用。
- 墙钟（监督器实测）：8 例全部在 1318.437 – 1323.031 s 区间（最长为 `B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-F2` 1323.031 s），与 ZFINE2 实测 1438.6 s / 1334.8 s 同量级；gate `wall_time_note` 照录："per-run walls 1318.4-1323.0 s, consistent with ZFINE2 measured 1438.6/1334.8 s minus compile cache reuse; total serial batch wall about 10569 s"。
- Job 峰值提交内存（Windows Job committed memory，非 RSS、非 GPU 显存）：15,246,221,312 – 15,248,601,088 B（约 14.2 GiB），全部低于 20 GiB（21,474,836,480 B）硬上限（规格书 §6-P5）。
- h5 元数据（8 例一致项）：dtype float64、`all_finite=true`、网格 1×2560×16000、dx_dy_dz [0.0125, 0.0125, 0.003125] m、dt 1.0112647039820337e-11 s、接收分量 Ex @ (16.65, 45.0) m、接收波形 (118665,) float64 全有限；时窗 1200 ns，尾窗诊断固定 200/400 ns（`tail_taper_fraction` 200ns=0.16666175082586127、400ns=0.33332560844063913）；分析频点网格 `linspace(20e6,170e6,501)`（501 频点）。
- **iterations 约定（如实说明）**：实测 `iterations=118665`、`transformed_samples=118664`，与契约 `time_steps=118663` 相差 2（`iterations_minus_contract_time_steps=2`）。求解器迭代数沿用与粗档 BASE（20352）及归档 DEP_BG_ZFINE2 参考运行（118665）相同的 ceil(window/dt)+1 约定，契约字段 118663 为 floor 计划值。`results.json` `time_steps_note` 照录："the solver iteration count follows the ceil(window/dt)+1 convention seen on the BASE tier (20352) and on the archived DEP ZFINE2 reference (118665); the frozen contract field is a floor-based planning figure"。

### vctip 处置（单列，零干预）

粗档批次曾以外部哨兵在两条件同时满足时（监督器处于 `waiting_descendants` 阶段、vctip 已成孤儿）清除 vctip.exe（MSVC 编译器遥测进程，[粗档报告](2026-09-26_batch2d_v1_results.md)单列记载）。本批哨兵随批运行但**零干预**：8 例 `supervision.json` 逐例核验均为 `waited_for_descendants=false`，无任何清除事件，无本批处置日志文件。事实原因：细档求解约 1320 s，超过 vctip 孤儿空转存活期（粗档观察约 900 s），求解结束时 vctip 已自行退出，哨兵无需触发。本条仅作事实记录，不构成关于 vctip 行为机制的新证据。

## 复算校验：与归档 DEP_BG_ZFINE2 的逐位零差分（单列）

**B2D-C3m-BG-F2 vs 归档 DEP_BG_ZFINE2（运行间可复现性证据，不是收敛或精度结论）**：两运行输入文件均为 16 行，`identical_except_title=true`、`differing_line_indices_excluding_title=[]`（仅 `#title` 行不同，同网格、同接收、同材料、同源）；因此任何非零差分都将是运行间不可复现的证据。实测结果：接收 Ex 波形与源波形 `max_abs_difference=0.0`、`relative_L2_difference=0.0`（各 118665 样本，接收相关系数 0.9999999999999998、源相关系数 1.0，dt 一致）。即两运行输出逐位相同。

`hard_limits.reproduction_pair_is` 照录："run-to-run reproducibility evidence on one grid tier, not a convergence result and not an accuracy statement"。本条只证明该网格档上运行间可复现性，不证明对精确解的收敛，也不构成任何精度声明；`reference_state` 保持 `numerically_unresolved`。

## P2 方向一致性（核心，按族分层，不跨族平均、不跨档合并）

粗档锚点 `B2D-C{h}m-D10m-W4m-T0.5m-E20-S0.02`（BASE）与细档同母模型锚点 `-F2`（FINE2）逐族配对比较。`results.json` `reporting` 照录："rows are family x grid-tier-stratified: no cross-family average, no cross-family differencing, and the BASE and FINE2 results are never merged, subtracted as error, or used to certify grid convergence (spec §6-P1/P2, §7.8); each tier keeps its own native dt, its own tail window samples and its own grid_tier label"。各 tier 各自使用原生 dt（BASE 5.896635841874211e-11 s、FINE2 1.0112647039820337e-11 s）、各自尾窗样本数；无通过阈值（`direction_consistency_threshold=null`，`threshold` 照录："none declared: the numbers are reported as directions/agreement rates only; no pass threshold exists in this analysis or in the frozen contract (spec §6-P2)"）。

| 族 | 差分谱形状相关（主 200ns / 稳健 400ns） | 逐频点符号一致率 real / imag | 包络峰到时 BASE→FINE2 (ns) | FINE2−BASE | 谱 L2 比 FINE2/BASE | 时域峰值比 FINE2/BASE | 谱峰 BASE→FINE2 (MHz) |
|---|---|---|---|---:|---:|---:|---|
| C1 | 0.999548930464859 / 0.9995489302578601 | 93.81% / 93.81% | 320.0104271385134 → 302.7726523722209 | −17.237774766292546 | 1.0228315342276035 | 1.6657024470013109 | 170.0 → 166.1 |
| C3 | 0.9994178181733713 / 0.9994178183239077 | 92.22% / 92.22% | 331.8626651806806 → 316.1921349940624 | −15.67053018661818 | 1.039002729692549 | 1.903289349007366 | 169.7 → 170.0 |
| C5 | 0.9990272372038344 / 0.9990272371914642 | 89.82% / 92.61% | 346.7221875022036 → 329.61161761590404 | −17.110569886299572 | 1.0558752689384794 | 1.918323952556588 | 168.2 → 168.8 |
| C8 | 0.9981254455630282 / 0.998125445771708 | 88.02% / 88.62% | 368.8345719092319 → 349.7155599310669 | −19.11901197816502 | 1.0803007766438417 | 2.1139523083257004 | 169.4 → 170.0 |

逐项说明：

- **差分谱形状相关**：四族主 200ns 窗与稳健 400ns 窗一致（同族两窗相关值之差 < 3e-10），范围 0.99813（C8）– 0.99955（C1）。符号一致率统计照录 `zero_masking="none: every requested frequency point is counted, including zero parts"`，且 `n_points_with_zero_real=0`、`n_points_with_zero_imag=0`（501 点全部计入）。
- **包络峰到时**：四族细档一致早于粗档 15.7 – 19.1 ns，**各族同向**（`direction="FINE2 envelope peak earlier"`）。此为系统性档间偏移而非机制改变：两档 dt 不同，各档在各自原生时间轴上读取（note 照录："each tier is read on its own native time axis; the two dt differ, so the comparison is of the peak index time, not of a common grid"）；跨档比较只取偏移方向，不取偏移量作误差。
- **幅值比方向**：四族谱 L2 比与 |时域峰值| 比均显示 FINE2 更大（`spectral_L2_direction="FINE2 larger"`、`time_peak_abs_direction="FINE2 larger"`）。note 照录（`amplitude_ratio_direction_is`）："cross-tier ratio reported as a direction indicator only; it mixes grid dispersion, source discretisation and taper sample count, so it is not an accuracy or grid-error measure"。即该比值混合了网格色散、源离散与尾窗样本数，只作方向指标，不作精度或网格误差度量。
- **谱峰频差**：如实列出（上表末列），不设"通过"判读。

## 细档自身锚点指标（按族）

配对规则（`pairing_policy`）：每个非 BG 细档算例只与同族同档 BG 差分（C1/C3/C5/C8 族 BG 分别为 `B2D-C1m-BG-F2`/`B2D-C3m-BG-F2`/`B2D-C5m-BG-F2`/`B2D-C8m-BG-F2`），`cross_family_differencing=false`、`cross_tier_differencing=false`。角色构成（`stratified_by_family`）：四族各 2 例（BG 1 / target 1），每族 `n_differences_computed=1`；BG 例 `difference.computed=false`（note 照录："family BG is the differencing reference itself (spec §6-P7)"）。各目标例 note 照录："time_diff_peak_V_m carries the pair-internal excitation scale; no absolute detectability claim"；`source_samples_identical_to_family_bg=true`（四族一致）。

| 族 | 目标 run_id（省略公共前缀 B2D-C{h}m-） | 带内能量比 DIFF/BG | (dB) | 包络峰到时 (ns) | 谱峰 (MHz) |
|---|---|---:|---:|---:|---:|
| C1 | D10m-W4m-T0.5m-E20-S0.02-F2 | 7.872880397772246e-06 | −51.03866346167797 | 302.7726523722209 | 166.1 |
| C3 | D10m-W4m-T0.5m-E20-S0.02-F2 | 3.092143334046171e-07 | −65.09740382883551 | 316.1921349940624 | 170.0 |
| C5 | D10m-W4m-T0.5m-E20-S0.02-F2 | 1.2121864557034453e-08 | −79.1643057286364 | 329.61161761590404 | 168.8 |
| C8 | D10m-W4m-T0.5m-E20-S0.02-F2 | 9.369628603300173e-11 | −100.2827762349152 | 349.7155599310669 | 170.0 |

`in_band_energy_ratio_is` 照录："diagnostic energy ratio on the 501-point analysis grid; not detectability, not SNR, not a physical threshold"。

- **族序随覆盖厚度单调递减**：C1 7.87e-06 → C3 3.09e-07 → C5 1.21e-08 → C8 9.37e-11，与粗档族序方向一致（粗档同锚点为 7.52e-06 / 2.86e-07 / 1.09e-08 / 8.03e-11）。
- **包络峰到时随覆盖厚度递增**：302.77 → 316.19 → 329.61 → 349.72 ns，与粗档同向（粗档 320.01 / 331.86 / 346.72 / 368.83 ns）。
- **C8 仍处数值分辨关注区**：细档 9.369628603300173e-11（−100.2827762349152 dB）与粗档 8.03e-11 同量级。本批不基于该量级声明任何可探测性结论（含"不可探测"），C8 仅作厚覆盖吸收压力下的机制记录（规格书 §6-P6 风险诊断定位）。

## 结论（按规格书 §6-P2）

粗档[结果报告](2026-09-26_batch2d_v1_results.md)的三项机制结论——带内能量比随覆盖厚度深度单调递减、包络峰到时随覆盖厚度递增、电性对比增强方向上差分响应上升——在细档锚点上获得**方向一致支持**：四族差分谱形状相关 0.99813–0.99955、逐频点符号一致率 88.0%–93.8%、包络峰到时四族同向（FINE2 一致早 15.7–19.1 ns，系统性档间偏移）、带内能量比族序与包络峰到时族序均与粗档同向。无机制需要按 §6-P2 降级为"待复核"。该支持保留为**机制/算子评价证据**，不是物理充分性认证：细档复核只判方向一致，不用于生成绝对幅度/相位结论（规格书 §6-P3），亦不得把方向一致性外推为粗档精度达标（§7.8 禁止事项）。

## 硬限制（`results.json` hard_limits 全字段，逐条照录）

| 字段 | 值 |
|---|---|
| `absolute_accuracy_claimed` | false |
| `amplitude_ratio_direction_is` | "cross-tier ratio reported as a direction indicator only; it mixes grid dispersion, source discretisation and taper sample count, so it is not an accuracy or grid-error measure" |
| `clean_truth_generated` | false |
| `conclusions_scope` | "mechanism and operator-evaluation evidence for the four batch2d_v1 fine2 scenario families only (spec §6-P2/P6); not field performance, not 3D, not absolute accuracy, no maximum-depth claim" |
| `cross_family_average` | false |
| `cross_family_differencing` | false |
| `cross_tier_differencing` | false |
| `direction_consistency_threshold` | null |
| `family_bg_required` | true |
| `fdtd_solver_invoked` | false |
| `grid_convergence_certified` | false |
| `grid_tier` | "FINE2" |
| `grid_tiers_merged` | false |
| `in_band_energy_ratio_is` | "diagnostic energy ratio on the 501-point analysis grid; not detectability, not SNR, not a physical threshold" |
| `physical_acceptance_threshold` | null |
| `physical_label_eligible` | false |
| `reference_state` | "numerically_unresolved" |
| `reproduction_pair_is` | "run-to-run reproducibility evidence on one grid tier, not a convergence result and not an accuracy statement" |
| `training_eligible` | false |
| `training_labels_generated` | false |

即：本批不产生绝对精度声明、不产生 clean 真值、不产生物理阈值、不产生物理/训练标签；结论限定于 batch2d_v1 细档四族场景的机制与算子评价证据，不是实测性能、不是 3D、不含最大深度声明；不跨族平均、不跨族差分、不跨档差分、不合并网格档；方向一致性无通过阈值；网格收敛未认证；参考状态保持 `numerically_unresolved`；分析脚本仅做后处理复算（`fdtd_solver_invoked=false`）。

## 其他限制

- 结论限于四族锚点场景与所分析网格（FINE2，dy12.5/dz3.125 mm）：本批每族仅 1 个配对例，无深度/尺寸/电性单因子细档扫描，族间结论不得用于未覆盖的档位。
- 逐位零差分复现对只证明运行间可复现性（同输入同机同链），不证明收敛或精度；跨档比值混合网格色散/源离散/尾窗样本数，只作方向指标。
- C8 低量级比值处于数值分辨关注区，不作可探测性声明；其零差分确定性由粗档同族负控证据支持，本批细档段未设 NC/OFF 负控（规格书 §1.6 细档段仅 4 配对例）。
- 材料参数（覆盖层 εr16/σ0.01、砂岩 εr9/σ0.001、目标 εr20/σ0.02）为研究假设，不是场地实测值；二维结果不外推三维或实测（规格书 §6-P10 三维校核子集未确认）。

## 归档与复现

- 运行归档：8 个 run 目录 `artifacts/simulations/2026-09-26_B2D-*-F2/`，各含原始 h5、输入、stdout/stderr、`supervision.json`、`live_status.json`；全部 h5 与输入的路径及 SHA-256 见 `analysis_r1/results.json` 的 `inputs` 字段（另含契约 `configs/research/batch2d_v1_fine2/cases.json` SHA-256 `1efe6e29f44950f1eb52ecc781b1c6caf378903e00b302d9480189d20f7a49c3`、`groups.json` SHA-256 `2a39f0f24d492a782af2a3728b8146c915739621c31fb74412fd778ae1cbb5fc`）。
- attempt 记录：`artifacts/research_checks/2026-09-26_B2D-{C1m,C3m,C5m,C8m}-BG-F2_attempt.json` 与 `2026-09-26_B2D-{C1m,C3m,C5m,C8m}-D10m-W4m-T0.5m-E20-S0.02-F2_attempt.json` 共 8 份（gate 契约 `attempt_record` 字段逐例指向，packet_id `BATCH2D-V1-FINE2`）；gate 为 `configs/research/gprmax_v4_execution_gate.json`（`batch_id` 已标记 `batch2d_v1_fine2_consumed`、gate 已关闭，8 例 exit/墙钟/Job 峰值与验收摘要已回填 `execution_outcome`）。
- 分析与复算：`artifacts/research_checks/2026-09-26_batch2d_fine2_analysis_r1/` 与 `_r2/`，两次连跑 `results.json` 逐字节一致（SHA-256 前 16 位 `dd2ba3328a6f24b1`）。
- 复现入口：`scripts/analyze_batch2d_fine2.py`（CPU 后处理复算，`fdtd_solver_invoked=false`，不调用求解器）。
- 8 次 attempt 各自一次，已消耗，不可复用；旧会话不重启。

## 下一步

- 粗档/细档两段均已完成后，按规格书 §6-P2 将分层机制结论并入算子契约评价的输入材料；不生成物理/训练标签，不做 3D/实测外推，不宣布网格收敛。
- P10 三维校核子集：决策已落盘（[决策文档](2026-09-26_3d_validation_subset_decision.md)，首选 C3 族锚点 A0 配对、时机在本批之后、契约另行冻结并经用户敲定）；本批细档方向一致性结果即该决策 P-A/P-B 待定项所依赖的复核输入，后续按决策文档依赖链推进，本批不预设。
