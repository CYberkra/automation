# batch2d_v1 粗档批量结果：29 例全部完成，四族机制扫描与零差分负控

按[规格书](2026-09-26_batch_2d_spec_v1.md)执行：29 次 FDTD 运行（粗档 BASE，dy=dz=25 mm）全部串行完成，每例 1 次 attempt、无重试、无 CPU 回退（`retries=0`；批量契约形态 `max_fdtd_runs=29`，gate `batch2d_v1`），`reason=completed`、exit 0。4 GiB Job 提交内存与 20 min 墙钟硬上限（规格书 §6-P5 的 10 min 经 gate `scope_expansion_basis` 留有依据地修正为 20 min：dep3d_gold 实测墙钟 910.6–920.2 s 表明求解后 Job 滞留约 900 s，10 min 会误杀健康运行）均未触及；全批实测合计墙钟 1400.62 s（约 23.4 min）。

以下区分来源：墙钟、exit、Job 峰值提交内存为监督器实测（各 run 目录 `supervision.json`，经分析脚本并入结果文件）；差分、谱与到时指标全部来自分析脚本对已归档 HDF5 输出的实测复算（`analysis_r1/results.json`），不是求解器声明。r1 与 r2 两次连跑的 `results.json` 已核验逐字节一致。冻结契约 `configs/research/batch2d_v1/cases.json` 共 29 例（`n_cases_contract=29`）：规格书 §1.6 规划的 32 例粗档中有 3 例（`B2D-C3m-D2m-W4m-T0.5m-E20-S0.02`、`B2D-C5m-D5m-W4m-T0.5m-E20-S0.02`、`B2D-C8m-D5m-W4m-T0.5m-E20-S0.02`）在冻结时移除，契约 `removed_cases` 照录理由："target not fully inside sandstone (inside cover layer or straddling cover interface); removed 2026-09-26 per spec 1.5, cover-embedded anomaly role out of scope"。`unavailable_cases` 为空，29/29 全部 `analyzed`。

## 执行记录

- 29/29 `exit_code=0`、`reason=completed`，attempt 各 1 次，无重试、无 CPU 回退。
- 墙钟（监督器实测）：首例 `B2D-C3m-BG` 532.156 s；其余 28 例 29.859–59.218 s（四舍五入 29.9–59.2 s；最长为 `B2D-C3m-D5m` 59.218 s）。首例偏高的逐项构成未被单独归因，本批不做声明。
- 全批合计墙钟 1400.62 s ≈ 23.4 分钟完成全批；规格书 §3.2 规划包络为 2.5 min/例（152 s）——除首例外单例墙钟均低于该包络值。包络值是保守规划口径、不是 ETA 承诺，估算与实测的差异本身是规格书 §3.1 已知限制所述情形。
- Job 峰值提交内存（Windows Job committed memory，非 RSS、非 GPU 显存）：2,914,930,688 – 2,918,010,880 B（约 2.72 GiB），全部低于 4 GiB（4,294,967,296 B）硬上限（规格书 §6-P5）。
- h5 元数据（29 例一致项）：dtype float64、`all_finite=true`、网格 1×1280×2000、dx_dy_dz 0.025 m、dt 5.896635841874211e-11 s、iterations 20352、transformed_samples 20351、接收分量 Ex @ (16.65, 45.0) m；时窗 1200 ns，尾窗诊断固定 200/400 ns（`tail_taper_fraction` 200ns=0.16666257186379047、400ns=0.33333742813620953）；分析频点网格 `linspace(20e6,170e6,501)`（501 频点）。

### vctip 处置（单列）

dep3d_gold 批次观察到的"求解完成后 Job 内成员滞留"现象（[该文档根因补记](2026-09-26_dep3d_gold_results.md)在同环境同工具链下鉴定为 `vctip.exe`，MSVC 编译器遥测进程）在本批升级为主动处置：外部哨兵仅在两个条件同时满足时清除该进程——监督器处于 `waiting_descendants` 阶段、且 vctip 已成孤儿（其父编译链进程已退出）；除此之外不做任何干预。全批 29 次清除事件逐条留痕于 `artifacts/research_checks/2026-09-26_batch2d_v1_vctip_intervention_log.jsonl`（29 行）。监督器完成语义（root exit 0 + Job 清空后判定 completed）未被修改；求解器进程、输入文件与输出 h5 未被触碰。该处置只影响求解完成后的空转滞留，不参与任何求解计算；各例输出 h5 均完整归档并通过分析前校验（`all_finite=true`）。

## 负控（单列，不并入目标行）

**零对比负控 NC（四族各 1，均 D10m，几何存在、εr9/σ0.001 与围岩完全相同）**：四族 NC 的配对差分全部 `zero_difference=true`，`in_band_energy_ratio_DIFF_over_BG=0.0`（ratio 0.0），在 200/400 ns 两个尾窗与时域逐样本恒为零。dB 记 null 的原因：差分逐样本恒为零时对数能量比（10·log10）无定义，分析脚本按其 note 记 null（原文："paired difference is identically zero on both tail windows and in time (expected for the null-contrast control); peak/phase/dB metrics not applicable, reported as null"）；包络峰到时、谱峰频率/相位同样为 null。各 NC 的 `source_samples_identical_to_family_bg=true`。四族 NC 恒零表明：同族同格同源配对差分的确定性成立，几何/网格本身不产生伪响应（规格书 §1.4 预期"预期近零，偏离指示几何/平滑残差"）。

| 族 | run_id | zero_difference | 带内能量比 | dB |
|---|---|---|---:|---:|
| C1 | B2D-C1m-D10m-NC | true | 0.0 | null |
| C3 | B2D-C3m-D10m-NC | true | 0.0 | null |
| C5 | B2D-C5m-D10m-NC | true | 0.0 | null |
| C8 | B2D-C8m-D10m-NC | true | 0.0 | null |

**远偏移负控 OFF（`B2D-C3m-D10m-OFF`，仅 C3 族 1 例，单列）**：目标（εr20/σ0.02，y19–23 m）保留电性对比但横向移至距收发基线中心 5 m 处，是无目标假阳性对照，不是可探测性证据。差分指标：`in_band_energy_ratio_DIFF_over_BG=3.0406928886248514e-08`（−75.17027441582566 dB）、包络峰到时 319.3617971959073 ns、谱峰 36.8 MHz。照录其例外标注：`exceptions=["side_pml_margin_below_13m"]`，静态核验 `side_pml_margin_m=8.0`（< 13 m 族参考值；域 y0–32、横向 PML 物理厚度 1 m、PML 内边 y1/y31）。按规格书 §1.4 冻结处置选项 (a)：保留并加注此例外；其结果不与目标行合并使用。

## 分层报告（按族 × 角色，不跨族平均、不跨族差分）

先照录本报告主指标的硬限制（`results.json` hard_limits 原文，完整字段见下文"硬限制"节）：`in_band_energy_ratio_is` = "diagnostic energy ratio on the 501-point analysis grid; not detectability, not SNR, not a physical threshold"。即下列**带内能量比**（`in_band_energy_ratio_DIFF_over_BG`）是 501 频点诊断网格上"TGT−BG 差分带内能量 / 背景带内能量"的诊断性能量比，**不是可探测性、不是 SNR、不是物理阈值**。各目标例 note 照录："time_diff_peak_V_m carries the pair-internal excitation scale; no absolute detectability claim"。

配对规则（`pairing_policy`）：每个非 BG 算例只与同族同格 BG 差分（C1/C3/C5/C8 的族 BG 分别为 `B2D-C1m-BG`/`B2D-C3m-BG`/`B2D-C5m-BG`/`B2D-C8m-BG`），`cross_family_differencing=false`。角色构成（`stratified_by_family`）：C1 族 5 例（BG 1 / NC 1 / target 3）、C3 族 14 例（BG 1 / NC 1 / OFF 1 / target 11）、C5 族 4 例（BG 1 / NC 1 / target 2）、C8 族 6 例（BG 1 / NC 1 / target 4）。

### C1 族（覆盖层 1 m，group B2D-C1m，BG=B2D-C1m-BG）

| run_id（省略公共前缀 B2D-C1m-） | 角色 | 带内能量比 | (dB) | 包络峰到时 (ns) | 谱峰 (MHz) |
|---|---|---:|---:|---:|---:|
| BG | 族 BG | 差分参照自身 | — | — | — |
| D10m-NC | 零对比负控 | 0.0（zero_difference=true） | null | null | null |
| D2m-W4m-T0.5m-E20-S0.02 | target（浅部对照） | 7.479575188794756e-05 | −41.26123067697229 | 144.2317126922432 | 158.6 |
| D5m-W4m-T0.5m-E20-S0.02 | target（深度轴） | 3.188037626649901e-05 | −44.96476561507421 | 205.02602822196633 | 170.0 |
| D10m-W4m-T0.5m-E20-S0.02 | target（族锚点） | 7.522928910775171e-06 | −51.23613042102644 | 320.0104271385134 | 170.0 |

### C3 族（覆盖层 3 m，group B2D-C3m，BG=B2D-C3m-BG；含全批锚点 A0）

| run_id（省略公共前缀 B2D-C3m-） | 角色 | 带内能量比 | (dB) | 包络峰到时 (ns) | 谱峰 (MHz) |
|---|---|---:|---:|---:|---:|
| BG | 族 BG | 差分参照自身 | — | — | — |
| D10m-NC | 零对比负控 | 0.0（zero_difference=true） | null | null | null |
| D10m-OFF | 远偏移负控（余量例外单列） | 3.0406928886248514e-08 | −75.17027441582566 | 319.3617971959073 | 36.8 |
| D10m-W1m-T0.5m-E20-S0.02 | target（宽度单因子） | 2.063928687363525e-08 | −76.85305312480656 | 335.6954784778988 | 170.0 |
| D10m-W2m-T0.5m-E20-S0.02 | target（宽度单因子） | 7.430583938942992e-08 | −71.28977055486301 | 331.7447324638431 | 170.0 |
| D10m-W4m-T0.25m-E20-S0.02 | target（厚度单因子，薄层） | 2.314087836392468e-07 | −66.35620160450573 | 331.68576610542436 | 170.0 |
| D10m-W4m-T1m-E20-S0.02 | target（厚度单因子，厚异常） | 2.4986913266274587e-07 | −66.0228739070118 | 314.5265558055704 | 152.6 |
| D10m-W4m-T0.5m-E12-S0.005 | target（电性弱档） | 4.972339721308896e-08 | −73.03439207051788 | 329.9757417112808 | 133.7 |
| D10m-W4m-T0.5m-E20-S0.02 | target（**全批锚点 A0**） | 2.8635252336237045e-07 | −65.43098985494515 | 331.8626651806806 | 169.7 |
| D10m-W4m-T0.5m-E28-S0.05 | target（电性强档） | 4.809910142057482e-07 | −63.17863036967773 | 319.18489812065104 | 147.2 |
| D10m-W4m-T0.5m-E20-S0.005 | target（电导单因子） | 3.8268347934505973e-07 | −64.17160286100942 | 335.6954784778988 | 169.7 |
| D10m-W4m-T0.5m-E20-S0.05 | target（电导单因子） | 3.1493307602392536e-07 | −65.01781724945178 | 318.59523453646364 | 169.7 |
| D5m-W4m-T0.5m-E20-S0.02 | target（深度轴） | 1.2187313192682688e-06 | −59.140920281146464 | 218.64725701669576 | 169.7 |
| D20m-W4m-T0.5m-E20-S0.02 | target（深部弱响应） | 1.651093944888008e-08 | −77.82228215298248 | 536.9476597610657 | 169.7 |

### C5 族（覆盖层 5 m，group B2D-C5m，BG=B2D-C5m-BG）

| run_id（省略公共前缀 B2D-C5m-） | 角色 | 带内能量比 | (dB) | 包络峰到时 (ns) | 谱峰 (MHz) |
|---|---|---:|---:|---:|---:|
| BG | 族 BG | 差分参照自身 | — | — | — |
| D10m-NC | 零对比负控 | 0.0（zero_difference=true） | null | null | null |
| D10m-W4m-T0.5m-E20-S0.02 | target（族锚点） | 1.0869843969359355e-08 | −79.63776689927744 | 346.7221875022036 | 168.2 |
| D20m-W4m-T0.5m-E20-S0.02 | target（厚覆盖+深目标） | 6.245832594764899e-10 | −92.04409660218558 | 546.1464116743895 | 169.1 |

### C8 族（覆盖层 8 m 吸收压力族，group B2D-C8m，BG=B2D-C8m-BG）

| run_id（省略公共前缀 B2D-C8m-） | 角色 | 带内能量比 | (dB) | 包络峰到时 (ns) | 谱峰 (MHz) |
|---|---|---:|---:|---:|---:|
| BG | 族 BG | 差分参照自身 | — | — | — |
| D10m-NC | 零对比负控 | 0.0（zero_difference=true） | null | null | null |
| D10m-W4m-T0.5m-E20-S0.02 | target（族锚点） | 8.026249736550938e-11 | −100.95487331330483 | 368.8345719092319 | 169.4 |
| D10m-W4m-T0.5m-E28-S0.05 | target（厚覆盖+强对比） | 1.3575229363982578e-10 | −98.67252823936502 | 354.9774776808275 | 142.7 |
| D20m-W4m-T0.5m-E20-S0.02 | target（厚覆盖+深目标，不可观测压力） | 4.598276710920589e-12 | −113.37404897701373 | 567.5511997803928 | 169.7 |
| D20m-W4m-T0.25m-E20-S0.02 | target（厚覆盖+深薄层） | 3.8167380710479025e-12 | −114.1830764309039 | 563.1287228989871 | 170.0 |

### 深度扫描：带内能量比随深度单调递减

族内同格差分、锚点尺寸（W4×T0.5）与锚点电性（E20-S0.02）下：

- C1（三档全）：D2m 7.479575188794756e-05 → D5m 3.188037626649901e-05 → D10m 7.522928910775171e-06，单调递减；
- C3（三档全）：D5m 1.2187313192682688e-06 → D10m 2.8635252336237045e-07 → D20m 1.651093944888008e-08，单调递减；
- C5（两档）：D10m 1.0869843969359355e-08 → D20m 6.245832594764899e-10，同向；
- C8（两档）：D10m 8.026249736550938e-11 → D20m-W4m-T0.5m 4.598276710920589e-12，同向（另有 D20m-W4m-T0.25m 3.8167380710479025e-12）。

C1/C3 完整三档严格单调递减，C5/C8 两档方向一致。按规格书 §6-P6，此处只报告机制方向，不生成可探测性或最大探测深度结论。

### 包络峰到时随深度递增

- C3：D5m 218.64725701669576 → D10m 331.8626651806806 → D20m 536.9476597610657 ns（218.6 / 331.9 / 536.9 ns）；
- C1：D2m 144.2317126922432 → D5m 205.02602822196633 → D10m 320.0104271385134 ns；
- C5：D10m 346.7221875022036 → D20m 546.1464116743895 ns；C8：D10m 368.8345719092319 → D20m-T0.5m 567.5511997803928 ns（T0.25m 为 563.1287228989871 ns）。

四族到时次序均随深度增加而变晚，与双程走时方向一致；本批不做时窗端点约定换算与物理深度声明。

### 电性档对比（C3 族 D10m，同尺寸 W4×T0.5）

介电—电导三档：E28-S0.05 4.809910142057482e-07（−63.17863036967773 dB）> E20-S0.02 2.8635252336237045e-07（−65.43098985494515 dB）> E12-S0.005 4.972339721308896e-08（−73.03439207051788 dB），即对比度增强方向上带内能量比上升。

电导单因子（固定 εr20，D10m）仅报数值、不声称单调序：S0.005 3.8268347934505973e-07、S0.02 2.8635252336237045e-07、S0.05 3.1493307602392536e-07（S0.005 与 S0.05 均高于 S0.02）。

### 与既有 DEP 2D 链的对照（C3 族，方向对照）

C3 族母模型几何与 2026-09-25 DEP 2D 链一致（同域 inf 32 50、同地表/覆盖层界面定义、同天线 Tx(15.35,45)/Rx(16.65,45)、同材料假设 εr16/σ0.01 覆盖层与 εr9/σ0.001 砂岩）；本批 C3 族 D5/D10/D20 包络峰到时（218.6 / 331.9 / 536.9 ns）随深度递增，与既有链 DEP_05/10/20 的方向一致。[dep3d_gold 结果](2026-09-26_dep3d_gold_results.md)中同链 2D 对照 B2D5CM（dy=dz=50 mm）的包络峰到时 539.8960 ns、谱峰 162.8 MHz，与本批 C3-D20m 的 536.9476597610657 ns 同量级。规格书 §0 限制照录：旧 DEP_05/10/20 单点在旧网格（2.5 cm/FINE/ZFINE 系列）上完成，与本批粗档结果**不可混用、不可跨格相减**；此处仅作方向对照，不作数值合并。

## C8 族低量级数值的注意事项

C8 族差分带内能量比整体低 8–12 个数量级：族锚点 D10m 8.026249736550938e-11（−100.95487331330483 dB），D20m-W4m-T0.5m 低至 4.598276710920589e-12（−113.37404897701373 dB）、D20m-W4m-T0.25m 3.8167380710479025e-12（−114.1830764309039 dB）。如实记录：这些比值的绝对量级已接近数值分辨关注区，单看数值不能排除数值底噪贡献。但同族零对比负控 `B2D-C8m-D10m-NC` 的配对差分逐样本恒为零（`zero_difference=true`），表明同格同源配对差分的确定性成立——C8 目标例的非零差分来自模型差异，而非配对运算引入的随机涨落。即便如此，本批不基于 C8 低量级数值声明任何可探测性结论（含"不可探测"）；C8 D20m 档仅作为厚覆盖吸收压力下的机制记录（规格书 §6-P6 的风险诊断定位）。

## 硬限制（`results.json` hard_limits 全字段，逐条照录）

| 字段 | 值 |
|---|---|
| `absolute_accuracy_claimed` | false |
| `clean_truth_generated` | false |
| `conclusions_scope` | "mechanism and operator-evaluation evidence for the batch2d_v1 BASE scenario families only (spec §6-P2/P6); not field performance, not 3D, not absolute accuracy, no maximum-depth claim" |
| `cross_family_differencing` | false |
| `family_bg_required` | true |
| `fdtd_solver_invoked` | false |
| `grid_convergence_certified` | false |
| `grid_tier` | "BASE" |
| `grid_tiers_merged` | false |
| `in_band_energy_ratio_is` | "diagnostic energy ratio on the 501-point analysis grid; not detectability, not SNR, not a physical threshold" |
| `negative_controls_merged_into_targets` | false |
| `physical_acceptance_threshold` | null |
| `physical_label_eligible` | false |
| `reference_state` | "numerically_unresolved" |
| `training_eligible` | false |
| `training_labels_generated` | false |

即：本批不产生绝对精度声明、不产生 clean 真值、不产生物理阈值、不产生物理/训练标签；结论限定于 batch2d_v1 BASE 场景族的机制与算子评价证据，不是实测性能、不是 3D、不含最大深度声明；不跨族差分、不合并网格档、负控不并入目标行；网格收敛未认证；参考状态保持 `numerically_unresolved`；分析脚本仅做后处理复算（`fdtd_solver_invoked=false`）。

## 其他限制

- 粗档 BASE 的物理充分性未认证（规格书 §6-P1/P2）：本批只给出粗档上的机制方向与差分诊断量，粗档结论须细档子集方向一致方可保留，本报告不预支细档结果。
- 不做 3D 与实测外推：规格书 §6-P10 三维校核子集未确认前，任何结论不得表述为三维校核通过或实测可用；二维幅度—深度趋势不得直接迁移到三维。
- 材料参数（覆盖层 εr16/σ0.01、砂岩 εr9/σ0.001、目标族 εr12/20/28 与 σ0.005/0.02/0.05）为研究假设，不是场地实测值。
- OFF 档侧向 PML 余量 8.0 m < 13 m 的例外已单列（见"负控"节），其结果不并入目标行使用。
- 电性档与电导单因子的比较为族内同格差分的方向观察，不构成跨族或跨档的绝对幅值结论。

## 归档与复现

- 运行归档：29 个 run 目录 `artifacts/research_checks/2026-09-26_B2D-*/`，各含原始 h5；全部 29 个 HDF5 的路径与 SHA-256 见 `analysis_r1/results.json` 的 `inputs` 字段（另含契约文件 `configs/research/batch2d_v1/cases.json` SHA-256 `7145d2c735b702591d1cd5b12395b6db42be8c61edac0a3a7cf256aeb917174e`、`groups.json` SHA-256 `a51c00dd36048c2269cd854fc3f758fb36973c54585b7e474f1210eb6524b8df`）。
- 分析与复算：`artifacts/research_checks/2026-09-26_batch2d_analysis_r1/` 与 `_r2/`，两次连跑 `results.json` 逐字节一致。
- vctip 处置日志：`artifacts/research_checks/2026-09-26_batch2d_v1_vctip_intervention_log.jsonl`（29 行，逐事件留痕）。
- 复现入口：`scripts/analyze_batch2d.py`（CPU 后处理复算，`fdtd_solver_invoked=false`，不调用求解器）。
- 29 次 attempt 各自一次，已消耗，不可复用；旧会话不重启。

## 下一步

- **细档子集复核（规格书 §6-P2/P3）**：粗档机制结论须细档方向一致方可保留；细档段（FINE2，每族锚点配对例）需另行冻结契约后执行，本文不预支其结果。
- **P10 三维校核子集**：按规格书 §6-P10 另行敲定场景选取、规模与时机，本批不预设。
- 不生成物理/训练标签，不做 3D/实测外推，不宣布网格收敛。
