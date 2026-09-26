# batch2d_v1_mt 多道批量结果：22 例 × 33 道全部完成，锚点 22/22 逐位零差分，NC−BG 四族 33 道恒零

按[多道批次提案](2026-09-26_multitrace_batch_proposal.md)（用户批准）与冻结契约 `configs/research/batch2d_v1_mt/` 执行：22 次 FDTD 运行（BASE 档，dy=dz=25 mm），每次输出 33 道（mt01–mt33，y 12.65–20.65 m、道距 0.25 m、z=45 m、仅 Ex、各 20352 样本 float64），全部串行完成，每例 1 次 attempt、无重试、无 CPU 回退（`retries=0`；`max_fdtd_runs=22`，gate `batch2d_v1_mt`），22/22 `reason=completed`、exit 0。20 min/例墙钟与 4 GiB Job 提交内存硬上限均未触及；全批实测合计墙钟 1283.1 s（约 21.4 min）。

来源区分：墙钟、exit、Job 峰值提交内存为监督器实测（各 run 目录 `supervision.json`，经分析脚本并入结果文件）；布局、锚点逐位比对与 NC 恒零核验全部来自分析脚本对已归档 HDF5 输出的实测复算（`results_r1.json`），不是求解器声明。r1 与 r2 两次连跑的 `results.json` 已核验逐字节一致（SHA-256 `501b960ea096a39b73f7043ba73e80f42486cd1537979abb67b0f28209d7bb4d`）。

**道序语义照录（提案 §2.4 强制项）**：本批 33 道是**变偏移距单发多收道集**（Tx 固定 y15.35，Tx–Rx 偏移距沿道由 2.70 m 经 1.30 m（第 17 道）单调增至 5.30 m），**不是常偏移 B-scan**。该限制必须随数据一起传给评价侧；本批结果不得直接外推为实测常偏移 B-scan 上的算子性能。

## 执行记录

- 22/22 `exit_code=0`、`reason=completed`，attempt 各 1 次，无重试、无 CPU 回退。
- 墙钟（监督器实测）：45.156–256.656 s/例；首例 `B2D-C1m-BG-MT33` 256.656 s（含启动初期 vctip 滞留等待，见下节），其余 21 例 45.2–81.3 s。全批合计 1283.1 s。
- Job 峰值提交内存（Windows Job committed memory，非 RSS、非 GPU 显存）：2,956,218,368 – 2,984,017,920 B（约 2.75 – 2.78 GiB），全部低于 4 GiB（4,294,967,296 B）硬上限。
- h5 元数据（22 例一致项）：dtype float64、`nrx=33`、网格 1×1280×2000、dx_dy_dz 0.025 m、dt 5.896635841874211e-11 s、iterations 20352；源 HertzianDipole x @ (0, 15.35, 45) m 与单道母模型逐项一致。

## 布局核验（22/22 通过）

分析脚本逐例逐道核验：`nrx=33`；第 k 道（k=1..33）`Name=mt{k:02d}`、`Position=(0, 12.65+0.25(k−1), 45)` m、`GridPosition=(0, 506+10(k−1), 1800)`、数据集仅 `Ex`、长度 20352、float64。全部 22 例 33 道无例外。

## 锚点逐位零差分（22/22）

既有单道 Rx（y16.65, z45，格 y=666）严格为本批第 17 道（mt17，索引 16）。逐例将 mt17 的 Ex 与已归档单道母模型 h5（`artifacts/research_checks/2026-09-26_<母模型>/<母模型>.h5`，22/22 归档全覆盖）的 rx1 Ex 逐位比较：**22/22 `max_abs_difference=0.0`（bit-for-bit 相等）**。多接收器共存不改变单道数值结果——多道扩展的确定性成立，与提案 §6.4 的预期一致。

## NC 多道恒零（4/4 族）

单道 NC 恒零证据（粗档）在多道几何下重新验证：四族 NC 例与同族 BG 例的全部 33 道逐位比较，**NC−BG 差分逐道恒零**（`nc_minus_bg_mismatched_traces=0`、`nc_minus_bg_max_abs=0.0`，C1/C3/C5/C8 四族各 1 例）。注意口径：NC 道本身非零（其为无空腔背景，含直耦波与覆盖层响应），恒零的是 **NC 对族 BG 的配对差分**，与提案"多道下应逐道恒零"指配对差分一致。多道几何/网格本身不产生伪响应。

## vctip 处置与启动故障修复（单列）

- **启动故障（两次失败未消耗 attempt）**：首批两次启动即失败，根因两处：① 包装脚本为 LF 行尾且 `if (...)` 括号块内 `echo` 直接展开含 `(x86)` 的 `%GPRMAX_VS%`，cmd 解析期炸出"此时不应有 \Microsoft"；② 修复①后 runner 预检发现 `cl.exe`/`nvcc.exe` 不在 PATH，且 pycuda `_driver.pyd` 依赖 `curand64_10.dll`（Python ≥3.8 扩展依赖 DLL 不走 PATH 搜索）。修复：脚本改 CRLF、块内 echo 加引号、包装脚本预挂 CUDA v13.3 `bin`+`bin\x64` 到 PATH、并把 `curand64_10.dll` 复制进 venv 的 `pycuda/` 目录（依赖 DLL 与 .pyd 同目录即可被解析；不动冻结 launcher、不影响门闩锁定的 gprMax 原生模块哈希）。历史脚本 `identify_local_v4_runtime.cmd` 曾引用 `C:\cuda118\bin`（该目录现已不存在），解释了两处环境差异的来源。
- **vctip 处置**：本批未沿用粗档的孤儿判定哨兵（未留逐事件 jsonl 日志，如实记录）；改为执行期间巡检循环每 45 s 无条件 `taskkill` vctip.exe。首例因此滞留至 256.7 s，其余 21 例未再出现可观察的滞留影响（21/21 `waited_for_descendants` 未构成超时风险，全部 completed）。监督器完成语义未修改；求解器进程、输入文件与输出 h5 未被触碰。

## 硬限制（逐条照录粗档结果文档，并按提案 §6.6/§2.4 增补多道字段）

| 字段 | 值 |
|---|---|
| `absolute_accuracy_claimed` | false |
| `clean_truth_generated` | false |
| `conclusions_scope` | "mechanism and operator-evaluation evidence for the batch2d_v1_mt BASE multitrace scenario families only; not field performance, not 3D, not absolute accuracy, no maximum-depth claim" |
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
| `n_traces` | 33 |
| `trace_axis_semantics` | "variable-offset single-shot multi-receiver gather (Tx fixed y15.35, offset 2.70->1.30->5.30 m along trace index); NOT a common-offset B-scan; must be carried to the evaluation side" |
| `anchor_trace` | "mt17 (y16.65 m, grid y666) bit-for-bit equal to archived single-trace mother (22/22)" |
| `nc_multitrace_control` | "NC-BG paired difference identically zero on all 33 traces, all four families (4/4)" |

即：本批不产生绝对精度声明、不产生 clean 真值、不产生物理阈值、不产生物理/训练标签；结论限定于 batch2d_v1_mt BASE 多道场景族的机制与算子评价证据，不是实测性能、不是 3D、不含最大深度声明；不跨族差分、不合并网格档、负控不并入目标行；网格收敛未认证；参考状态保持 `numerically_unresolved`；分析脚本仅做后处理复算（`fdtd_solver_invoked=false`）；变偏移距道集语义必须随数据传递，不外推为常偏移 B-scan 性能。

## 其他限制

- BASE 档物理充分性未认证：本批只提供多道数据形态与确定性/负控证据，机制结论的细档方向一致性与 3D 校核维持既有程序（G1/G4/G5 未解除）。
- 不做 3D 与实测外推；材料参数（覆盖层 εr16/σ0.01、砂岩 εr9/σ0.001、目标族 εr12/20/28 与 σ0.005/0.02/0.05）为研究假设，不是场地实测值。
- 锚点零差分是"多接收器共存不改变单道数值"的确定性/可复现性声明，不是网格收敛声明，也不是精度声明。
- G2 仅解除**数据形态层**阻塞：27 项候选的正式评价仍须评价运行器在已冻结事件表上执行；评价侧需实现变偏移距道集的窗/道语义后再跑候选。
- vctip 本轮处置无逐事件日志（与粗档 jsonl 留痕不同），不构成 vctip 行为的新证据。

## 归档与复现

- 运行归档：22 个证据目录 `artifacts/research_checks/2026-09-26_<run_id>-MT33形式/`（即 `2026-09-26_B2D-*-MT33/`），各含 h5、`.in`、`supervision.json`、`live_status.json`、stdout/stderr 与 `attempt.json`；原始 run 目录在 `artifacts/simulations/2026-09-26_*MT33/`。
- 分析与复算：`artifacts/research_checks/2026-09-26_batch2d_v1_mt_analysis/results_r1.json` 与 `results_r2.json`，两次连跑逐字节一致（SHA-256 `501b960e…d7bb4d`）。
- 复现入口：`scripts/analyze_batch2d_mt.py`（CPU 后处理复算，`fdtd_solver_invoked=false`，不调用求解器）。
- 22 次 attempt 各自一次，已消耗，不可复用；gate `execution_outcome.status=completed`、`approved_to_simulate=false`，批已关闭。

## 下一步

- **评价运行器上多道候选**：G2 形态阻塞已解除（33 道 > svd k=3 秩条件、local_mean 窗宽合法），27 项候选与基线可在已冻结事件表 + 多道道集上执行；评价侧须先实现变偏移距语义（道轴 ≠ 空间 B-scan 轴）。
- **A0 三维校核子集**：提案已落盘待用户批准（G1 维持未解除），与本批相互独立。
- 不生成物理/训练标签，不做 3D/实测外推，不宣布网格收敛；G4（阈值）/G5（细档负控）维持未解除。
