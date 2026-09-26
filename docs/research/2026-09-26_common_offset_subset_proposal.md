# 常偏移 B-scan 稀疏配对子集提案（batch2d_v1_co，§2.4-B 落实稿）

状态：**待用户敲定的契约提案**（未冻结、未执行）。起草 2026-09-26，依据 [多道批提案](2026-09-26_multitrace_batch_proposal.md) §2.4-B 的预留条款与 [加速选项实测](2026-09-26_speedup_probe_results.md) 结论。

## 1. 目的与唯一研究问题

MT33 变偏移距道集上，直耦波与平层界面存在道间正常时差（直耦波道间时差约 147 样本），可能污染 mean/SVD 类沿道算子的残差结构。本子集只回答一个问题：

> **去掉道间时差（真实常偏移 B-scan 语义）后，同一配对差分上的算子诊断结论是否改变？**

- 规模刻意最小：仅 C3 族两个母模型（`B2D-C3m-BG`、`B2D-C3m-D10m-W4m-T0.5m-E20-S0.02` 即 A0），每例 11 道常偏移剖面，共 **22 次 FDTD 运行**。
- 结论只作"几何敏感性对照"：与 MT33 批**分层报告、不得合并、不得跨几何差分**（常偏移道与变偏移距道永远不构成配对）。

## 2. 采集几何（冻结值）

沿用母模型几何关系：Rx 在 Tx 前方 **1.30 m**（母模型 Tx 15.35 / Rx 16.65），收发对沿 y 同步平移，偏移距恒定——这是真实常偏移 B-scan 语义，与 MT33 的单发多收有本质区别（此处 Tx 每挪一道必须重跑一次 FDTD，22 次运行不可省）。

| 参数 | 值 | 依据 |
|---|---|---|
| 道数 | 11 | 覆盖目标区所需的最小稀疏剖面 |
| 道距 | 1.00 m | Rx y ∈ [12.65, 22.65]，格胞 506+40k（40 道/米，全部整格） |
| 偏移距 | 1.30 m（恒定） | 与母模型一致；Tx y ∈ [11.35, 21.35]，格胞 454+40k |
| 锚点道 | **t05**（Rx 16.65 / Tx 15.35） | 与母模型逐字节同输入 ⇒ 输出应与已归档单道 h5 逐位一致（CUDA double 确定性已由 [加速实测](2026-09-26_speedup_probe_results.md) 附带证实） |
| 目标覆盖 | A0 目标 y19–23 | Rx 19.65–22.65 五道位于目标正上方；t01–t04 为目标外对照道 |
| 分量 / 高度 | Ex @ z=45 m | 与母模型一致 |
| 网格 / 时窗 / PML | BASE 25 mm、1.2 μs（20352 步）、40 格 | 与 batch2d_v1 全项一致，保证可与本族档案同精度配对 |
| 精度 | CUDA **double** | 维持契约；FP32 档未启用 |

Tx/Rx 均在 PML 内边（y=1 m / 31 m）之内：11.35 > 1、22.65 < 31【自算】。

## 3. 输入生成规则

对每个母模型的 11 个道位，仅改两行、其余字节不动：

- `#hertzian_dipole: x inf <tx_y> 45 impulse`
- `#rx: inf <rx_y> 45 measurement Ex`

**锚点例外**：t05 的 .in 与母模型 .in **逐字节相同**（含 #title），冻结时断言 `sha256(CO输入) == sha256(母模型输入)`；其余 10 道 #title 改为 run_id。静态核验断言：除这两/三行外与母模型逐行一致（复用 MT 冻结脚本的 strip 比对法）。

run_id 规则：`<母模型>-CO11-t<kk>`（kk=01..11），如 `B2D-C3m-BG-CO11-t01`。批次目录 `configs/research/batch2d_v1_co/`。

## 4. 启动模式：单进程多例（选项 3，已实测）

- 不再每例开进程。按**块**执行：块 = 同一母模型的连续若干道，块内一个 Python 进程经 `gprMax.run()` 串行执行（官方 CLI `-n` 同一内部路径）。
- 块大小 6/6/6/4（BG 前 6、BG 后 5 + TGT 前 6……具体以冻结脚本为准），单块墙钟预算 ~3.5 min（实测 21 s/例 × 6 + 进程费 ~35 s）。
- 收益（[加速实测](2026-09-26_speedup_probe_results.md)）：每例 45–58 s → ~21 s，本批 22 例预计 **~10 min**（旧模式 ~17 min）；本批同时是选项 3 的首次契约级试点。
- 代价：块内一例崩溃则该例之后同块各例中断（后续块不受影响，attempt 不消耗、可重冻结重跑）。无重试。

## 5. 契约新增字段（只增不改）

- `common_offset`: `{n_traces: 11, trace_spacing_m: 1.0, offset_m: 1.3, rx_y_range_m: [12.65, 22.65], tx_y_range_m: [11.35, 21.35], anchor_trace_id: "t05", anchor_trace_y_m: 16.65, trace_ids: [t01..t11], trace_order: "construction", acquisition_geometry: "common-offset b-scan (Tx-Rx pair translated along y, fixed 1.30 m offset)", rx_output_components: ["Ex"]}`
- 评价侧语义字段（随数据传递，不得省略）：`common_offset_bscan_v1=true`；与 MT33 的 `variable_offset_csg_not_common_offset_bscan_v1` **互斥标注**，跨几何配对在静态核验中直接拒绝。
- 逐字段照录 batch2d_v1_mt 全部 hard_limits（`physical_acceptance_threshold=null`、`reference_state="numerically_unresolved"`、`training_eligible=false` 等），本批不升级任何资格。

## 6. 验收项（执行后必须实测确认）

1. **锚点逐位一致**：t05 两个母模型（BG、A0）的 h5 与已归档单道 h5（`artifacts/research_checks/2026-09-26_B2D-C3m-BG/` 与 A0 对应目录）逐字节零差分；输入 sha 与母模型输入 sha 相等。
2. **配对差分可用**：BG↔A0 同道位配对差分在 ~332 ns 附近呈现与 MT33 批同量级的空腔响应（诊断量比，不作可探测性声明）。
3. **几何正确性**：t01–t04（目标外）与 t08–t11 的差分能量显著低于目标正上方道（分层报告，不设阈值声明）。
4. 每例 exit 0、attempt 各 1 次、Job 峰值 < 4 GiB、单块墙钟 < 20 min。

## 7. 预算

| 项 | 值 | 依据 |
|---|---|---|
| 运行数 | 22（max_fdtd_runs=22） | 2 × 11 |
| 单例输出 | 20352 × 8 B ≈ 0.16 MB | 单接收道 |
| 单例 Job 峰值 | ~2.9 GiB（与 batch2d_v1 同网格同口径） | 实测包络 |
| 批次墙钟 | 规划包络 22 × 150 s；实测预期 ~10 min | 选项 3 实测 21 s/例 |
| 硬上限 | 单块 20 min、Job 4 GiB、输出 2 GiB、retries 0 | 沿用 |

## 8. 本批明确不做

1. 不扩母模型（C1/C5/C8 不入此批）、不加 NC/OFF、不加 FP32 档。
2. 不与 MT33 批跨几何差分/合并；结论不得外推其他族。
3. 不在本批内评价算子（算子评价仍属评价运行器 + G3/G4 范畴，本批只供输入）。
4. 不产生物理阈值、不升级 reference_state、不产生训练标签。
5. 不改动 scripts/research_operator_contract.py 与既有算子契约。
6. 常偏移剖面只用于回答 §1 的几何敏感性问题，不替代 MT33 主数据形态。

## 9. 冻结与执行流程

1. 用户敲定本提案 → 2. `freeze_batch2d_v1_co.py` 生成 22 个 .in + cases/groups/budget/static_check + 新启动器（gate `approved_to_simulate=false`）→ 3. 静态核验与用户批准执行 → 4. 执行（4 块顺序：BG 全部先于 A0，块内 t01→t11）→ 5. 验收项实测回填 → 6. gate 关闭（`_consumed`）并归档。

## 10. 风险与回退

- **锚点不一致**：若 t05 输出与归档单道出现任何字节差，立即停批，查明（预期不可能：同输入同精度同求解器已证确定性）。
- **块内崩溃**：后续块不动；已崩块内已完成各例输出保留但整批标记 partial，由用户决定补跑（新 attempt 记录）。
- **结论若显示几何敏感**（常偏移 vs 变偏移距算子结论不同）：如实分层报告，并在后续批次设计时把采集几何列为独立因素——这正是本批的价值。
