# S4 基线组（I/F/L/R/O）冻结提案 v0.1

日期：2026-09-26。入口：用户 19:10 "继续"（批准按 [eval_batch2d_co 报告](2026-09-26_eval_batch2d_co_eval_v0.1.md) 的下一步推进 S4 基线组开发冻结流程）；18:48 仿真启动自主权与 11:28 推送自主权登记于 [decision_log](decision_log.md) 置顶条目。本文档为纯设计+冻结提案，**不运行任何基线、不产生任何分数**；批准后由 `scripts/freeze_s4_baseline_group.py` 程序化冻结。

## 1. 冻结什么、不冻结什么（G4 未解除的诚实口径）

依据 [评价与标签协议 v0.2](2026-09-24_evaluation_and_labels_v0.2.md) 与 [运行器设计 v0.1](2026-09-26_eval_runner_and_selector_design_v0.1.md) §B：

- **参考状态 `numerically_unresolved` 下，`D/A/H/ρ` 不可计算、可行性全 null 界限**——因此 F 的单配置选择、R 的门槛选择、L 的 λ/w/q 选择在**质量指标上全部 undetermined**，本轮**不冻结任何配置选择**。
- 本轮冻结的是：**搜索空间、选择协议、开发/测试组划分、L/R 的评价口径**。任何"哪一个是 F"的答案在 G4 解除前一律不得产出（防止把机制诊断读成选型结论）。
- O（候选内离线最优）需要 clean 真值参考（`clean_truth_generated=false`）⇒ **O 声明不可用、不得以任何方式近似**（设计 §B.2 精神）。

## 2. 五类对象的具体口径

| 对象 | v0.1 冻结口径 | 本轮产物 |
|---|---|---|
| I identity | 27 项目录中的 `B0_G1_BG`；已有 [MT 评价](2026-09-26_eval_batch2d_mt_eval_v0.1.md) 与 [CO 评价](2026-09-26_eval_batch2d_co_eval_v0.1.md) 记录，不重跑 | 引用既有记录 |
| F 最佳固定流程 | 选择协议冻结（§3）；候选集 = 27 项目录；**配置选择 undetermined（G4）** | 开发视图 = 既有评价记录在开发组的投影（无新计算） |
| L 局部均值 | 外部参照网格：λ∈{0.25,0.5,1} × 奇数窗宽 w∈{5,11,21} × q∈{1,2,4}（沿用[传统基线 §4](2026-09-24_cycle03_traditional_baselines.md) 提案值，开发起点非最优）；**合法性规则：`w ≤ n_traces`，道数不够则该组合对该窗不可用，不自动缩短**（MT 33 道全合法；CO 11 道 w=21 不可用）；顺序固定 = 局部均值 → 固定共享增益（BG 类序，无 GB 变体）；评价走 §A 同一套指标函数 | 实跑（两几何全部可评价行） |
| R 输入谱规则 | 规则冻结：对输入窗计算 `d_k=(σ_k−σ_{k+1})/σ_1`（k∈{1,2,3}，与算子契约同一 SVD）；取最大 d_k，若 `max d_k < τ` 则 identity；**门槛候选 {0.05,0.1,0.2} 为未校准开发起点，τ 选择 undetermined（G4）**；保留数值截断拒绝规则（数值秩底限以下并列不计为已解析分量，与 `gap_rtol=1e-8` 同口径）；先 q=1 消融（完整配置 q 选择随 τ 一并 undetermined） | 实跑：逐窗谱诊断 + 三个 τ 候选各自的选中结果（机制诊断，非选型） |
| O 离线最优 | **unavailable_no_computable_reference**：`clean_truth_generated=false`、参考未校准；禁止近似 | 逐行声明不可用 |

## 3. F 选择协议（预注册，G4 解除前不执行选择）

当且仅当 G4 解除（合法效用 U 与容差冻结）后，按以下顺序在**开发组**内选取唯一 F 配置：

1. **保护约束门**：可用性（无 ConfigUnavailable）、增益安全（零截幅/零溢出）、负控窗 `N_b` 不劣于 identity。
2. **受保护事件保真门**：弱对比 E12-S0.005、深部 D20m、薄层 T0.25m、小尺度 W1m/W2m 逐类单列，任一项不通过即该候选出局（v0.2 §4"不按事件能量加权掩盖弱层"）。
3. **排序**：在通过门控的候选上按冻结效用 U 的区间偏序排序；并列允许，证据充分才唯一（v0.2 §7）。
4. 选择结果连同开发证据冻结为 F_v0.2，在**测试组**上只做一次确认性评价。

同等动作集合与调参机会纪律（设计 §B.4）：L/R 若显示额外收益，先扩充动作集合并对全部参照同等的重新开发，不得让任何一方独享更强算子。

## 4. 开发组 / 测试组划分（group_id 单位）

- 划分单位 = `group_id`（母模型族），同族全部变体同组（设计 §C.5）。
- **开发组 = {B2D-C1m, B2D-C3m}**：覆盖受保护事件类别最全的两族（E12 弱对比、T0.25 薄层、W1m/W2m 小尺度在 C3；D20m 在 C3），且含两族满足"稳定差异"判定的最小族数。
- **测试组 = {B2D-C5m, B2D-C8m}**：保留深部/吸收压力族（C5/C8 均含 D20m；C8 含 T0.25m），受保护类别在测试侧仍有覆盖。
- 本划分在 G4 解除前**不参与任何选择**（无选择发生），仅约束报告视图；G4 解除后若要改动划分，须重新冻结并留痕。
- 4 个族的总样本量不支持任何泛化声明（设计 §C.5 警告照录）；CO 数据仅 C3 族，属几何语义层，与开发/测试划分正交。

## 5. 评价口径（与 §A 运行器同一套函数）

- 窗口：同一冻结事件表、同一"原样施加不加道间平移"映射、同一时间轴约定。
- 指标：`N_b`（`energy_metrics`）、增益风险诊断（同曲线定义、BG 序 pre-gain 持久化 + SHA-256）、SVD 谱诊断、资源。**`D/A/H/ρ` 与 `R_c` 维持不可计算/不计算**。
- 语义字段：每条记录携带所在几何的 `gather_semantics.version`（MT：`variable_offset_csg_not_common_offset_bscan_v1`；CO：`common_offset_bscan_v1`）；CO 与 MT **分层报告、不合并、不跨几何配对**。
- 大数组持久化到 git 忽略目录带 SHA-256，不内联 JSON。
- 可行性：`configuration_labels` 全 null 界限 ⇒ `partial_only`，不排名不选优。

## 6. 纪律照录

`training_eligible=false`、`training_labels_generated=false`、`physical_label_eligible=false`、`clean_truth_generated=false`、`reference_state="numerically_unresolved"`、物理阈值 null、不跨族/跨档差分、负控不并入目标行、网格收敛未认证、BASE 档。G1（3D/A0 待用户敲定）/G4（阈值）/G5（细档多道）维持未解除。

## 7. 产物清单（批准后执行）

1. `scripts/freeze_s4_baseline_group.py` + `configs/research/s4_baseline_group_v0.1.json`（程序化冻结、整文件哈希、重跑字节一致）。
2. `scripts/run_baseline_s4.py`：L 实跑 + R 规则诊断 + F 开发视图引用 + O 不可用声明；`--output-dir` 必须不存在；r1/r2 两遍验收。
3. 报告 + registry + START_HERE + decision_log 回填；`check_handoff.py` 与 `verify_workspace.py` 全过后提交推送。
