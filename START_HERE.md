# 接续入口：无需聊天上下文

更新：2026-09-26。当前接续工作区为 `E:/automation_djh/automation_repo`（本机）；另一台电脑的历史工作区为 `D:/自动处理`，远端 `CYberkra/automation`。以本页当前状态为准；历史报告中的旧“下一步”不代表待办，本夜（2026-09-26）已作废的旧下一步见“当前状态”段末的作废声明。

## 当前状态（优先于下方历史记录）

**2026-09-26（最新，本机 `E:/automation_djh/automation_repo`）：粗档与细档批量均已完成，P10 三维校核子集决策已落盘、契约待冻结并经用户敲定。**

|单元|状态|报告/证据|
|---|---|---|
|ZFINE2 垂向细化配对|完成（2/2，各 1 次 attempt）|[结果](docs/research/2026-09-25_deep_zfine2_results.md)；`artifacts/research_checks/2026-09-25_DEP_BG_ZFINE2/`、`2026-09-25_DEP_20_ZFINE2/`、`2026-09-25_deep_zfine2_analysis_r1`|
|dep3d_gold_v1|完成（8/8）|[设计](docs/research/2026-09-26_dep3d_gold_design.md)、[结果](docs/research/2026-09-26_dep3d_gold_results.md)|
|batch2d_v1 粗档 BASE|完成（29/29，全批约 23.4 min）|[规格书](docs/research/2026-09-26_batch_2d_spec_v1.md)、[结果](docs/research/2026-09-26_batch2d_v1_results.md)|
|batch2d_v1_fine2 细档|完成（8/8）|[结果](docs/research/2026-09-26_batch2d_v1_fine2_results.md)；分析 `artifacts/research_checks/2026-09-26_batch2d_fine2_analysis_r1`（r1/r2 字节一致）；运行归档 `artifacts/simulations/2026-09-26_*-F2/`；gate `configs/research/gprmax_v4_execution_gate.json`（已消耗并关闭）|
|P10 三维校核子集|决策落盘，**本批零 3D**，未执行|[决策](docs/research/2026-09-26_3d_validation_subset_decision.md)|
|算子评价仿真证据材料 v0.1|完成（纯整理文档）|[证据材料](docs/research/2026-09-26_operator_evaluation_sim_evidence_v0.1.md)；registry 条目 `operator_evaluation_sim_evidence_v0_1`|
|A0 三维校核契约提案|已验收，**待用户敲定**（未执行）|[提案](docs/research/2026-09-26_a0_3d_validation_contract_proposal.md)；registry 条目 `a0_3d_validation_contract_proposal_v0_1`|
|事件表提案 v0.1（G3 前置）|已验收，**待用户敲定**（未冻结）|[提案](docs/research/2026-09-26_event_table_proposal.md)；registry 条目 `batch2d_v1_event_table_proposal_v0_1`|
|评价运行器与选择器设计 v0.1|完成（纯设计文档）|[设计](docs/research/2026-09-26_eval_runner_and_selector_design_v0.1.md)；registry 条目 `eval_runner_and_selector_design_v0_1`；关键结论：单道输入下 27/27 候选全部被门禁挡下，G2 多道为总阻塞点|
|多道场景族批次设计提案 v0.1（解除 G2）|已验收，**待用户敲定**（未执行）|[提案](docs/research/2026-09-26_multitrace_batch_proposal.md)；registry 条目 `multitrace_batch_proposal_v0_1`；22 例 33 道、预算 ≤65 min|

**2026-09-26 batch2d_v1_fine2 细档完成（8/8，结果已归档可复算）。** FINE2 档（dy12.5/dz3.125 mm，与 ZFINE2 同格）8 例 = 四族锚点配对（C3 族先跑，再 C1/C5/C8），串行、失败即中止；gate `batch2d_v1_fine2`：墙钟 40 min/例、Job 20 GiB、CUDA double、retries=0、无 CPU 回退。用途是规格书 §6-P2/P3 的粗档机制方向一致性复核。分析脚本 `scripts/analyze_batch2d_fine2.py` 已完成真实分析与 r1/r2 复算（`results.json` 字节一致）。实测结果（全部来自已归档输出与监督器记录，可复算）：8/8 `exit_code=0`、`reason=completed`，各 1 次 attempt、无重试、无 CPU 回退；墙钟 1318.437–1323.031 s（未触及 40 min/例上限），Job 峰值提交内存 15,246,221,312–15,248,601,088 B（约 14.2 GiB，未触及 20 GiB 上限）。复算校验：`B2D-C3m-BG-F2` 与归档 `DEP_BG_ZFINE2` 仅 `#title` 行不同，接收与源波形 118665 样本逐位零差（`max_abs_difference=0.0`），只证明该档运行间可复现性，不是收敛或精度结论。P2 方向一致性（按族分层，无通过阈值）：四族差分谱形状相关 0.99813–0.99955、逐频点符号一致率 88.0%–93.8%、包络峰细档一致早 15.7–19.1 ns（系统性档间偏移）；细档带内能量比 C1 7.87e-06 / C3 3.09e-07 / C5 1.21e-08 / C8 9.37e-11，族序与粗档同向，C8 仍处数值分辨关注区、不作可探测性声明。结论：粗档三项机制结论（带内能量比随覆盖厚度递减、包络峰到时递增、电性对比增强方向差分上升）获细档方向一致支持，无机制需按 §6-P2 降级为“待复核”，保留为机制/算子评价证据，不是物理充分性认证，不得外推为粗档精度达标。vctip 哨兵随批运行但零干预（8 例 `waited_for_descendants=false`、无处置日志），仅为事实记录，不构成 vctip 行为机制的新证据。分析 r1/r2 `results.json` 字节一致（SHA-256 前 16 位 `dd2ba3328a6f24b1`）。`reference_state` 保持 `numerically_unresolved`，网格收敛未认证，无物理/训练标签、无 clean 真值。[结果](docs/research/2026-09-26_batch2d_v1_fine2_results.md)，证据 `artifacts/research_checks/2026-09-26_batch2d_fine2_analysis_r1`（`_r2` 同）。

**2026-09-26 batch2d_v1 粗档完成（29/29）。** 粗档 BASE（dy=dz=25 mm）29 例全部 exit 0，每例 1 次 attempt、无重试、无 CPU 回退；全批合计墙钟 1400.62 s（约 23.4 min），Job 峰值提交内存 2,914,930,688–2,918,010,880 B（约 2.72 GiB），未触及 20 min/例与 4 GiB 上限。四族零对比负控 NC 的配对差分逐样本恒零（`zero_difference=true`）；带内能量比随深度单调递减（C1/C3 三档严格单调、C5/C8 两档同向）、包络峰到时随深度递增；带内能量比是 501 频点诊断网格上的诊断性能量比，**不是可探测性、不是 SNR、不是物理阈值**。分析 r1/r2 字节一致。粗档物理充分性未认证，结论须细档方向一致方可保留（规格书 §6-P2；该条件已于 2026-09-26 由细档满足，见上段，保留为机制/算子评价证据）。求解后约 900 s 滞留的根因在同环境下鉴定为 `vctip.exe`（MSVC 编译器遥测孤儿进程）；处置属**操作层干预**：外部哨兵仅在“监督器处于 `waiting_descendants` 阶段且 vctip 已成孤儿”两条件同时满足时清除该进程，29 次事件逐条留痕于 `artifacts/research_checks/2026-09-26_batch2d_v1_vctip_intervention_log.jsonl`，不改监督器完成语义、不触碰求解器进程/输入/输出 h5。`reference_state` 保持 `numerically_unresolved`，无物理/训练标签、无 clean 真值。[结果](docs/research/2026-09-26_batch2d_v1_results.md)

**2026-09-26 dep3d_gold_v1 完成（8/8）。** 4 个 2D 同格距对照 + 2 组 3D 配对，各 1 次 attempt、exit 0、CUDA double、无重试、无 CPU 回退；墙钟 910.6–920.2 s，Job 峰值提交内存约 2.32–14.98 GB（<20 GiB 上限）。关键比较值（分析脚本对归档 h5 的复算，r1/r2 字节一致）：3D 纯 dz 细化（5 cm→2.5 cm）谱形状相关 **0.9951**；2D/3D 同格距归一化谱形状相关 **0.935**（B2D5CM vs DEP3D_5CM）/ **0.958**（B2DANISO vs DEP3D_ANISO）；远场 3D→2D 变换后 **0.967 / 0.979**（探索性证据单列，不作结论依据）。限制：3D 网格收敛未认证、不跨维比较绝对幅值、相位只以线性拟合形式量报告且不跨维迁移、维度/侧向域宽/目标 y 位置三类混淆因素未排除；无物理阈值、无训练标签、无 clean 真值，结论限定于 dep3d_gold_v1 场景族与所分析网格。[设计](docs/research/2026-09-26_dep3d_gold_design.md)、[结果](docs/research/2026-09-26_dep3d_gold_results.md)

**2026-09-26 P10 三维校核子集决策已落盘（未执行，本批零 3D）。** 本批（含 fine2）不含任何 3D；首选 3D 校核对象为 **C3 族锚点 A0 配对**（`B2D-C3m-D10m-W4m-T0.5m-E20-S0.02`，5 cm 各向同性档，2 次运行），前置条件（fine2 完成且方向一致性复核落盘）已于 2026-09-26 满足，契约仍须另行冻结并经用户敲定后方可执行；C3-D20m 类不新增（dep3d_gold 已实测其 3D 对应物），NC/OFF 负控不列入；C8 锚点（P-A）、C5 锚点（P-B）现已可依 fine2 复核结果评估，电性/电导档（P-C）还依赖 A0 3D 结果。预算只引 dep3d 实测（910.6–920.2 s/run、7.98 GiB / 13.95 GiB），不是 ETA 承诺；vctip 处置对 3D 墙钟无实测记录，不下调 3D 预算。[决策](docs/research/2026-09-26_3d_validation_subset_decision.md)

**ZFINE2 垂向细化配对已完成。** `DEP_BG_ZFINE2`/`DEP_20_ZFINE2` 各 1 次 attempt、exit 0（墙钟 1438.6 s / 1334.8 s，Job 峰值 15.28 GB，网格 1×2560×16000 = 40,960,000 单元）。2D 网格收敛链全带相对 L2 变化 **66.6647% → 17.2452% → 4.1216%**，最大相位变化 64.8431° → 16.2741° → 3.8847°；最近两步仅细化垂向，比值约 4.18 / 4.19，是垂向二阶行为的直接指示。未做 Richardson 外推，相邻差不等于相对精确解误差；分析 r1/r2 字节一致。`reference_state` 保持 `numerically_unresolved`，`physical_label_eligible=false`、`training_eligible=false`，物理阈值仍为 null。[结果](docs/research/2026-09-25_deep_zfine2_results.md)

**下一步（按序，取代本页所有旧下一步）：** 已完成项（2026-09-26）：① fine2 验收（8/8 exit 0，墙钟/Job/attempt 记录核对，无重试、无 CPU 回退）[已完成]；② 原始输出归档（8 个 run 目录含 h5、输入、stdout/stderr、`supervision.json`，哈希登记于 `results.json` 的 `inputs`）[已完成]；③ `scripts/analyze_batch2d_fine2.py` 分析 r1/r2（`results.json` 字节一致）[已完成]；④ 跨档方向一致性结论（四族方向一致，无机制按规格书 §6-P2 降级）[已完成]；⑤ batch2d_v1 细档结果报告落盘 [已完成]。当前待办：① 将粗档/细档分层机制结论并入算子契约评价的输入材料 [已完成，见上表证据材料 v0.1]；② 按 P10 依赖链起草 3D 校核契约 [已完成提案并验收，**待用户敲定**]；③ 事件表提案 v0.1 [已完成提案并验收，**待用户敲定**]；④ 评价运行器与配置选择器原型设计 [已完成，见上表设计 v0.1；关键发现：单道输入下 27/27 候选全部被 `min(shape)>=2` 门禁挡下，G2 多道场景族为一切下游的总阻塞点]；⑤ G2 多道场景族批次设计 [已完成提案并验收，**待用户敲定**]。**当前全部待办均等待用户决定**：② A0 3D 校核、③ 事件表冻结、⑤ 多道批次三项提案均已备妥，用户敲定任一项后即可按各自冻结流程推进（首选 C3 族锚点 A0 配对 `B2D-C3m-D10m-W4m-T0.5m-E20-S0.02`，5 cm 各向同性档、2 次运行；契约另行冻结并**须用户敲定后方可执行**；P-A（C8 锚点）/P-B（C5 锚点）现已可依 fine2 复核结果评估，P-C（电性/电导档）仍依赖 A0 3D 结果）。

**旧下一步作废声明：** 本页 2026-09-25 及更早段落中的“下一步”（含“继续准备垂向细化资源方案”、“下一步建立算子误差预算并继续准备垂向细化资源方案”、“预算 dy12.5mm、dz3.125mm 配对”等）已由上一步列表取代，不再作为待办——垂向细化已由 ZFINE2 完成、算子误差预算已由 2026-09-25 单元完成。历史结果本身保留，不改写。


**2026-09-25 算子误差预算证书完成（纯数组，最新）。** 新契约 `configs/research/error_budget_contract_v0.1.json` 与实现 `scripts/research_error_budget.py`：方向诊断预算只用精确构造、`math.fsum` 精确重算加计数舍入 slack、Wedin/Davis–Kahan sin Θ 标准摄动界，禁止 eps×常数。损伤试验 91 个余弦评价恢复 51 个（40 个因输出范数≤预算合规缺失，全均值/SVD 类为主）；三网格 84 个 FDTD 模板余弦保持缺失，原因细化为 `fdtd_numerically_unresolved`。回归 157→182 项全绿。预算只认证数值可辨识性，物理/训练标签仍关闭。[报告](docs/research/2026-09-25_error_budget.md)，证据 `artifacts/research_checks/2026-09-25_error_budget/`。下一步：继续准备垂向细化资源方案；FDTD 来源方向诊断在参考收敛前保持关闭。

**2026-09-25 2D/3D 配对批次 dim23_pair_v1 已完成（路径 A）。** 本机 V4/CUDA 环境当日重建并核验（Python 3.12 + gprMax 4.0.0 + pycuda 2026.1，RTX 3060 Laptop，CUDA double kernel 通过）；中文路径导致的 nvcc/cl 编译失败用 junction `D:\gprmax_v4_gpu_env`、`C:\cuda118` 解决。4 算例（2D_BG/2D_TGT/3D_BG/3D_TGT，5cm 同格距，3D 3,120 万单元）各 1 次 attempt 全部 exit 0，3D 每次约 8 分钟墙钟。配对差分：到时差 9.75ns，波形相关 0.716，归一化谱形状相关 0.902（远场 3D→2D 变换后 0.954），3D 差分谱显著高频化；未做跨维度绝对幅度对比。分析复算字节一致。gate 已关闭，attempt 已消耗。[结果](docs/research/2026-09-25_dim23_pair_results.md)，证据 `artifacts/research_checks/2026-09-25_dim23_pair/`。

**本机已接续并完成三项审查修复。** 比较脚本另存且拒绝覆盖；监督器等待进程组清空；余弦诊断缺少数值误差预算时标为缺失。157项数组、10项监督器检查通过，42行共同层/目标保真结论保持。见[修复及评价增补](docs/research/2026-09-25_direction_diagnostic_update.md)与[验证记录](artifacts/research_checks/2026-09-25_review_fixes/record.json)。下一步建立算子误差预算并继续准备垂向细化资源方案；没有新正演或训练。

**执行契约状态（2026-09-26，本句取代此前所有“当前无待执行仿真契约”及“唯一在执行的是 batch2d_v1_fine2”表述）**：当前**无**待执行且无在执行契约。`batch2d_v1_fine2` 已 8/8 完成、8 次 attempt 各 1 次已消耗，gate 关闭（`batch_id=batch2d_v1_fine2_consumed`、`approved_to_simulate=false`、`execution_outcome.status=completed`，墙钟/Job/验收摘要已回填）；`dep3d_gold_v1` 与 `batch2d_v1` 粗档契约同样已消耗并关闭。P10 三维校核子集尚无冻结契约，须另行冻结并经用户敲定后方可执行。YFINE/ZFINE2 等更早契约保留为历史，不因后续代码修复而改写原授权快照。用户自主 GPU 研究授权保留；本机 V4/CUDA 环境已按上述两个批次的 CUDA double 执行核验。旧路径、PID、会话号与硬件耗时属于原机器。

**2026-09-25 2D/3D 配对调研完成（未执行）。** 仓库从未做过 2D/3D 配对数值对比；本机实测 RTX 3060 Laptop（约6GiB显存）/15.8GiB 内存/Python3.10，无 V4 环境与源码，当前跑不了求解，但 ≤5,000万单元的小 3D 配对经估算可行。见[调研与设计草案](docs/research/2026-09-25_2d3d_paired_plan.md)：路径 A/B/C 待用户选择；契约未冻结，gate 仍为 false。算子误差预算（纯数组）不依赖环境，可立即推进。

**YFINE两组已完成，无本批求解任务运行。** GPU共36分7秒；相对ZFINE，目标全带复数变化0.2103%、最大相位变化0.203°，耗时约2.49倍。原始输出已归档，分析复算一致。[结果及下一步](docs/research/2026-09-25_deep_yfine_results.md)。当前平层场景优先保留横向12.5mm，下一步预算垂向3.125mm配对；尚未构建或启动。旧会话68122与21056均已结束，不重启attempt。

三网格CPU复核已完成42项处理机制与33项不匹配负控；旧默认和新结果复算一致。[最新报告](docs/research/2026-09-25_three_grid_mechanisms.md)。原[垂向细化](docs/research/2026-09-25_deep_zfine_results.md)方向有效，但20m全频带物理真值仍未认证；禁止物理/训练标签，不将模板当真实航线B扫。

共同层界联合检查完成42项，复算一致：干扰能量下降可伴随地层删除，不能仅凭压背景选算法。[结果](docs/research/2026-09-25_joint_roles_results.md)。

背景独立参考检查完成：FINE→YFINE均匀细化使总场L2误差约降低四倍；这不代表20m目标差分已收敛。[背景检查](docs/research/2026-09-25_deep_background_reference.md)。

## 已完成的二维研究（原报告中的“下一步”是历史记录）

|研究单元|已确认的结果|证据|
|---|---|---|
|维度策略|YZ保留高度及横向基线，沿航线不变；2D初筛配合代表性3D验证|[维度策略](docs/research/2026-09-25_dimension_strategy.md)|
|空气/PEC校准|二维线源与官方SFCW链已核验|[校准](docs/research/2026-09-25_yz2d_results.md)|
|有损平层|独立分层参考；细网格复数误差由2.21%降至0.558%|[平层](docs/research/2026-09-25_yz_layers_results.md)、[细化](docs/research/2026-09-25_yz_fine_results.md)|
|5/10/20m异常|深度差分初筛完成，不代表现场可探测性|[深度](docs/research/2026-09-25_yz_depth_results.md)|
|扩大域及均匀细化|扩大域影响极小；20m全带相位对格距敏感|[对照](docs/research/2026-09-25_deep_controls_results.md)|
|参考资格|数值未收敛，物理/训练标签关闭|[资格](docs/research/2026-09-25_deep_eligibility.md)|
|垂向细化ZFINE|背景及目标均完成；相对FINE最大相位变化16.27°，仍未认证收敛|[结果](docs/research/2026-09-25_deep_zfine_results.md)|
|构造机制与负控|共同目标可被误删；配对保真不保证无伪响应|[三网格](docs/research/2026-09-25_three_grid_mechanisms.md)、[联合角色](docs/research/2026-09-25_joint_roles_results.md)|

旧HSG校验已完成：[可行性](docs/research/2026-09-24_hsg_feasibility.md)、[小型对照](docs/research/2026-09-25_hsg_smoke_results.md)、[时钟对照](docs/research/2026-09-25_hsg_clock_results.md)。内部源/耦合时钟尚有专项问题，不阻塞当前普通二维研究。历史会话21056已结束，不应重启。

## 已完成的仿真与参考链

- [垂向4cm配对](docs/research/2026-09-24_dz4_results.md)：资源预算、两次GPU求解、独立参考对比和复算完成；无需重复。

- [GPU环境](docs/research/2026-09-24_gpu_environment.md)：V4.0.0独立CUDA double环境与复建说明。FDTD必须GPU，禁止CPU回退；官方后处理和独立频域积分可用CPU。当前求解器Python为`artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe`。
- [全空气与官方SFCW](docs/research/2026-09-24_M00_sfcw_results.md)、[尾窗/PML对照](docs/research/2026-09-24_tail_control_results.md)：有限窗静电残留机制及渐消敏感性已检查。渐消参数不能直接移植到20m晚到模型。
- [初始介质界面](docs/research/2026-09-24_interface_results.md)、[PEC严格镜像参考](docs/research/2026-09-24_PEC_results.md)：发现约5.4°高频长路径相位误差，原始介质提前分量当时未归因。
- [固定PML干预](docs/research/2026-09-24_fixedPML_results.md)：0–80ns差分从0.0482V/m变为逐样本零，支持PML自动参数差异机制；既定渐消后带内新旧最大变化约0.00037%。这一项已完成。
- [独立半空间参考](docs/research/2026-09-24_halfspace_reference_results.md)：empymod1.10.6官方全波内核加本项目分段积分；只做校验，不替代官方SFCW。2.6导入失败和默认积分失败保留。

所有已运行attempt均已消耗，不重复执行。原CPU r1/r2是主动停止的历史试跑，无完整结果；历史“正在运行”文字不可用来重启它们。[当前执行记录](configs/research/gprmax_v4_execution_gate.json)随当前YFINE批次阶段更新；每个attempt仅执行一次，已消耗记录不能复用。

## 已确认的项目事实

- 雷达：**SFCW 20–170 MHz，步长0.3 MHz**。等间隔且含两端点时推导为501个频点、周期重构时窗约3.333 μs；尚未核对设备频率向量，不能因此认定旧CSV的501个时域样点与频点逐项对应。
- 收发几何：用户补充工作高度约15m、天线间距1.3m；暂按离当地地面高度理解。基线已确认相对飞行方向左右排列（横跨航线）；极化与高度起伏未知。设备图片已补充并归档，标注收发距1300mm，但未标极化。旧0.5m/0.4m校准输入已失效，须重设计空气域、横向域及时窗。
- 地层背景：地表覆盖层**可以是粉质粘土**，地下基岩**一般为砂岩**；不是所有场景都固定如此。厚度、含水状态、风化/裂隙及复介电参数待确定。项目是浅层非显性滑坡，“20m浅”暂解释为地表以下约0–20m，精确范围待确认。
- 首版研究范围：背景与相干杂波抑制、增益与衰减补偿；不扩大为自动滑坡识别或成像。
- 现有实测测线留给后期，当前不读取用于训练/调参/ROI/实测杂波混合。
- **用户已明确允许仿真并授权自行设计；FDTD必须GPU。** 每次执行仍须先落盘具体参数、资源上限和单次attempt，当前范围以全局gate与最新设计包为准；旧M00单次授权不能覆盖后续执行记录，也不能否定后来的自主设计授权。
- 用户要求主动归档、提交和推送，确保跨地点、跨会话接续；不是只在聊天中汇报。

## 已经做完，不要重复从头开始

1. [评价与标签 v0.2](docs/research/2026-09-24_evaluation_and_labels_v0.2.md)：指标/缺失/偏序集合已有代码，物理阈值仍为 null。
2. 27 项构造数组算子；本轮候选能力只比较其中 7 项 q=1 背景配置。
3. [损伤试验](docs/research/2026-09-24_damage_pilot_findings.md)：63 组损伤、12 组弱事件删除、42 行候选；72 项机制检查通过。完整有效交付是 `artifacts/research_checks/2026-09-24_damage_pilot_r2/`，不是中间绘图失败目录。
4. [V4 阅读报告](docs/research/2026-09-24_gprmax_v4_review.md)：官方有关章节与本地源码已对照；19 份文件指纹、10 个官方页面来源/缓存哈希已归档。

上述四项是早期代码/数组/来源阅读证据；后续FDTD结果见本页研究归档。尚无训练集或实测性能结论。

## 关键判断与容易重复犯的错误

- 全事件平均可掩盖局部删除，需要固定局部诊断；D²=A²+H²，不应当成独立三份证据。
- 参考不确定性会改变可行性；例子中的 0.1 不是物理阈值。
- 现有整窗均值/SVD 对平层可能一起删除有用层。尚未证明选择器有超出固定配置/简单规则的稳定收益。
- V4 源码要求 Python 3.11–3.13；当前本机名为 `gprmax` 的旧环境包元数据为 V3.1.6/Python3.10。上轮借其绘图库不等于使用过 V4。
- 本地源码实际在 `E:\gprMax-v.4.0.0\gprMax-v.4.0.0`，不随本仓库上传。新电脑不必复建盘符，提供自己的源码路径并核对台账哈希。
- 内置 80 MHz Ricker 中心约 17.68 ns，旧契约写 37.5 ns；使用前须明确选择并修订时序，不偷偷替换。
- 旧 P1 只有约 7 m 地下域和 240 ns 请求窗，仅用于小例；不能说已覆盖 20 m 项目。早期 45 次正演/18 个母模型/1476 次 A-scan 都不是授权。
- V4 `.h5`、接收器身份、源+dt/2 与场0/−dt/2、SFCW 频域处理须按当前元数据；不能套 V3 路径和默认时间轴。

## 新地点的最短接续步骤

1. 获取最新main，检查git status，保留已有未提交改动。原电脑曾有设备图片缺失的本地状态；本次接续图片存在且交接检查通过，不把旧缺失记录当当前状态。
2. 先读本页当前状态、[YFINE设计](docs/research/2026-09-25_deep_yfine_design.md)与[ZFINE结果](docs/research/2026-09-25_deep_zfine_results.md)；需要背景再查[V4阅读](docs/research/2026-09-24_gprmax_v4_review.md)、[粉质粘土—砂岩简报](docs/research/2026-09-24_cover_sandstone_brief.md)、[损伤研究](docs/research/2026-09-24_damage_pilot_findings.md)。保留实测数据不用于当前开发。
3. `python scripts/verify_workspace.py`为182项数组检查（八组：10/23/5/26/72/21/20/5），不调用FDTD/实测/网络。无相关改动或疑点不重复跑。`python scripts/check_handoff.py`检查文件、链接、历史来源提交和当前执行契约。
4. 环境不随Git上传。GPU环境按上述说明复建；独立参考环境按最新报告及锁定依赖复建。不要把本机成功当成新机器核验，也不要调用可能启动FDTD的示例作导入检查。
5. 每个研究单元保留设计、数据、代码、日志、来源/hash、失败与限制，更新本页/进度/运行索引，再显式暂存、提交、推送。跨地点并行用codex/前缀分支，不强推或重写共享历史。

早期[远端克隆验证](artifacts/research_checks/2026-09-24_remote_handoff_verification.json)对应提交251986e，只证明当时86文件与136项数组检查，不是当前全部产物的新机认证。

可交给新会话的任务：先读本页“当前状态”段的 2026-09-26 内容。**旧句中的“下一步依据实测资源预算 dy12.5mm、dz3.125mm 配对”已作废**——该配对即 ZFINE2，已完成；当前待办是：将粗档/细档分层机制结论并入算子契约评价的输入材料；按 P10 依赖链起草 3D 校核契约（须用户敲定后方可执行）。fine2 已 8/8 完成，不重复启动，已消耗 attempt 不复用；同时保留共同层界/目标/无目标负控的处理评价。物理真值与训练标签仍关闭，保留实测数据不参与开发。

## 记录索引

- [研究产物/运行索引](docs/research/research_artifact_registry.json)
- [进度](docs/research/continuous_research_status.md)、[决策](docs/research/decision_log.md)、[Git/环境指南](docs/WORKFLOW.md)
- [V4源码台账](artifacts/research_checks/2026-09-24_gprmax_v4_source_audit.json)、[早期资源预算](docs/research/2026-09-24_airborne_3d_budget.md)
- 凭据只由系统凭据管理器管理，不写入仓库、日志或远端URL。
