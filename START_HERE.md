# 接续入口：无需聊天上下文

更新：2026-09-25。当前接续工作区为 `D:/自动处理`；另一台电脑历史工作区为 `E:/automation_djh/automation_repo`，远端 `CYberkra/automation`。以本页当前状态为准；历史报告中的旧“下一步”不代表待办。

## 当前状态（优先于下方历史记录）

**本机已接续并完成三项审查修复。** 比较脚本另存且拒绝覆盖；监督器等待进程组清空；余弦诊断缺少数值误差预算时标为缺失。157项数组、10项监督器检查通过，42行共同层/目标保真结论保持。见[修复及评价增补](docs/research/2026-09-25_direction_diagnostic_update.md)与[验证记录](artifacts/research_checks/2026-09-25_review_fixes/record.json)。下一步建立算子误差预算并继续准备垂向细化资源方案；没有新正演或训练。

当前无待执行仿真契约；YFINE已消耗契约保留为历史，不因修复代码而改写原授权快照。用户自主GPU研究授权保留，本机V4/CUDA环境需独立核验。旧路径、PID、会话号与硬件耗时属于原机器。

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
3. `python scripts/verify_workspace.py`为157项数组检查，不调用FDTD/实测/网络。无相关改动或疑点不重复跑。`python scripts/check_handoff.py`检查文件、链接、历史来源提交和当前执行契约。
4. 环境不随Git上传。GPU环境按上述说明复建；独立参考环境按最新报告及锁定依赖复建。不要把本机成功当成新机器核验，也不要调用可能启动FDTD的示例作导入检查。
5. 每个研究单元保留设计、数据、代码、日志、来源/hash、失败与限制，更新本页/进度/运行索引，再显式暂存、提交、推送。跨地点并行用codex/前缀分支，不强推或重写共享历史。

早期[远端克隆验证](artifacts/research_checks/2026-09-24_remote_handoff_verification.json)对应提交251986e，只证明当时86文件与136项数组检查，不是当前全部产物的新机认证。

可交给新会话的任务：先读本页及YFINE结果。两组已归档并复算，不重复运行旧attempt。下一步依据实测资源预算dy12.5mm、dz3.125mm配对，冻结输入与上限后才启动新批次；同时保留共同层界/目标/无目标负控的处理评价。物理真值与训练标签仍关闭，保留实测数据不参与开发。

## 记录索引

- [研究产物/运行索引](docs/research/research_artifact_registry.json)
- [进度](docs/research/continuous_research_status.md)、[决策](docs/research/decision_log.md)、[Git/环境指南](docs/WORKFLOW.md)
- [V4源码台账](artifacts/research_checks/2026-09-24_gprmax_v4_source_audit.json)、[早期资源预算](docs/research/2026-09-24_airborne_3d_budget.md)
- 凭据只由系统凭据管理器管理，不写入仓库、日志或远端URL。
