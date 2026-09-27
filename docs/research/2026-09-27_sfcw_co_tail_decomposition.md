# C3 CO 背景记录尾部处理分解

本诊断仅使用既有 C3 CO 背景 11 道 H5 与官方 SFCW adapter 归档，不运行 FDTD、不涉及测试组或训练。固定方案见 [诊断契约](../../configs/research/sfcw_co_tail_decomposition_v1.json)，执行脚本为 [diagnose_sfcw_co_tail.py](../../scripts/diagnose_sfcw_co_tail.py)，数值与审计数组见 [results.json](../../artifacts/research_checks/2026-09-27_sfcw_co_tail_decomposition/results.json) 和同目录 `diagnostics.npz`。

输入是 11 道、每道 20,352 个 Ex 样本，采样间隔约 0.0589664 ns；最后采样时刻为 1200.02436018 ns。全局中心化平方范数分母包含全部 20,352 点。冻结的前五个半开物理窗合计覆盖 20,351 点；名义 1200 ns 后的一个样本单独记为未分配，其中心化平方范数占全局分母约 8.80×10⁻⁶，未移动窗口或从分母中删除。五窗占全局中心化平方范数的比例依次为：0–60 ns 约 1.30×10⁻¹³、60–240 ns 53.5291%、240–420 ns 23.0838%、420–1000 ns 20.8706%、1000–1200 ns 2.51569%。因此观察到的跨道原始差异分布在多个时段，不能归结为只发生在尾部。

按官方 `apply_tail_taper` 的 200 ns 接收端权重定义实际改变支持（W<1），该支持含 3,391 个样本，占全记录中心化平方范数 2.51482%。令原始接收道集为 X、加窗记录为 WX、被移除部分 U=X−WX；对 U 使用相同源记录及官方 `direct_frequency_response` 和矩形逆变换。官方复包络变换下的 F(U) 与归档无窗/加窗结果之差一致：相对 L2 闭合误差约 2.50×10⁻¹³。对应恒等式为 C(E_no)=C(E_tail)+C(F(U))，其中 C 是逐时刻跨道复数去均值；范数平方含交叉项，不能将三者平方范数当作独立份额。

固定 D10 窗（240–420 ns）中心化复包络范数从无窗的 0.0521243892 降至加窗后的 1.98693×10⁻⁸，后者/前者范数比为 3.81190×10⁻⁷。该数值支持一个有限结论：既有 200 ns taper 移除的原始尾部，经同一官方 SFCW 转换后，足以解释此窗复包络跨道残差的大幅变化。它只归因于确定的处理分解；不能说明原始跨道差异由 PML、有限计算域或其他物理/数值机制造成，也不提供 clean truth、误差预算或收敛认证。既有空气模型尾部/PML观察不外推为本层状 CO 数据的物理解释。

合成自检验证复数中心化与线性分解、零分母返回 null 且不加 epsilon。执行命令：`artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe scripts/diagnose_sfcw_co_tail.py`；正式输出目录运行前不存在，脚本拒绝覆盖既有目录。结果 SHA256 为 `2f47380bce151a19c0aa2a0c01da742607d3886ca468822f397fd3bf7e747db5`，NPZ SHA256 为 `b5e57dc2129acc3290987539d9fdb5f8b8b6a096e8e7fc569fe94794fa8506f0`。冻结契约 SHA256 为 `270912bf95392fa7b80e5af38849802fd049544fb105d9b2e09dd3e98fd14c74`。

独立验收核对92个标量，进一步逐位核对 raw time、中心化原始数组与 taper 权重（分别 20,352、223,872、20,352 个元素），并对 5,511 个复频响值和 22,044 个复包络值做显式 DFT及独立逆变换对照，相对 L2 差约 1.29334×10⁻¹⁰；该容差只衡量实现一致性。11 个对应 `.in` 文件除 title/Tx/Rx 外相同，域宽 32 m、网格 0.025 m、侧边 PML 40 格，故本报告使用的 y=1/31 m 内边界坐标有输入依据。独立验收确认未运行求解器或训练，且未认证 PML 物理成因。正式产物生成后再次调用脚本两次，均按预期在读取数据前拒绝覆盖现有目录，产物哈希保持不变。

接续可运行[独立计算](../../artifacts/research_checks/2026-09-27_sfcw_co_tail_decomposition/check_co_tail_independent.py)与[验收脚本](../../artifacts/research_checks/2026-09-27_sfcw_co_tail_decomposition/accept_co_tail.py)；已保存[验收记录](../../artifacts/research_checks/2026-09-27_sfcw_co_tail_decomposition/root_acceptance.json)。下一步起草能区分边界与记录截断影响的最小C3验证方案，附网格精度检查、固定处理规则及GPU预算，先设计不启动；不重跑S5、不训练、不解除G4。
