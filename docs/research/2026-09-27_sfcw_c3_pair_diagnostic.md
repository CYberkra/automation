# C3官方SFCW配对变化诊断

本诊断只读取正式开发侧C3复包络产物：CO11与MT33各自匹配BG/D10 TGT，使用无taper及200 ns taper版本。没有读取C5/C8、没有重跑S5、没有求解器调用。固定窗口取自先于诊断输出冻结的[窗口契约](../../configs/research/sfcw_c3_pair_diagnostic_v1.json)，依据平面法向传播与名义层深，只是物理时间上下文，不由信号峰值选择。算子只有identity和逐时刻跨道复数均值移除；在完整道集上执行后才切窗。

输出[结果JSON](../../artifacts/research_checks/2026-09-27_sfcw_c3_pair_diagnostic/results.json)含48条窗口/几何/taper/算子指标、12条背景与配对差异taper敏感性，以及264条逐道差异范数记录。六个半开窗在保存的2004个时间样本上恰好覆盖一次；样本数依次为37、108、108、349、120、1282。每窗实际首末采样时刻均在JSON中。taper敏感性同时给出同窗无taper范数分母和全周期无taper范数分母，字段名明确区分；分母恰为零时记录null与原因，不加epsilon。

D10窗口的无taper复数跨道均值移除结果：CO背景残留范数比0.00319614、差异保留范数比0.705002、配对变化误差范数比0.709205；MT对应为0.987785、0.390818、0.920468。retention和error都是范数比，不是能量百分比；平方后才对应离散复样本平方和之比。差异是TGT−BG配对变化，不是clean真值；小的背景残留不认证零背景。D10窗口中，200 ns taper相对无窗的同窗背景变化比CO为1.3224%、MT为1.3534%，配对差异变化比分别为7.67×10⁻⁷与5.72×10⁻⁷；这只描述本批数据与此窗，不能外推到弱差异的可信误差界。

CO在200 ns taper下的D10背景残留范数比进一步降至1.21×10⁻⁹，而差异保留仍约0.705；如此强的残留变化说明不能把渐消后的近零背景作为物理质量保证。MT同条件背景残留仍约0.989。CO是横航迹平移的2D代理，MT是固定发射源的变偏移道集，二者都不是真实沿航迹B扫；不能据此宣布某算法适用于现场。

合成自检通过：完全跨道相同的复数目标被均值算子去除；零分母返回带原因的null；复数数据未取模；合成样例线性算子恒等式的最大误差为0（实际数据最大相对线性残差约2.50×10⁻¹¹）。覆盖测试用独立样例，不是clean认证。主控独立脚本对456个核心数值完成交叉验收（实现比较容差`rtol=1e-11, atol=1e-24`仅用于浮点实现对照），结果通过。再次执行脚本时，输出目录存在保护立即拒绝覆盖；没有重写结果。

复核入口：[实现](../../scripts/diagnose_sfcw_c3_pair.py)（`--self-test-only`不读取数据）、[独立计算](../../artifacts/research_checks/2026-09-27_sfcw_c3_pair_diagnostic/sfcw_pair_root_check.py)、[456值验收](../../artifacts/research_checks/2026-09-27_sfcw_c3_pair_diagnostic/accept_sfcw_pair.py)、[验收记录](../../artifacts/research_checks/2026-09-27_sfcw_c3_pair_diagnostic/root_acceptance.json)。下一步仅分析已有C3背景跨道变化与记录末端的关系，当前单一尾窗对照不足以确定边界/截断原因。无新FDTD、S5、训练，G4不解除。

本批来源SHA：契约`abf641d21602eb6f013b5fa78779664ba8732a90bd409e86e1c14328bf25a170`，输入manifest`00ecdf7c62fabd62b45af7c43aaed125161f9456ea95c20d077c933bb3c6125a`，输入NPZ`52f544f8516dceb9a8573e4f1f35efbb78f6b971d7c315526fb95df4f6a30984`，本诊断脚本`5d40dd215b6ba6942e44aeccfcdf2622cc9abfbd30880c61f199cab5f866b204`，结果JSON`57df57521355bf82cdc4a7124ac3e1ba1996a9e0dbb0ff7bd8734d3dbd4608ba`。窗口契约原始CP936副本保存在ignored local_checks中；编码修复发生在任何诊断结果生成前，窗口数值及算子语义保持不变。另一次初始预检因相邻窗口浮点值不能逐位相等而停止，未创建目录或读取结果；随后按ps边界容差并验证每个保存样本恰落入一个窗口后才执行正式诊断。
