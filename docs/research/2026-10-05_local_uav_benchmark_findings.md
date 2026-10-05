# 本机轻量前向自检：资源、结果与未完成部分

日期：2026-10-05。源码起点 `8d72ace`。用户授权“这三轮需要很多资源吗？若是本地方便快速跑完就本地做完吧”。本单元实际执行 V4.0.0 双精度 FDTD，并保留失败与资源拒绝；没有训练、拟合营山材料、修改旧 gate/G4 或重用已消耗 attempt。

## 资源结论与交付范围

基础自检适合本机，完整高航高三维验证不属于同一个资源量级。当前机器为 RTX 3060 Laptop 6GiB、约16GiB RAM，运行期间可用内存约1–2GiB且波动；不是历史4090/63GiB机器。默认 Python 的 gprMax 是 V3.1.6，本轮仅使用 `D:\gprmax_v4_gpu_env\Scripts\python.exe` 的 V4.0.0。该环境为链接目录，契约记录解析后的真实 Python 路径，运行源码/编译扩展另有哈希。

| 工作 | 实际执行 | 测得开销 | 状态 |
|---|---|---|---|
| 第一轮空气三维偶极 | 10cm，x/y极化，2道×4接收方向 | GPU监督墙钟合计209.75s；峰值自有进程RSS约0.387GiB | 基础复幅相检查通过；5cm细化和互易性未执行 |
| 第二轮正入射平层辅助 | 5cm宽域5道、1.25cm窄域5道、5cm窄域5道 | CPU监督墙钟合计595.82s；最大Job提交内存约0.224GiB | 细网格8项解析诊断全部通过；粗网格失败保留 |
| 接收/频谱口径辅助 | 复用第一轮复谱，等电流元/等自由空间辐射功率对照 | CPU数秒，无新求解 | 假设敏感性完成，未建立真实设备端口模型 |
| 悬空点源平层、横向布局、起伏与深度 | 未执行 | 示例16×16×24m@5cm约4915万单元，基础主机数组下界约4.30GiB，另需PML/色散/上下文/构建内存 | 不属于本机当时可快速完成的轻量检查；尺寸本身也未通过域收敛验证 |

监督墙钟含进程启动，GPU还含CUDA编译；不等于本会话总耗时。GPU两道日志中的实际求解约0.35s/道，编译占主要时间。CPU峰值是**Job提交内存**，GPU表格为**进程树RSS**，两者定义不同。

共归档22份新native H5，其中5份为失败的辅助三维平面波配置，**17份属于稳定的空气/平层检查**。不能将这个数量描述成22份有效训练样本或“完整三轮通过”。三维空气细网格两份准备契约在RAM/VRAM预检处拒绝，没有STARTED、H5或attempt消费；拒绝只描述当次可用资源，不证明机器永远不能运行该空气小模型。

总索引：`artifacts/research_checks/2026-10-05_local_uav_benchmark_archive_r1/summary.json`，由 `scripts/summarize_uav_local_benchmarks.py` 独立核验输入、native输出与冻结源码身份后生成。

## 第一轮：源、坐标和复相位没有显示明显基础错误

空气8×8×8m、10cm网格、六面1m HORIPML，100MHz Ricker；x/y水平偶极，收发1.3m，四接收点为±x/±y。Hertzian源实际电流元长度随极化方向网格变化，源幅度设为1/dl，保持电流元峰值1A·m，频响报告 **E/(Iℓ)**。不是端口S21，也不是恒设备发射功率。

包含1/r、1/r²、1/r³全部项的连续偶极解析场直接比较20–170MHz的501个复频点，没有经验幅度/相位拟合。实际源与接收器半步时钟用于官方direct链；另用7个频点直接求和复算。

- 8个共极化响应的相对复L2误差：0.6449%–1.1103%。
- 最大相位误差：0.5635°；幅值比误差约−0.0564至+0.1497dB。
- x↔y旋转复频响相对差：约4.65e−16/6.46e−16。
- 独立实际时钟DFT最大相对差：3.86e−13。
- native源/接收输出为float64；版本、网格、dt和所有实际站位核验通过。

通过的是预声明5%复L2/5°有限网格工程诊断门限，不是严格收敛、互易性或现场物理认证。空气5cm细化/反向交换源接收器尚缺native证据。

图：`artifacts/research_checks/2026-10-05_local_uav_free_space_coarse_results_r2/free_space_comparison.png`。右下是**单站四接收方向道集，非空间完整B-scan**，图中如实标注。

## 第二轮：粗网格的相位误差可由离散传播与半格界面解释

这里新增的是**正入射平面波辅助标杆**，不能替代8m悬空偶极的分层格林函数验证。使用本地 `#plane_wave_axial`、内部TFSF、全横向均匀介质、非平均阶梯界面。地表z=12m，覆盖层3m、ε∞=9，基岩εr=4；分别检查空气、εr9无损半空间、无损平层、σ=0.001S/m平层、再增加Δε=1/τ=6ns的Debye平层。参数是机制假设，不是营山电性标定。

上方z=8m接收场减去匹配空气入射并按入射复谱归一；下方z=16m总场除以同点空气场。采用连续介质Fresnel和含多次反射的有限厚度平层公式，ε*=ε∞+Δε/(1+jωτ)−jσ/(ωε0)。没有将平面波归一称作天线端口S参数。

| 无损半空间诊断 | 5cm | 1.25cm |
|---|---:|---:|
| 反射复L2误差 | 10.39% | 2.69% |
| 透射复L2误差 | 29.76% | 3.86% |
| 最大反射相位误差 | 9.13° | 2.48° |
| 最大透射相位误差 | 38.29° | 4.25° |

四类介质×反射/透射共8项：5cm没有通过原先5%/5°门限；1.25cm全部通过，复L2范围2.25%–3.86%，最大相位4.25°。**没有放宽门限来通过**。

为区分网格细化与缩域，另跑5cm/2m横向域，与原5cm/12m域比较：8个复响应最大相对差仅1.999e−15。因此本辅助TFSF标杆中缩域没有改变该接收响应；随后固定2m域做5cm→1.25cm对照才用于网格判断。不能把该结果当成悬空点源的侧向PML认证。

进一步用**无拟合预测**核对相位：沿传播轴，Yee色散关系为

`k_num = 2/dl * asin(sqrt(er)/S * sin(omega*dt/2)), S=c*dt/dl=1/sqrt(2)`。

非平均的介电节点跳变等效界面比声明箱体边缘提前半个z格。4m εr9路径相对空气的离散传播相位误差，在170MHz为5cm时−27.91°、1.25cm时−1.69°；透射再计入半格位置，分别约−38.12°和−4.24°。预测与native半空间透射相位曲线的最大残差分别0.1744°、0.00264°；反射预测残差分别0.00678°、0.000106°。这是当前平层辅助配置的具体数值机制，**不是历史UAV B-scan唯一根因证明**。

参考公式另通过无损能流守恒、半空间已知反射负极性/传播时延、导电及Debye被动性检查。独立实际时钟DFT最大差在粗网格约1.08e−11、细网格1.71e−12；均远小于这里的连续解析偏差。

图与结果：

- 主对照：`artifacts/research_checks/2026-10-05_local_uav_grid_comparison_r1/grid_comparison.png`。
- 5cm宽域：`artifacts/research_checks/2026-10-05_local_uav_plane_results_r5/plane_comparison.png`。
- 1.25cm细化：`artifacts/research_checks/2026-10-05_local_uav_plane_results_fine_r6/plane_comparison.png`。
- 5cm窄域：`artifacts/research_checks/2026-10-05_local_uav_plane_results_domain_r7/plane_comparison.png`。

每个图有中文材料、处理和物理轴标注。多场景同一站位A-scan不被伪装成空间B-scan；这轮没有移动天线测线，也没有新波场快照。

## 接收/频谱口径：可本地快速检验，但不是实测补偿

复用已核验三维空气频响、相同Hann501频带，分别权重1与95MHz/f，二者在95MHz固定相同电流元锚点。后者表示理想空气偶极的等**自由空间辐射功率**假设，不等于设备恒端口接受功率。所有曲线使用同一个等电流参考尺度，未逐道归一。

该假设下轴向/侧向带内能量幅比从−1.189dB变为+0.125dB。证明总体耦合比较依赖频率权重；没有证明左右布局一定改善，也没有将1/f直接写入生产处理链。空气近场随频率的轴向/侧向幅比可反转，与此前独立偶极解析判断一致。

图：`artifacts/research_checks/2026-10-05_local_uav_observation_r1/observation_sensitivity.png`。该图不含地层、真实天线、电路端口或现场CSV。

## 失败证据与源码保存

所有历史输入/输出保留，重新设计使用新目录/契约，没有在原attempt内重试：

1. `free_space_r1`与`free_space_fine_r3`：资源预检拒绝，无求解。
2. `plane_r1`：OpenBLAS在Job提交内存上限内初始化失败；后续声明并限制BLAS/OMP线程。
3. `plane_r2`：本地V4明确拒绝二维对称边界，无H5。
4. `plane_r3`：5道native完成且float64有限，但错误极化/对称面与贴域面的TFSF配置使场增长至1e70量级；**物理拒绝**，不能用native身份PASS冒充有效结果。未改求解器，也不把它归为旧UAV偶极模型缺陷。失败图：`artifacts/research_checks/2026-10-05_local_uav_benchmark_archive_r1/invalid_plane_configuration.png`。
5. `plane_r4`：TMy中的ψ90投影到结构性不更新的Ex/Hy，本地V4拒绝；源码投影公式确认ψ0对应本检查所需Ey/Hx。
6. `plane_r5/r6/r7`：内部TFSF、正确极化、native稳定；粗/细网格物理诊断状态分别独立保存。

旧源码在胶囊内以 `frozen_runner.py`/`frozen_launcher.cmd` 保留，与当时契约哈希匹配；原始粗网格分析也保留 `frozen_analysis.py`。总索引脚本逐条检查匹配，不用当前文件哈希覆盖历史记录。

## 复现与下一项

新CPU平层执行入口（输出目录必须不存在；已消费attempt不可重用）：

```powershell
& 'D:\gprmax_v4_gpu_env\Scripts\python.exe' scripts/uav_local_plane_benchmark.py freeze --dl 0.0125 --width 2 --out artifacts/research_checks/NEW_PLANE_CAPSULE
& 'D:\gprmax_v4_gpu_env\Scripts\python.exe' scripts/uav_local_plane_benchmark.py run --out artifacts/research_checks/NEW_PLANE_CAPSULE
& 'D:\gprmax_v4_gpu_env\Scripts\python.exe' scripts/analyze_uav_local_plane.py --study artifacts/research_checks/NEW_PLANE_CAPSULE --out artifacts/research_checks/NEW_PLANE_ANALYSIS
```

冻结前阅读资源上限与预声明诊断范围；V4 Python源码和编译扩展身份变化时另冻目标机器契约。CPU监督器使用Windows Job Object限制整个子进程组；磁盘上限是轮询日志目录，H5尺寸还受固定接收器/时窗约束，不宣称硬磁盘配额。

下一项仍是**有限高度点源的半空间/分层格林函数独立参考与单道受控验证**，随后才做12道布局/极化/场景及深度对照。正常入射Fresnel不能签认此步骤。若换大内存机器，先重新冻结本机V4/精度/材料/网格/物理PML与单道资源契约，实际测峰值，再决定整批；示例三维尺寸只是资源下界，不是已准备、已验证或已执行的模型。

本单元证明了基础空气链自洽，及平层粗网格相位误差的具体来源；**没有彻底解释营山与原仿真的B-scan差距，也没有完成完整三轮物理验收**。场地/设备未知仍限制定量解释，G4、训练资格与保留实测边界不变。

## 本轮读取的实现依据

读取本地V4文档 `docs/source/input_hash_cmds.rst` 的plane_wave_axial与symmetry_boundary节，以及plane-wave outputs说明；实际安装的 `sources.py` 中HertzianDipole电流元长度、DPW投影/二维活跃分量验证，`user_objects/cmds_multiuse.py` 的TFSF离散位置/本机二维对称边界限制，官方SFCW processing实际时钟DFT/direct实现。文档与本机实现能力不同之处以native失败记录和冻结源码为准。

连续空气偶极公式沿用前一轮MIT6.013相关章节依据；平层公式经过独立能流、被动性与已知时延反例检查。以上分别是解析参考、实现身份、native仿真和有限误差诊断证据，均未被写成实测性能。
