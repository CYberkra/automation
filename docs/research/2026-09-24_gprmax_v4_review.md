# gprMax 4.0.0：官方资料、源码核对与本项目接入

日期：2026-09-24。范围：官方 V4 迁移、安装、输入、输出、源/端口、studies、精度、SFCW 和网格/边界相关章节；不是通读全部手册。仅阅读和静态核对，未安装/升级求解器、未导入本地 V4 包、未执行官方示例/测试、未启动仿真。

## 版本与来源身份

- 用户给定目录 `E:\gprMax-v.4.0.0` 内还有一层同名目录；实际源码根为 `E:\gprMax-v.4.0.0\gprMax-v.4.0.0`。
- `_version.py` 声明 4.0.0；下载目录没有 `.git`，没有据此伪造上游提交号，也未证明这个目录逐字节等同某个官方发布压缩包。
- [源码与官网台账](../../artifacts/research_checks/2026-09-24_gprmax_v4_source_audit.json)保存 19 份选定文件的路径、阅读范围、字节数、SHA-256 和选定证据集合指纹，另有 10 个官方页面的 URL、缓存哈希。指纹仅覆盖所列文件，不是整个分发包的校验。
- 官网 `latest` 当前显示 4.0.0，但 URL 可变。缓存放在忽略的 `artifacts/research_sources/2026-09-24_gprmax_v4_official/`；Git 保留来源台账和本报告。缓存成功不等于页面每一节均已阅读。
- 既有 `D:\LongZhiBei\miniconda3\envs\gprmax` 环境的静态包元数据是 **gprMax 3.1.6 / Python 3.10.21**。上轮只借用其中 NumPy/Matplotlib 作图，没有调用其求解器。官方迁移指南的基线为 3.1.7，不能声称涵盖本机 3.1.6 的每一项差异。
- 本地 V4 `setup.py` 要求 **Python >=3.11,<3.14**；研究数组环境与求解器环境应分开。目录存在、版本文件正确，不代表 V4 已安装可运行。

核对方式：[静态来源审计脚本](../../scripts/audit_gprmax_v4_sources.py)，不导入或执行 gprMax。未来换机器重新指定源码根并生成新台账，比较所列哈希。安装准备按[官方安装章节](https://docs.gprmax.com/en/latest/inc_README.html#installation)，不要直接升级旧 V3 环境。

## V3→V4 对项目有实际影响的变化

| 项目 | V4 契约/本地核对 | 本项目处理 |
|---|---|---|
| Python API | `gprMax.run(inputfile=...)`；首个位置参数是 `scenes` | 使用关键字，不照搬旧 `gprMax.gprMax.api` |
| toolbox | `gprMax.toolboxes.*` | 后处理和求解器同版本，不混用旧工具 |
| 2D | 显式 `#domain_mode: TM` 与 z=`inf`；TM/TE 内部厚度不同 | 本项目 z 不变、Ez 的 TM 模型；不把 TE 当兼容选项 |
| 文件 | 普通结果 `.h5`；V3 常用 `.out` | 读 schema/元数据，改扩展名不是迁移 |
| 接收器 | `/rxs/rxN` 是文件内编号；当前有构建顺序/身份元数据 | 命名 measurement，配对/合并核对 Name/StudyID、分量及坐标；不固定相信 rx1 |
| 时间轴 | 接收 Ez 偏移 0；磁场/环流偏移 −dt/2；电流源通常 +dt/2 | 每个量读自身 SampleInterval/TimeSampleOffset，不统一平移半步 |
| 源历史 | `/srcs/srcN/excitation/samples` 保存实际使用的标量采样及量纲/空间尺度 | 以存储源历史做时序与频谱检查，不能只重建同名波形 |
| 精度 | CPU 与加速器默认 single；CPU/CUDA/OpenCL 可选 double | 精度必须显式记录，NumPy 后处理 float64 不提升求解器原始精度 |
| studies | 固定几何/网格的重复采集；材料、网格、负载变更需新建模型 | 网格收敛、M01/M02 开关不能复用同一几何；不混入普通 src/rx_steps；study 不支持 MPI 分解/taskfarm |
| 启动参数 | 普通续跑用 `-i`；`n` 是数量；taskfarm 与空间 MPI 不同 | 后续执行包明确输入/输出 basename、病例数及并行方式 |
| 其他数值变化 | 磁导平均、Drude 系数和 snapshot 时间/格式有变 | 首批无损非磁机制例之外，使用相关功能时另行验证；不假设 V3 位级一致 |

依据：[官方迁移指南](https://docs.gprmax.com/en/latest/migration_v3_v4.html)、[输入命令](https://docs.gprmax.com/en/latest/input_hash_cmds.html)、[输出](https://docs.gprmax.com/en/latest/output.html)、[studies](https://docs.gprmax.com/en/latest/studies.html)、[精度和后端](https://docs.gprmax.com/en/latest/accelerators.html)。本地对应文件及阅读范围见台账；表中 project action 是本项目安排，不是官方性能结论。

## 已发现的旧方案不一致：Ricker 的中心

旧 P1 契约定义 fc=80 MHz、中心 tc=37.5 ns。本地 `gprMax/waveforms.py` 的 `Waveform.calculate_coefficients` 对 Ricker 取 chi=√2/fc，`calculate_value` 以 time−chi 计算；因此无额外起始延时的内置中心约 **17.67767 ns**，相差约 **19.82233 ns**。

这是源码公式直接计算，不是运行结果，也未声称是 V3→V4 新增变化。仅把输入写成 `#waveform: ricker 1 80e6 ...` 不会复现旧契约。后续须在“保留原解析波形”与“修订契约使用内置波形”之间明确选定，并更新窗口。单纯设置 source start 也要核对启停截断和半步采样，不能宣称整个有限记录自动等价。

`HertzianDipole.calculate_waveform_values` 在 `(n+1/2)dt−start` 评价激励；writer 标注 +dt/2。对 z 向源，定义 Js=I·dl/(dx·dy·dz)、dl=dz，代数上成为 I/(dx·dy)。这说明二维线源与三维偶极缩放不能混用；实际细网格的离散位置和归一化仍需以后经批准的小例验证。

`FDTDGrid.calculate_dt` 的 TMz 分支使用 x/y 两个活动方向。通用文档中的三维 CFL 公式不应直接套到此二维实现。实际输出 dt 与样本数仍以文件为准，不能用请求的 time_window 创建 `linspace` 冒充实际时间轴。

## SFCW 20–170 MHz 的接入边界

用户已确认设备类型 SFCW，标称频段 20–170 MHz。频点数、步长、是否均匀、驻留时间、天线/接收链和输出定义尚未确认，不能把旧 CSV 的 501 个时域样点当成频点数。[项目背景记录](../../configs/research/project_context_v1.json)

官方 SFCW 工具以宽带 FDTD 响应进行频域合成：direct 使用实际源谱和源/接收物理时刻；homodyne 是单样本冲激、单 A-scan 的独立处理途径。它们描述理想离散线性系统的处理，不自动包含本设备的天线端口、扫频驻留或硬件校准。[官方 SFCW 说明](https://docs.gprmax.com/en/latest/inc_SFCW.html)

本项目后续应保存频率、复响应、源有效掩码和处理窗，直到输出定义对齐后才转为实数 B-scan。不能将 Ez（V/m）或场/电流比直接称为实测 S21；Hertzian 源没有电路端口。[官方源与端口](https://docs.gprmax.com/en/latest/sources_ports.html)

20–170 MHz 的总带宽是 150 MHz。若以后确认均匀频率步长 Δf，则周期重构时窗为 1/Δf；目前不能冻结该时窗。FDTD 记录长度应覆盖所需传播与衰减，和 1/Δf 是两回事；零填充不会增加物理带宽或补回晚到波。20 MHz 端的实际源谱是否足够，须用源谱及数值有效性检查，不能仅凭 80 MHz 的中心频率断言覆盖全部设备频段。

## “约 20 m 浅层”对模型的影响

用户原话“20m浅的非显性滑坡”暂解释为地下约 0–20 m 的目标范围，此解释尚待精确确认；它不是设备已证明能穿透 20 m 的性能声明。

旧 P1 的地表 y=7 m、底界 y=0 m，地下总范围仅 7 m（还未扣 PML）；层界面为 2 m。它仍可用于源/算子机制检查，不能成为该项目全深度代表模型。

只作时间尺度估算：在均匀、无损、非磁介质，垂直两程 t≈2d√εr/c。d=20 m、εr=4/9/16 时约为 **267/400/534 ns**。这些数值尚未加源延时、空气路程、斜路径和晚到波，不能当正式采集时窗；旧 240 ns 显然不足以覆盖这些假设下的 20 m 两程响应。

模型深度、收发几何、边界、材料损耗/色散与完整时窗必须重新形成一版项目场景方案。用现有理想介电常数作时间尺度算例，不等于给现场材料赋值。首个小模型无需立刻扩到 20 m，但必须标为机制校准，避免混淆用途。

## 仿真开始前需要提交用户敲定的执行包

遵照用户明确要求，当前 [执行记录](../../configs/research/gprmax_v4_execution_gate.json)为 `approved_to_simulate=false`。

下一工作可继续准备：独立 V4 环境方案、静态输入、源时序选择、M00/M01 校核设计及 20 m 场景尺寸/时窗草案。运行前给用户一份具体且可审查的包，列清：病例与次数、精确源码/构建版本、CPU/GPU 与精度、时间/内存/存储预算、输入哈希、参考与验收、失败停止规则。未敲定前不运行任何包含求解的示例、pytest、benchmark、study 或自动重试。此次也未运行 geometry-only。

目前没有可靠运行时估计，不虚报几秒完成，也不把先前 45 次 A-scan 当已批准预算。首次批准的小试结束后记录实测成本，再讨论扩展。
