# 九号线首次二维仿真任务书：5试点 → 99站粗采样SFCW

用户授权：“提交的推送把任务要求等等（包括sfcw要求等等）详细列一下以免跑偏，就先把2D这个模型跑一下吧？”；随后“继续，不过391道预计4090要多少时间？首次能不能适度降采样让它别那么慢跑完？”。

本任务先交付一个完整区段的粗采样B-scan。**只降低沿线采样密度：0.5m/391站改为约2m/99站；原地质、网格、材料、航高、求解域、时窗、精度和501个SFCW频点不改。先5个全域试点，原生与SFCW检查及结果审查通过后再跑剩余97站。** 三维、完整391站、处理算法优化/训练均不在首轮执行任务内。

机器无关任务：[line9_2d_first_run_task_v0_1.json](../../configs/research/line9_2d_first_run_task_v0_1.json)；站位和预算由脚本生成，不能手工更改事件表/契约数值。这个JSON是批准任务，目标机执行契约由 `scripts/run_line9_2d_first.py freeze` 独立生成；不能复用其他批的已消费attempt。

## 1. 时间预算与降采样

证据：[planning_evidence.json](../../artifacts/research_checks/2026-10-06_line9_2d_first_run_planning_r1/planning_evidence.json)。参考为仓库 `2026-10-05_slope_bscan_2d_fine2m_r1` 原始执行日志：**4090 Laptop GPU**，726万单元、6785步，26道中位约23.94s。当前4800万单元、CFL步长约58.97ps、1200ns约20352步；单元×迭代量约19.85倍，线性外推约7.92分钟/道。

| 方案 | 站数/独立求解数 | 按旧4090 Laptop批线性外推 |
|---|---:|---:|
| 原完整0.5m扫描 | 391 | 约51.6小时 |
| 首次约2m扫描 | 99 | 约13.1小时 |
| 首次扫描加5试点并复用2站 | 102 | 约13.5小时 |

这是工作量外推，**不是当前模型实测性能，也不是桌面4090/4090D的准确计时**；不同材料数、PML、几何导入、编译和运行时频率会改变耗时。桌面卡可能更快，但目前不乘一个臆造加速系数。完整扫描先按天级（约1–3天）、粗扫描按半天至一天量级规划，不能当保证上限；目标机完成前三道后以中位/最大真实耗时重新给出预算。冻结的60分钟/道、5小时试点批、48小时粗扫批仅为停止线，不能当ETA。

采集X220→25m，`s=220−X`。从既有391站中选索引0/4/8/…/388/390：X220、218、…、26、25，**最后间隔1m**，含两端共99站。99站正演量约为391的25.3%，在同机同设置下约快3.95倍；单道显存/内存需求不减少。2m首版只看整体回波形态，不能认证空间Nyquist、细层横向变化或成像采样充分；图按真实站距绘制，不补造未计算道、不把插值当实测/正演道。

全域试点X220、196.75、166.5、122.5、25。X220/25与99站重合，复用两份经核验原生输出；其余三份仅为钻孔/尖灭/中部诊断。最多5+97=102次独立求解，不重复计算重合站；原391站输入保留可后续补密，未经新任务不得自动启动。

## 2. 必须取得的私有模型包

Git包含任务、代码、聚合证据；**不含场地H5/接触线/站位高程**。目标电脑还需用户在本地复制：

- `D:\自动处理\artifacts\local_checks\2026-10-06_line9_material_model_r1_portable.zip`
- ZIP SHA-256：`c0ca87b431b48538cf8baa728483a8d6cc69de1a5165c8f5c70cf950aadd925b`
- 解压后manifest SHA-256：`ef1b49955ff74df341c744b05bbda698af72fe4a2a3f02ef59a6445fa54e9cde`

ZIP是上一个模型制作单元的不可变包，里面的旧任务说明未包含本轮99站安排；**以本Git任务及目标机新执行契约为准**。ZIP不覆盖/不重新上传；相对输入路径保持 `cases/<id>/profile.in`→`../../geometries/`。原始PDF/旧H5/材料JSON禁止改写或强制加入Git。

## 3. 不变量（不要为了速度改变）

| 项目 | 固定要求 |
|---|---|
| gprMax | 用户指定V4.0.0，核验真实安装/构建；不是旧V3环境 |
| 几何 | 原X−50–350m，高程405–480m，域400×75×0.025m；XYZ=横距/高程/横线 |
| 网格 | XYZ均0.025m，16000×3000×1，严格二维XY剖面 |
| 本构 | 空气自由空间；粉质黏土ε∞11/σ.001/Δ.5/τ6.4567ns；泥岩ε∞11.739510249300162/σ.003/Δ6.200368247790598/τ8ns；砂岩ε∞8.837547784813761/σ.001/Δ2.575432653759346/τ6.4567ns |
| 航高 | 地形跟随收发中点AGL15m，沿用冻结的网格量化位置 |
| 收发 | 沿线1.3m二维代理；实际横线1.3m只能由三维检查，本轮不假称实际实体天线 |
| 源/接收 | z向理想Hertzian，Ricker100MHz/40A；0.025m电流元、名义1Am；接收Ez |
| 边界 | HORIPML，80/80/0/80/80/0；活动轴2m，地质延续穿过PML，平均n |
| 时窗 | 1200ns；不能直接截短以省时，完整性待试点检验 |
| 精度 | CUDA `-gpu_precision double`；原生接收Ez和源samples必须float64，转型不算FP64 |

材料是研究初值，非营山实测拟合；查看模型的横线不变外延及源剖面以外的域填充假设沿用制作报告。不能改变损耗去追求更漂亮的图，也不能把物理不可见改称求解失败。

## 4. 目标机器检查与停止条件

- 当前聊天机器实际RTX3060 Laptop6GiB/16GB RAM；本次执行冻结前容量拒绝，无求解/新B-scan。新站位稀疏化不会使本机显存容得下单道。
- 二维每道检查线：**可用RAM≥13.94GiB、空闲GPU显存≥11.66GiB**。为主要数组下界×1.3+预留，不是实测峰值；系统总RAM和显卡标称容量不能代替空闲值。
- Windows既有监督器：单个可见GPU（GPU0）、同机独占锁、一次1个求解进程、no_retry；记录GPU UUID，freeze/run间不得换卡。Linux/多GPU须另配相应监督器并冻契约，不照抄Windows执行器。
- 每道实际运行前重查RAM/VRAM；本次自有求解进程树RSS上限24GiB、系统剩余RAM下限1.5GiB、每道60min、试点全批5h、后续97道48h。越界只停止本任务自有求解树，保留stdout/stderr/执行事件/失败输出；不杀用户其他进程、不清内存、不转CPU、不降精度、不缩域。
- 非零退出、缺H5、dtype/形状/dt/站位/源参数/哈希不符、频点无源支撑、DFT/逆变换不一致均停止后续。旧attempt不覆盖、不原地重跑；后续恢复仅能在新契约中显式复用已验证结果、列出未完成项。
- 源/尾端/弱事件/PML/网格疑问单列。没有预设“地下必须清晰”的成功门槛；数值链通过不等于物理签认，合理保留风险后可以授权继续粗扫描诊断。

## 5. SFCW要求（不能只交Ricker图）

1. 原始`.h5`、源`srcs/src1/excitation/samples`、接收`rxs/rx1/Ez`只读留存并记录SHA、dt、样本数、版本、原始dtype、网格站位。不得从截图或实数CSV重建虚构复谱。
2. 频率精确 `20MHz + k×0.3MHz, k=0…500`；两端20/170MHz、501点。降低沿线站数**不降低频点数/频带**，也不做501次独立CW正演。
3. 使用已核验V4官方 `toolboxes.SFCW.processing.direct_frequency_response`，processing源码SHA `adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b`。保留复数 `Y(f)/X(f)/source.spatial_scale`，实际SpatialScale必须0.025m；按原生source/receiver各自TimeSampleOffset作精确频率DTFT，不取最近FFT格点、不漏传播相位。
4. 检查全部501个源频点支撑（库默认source_floor_db=−100，是除法条件数检查，不是FDTD误差底），报告带内最小/最大源谱比。无支撑不得填零/盲目外推。
5. 源归一化已处理源相位/已知子波延迟；额外time_shift=0。禁止再次整体平移去掉源延迟，也不擅自清零15m空气双程约100ns。时间轴从原生时基/逆变换推导；不是地下深度轴。
6. 矩形与Hann两版均保留，Hann不是已确认现场窗；官方窗均值归一化保持一致。zero_pad_factor=8，输出完整复轮廓及带符号 `2*real(complex_bandpass)`；取模/包络另域，不作为带符号数据替代。501×8点仅细化显示采样，不能称带宽/分辨率提升。
7. Δf对应3.333μs周期，零填充后采样约0.832ns；图显示0–1200ns，保存完整逆变换时轴/复数组。均匀频率支撑与原生1200ns记录是两个不同约束。
8. 主结果不做尾渐消、增益、AGC、逐道归一化、背景扣除/SVD、图像美化或成像。另做最后20ns尾渐消敏感性，报告最后20ns峰/RMS、末点/全记录峰及复谱变化；不把压低尾部当完整性证据。没有匹配H0，不新造“clean”或理想差分。
9. 各道在7个指定频率独立DTFT复算，矩形/Hann逆变换另按复指数求和复算；相对L2≤1e−9只是实现一致性检查，**不认证弱回波FDTD误差/物理正确性**。
10. 单位是场响应/电流矩，绝不写端口S21；二维线源结果不与三维点源绝对幅值直接等同。

## 6. 交付与试点审查

运行归档包括machine-specific execution_contract、原输入/几何/材料哈希、原生H5、stdout/stderr、execution.jsonl、completed_verification；耗时必须区分首道编译和后续中位/最大值，并报告本任务进程树RSS。容量检查与冻结输入位于私有忽略目录，不能强制加入Git。

SFCW交付 `sfcw_arrays.npz`、analysis_report.json、**中文标注三行灰度图**：原始Ricker总场 / SFCW矩形总场 / SFCW Hann总场。总场、窗/处理、15m、二维沿线代理、材料假设、轴单位逐一标注。两SFCW行共享色标，原始行单位/量纲不同使用独立色标；不逐道归一化。5试点之间画白色未计算区；99站按真实s绘制，最后1m不能误当2m。

图中约100ns线仅是平地镜面地表参考，不是假定地下界面真值；若叠地下到时线须由脚本计算并明确直射/色散/双基地近似。当前未设计波场快照，不生成虚构GIF；若发现到时/反弹身份不清，另冻单道有界快照任务。

试点SFCW分析会生成 `pilot_review.template.json`，状态REQUIRES_REVIEW。接续agent需实际查看图与尾端、检查运行峰值和耗时，填写四项理由：原生/处理检查、尾端与弱事件风险、边界/网格限制、资源与新ETA。写入试点胶囊 `pilot_review.json`、status=APPROVE_PREVIEW，并绑定completed_verification及analysis_report哈希；**这是已授权任务内的agent验收，不要求再次询问用户批准**。尚有必须修复的输入/时间截断等问题时停止预览并报告，不能仅改状态放行。

每个求解批结束先交付该批中文图和限制，再接续；图及聚合报告按用户要求入库推送，原始PDF/地质H5/接触线/现场CSV留本机。提交/推送核验远端HEAD；任务结束更新START_HERE与产物索引，不把这一轮变成材料拟合、训练或G4签认。

## 7. Windows4090接续命令

先拉取分支 `codex/measured-material-spectra-20261006`，复制上述ZIP并核验哈希、解压至忽略目录。示例变量需按实际目标机填写；不安装/迁移本机环境，也不重新配置调度。

```powershell
$env:GPRMAX_PYTHON = 'E:\automation_djh\artifacts\local_checks\gprmax_v4_gpu_env2\Scripts\python.exe'
$env:GPRMAX_VCVARS = 'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
$env:GPRMAX_CUDA_BIN = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.3\bin'
# 独立worktree接续同一GPU时，冻结到既有监督器的同一个锁；普通单目录可不设。
$env:GPRMAX_GPU_LOCK = 'E:\automation_djh\artifacts\local_checks\hs4_gpu_exclusive.lock'
$line9Package = 'E:\automation_djh\artifacts\local_checks\2026-10-06_line9_material_model_r1'
$line9Pilots = 'E:\automation_djh\artifacts\local_checks\line9_2d_pilots_target_r1'
$line9Preview = 'E:\automation_djh\artifacts\local_checks\line9_2d_preview_target_r1'

# 以下包装器仅配置已有MSVC/CUDA，实际版本与哈希在freeze中绑定。
scripts/run_line9_2d_v4.cmd scripts/run_line9_2d_first.py capacity
scripts/run_line9_2d_v4.cmd scripts/run_line9_2d_first.py freeze --package "$line9Package" --out "$line9Pilots" --stage pilots
scripts/run_line9_2d_v4.cmd scripts/run_line9_2d_first.py run --out "$line9Pilots"
scripts/run_line9_2d_v4.cmd scripts/analyze_line9_2d_sfcw.py --study "$line9Pilots" --out artifacts/research_checks/line9_2d_pilots_target_r1_sfcw

# 接续agent完成真实试点审查并写入 $line9Pilots/pilot_review.json 后：
scripts/run_line9_2d_v4.cmd scripts/run_line9_2d_first.py freeze --package "$line9Package" --out "$line9Preview" --stage preview --pilots "$line9Pilots"
scripts/run_line9_2d_v4.cmd scripts/run_line9_2d_first.py run --out "$line9Preview"
scripts/run_line9_2d_v4.cmd scripts/analyze_line9_2d_sfcw.py --study "$line9Preview" --out artifacts/research_checks/line9_2d_preview_target_r1_sfcw
```

脚本拒绝覆盖已有胶囊/分析结果，监督器拒绝已消费attempt。当前对新执行器/处理器完成的是数组/HDF5合约检查（包含FP32拒绝、错误时钟/极化拒绝、半步相位及电流元归一化、独立逆变换），不是GPU求解复现。实际第一次原生求解/性能/尾端验收由目标机完成。

## 8. ROG远端启动记录

2026-10-06用户完成SSH服务、密钥授权与防火墙配置后，客户端核对服务器指纹并实际登录成功。远端为4090 Laptop 16GiB、63.21GiB RAM；启动前空闲显存约13.84GiB、可用RAM约49GiB，容量线通过。模型ZIP传输后SHA-256复验一致。既有`E:\automation_djh`保留原分支和未跟踪工作，在`E:\line9_rog_first_20261006`建立独立worktree；运行时仍使用既有V4环境。

旧`E:\msvc2022bt`和`C:\cuda133`为junction，SSH下拒绝遍历，改用上方真实安装路径。首次freeze因`cl`中文版本输出按GBK解码失败，在胶囊创建和求解前退出；版本日志读取改为显式UTF-8并容许替换无法解码的展示字符，编译器二进制哈希仍逐字节记录，科学数组与源码身份不改变。独立worktree新增`GPRMAX_GPU_LOCK`以冻结到同机既有共享锁。

目标机修复后复跑[数组/HDF5检查](../../artifacts/research_checks/2026-10-06_line9_2d_first_run_planning_r1/rog_contract_checks_r1.json)通过，半步相位相对L2约9.82e-13；5站胶囊freeze通过，`line9_2d_pilots_rog_r1`已启动。这里只确认启动，不声称5站完成或物理验收通过；最终结果应另附当批中文可视化报告。
