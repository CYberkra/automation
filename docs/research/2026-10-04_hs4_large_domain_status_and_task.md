# 108m大域：已准备，RAM预检阻止启动；完整跨机任务

> **2026-10-04 后续：本机（RTX 4090 Laptop）已按下方"另一台电脑直接接续"流程新冻 target_r1 并完成六道，带限固定窗与 36m 归档一致到 ≤2.2e-8，原始宽带晚时 O(1) 不同。结果见 [108m 大域六道结果](2026-10-04_hs4_large_domain_findings.md)。下方为执行前的历史状态。**

用户要求“能不能跑一个很大的域看看？”。已按[执行前方案](2026-10-04_hs4_large_domain_plan.md)生成中心/左/右三站各起伏/全覆盖层六份V4输入，横向36→108m、570.24万单元，保持2.5cm X/Z网格、33m高、15m航高、原Debye、PML物理厚度与600ns。原ROI及收发刚性+36m，新增边缘地层水平延续；全材料解析与原图edge-pad独立比较，六道PASS。

**本轮新增正演0道，没有大域接收结果、B-scan或边界归因结论。** 通过完整vcvars64/CUDA环境启动监督器后，实时资源预检抛`live RAM/VRAM preflight: no solve authorized by capacity`，发生在创建execution.jsonl与调用求解器之前，因此未消耗求解attempt。没有缩域、粗化、改材料、降精度、停止其他应用或重跑旧仿真。复查可用RAM约1.07GiB、VRAM4.502GiB，冻结门槛RAM1.8/VRAM4.4GiB；复查数值不能倒填为启动时刻测量。

## 证据与准备状态

- 初始冻结/拒绝：[r1](../../artifacts/research_checks/2026-10-04_hs4_large_domain_r1/execution_contract.json)、[拒绝记录](../../artifacts/research_checks/2026-10-04_hs4_large_domain_r1/capacity_rejection.json)。包含当时脚本逐字节副本`frozen_code/`，主脚本后来增加跨机运行时冻结，r1保留历史，不用当前代码执行这个旧合同。
- 当前可执行准备包：[ready_r2](../../artifacts/research_checks/2026-10-04_hs4_large_domain_ready_r2/execution_contract.json)、[六道几何预检](../../artifacts/research_checks/2026-10-04_hs4_large_domain_ready_r2/preflight_verification.json)。新冻结未运行，输入内容与r1相同。319项Python源码须与已审计V4.0.0一致，目标机器本地编译扩展逐个重新哈希；构建不同不宣称跨机数值相同。
- 新[执行脚本](../../scripts/hs4_large_domain_controls.py)和[分析/独立DFT脚本](../../scripts/analyze_hs4_large_domain_controls.py)语法编译通过；冻结及全材料变换检查实际通过。新分析入口因缺正演输入未做端到端运行，不能把复用已测量DFT和比较函数说成新大域验证完成。

## 当前电脑资源恢复后的命令

从项目根目录执行。ready_r2的运行环境/代码已冻结，若环境/代码变化须新freeze到另一不存在目录，不覆盖旧合同。

```powershell
& scripts/run_hs4_large_domain_controls.cmd run --execute --out artifacts/research_checks/2026-10-04_hs4_large_domain_ready_r2
& 'D:\gprmax_v4_gpu_env\Scripts\python.exe' scripts/analyze_hs4_large_domain_controls.py --capsule artifacts/research_checks/2026-10-04_hs4_large_domain_ready_r2 --out artifacts/research_checks/2026-10-04_hs4_large_domain_results_r1
```

启动仍须实时RAM1.8/VRAM4.4GiB，逐道RSS2.35GiB/系统余量0.2GiB/900s，5400s整批。这些是有旧同版本运行估计支撑的预算，**不是已量到大域峰值**。显存很接近预算；若初始化真实不足，留失败attempt并停止，不降精度/缩域/自动重试。无执行事件的资源预检拒绝不消耗attempt，可在资源恢复后按原合同启动。已有execution.jsonl的失败/完成attempt不得重复run。

## 另一台电脑直接接续

先拉main，进入实际项目根目录。复用用户V4.0.0/CUDA原生双精度环境；不得安装latest替换指定版本。至少满足上述实时余量，建议空闲RAM>=4GiB及GPU>=8GiB以留出初始化余量，仍以实际日志与运行保护为准。设置目标机已有路径，随后新冻目标机器胶囊。以下路径为示例，必须换成该机真实路径；无须复制本机venv、绝对路径合同或CUDA缓存。

```powershell
$env:HS4_V4_PYTHON = 'D:\gprmax_v4_gpu_env\Scripts\python.exe'
$env:HS4_VCVARS64 = 'E:\sisual stdio 2022\VC\Auxiliary\Build\vcvars64.bat'
$env:HS4_CUDA_BIN = 'C:\cuda118\bin'
$hs4TargetCapsule = 'artifacts/research_checks/2026-10-04_hs4_large_domain_target_r1'
$hs4TargetResults = 'artifacts/research_checks/2026-10-04_hs4_large_domain_target_results_r1'
& $env:HS4_V4_PYTHON scripts/hs4_large_domain_controls.py freeze --out $hs4TargetCapsule
& scripts/run_hs4_large_domain_controls.cmd run --execute --out $hs4TargetCapsule
& $env:HS4_V4_PYTHON scripts/analyze_hs4_large_domain_controls.py --capsule $hs4TargetCapsule --out $hs4TargetResults
```

所有输出目录必须不存在。冻结沿用Git中本机已审计胶囊的Python源码身份，仅允许目标机扩展重新编译且全部记录哈希。运行逐道实际审计完整cell-map、native4.0.0/float64/有限值、CFL/完整时窗、源样本/时序和实际收发格点；未通过不得处理。成功结果输出summary.json、comparison_arrays.npz、同物理尺度large_domain_ascans.png；分析需12份原始接收（六份新输出及六份已归档参照），七个频点独立直接DFT及标量指标复算。全材料图保留目标机，Git归档审计哈希/可复现输入、原生接收、日志、结果与报告；不得上传保留实测资料。

重点看固定160–220ns的界面对比复数/有符号/包络变化及峰时，同时报告两个子窗。总响应与原始时域另列，不能把宽带晚时变化当地下带限窗同等变化。108m最短侧边空气返回约325ns，但有限带宽/窗函数会把晚时响应混入别的窗，不能单凭传播估计宣布排除边界。

若差异很小，结论限定为扩大X/外延在当前网格/材料/二维三站与该窗下影响小；若大，进一步检查边界及新增外延路径，不能预设侧反射主因。顶部/底部未远移，三站不是完整B-scan，二维线源不等于三维有限天线，网格收敛、材料真实性、实测性能均不由本对照认证。旧15道因素结论不被未执行的大域提案替代，G4/训练/首版范围不变。

给另一台Codex的任务：读取START_HERE本条和本文件，按该机V4.0.0环境新freeze六道108m配对；先完成中心，再左右，不重跑旧attempt。运行前守预算，完成后独立审计/分析，记录失败及限制、更新入口/台账并提交推送；RAM不足时保留任务，不清理用户进程或改变数值参数。
