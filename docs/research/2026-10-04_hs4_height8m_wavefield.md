# 2026-10-04 8 m 航高密集波场（~1 ns 帧）与界面回波成因分析

> 用户当日需求：做一个 2D 模型，无人机（天线）距地面 8 m，波场快照约 1 ns 一帧，分析"3D 与对应 2D 仿真 B-scan 深部看不到基覆界面形状"的成因。此前多轮尝试已确认地下介质传输正常，"随后就不清楚为什么了"。本单元在 8 m 中间航高上用密集快照把回波生成—上行—到达天线的全过程定量化，并与既有 15 m/2 m 证据衔接。

## 环境重建（本单元前置工作）

2026-10-03  reconcile 清空了旧 venv（`D:\gprmax_v4_gpu_env` 与 `artifacts/local_checks/gprmax_v4_gpu_env` 均不存在）。本次在本机（RTX 4090 Laptop 16GB，驱动 610.88，CUDA Toolkit v13.3，MSVC Build Tools 2022 14.44.35207）重建：

- Python 3.12.15（uv 管理的本机既有解释器，未新下载）；venv 于 `artifacts/local_checks/gprmax_v4_gpu_env/`（gitignored），锁定依赖按 `configs/research/gprmax_v4_gpu_win_py312_dependencies.txt`，pycuda 2026.1。
- pycuda `_driver` 加载需 CUDA DLL 近旁：已复制 `cudart64_13.dll`、`curand64_10.dll`（来自 `CUDA\v13.3\bin\x64`，注意 13.3 的 DLL 在 bin\x64 子目录）到 venv 的 pycuda 包目录。
- junction：`C:\cuda133`→CUDA v13.3，`E:\msvc2022bt`→BuildTools（用 .bat 内 mklink /J 创建；直接在 `cmd //c` 行内带引号传参会被 MSYS 转义破坏）。
- 全部环境变量必须写在 .bat 内部（Git Bash 向 cmd 行内传 `VAR=...` 会被路径转换破坏）；`call vcvars64.bat >nul` 的重定向会阻断环境传播，不能加。
- **NVML 怪癖（两次失败attempt的根因）**：在剥离环境（env -i）下，即使设了 SYSTEMROOT/PATH/ProgramData，nvidia-smi 仍报 "Failed to initialize NVML: Unknown Error"（exit 255）；继承完整用户环境启动同一个 .cmd 即正常。运行器必须以完整继承环境启动。
- wheel 本地构建：`artifacts/local_checks/wheelhouse/gprmax-4.0.0-cp312-cp312-win_amd64.whl`（sha256 f16053e1…），安装后 `scripts/inspect_cuda_runtime.py`（双精度核函数 4096 元 max_abs_error=0.0）与 `scripts/inspect_v4_runtime.py` 均 PASS，结果存 `artifacts/local_checks/2026-10-04_cuda_double_check/`。**与旧冻结契约对账：344 项 source_identities 中 319 个 .py 全部逐字节一致；25 个 .pyd 为未变源码的新本地编译产物，哈希已在新契约中重新记录并注明。**

## 仿真胶囊

- 脚本 `scripts/hs4_height8m_wavefield_v0_1.py`（设计/验收/执行器），`v0_2.py`、`v0_3.py` 为两次**仅启动器**失败的延续（第一次 cmd 未播基础环境变量，第二次仍用剥离环境启动；两次均未进入求解器、无任何 profile.h5，失败胶囊 `2026-10-04_hs4_height8m_wavefield`、`_b` 原样保留）。
- 执行胶囊 `artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c/`：与 HS4T2D 联合网格中心道完全相同的 0.8 m 起伏模型、材料（rock εr9/σ1mS；cover εr18.017/σ3mS+Debye）、2.5 cm 网格、PML（侧 2 m、上下 1 m）、收发偏移 1.3 m、Ricker 95 MHz；仅天线改为 z=20 m（地面 z=12 以上 8 m），快照步 17 迭代=**1.0024 ns**、0–600 ns 共 599 帧（每快照道 0.215 GiB），ROI X13.5–22.5、Z8–28.1 m、输出 0.15 m。
- 三道：mid8_rough、mid8_halfspace（各 599 帧快照）+ mid8_rough_passive（无快照对照）。
- 验收 PASS（`completed_verification.json`）：快照/原生探针空间闭合逐点一致（max_abs_difference 全 0.0，含 E/H 半步插值）；被动对照接收器与带快照道**逐比特一致**；材料几何与参考逐元素一致；dt/迭代数/网格属性一致。每道约 14–80 s。

## 分析（`scripts/analyze_hs4_height8m_wavefield.py` → `2026-10-04_hs4_height8m_wavefield_analysis/`）

差分定义：rough − halfspace，其中 halfspace 为全覆盖层半空间（`centre_halfspace`：z<12 全为 cover），因此差分是**纯基覆界面对比响应**——直达波与地表反射在两模型中完全相同而消去。

### 回波传播链（8 m，1 ns 帧直接可见）

| 时刻 | 事件 |
|---|---|
| 0–53 ns | 直达波与地面反射；差分场为数值零（对数图黑，<1e-12 全局峰） |
| ~60 ns | 差分场首次在界面深度出现（界面带能量过底噪 60.1 ns；垂直预测下行到时 8/0.3+3/0.07≈69 ns，波前有限宽度使首部更早） |
| 70–80 ns | 界面差分场先从天线正下方（x≈17–18.5 m）亮起，随后沿起伏横向铺开 |
| 87.2 ns | 界面带差分能量峰（生成峰） |
| 98.2 ns | 覆盖层带峰（上行中） |
| 131.3 ns | 低空带峰（穿出地表） |
| 140.3 ns | 天线面带峰（回波到达）；接收器差分包络峰 138.7 ns，与垂直双程预测 138.3 ns 吻合 |
| 120–200 ns | 覆盖层内同心弧状多次混响持续并逐渐衰减；岩石内可见下传透射差分弧 |

### 跨航高接收器对比（原生 H5；15m/2m 原始快照帧仅前机留档，本机已无）

| 航高 | 回波包络峰时 | 峰幅 | 相对全记录峰（直达波 91.0） |
|---|---|---|---|
| 2 m | 101.5 ns | 0.1448 | −56.0 dB |
| 8 m | 138.7 ns | 0.0481 | −65.5 dB |
| 15 m | 185.7 ns | 0.0271 | −70.5 dB |

峰时位移精确符合双程空气路径：15→8 m 差 47.1 ns（理论 46.7），8→2 m 差 37.2 ns（每米 6.2 ns，理论 6.67，斜距几何）。幅度随航高按约 1/r 空气路径扩散衰减（8/2 比 0.33 vs 路径比 0.25；15/8 比 0.56 vs 0.53）。

### 对"深部看不到基覆界面形状"的归因

1. **传播链完整，不是传输问题**：界面回波按时生成、按时返回天线（本单元 1 ns 帧逐步可见，且与双程预测闭合）。
2. **是动态范围问题**：回波相对直达波 −56 ~ −70.5 dB，任何共享色标/增益的常规显示下必然淹没；与上一单元界面回波预算（窗内带符号峰为全记录 0.015%–0.047%）一致。
3. **单站信息不构成形状**：8 m 单个共偏移站对界面的照射足迹（≥20% 峰）宽 8.85 m，铺满整个 ROI——界面在回波中表现为近镜面反射的宽足迹叠加，起伏形状只能靠多站三角定位（成像/偏移）恢复，不能从单道到时读取。这与既有"最强波瓣锁定共同浅拱顶、包络峰时与局部覆盖厚度相关仅 0.613"结论一致。
4. **覆盖层多次混响叠加**：120–200 ns 的同心弧图案表明界面响应与地表—界面间多次波在同一时间窗叠加，单道时间剖面上界面回波不是一条干净同相轴。

### 限制

- 单中心站位、TM 不变 y 二维、Ricker 95 MHz；不代表三维有限天线、多站位 B-scan 或实测性能。
- 15m/2m 原始快照帧不在本机（前机留档，哈希在案），跨航高比较仅用原生接收器 H5；如需逐帧跨航高对比须另冻有界契约重取。
- 快照 ROI X13.5–22.5 不含远侧坡面贡献；footprint 边缘受 ROI 截断。
- 两次失败 attempt 为启动器环境问题，未消耗求解器 attempt；契约链与修复说明在各自胶囊内。

## 产物索引

- 执行：`artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c/`（契约、输入、原生 H5、599×2 快照、几何、验收）
- 失败保留：`2026-10-04_hs4_height8m_wavefield`、`_b`
- 分析：`artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_analysis/`（summary.json + 4 图：band_energy、footprint、diff_atlas 固定色标、diff_atlas_log 对数包络、receiver_traces_by_height）
- 环境核验：`artifacts/local_checks/2026-10-04_cuda_double_check/`
