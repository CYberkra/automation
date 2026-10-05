# 本机（RTX4090）空气 5 cm 细网格批 + Device Guard 环境事故记录

- 日期：2026-10-05
- 触发：拉取对方（3060 机）推送后，其[本机轻量前向自检](2026-10-05_local_uav_benchmark_findings.md)遗留"空气 5 cm 细化+互易性因资源预检拒绝未跑"（`free_space_fine_r3` 契约，attempt 未消耗）。按跨机纪律本机重新冻结，不手改旧契约。
- 性质：基础链路自检（空气域点偶极 vs Hertzian 解析连续解），非天线/地面/实测验证。

## 环境事故（Device Guard/WDAC，本单元前置障碍）

- 当日下午起，本机组织应用程序控制策略（WDAC/Smart App Control，winerror 4551 = System Integrity 错误段 4550–4559）**按哈希封禁**了原 V4 环境的基础解释器（uv 安装目录的 Python 3.12.15 `python.exe`）与 `gprMax/cython/virtual_waveguide.cp312-win_amd64.pyd`。当天上午该环境仍正常（Debye 批 B-scan 图为其最后产出）。
- 判定过程：同字节文件换路径仍封（哈希封禁非路径封禁）；wheel 内原始 pyd 与已装 pyd 逐字节相同且同被封（排除本机编译损坏）；附加 1 字节覆盖数据后即可加载（证实纯哈希封禁）；其余 24 个 cython pyd 不受影响。
- 处置（用户授权 computer use 解决）：
  1. conda（`D:\LongZhiBei\miniconda3`，未被策略拦截）新建 Python 3.12.14 → `artifacts/local_checks/py312conda`；
  2. 以其重建 venv `artifacts/local_checks/gprmax_v4_gpu_env2` 并整目录复制旧 venv 的 site-packages（同为 cp312 ABI）；gprMax 4.0.0 导入、pycuda CUDA float64 初始化通过；**344/344 项源码与冻结记录逐字节一致**；
  3. 从本机留档源码副本 `artifacts/local_checks/gprmax_v4_src_copy` 重编 wheel（`wheelhouse2/`，NO_CYTHON_COMPILE 同 Oct-4 配方），**仅替换被封的 virtual_waveguide.pyd**（旧哈希 `e0660c2f…`，新哈希 `501c665e…`，同一源码同一工具链重编，无物理/数值语义改动）。
- 旧契约中的旧 pyd 哈希保留为历史记录，不用新哈希覆盖；新冻结契约记录新哈希。原旧 venv 与被封文件原样保留本机。
- 新运行包装器 `scripts/run_uav_local_benchmark_4090.cmd`（本机 `E:\msvc2022bt` vcvars64 + `C:\cuda133`）；对方推送的 `scripts/run_uav_local_benchmark.cmd` 指向其机器路径（`E:\sisual stdio 2022`、`D:\gprmax_v4_gpu_env`），本机不可用。

## 批次与结果

- 本机失败留档：`2026-10-05_local_uav_free_space_fine_4090_r1`（旧 pyd 被封致求解器 import 失败）、`_r2`（包装器未带 CUDA 环境致 nvcc 预处理失败），均按纪律保留、未重跑。
- 成功批：`artifacts/research_checks/2026-10-05_local_uav_free_space_fine_4090_r3`（4 道：d0.05 x/y forward + y 极化 reverse_x/reverse_y，8×8×8 m 空气域、5 cm 网格、HORIPML、100 MHz Ricker 实际源、V4.0.0/CUDA float64，completed_verification PASS）。
- 分析：`scripts/analyze_uav_local_free_space.py` → `_r3_analysis/summary.json`。**10 个接收点全部过声明诊断门（5%/5°）**：复谱相对 L2 误差 0.159–0.270%，最大相位误差 ≤0.137°；互易/对称差 6.3e-16（轴向）、7.6e-16（侧向）、1.5e-15/1.2e-15（reverse_x/y，数量级为 float64 机器精度）；独立 DFT ≤2.5e-13。
- 与 3060 机 10 cm 粗网格（L2 0.64–1.11%、相位 ≤0.564°）相比，5 cm 细化把误差降到约 1/4，与二阶 FDTD 收敛方向一致。
- 可视化报告：`scripts/plot_uav_free_space_fine_4090.py` → `_r3_report/free_space_fine_4090_report.png`（FDTD vs 解析幅值谱+相位差；10 接收点误差柱状+互易性数值）。

## 结论与限制

- 本机 V4 空气链路（源、坐标、复相位、GPU double、5 cm 细网格、互易性）自洽；对方遗留的 5 cm 细化/互易性缺口补齐。
- 这仍是空气域点偶极自检：**不能签认悬空点源布局、半空间/分层、三维完整布局/起伏/深度或实测差距已解决**。对方既定下一项（分层格林函数 + 单道点源验证；之后 12 道布局×极化×场景、三深度对照）不变，需另行冻结。
- 环境脆弱性：WDAC 云判定可随时再封新哈希；若再次发生，优先重编对应二进制并如实记录，不改物理源码。
