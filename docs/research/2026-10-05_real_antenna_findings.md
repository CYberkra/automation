# 真实结构天线 ×4 缩放批：P0/P1/P2 与 5W/20W 功率后处理（2026-10-05）

状态：机制诊断级。天线为 **gprMax 内置 GSSI-400 按电磁相似律 ×4 缩放的假定设计**（非已验证型号）；粗网格（4 cm）下天线薄片细节被量化/加厚（逐条记录于各 `input_manifest.json` 的 `coarse_adapt_log`)。噪声底为声明模型，非硬件实测。不做物理验收与实测外推。

## 设计

- 天线封装 `scripts/antenna_gssi400_x4.py`：几何坐标 ×4（外壳 1.2×1.2×0.712 m)、吸波 σ÷4=0.0095996 S/m、馈源/接收阻抗不变（257.97/288.93 Ω)、激励换 impulse 电压源（官方 SFCW 真实源链路）。收发沿 x 间距 0.648 m（与 2D 的 1.3 m 不同，已声明）。V4 subgrid 仅 CPU，故主网格 4 cm 下用 `coarse_adapt`：丢 PCB 薄板与馈针边、hdpe 滑板与 PEC 前壁加厚到一格。
- 3D 场景（P1/P2 共用）：域 12.0×3.2×24.0 m @4 cm(300×80×600,14.4M 单元），坐标平移 x'=x−11.6、z'=z−6；材料与 2D 档案一致（cover εr=18.017+Debye 7.878/6.4567ns、σ=0.003;rock εr=9、σ=0.001)；起伏界面取自 2D centre_rough 的 #box 剖面沿 y 不变拉伸（z 量化到 4 cm);halfspace=均匀 cover（无基岩，同 2D 语义）；地表 z=12 m；天线滑板底 z=27 m(15 m 航高）;PML x 50 格/z 25 格（对齐 2D 的 2 m/1 m),y 10 格（声明选择）;600 ns;GPU double。
- P2:13 站 x=15.25+0.5k(k=0..12)×{rough, halfspace}=26 道，冻结式（目录已存在即拒跑、每道 manifest 先行、audit FAIL 即停），全部 PASS。求解实测 ~91 s/道（RTX 4090 Laptop)。

## 结果

- **P0 自由空间链路**：收发非零，SFCW direct_frequency_response 在 20/105/170 MHz 有效；端口 S11 在 20–170 MHz 为 −3.3~−4.1 dB（失配大，振铃强）；网格频率上限 289.8 MHz;P2 域内因 cover 材料为 175.0 MHz。
- **总场 B-scan**(`antenna_bscan_power_panels.png` 第一行）：直达+地表波主导；中心站道包络（`trace_envelopes.png`）示直达 ~0.1(0–50 ns)、地表反射 ~5e-3(100–120 ns)、另在 ~205 ns 有 ~1e-4 鼓包——时间与"天线↔地表二次弹跳"（天线 S11 −3.5 dB 的再辐射，100+100 ns）一致，不是界面回波。
- **界面回波水平**:gated contrast（声明距离门 120–140 ns，理由见下）界面窗内回波/直达 = **−72.4~−76.6 dB**(13 站），比 2D 线源 15 m 档案（−65.9 dB）再低约 8–10 dB——3D 偶极子 spreading + 失配，比 2D 更不利。
- **差分数值底（重要）**:rough−halfspace 在直达/地表时段有非因果残差（−44 dB 级），机理定位为 HORIPML 板系数被远处材料污染 + 网格速度前驱双通道，详见 [跨场景差分底](2026-10-05_cross_scene_subtraction_floor.md)。同场景重跑逐位为零，排除求解器随机性。对比分析一律过 120–140 ns 距离门；总场不受影响。
- **5W/20W 功率后处理**（声明噪声模型：逐频点复高斯，D_ref=100 dB@36dBm,+1 dB 功率=+1 dB 动态范围，种子 20261005,K=8):37 dBm(5W）余量 **+38.9 dB**,43 dBm(20W)**+45.6 dB**（门限 6 dB)。图板第 4/5 行与无噪行几乎无差别。**结论与前批一致且更强：噪声底不是瓶颈，功率（5W→20W）解决不了相干直达/地表波/天线多次反射压制问题。**
- 门后 contrast 在 165–225 ns 有 ~1.5e-5 响应，但其峰位随站位右移（k0→k12 变晚）与起伏界面（k0→k12 变浅应变早）**相反**，且落在天线多次反射时段附近——不能归因为界面形状；标注为"疑似多次反射差分+数值底污染"。

## 产物

- `artifacts/research_checks/2026-10-05_antenna_gssi400x4_p0/`(h5 结构/端口 S11/粗化日志）
- `artifacts/research_checks/2026-10-05_antenna_gssi400x4_p1/`(rough/halfspace:manifest+audit+h5;57 MB vtkhdf 按用户体量约束留本机）
- `artifacts/research_checks/2026-10-05_antenna_gssi400x4_p1b/`（数值底诊断：重跑 diff 恰为 0)
- `artifacts/research_checks/2026-10-05_antenna_gssi400x4_p2/`(26 道 manifest+audit+h5+batch_index;k06  vtkhdf 留本机）
- `artifacts/research_checks/2026-10-05_antenna_bscan_power_r1/`（图板+包络图+summary.json)
- 脚本：`scripts/antenna_gssi400_x4.py`、`p0_…`、`p1_…`、`p1b_…`、`p2_…`、`analyze_antenna_bscan_power.py` 及对应 cmd 包装器。

## 背景抑制（2026-10-05 追加，用户"做个背景抑制呢")

脚本 `scripts/analyze_antenna_bscan_background.py`，结果 `artifacts/research_checks/2026-10-05_antenna_bscan_background_r1/`（图 `background_suppression.png`/`_gray.png` + summary.json)。对 P2 总场 rough B-scan（声明窗 0–300 ns,[sample, trace],SFCW real bandpass）应用冻结算子目录的背景算子，对照参考为门后 rough−halfspace contrast（带数值底警告），负控为 halfspace 上加 SVD rank-2。

- **直达+地表波在该窗内几乎是秩 2 结构**：相对奇异值 1 / 9.2e-3 / 1.0e-4。mean×1.0 与 SVD rank-1 把 0–140ns 能量压到 0.92–0.94%(−40dB);**SVD rank-2 压到 1.5e-5(−96dB)** 且秩 2/3 截断间隙（9.0e-3）可分辨；rank-3 的截断间隙（5.3e-5）接近数值并列，且界面窗相关降到 0.248——**rank-3 开始吃信号，不要用**。
- 界面窗（±25ns 沿真实起伏）与门后 contrast 参考的逐站相关：median 0.54–0.56(mean1.0 / svd1 / svd2 相近），逐站 −0.09~0.80，相对 L2 0.62–1.22——**部分对齐**，剩余偏差来自参考数值底、压制伪影及多次反射残留，机制级不裁决。
- RPCA(λ=λ0）残留 54% 早窗能量（rank_L=6)，对本数据不如直接 SVD rank-2；不作首选。
- 负控：halfspace 上 SVD rank-2 后界面窗残余 RMS 3.4e-8——**没有界面时不会凭空造出界面**，压制不引入假阳性。
- 声明：去除分量是"相干水平分量（直达+地表波）"，不是已认证噪声；若真实界面完全水平等时，这些算子同样会把它去掉——本场景界面起伏随站位移动所以幸存。

## 下一步建议（待用户决定）

1. 天线失配（S11≈−3.5 dB）主导振铃与多次反射：若要做"真实天线"结论，需要按实际设备的天线模型或至少加匹配段；当前 ×4 缩放假定设计只能给机制级结论。
2. 基覆界面成像的可行路线仍是 2D 已验证的对比+成像（RTM 类）;3D 真实结构天线批显示原始 B-scan 上界面回波比直达低 ~76 dB,**任何功率提升都无法改变这个相干比值**。
3. 若要继续用差分：距离门必须保留；或研究同源双场景共享 PML 系数的定制构型以压低数值底。
