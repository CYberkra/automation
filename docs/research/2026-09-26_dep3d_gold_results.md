# dep3d_gold_v1 批次结果：8 算例全部完成，3D 纯 dz 对照与同格距跨维比较

按[设计书](2026-09-26_dep3d_gold_design.md)执行：8 次 FDTD 运行（4 个 2D 同格距对照 + 2 组 3D 配对）全部串行完成，每例 1 次 attempt、CUDA double、无重试、无 CPU 回退，`reason=completed`、exit 0。20 GiB Job 提交内存与 40 分钟墙钟上限均未触及。

以下区分来源：墙钟、exit、Job 峰值提交内存、进程计数为监督器实测（`supervision.json`）；求解器"FDTD 求解耗时（Time taken）"与"gprMax 进程总时长（Simulation completed in）"为 stdout 自报；收敛链与跨维比较数值全部来自分析脚本对已归档 HDF5 输出的实测复算（`analysis_r1/results.json`），不是求解器声明。r1 与 r2 两次连跑的 `results.json` 已核验逐字节一致。

## 执行记录

|运行|attempt|exit|墙钟（实测）|求解耗时（自报）|gprMax 进程总时长（自报）|Job 峰值提交内存（实测）|
|---|---:|---:|---:|---:|---:|---:|
|B2D5CM_BG|1|0|910.813 s|2.2555 s|10.8488 s|2,324,672,512 B（约 2.16 GiB）|
|B2D5CM_TGT|1|0|910.562 s|2.2289 s|10.8537 s|2,322,751,488 B（约 2.16 GiB）|
|B2DANISO_BG|1|0|910.843 s|4.4439 s|13.2596 s|2,532,499,456 B（约 2.36 GiB）|
|B2DANISO_TGT|1|0|910.609 s|4.3576 s|13.0113 s|2,530,480,128 B（约 2.36 GiB）|
|DEP3D_5CM_BG|1|0|916.078 s|3 min 21.5987 s|3 min 35.8739 s|8,567,218,176 B（约 7.98 GiB）|
|DEP3D_5CM_TGT|1|0|916.156 s|3 min 21.6987 s|3 min 36.0577 s|8,565,354,496 B（约 7.98 GiB）|
|DEP3D_ANISO_BG|1|0|920.218 s|9 min 23.3906 s|9 min 41.6521 s|14,976,286,720 B（约 13.95 GiB）|
|DEP3D_ANISO_TGT|1|0|919.984 s|9 min 23.4664 s|9 min 41.8386 s|14,979,588,096 B（约 13.95 GiB）|

- 峰值 Job 提交内存范围 2,322,751,488–14,979,588,096 B（约 2.32–14.98 GB 十进制），全部低于 20 GiB（21,474,836,480 B）上限。实测"Job 提交/单元"比值随网格变化，设计书 §7.2 预估的悲观口径（22.91 GiB）未出现。
- 墙钟范围 910.6–920.2 s，远低于 2400 s 上限。

## 墙钟构成（按事实记录）

全部 8 例的墙钟都远超 gprMax 进程自身时长（FDTD 求解 2.2–563.5 s，含构建与写盘的进程总时长 10.8–581.8 s）。以监督器实测墙钟减去 stdout 自报的进程总时长，差值（派生量）为：

- 4 个 2D 对照：897.6–900.0 s；
- 3D 5 cm 配对：约 700.1–700.2 s；
- 3D ANISO 配对：约 338.1–338.6 s。

即求解完成后 Job 内仍有成员进程存活（滞留最长约 900 s）。监督器记录 `active_job_processes_at_stop=0`、`waited_for_descendants=true`，即停止时 Job 已清空、全部后代进程已被等待；滞留期间 Job 内成员进程的身份与数量未被监督器记录。该现象对结果无影响（输出 h5 均完整归档并通过分析前校验）。

**根因补记（2026-09-26，batch2d_v1 执行期间鉴定）**：在相同监督器/启动器/工具链环境下观察 batch2d_v1 首例的滞留窗口，Job 内唯一存活成员鉴定为 `vctip.exe`（Visual Studio 2022 BuildTools MSVC 14.44 的编译器遥测进程），由 pycuda 内核编译链中的 `cl.exe` 在算例开始时派生，编译链退出后成为孤儿并继续空转。其空转时长与"自最后一次编译器活动起约 900 s"的行为一致：本批 2D 例编译极早结束（滞留 897.6–900.0 s）、3D 5cm 例约 700 s、ANISO 例编译持续更久（滞留约 338 s），差异方向全部吻合。本批 8 例的滞留成员身份未被当时记录，以上为同环境同工具链下的鉴定推断，不作为对本批 8 例的直接观测。

## 各配对主结果（TGT−BG 差分，200 ns 尾窗）

|配对|维度|网格|包络峰到时 (ns)|谱峰频率 (MHz)|
|---|---|---|---:|---:|
|B2D5CM|2D TMx|1×640×1000，dy=dz=0.05|539.8960|162.8|
|B2DANISO|2D TMx|1×640×2000，dy=0.05/dz=0.025|536.5057|169.7|
|DEP3D_5CM|3D|160×200×1000，0.05 各向同性|541.2554|168.5|
|DEP3D_ANISO|3D|160×200×2000，dz=0.025|543.4823|170.0|

- 四例包络峰到时均落在设计 §9.3 的诊断窗 440–600 ns 内，窗内峰与全局峰比值均为 1.0；理论双程到时约 520 ns（设计 §9.3 自算）。
- 400 ns 稳健性对照：各配对相对主窗口的相对 L2 差为 0.000129–0.000295，尾窗拖尾不是主导差异。

## 3D 纯 dz 细化对照（设计 §9.2）

DEP3D_ANISO（dz 2.5 cm）相对 DEP3D_5CM（dz 5 cm）的 TGT−BG 差分比较：

|量|值|
|---|---:|
|谱形状相关系数|0.9951|
|谱峰频率差|1.5 MHz（170.0 − 168.5，向带边移动）|
|包络峰到时差|160.76 ns（细网格更晚）|
|时域波形最大归一化互相关|0.18（最优时移 21.86 ns）|

方向解读（仅报差异与方向）：谱形状高度一致（0.9951），谱峰频差仅 1.5 MHz，但时域差分波形相关很低（0.18）、包络峰到时差达 160.76 ns——dz 细化对时域差分波形的影响远大于对谱形状的影响。设计 §9.2 明确不给阈值承诺、无预声明收敛判据，本文只报告差异与方向，**不宣布收敛通过**，也不据此声明 5 cm 网格可用于或不可用于批量仿真。

## 2D/3D 同格距跨维比较（设计 §9.3）

两组格距对齐的配对（各维内部先做 TGT−BG 差分，只比形状，禁止跨维绝对幅值比较）：

|量|B2D5CM vs DEP3D_5CM|B2DANISO vs DEP3D_ANISO|
|---|---:|---:|
|归一化谱形状相关|0.935|0.958|
|包络峰到时差 (ns)|−71.02|−53.72|
|谱峰频率差 (MHz)|−5.70|−0.30|
|时域波形最大归一化互相关|0.738（时移 0.94 ns）|0.732（时移 0.37 ns）|
|相位差线性拟合：等效时移 (ns)/残差 std (°)|−0.137 / 5.51|−0.046 / 5.45|
|选频点归一化幅度趋势差最大值|17.20 dB @ 20 MHz|17.43 dB @ 20 MHz|

- 幅度趋势差在 20 MHz 处最大（17.20 / 17.43 dB），随频率上升收窄（170 MHz 处 −0.19 / −0.00 dB）；低频段 2D 与 3D 差分幅值趋势偏离最大。跨维绝对幅值不比较（`cross_dimension_absolute_amplitude_compared: false`）。
- 相位仅以"相位差对频率线性拟合"的形式报告：等效时移为形式量（−0.137 / −0.046 ns），拟合残差 std 5.51° / 5.45°，不当作物理时延，不跨维迁移。
- **混淆因素（按设计 §9.3 原文记录）**：即便格距对齐，每对比较仍同时混入 (i) 被测量的维度差异（2D 目标沿 x 无限/线源/柱面扩散 vs 3D 目标沿 x 有限 2 m/点源/球面扩散）；(ii) 侧向域宽差异（2D 32 m vs 3D 10 m）；(iii) 目标 y 位置差异（2D y 14–18 vs 3D y 3–7，域内相对位置相同、离侧向 PML 的绝对距离不同）。这些不是干净的维度因子；两对差值之间不得再相减当作可加分解。

### 远场 3D→2D 变换（探索性证据，单列）

`D_2D_equiv(f) ∝ D_3D(f) / √(i·2πf)`（常数省略，只比形状）：

|比较对|变换前谱形相关|变换后谱形相关|
|---|---:|---:|
|B2D5CM vs DEP3D_5CM|0.935|0.967|
|B2DANISO vs DEP3D_ANISO|0.958|0.979|

两个格距下变换均使谱形相关上升（+0.033 / +0.020），方向一致。限制随结果写出：近场区不适用、常数未定、只比形状；仅作探索性证据单列，不作为结论依据（`far_field_transform_is_conclusion_basis: false`）。

## 粗 2D 网格 vs DEP 2D ZFINE2 链（仅形状比较，格距不对齐）

参照：ZFINE2_200（DEP 2D 链 TGT−BG，200 ns 尾窗，dy 12.5 mm / dz 3.125 mm），来源 `2026-09-25_deep_zfine2_analysis_r1/arrays.npz`；谱在同一 501 频点网格上直接可比。

|比较对|相对 L2 差|归一化谱形状相关|
|---|---:|---:|
|B2D5CM vs ZFINE2|1.1379|0.982|
|B2DANISO vs ZFINE2|0.7927|0.999|

同为本批 2D 对照，垂向细化后的 B2DANISO 与 ZFINE2 链的谱形状几乎一致（0.999），粗 dz 的 B2D5CM 略低（0.982）；相对 L2 差 1.14 / 0.79 为格距不对齐下的整体差值，只报数值，不作收敛或精度声明。

## 硬限制（`results.json` hard_limits 全字段，逐条照录）

|字段|值|
|---|---|
|`cross_dimension_absolute_amplitude_compared`|false|
|`phase_not_transferable_across_dimensions`|true|
|`phase_reported_only_as`|"linear fit of the phase difference vs frequency; the equivalent time offset is a formal quantity, not a physical delay"|
|`physical_acceptance_threshold`|null|
|`training_labels_generated`|false|
|`clean_truth_generated`|false|
|`grid_convergence_certified`|false|
|`far_field_transform_is_conclusion_basis`|false|
|`conclusions_scope`|"limited to the dep3d_gold_v1 scenario family and the grids analysed; not extrapolated to field data"|
|`solver_invoked`|false|

即：本批不产生物理阈值、不产生训练标签、不产生 clean 真值；3D 网格收敛未认证；远场变换不作为结论依据；分析结论限定于 dep3d_gold_v1 场景族与所分析网格，不外推到实测数据；分析脚本仅做后处理复算（`solver_invoked=false`）。

## 其他限制

- 5 cm 3D 网格在 170 MHz 仅达官方起步准则（覆盖层 8.8 格/波长、砂岩 11.8 格/波长）；ANISO 例仅细化 z 向，不构成各向同性收敛认证。本批只给出"两个网格之间的差"，不给"距离真值的误差"。
- 2D 无限目标与 3D 有限目标的差异是被测量，不当作误差去"修正"。
- 材料参数（εr16/σ0.01、εr9/σ0.001、εr20/σ0.02）是研究假设，不是场地实测值。
- 墙钟滞留根因已在 batch2d_v1 执行期间鉴定为 `vctip.exe` 遥测进程（见"墙钟构成"根因补记），不影响输出完整性结论。
- 结论限定于本平层场景、既定几何与两套网格；不外推到实测性能、其他几何/材料/倾斜界面。

## 归档与复现

原始运行归档（每个目录含 h5、`.in`、stdout/stderr、supervision.json、record.json、authorization_at_start.json、显存采样）：

- `artifacts/research_checks/2026-09-26_B2D5CM_BG/`、`2026-09-26_B2D5CM_TGT/`
- `artifacts/research_checks/2026-09-26_B2DANISO_BG/`、`2026-09-26_B2DANISO_TGT/`
- `artifacts/research_checks/2026-09-26_DEP3D_5CM_BG/`、`2026-09-26_DEP3D_5CM_TGT/`
- `artifacts/research_checks/2026-09-26_DEP3D_ANISO_BG/`、`2026-09-26_DEP3D_ANISO_TGT/`

分析与复算：

- `artifacts/research_checks/2026-09-26_dep3d_gold_analysis_r1/` 与 `_r2/`：两次连跑 `results.json` 逐字节一致；`inputs` 字段覆盖全部 8 个 HDF5 及 ZFINE2 参照 `arrays.npz` 的 SHA-256。
- 输入契约与预算：`configs/research/dep3d_gold_v1/`（8 份 `.in` + `budget.json`）；静态核验脚本 `scripts/check_dep3d_gold_static.py`。
- 复现入口：`scripts/analyze_dep3d_gold.py`（CPU 后处理，仅复算比较，不调用求解器）。
- 8 次 attempt 各自一次，不可复用；旧会话不重启。

## 下一步

- 本批证据（3D 纯 dz 对照、同格距跨维形状差、墙钟/内存实测校准值）交由后续批量仿真规格书的网格与预算决策参考；是否升级任何 `reference_state` 由资格评估决定，本文不做声明。
- 墙钟滞留根因已鉴定为 `vctip.exe`（MSVC 编译器遥测）；batch2d_v1 采用"求解完成后清除孤立 vctip"的处置并全程留痕，处置效果与合规性见批量结果文档。
- 不生成物理/训练标签，不做超出 dep3d_gold_v1 场景族的外推。
