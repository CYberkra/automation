# 任务书：A0 3D 校核配对分析脚本（草稿，kimi 验收后入库）

日期：2026-09-27。委派：CodeBuddy（hy4-preview）。性质：CPU 只读分析脚本实现；不跑求解器、不 git、不联网、不改 configs/。

## 背景

仓库 `E:\automation_djh\automation_repo`（Python 一律 `artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe`，仓库根为工作目录）。A0 3D 校核批次 `a0_3d_v1`（契约提案 `docs/research/2026-09-26_a0_3d_validation_contract_proposal.md`，**先通读 §1–§3 与 §8**）正在执行：2 例 3D 5 cm 各向同性（`B2D-C3m-BG-3D5CM` / `B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM`，C3 族 D10m 锚点），输出 h5 将在 `artifacts/simulations/2026-09-27_<RUN_ID>/`（h5 结构与 dep3d 相同：根 attrs 含 `Iterations/dt/dx_dy_dz/nx_ny_nz`，数据在 `/rxs/rx1/Ex`，float64）。

## 任务

实现 `scripts/analyze_a0_3d.py`：**以 `scripts/analyze_dep3d_gold.py` 为模板改写**（先通读它），口径严格按 A0 提案 §2.1 的 7 项比较指标与 §2.2 的"明确不做"清单。

比较结构：
1. **3D 内部配对差分**：`TGT − BG`（同域同格同源）。
2. **跨维只比形状**：3D 差分 vs 2D 侧 A0 锚点差分，分别对两档（**不跨格相减、不跨档合并**，各自成对照）：
   - 2D 粗档 BASE：`artifacts/research_checks/2026-09-26_B2D-C3m-BG/B2D-C3m-BG.h5` 与 `artifacts/research_checks/2026-09-26_B2D-C3m-D10m-W4m-T0.5m-E20-S0.02/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02.h5`（dy=dz=25 mm；若文件不存在则查 `artifacts/simulations/` 同名目录，并把实际路径写进 provenance；h5 根 attrs 有 dt，以 h5 实际元数据为准，不要硬编码 dt）
   - 2D 细档 FINE2：`artifacts/simulations/2026-09-26_B2D-C3m-BG-F2/` 与 `...-D10m-...-S0.02-F2/` 下的 h5（dy 12.5/dz 3.125 mm）
3. **官方 SFCW 后处理**：501 频点 20–170 MHz，`gprMax.toolboxes.SFCW.processing` 的 `load_source/load_receiver/direct_frequency_response`（模板里已有用法）；尾窗主 200 ns / 稳健 400 ns；**实波形诊断窗 240–400 ns**（A0 提案 §1.3 差异点 4，D10m 自算双程到时 320.2 ns ± 80 ns）。
4. 指标（提案 §2.1 全 7 项）：归一化谱形状相关（主量）、包络峰到时方向（只报方向与差值 ns）、谱峰频率差、时域差分波形最大归一化互相关+最优时移、选频点（20/50/80/110/140/170 MHz）归一化幅度趋势差 dB、相位差线性拟合等效时移（形式量）+残差 std、远场 3D→2D 变换后谱形相关（**探索性单列**，`far_field_transform_is_conclusion_basis: false`，模板已有实现可复用）。
5. **逐字照录**提案 §2.2 的六条"明确不做"到结果文件 hard_limits：`cross_dimension_absolute_amplitude_compared: false`、`grid_convergence_certified: false`、`direction_consistency_threshold: null`、`physical_acceptance_threshold: null`、`training_labels_generated: false` 等；并加 `solver_invoked: false`、与 dep3d D20m 类结果"可互引量级参照但非预期值"的声明。
6. 参照量级（只读引用、写入 results.json 的 reference_context 字段，不参与判定）：dep3d D20m 类实测跨维谱形状相关 0.935/0.958、远场 0.967/0.979、包络峰到时差 −71.02/−53.72 ns；2D 粗档 A0 包络峰 331.8626651806806 ns、细档 316.1921349940624 ns（来源：提案 §2.3 表）。
7. 输出完整性 assert（分析前）：h5 `nx_ny_nz=[160,200,1000]`、`dx_dy_dz=[0.05]*3`、Iterations=12464、dt≈9.629166e-11、Ex float64、全部有限值；BG/TGT 源样本一致（模板有同款检查可复用）。

## 输出

- 两次运行到 `artifacts/research_checks/2026-09-27_a0_3d_analysis_r1/` 与 `_r2/`（`results.json` + 必要 `arrays.npz`），results.json 两遍**逐字节一致**（无随机数、无时间戳；输入以 SHA-256 标识）。
- **但你本次不执行**：3D 仿真尚未完成。只交付脚本本身 + 一段 `--selftest` 模式：用 dep3d 的 5CM h5（`artifacts/simulations/2026-09-26_DEP3D_5CM_BG/DEP3D_5CM_BG.h5` 与 TGT 同名）作为伪输入跑通全流程（结构断言放宽到从 h5 实际读），证明脚本可运行；selftest 输出目录用 `artifacts/local_checks/a0_3d_analysis_selftest/`。kimi 随后验收 selftest 结果合理（例如与 dep3d 分析同输入时指标应一致或可解释），仿真完成后再由 kimi 正式跑 r1/r2。
- 脚本头部 docstring 写明口径来源（A0 提案 SHA 可用文件哈希自检但不断言——契约文本可能修订；结构断言必须断言）。

## 禁止

不碰 gate/契约 JSON；不跑求解器；不 git；不联网；不改 `analyze_dep3d_gold.py` 或任何既有脚本；不把 3D 与 2D 直接相减（各自内部差分后再比形状）。
