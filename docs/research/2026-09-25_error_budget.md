# 算子链输出误差预算：证书、契约与方向诊断恢复

本单元为[方向诊断契约 v0.1](../../configs/research/direction_diagnostic_contract_v0.1.json) 的 `direction_error_bound` 建立**可计算的预算依据**，按其"后续需要分别建立输入形成、算子运算、差分相消及参考误差的依据"的待办执行。纯 CPU/NumPy 数组工作：无 FDTD、无训练、无实测读取。新契约 [error_budget_contract_v0.1](../../configs/research/error_budget_contract_v0.1.json)，实现 `scripts/research_error_budget.py`，检查 `scripts/check_error_budget.py`，应用 `scripts/study_error_budget_restoration.py` 与 `scripts/attest_three_grid_direction.py`。证据目录 `artifacts/research_checks/2026-09-25_error_budget/`，旧归档与哈希不变。

## 1. 预算组成与依据类型（契约要点）

预算 = 固定评价窗内处理输出的 L2 绝对误差上界（输入幅值单位），由四部分组成：

- **e_in 输入形成误差**：仅当输入与真值是同一构造数组对象时给零（`exact_construction` 证书）；FDTD 来源一律拒绝给预算（`fdtd_numerically_unresolved`）；其他来历不明拒绝（`unknown_input_provenance`）。
- **e_op 算子运算误差**：逐元素 IEEE 计数舍入界。每个浮点加/减/乘/除贡献 ≤ ½ ulp（`np.spacing/2`）；`math.fsum` 是正确舍入求和，贡献 ≤ ½ ulp。实测 `|float64 − fsum 重算路径|` 加上逐操作计数的 slack 即为证书。**没有任何 eps×常数或经验因子**；所有数字都由实际数组算出。
- **e_svd 子空间摄动（仅 SVD 截断算子）**：对每个保留的奇异三元组，用 fsum 路径算残差 `r_i = max(||Xv_i−σ_i u_i||, ||Xᵀu_i−σ_i v_i||)`，取与相邻奇异值的谱间隙 `gap_i`，按 Wedin/Davis–Kahan sin Θ 定理（Golub & Van Loan 型标准矩阵摄动理论）`sinΘ_i ≤ √2·r_i/gap_i`，分量输出摄动界 `b_i = r_i + 2√2·σ_i·sinΘ_i`（含 Weyl 的 |Δσ_i| ≤ r_i）。`gap_i ≤ 0` 或 `sinΘ_i ≥ 1` 时界失效，该行预算缺失并记 `svd_subspace_bound_vacuous`——允许缺失，不许硬凑。
- **e_diff 差分相消**：`e_diff ≤ e_1 + e_2`（三角不等式），另输出相消放大系数 `(||a||+||b||)/||a−b||` 作为风险诊断。
- **e_ref 参考误差**：单列，永不在输出预算里折叠。

**传播规则**：线性算子 `||Ae|| ≤ σmax(A)||e||`；恒等、行均值（λ∈[0,1] 的 I−λP 为正交投影类，σmax=1）、SVD 分量移除（正交投影移除，σmax=1；输入误差非零时不认证，记 `svd_input_error_not_certified`）、固定共享增益（σmax=end_gain，增益曲线按存储的 float64 数组定义为算子语义）。逐算子链式合成；当前步的逐元素算术界可按窗收紧，更早步与摄动项按全数组 L2 保守携带（任何窗内 `||E[mask]|| ≤ ||E||`，合法但偏保守）。

**边界**：预算只认证**数值可辨识性**（输出是否显著大于其数值误差），不认证物理精度；不是 GPR 探测阈值、训练标签或行业验收标准。

## 2. 证书与实测误差的一致性（20 项新检查）

`check_error_budget.py` 覆盖：恒等链预算恰为 0；均值/增益/链式证书 ≥ 对 fsum 参考及 Dekker TwoProd 精确乘积误差实测值（保守）；SVD 大间隙给预算且 ≥ fsum 重构实测算术误差；摄动公式数值核对；间隙失效（含算子层 `svd_cutoff_gap_unresolved`）与人为小间隙两情形均缺失；差分三角和与放大系数；构造数组恢复余弦与直接手算一致；相消后小输出仍缺失（`output_direction_numerically_unresolved`）；FDTD/未知来历拒绝；整体 ×1e6 换单位判定不变；窗预算 ≤ 全窗预算；非法输入界拒绝。典型实测：均值预算 4.3e-15 vs 实测 1.3e-15；链式预算 1.1e-14 vs 实测 1.3e-15（48×24 随机数组，O(1) 幅值）。

## 3. 损伤试验方向诊断恢复

损伤试验 63 组/42 行全部是构造数组（成分已知、真值即输入对象），预算来源显式声明为 `exact_construction + fsum 重算 + Wedin/Davis–Kahan`。重跑其方向诊断（42 行 × 整窗 + 49 个事件窗 = 91 个余弦评价）：

| 类别 | 恢复余弦 | 仍缺失 | 预算量级（O(1) 幅值数组） |
|---|---|---|---|
| identity | 13/13 | 0 | 恰为 0 |
| mean | 30/39 | 9 | ≤ 1.5e-14 |
| svd | 8/39 | 31 | 3.9e-14 ~ 1.2e-11（含摄动项） |
| **合计** | **51/91** | **40/91** | 最大 1.17e-11 |

整窗余弦恢复 24/42 行；事件窗 27/49。仍缺失的 40 个全部是输出范数 ≤ 预算（`output_direction_numerically_unresolved`）——主要是全均值/SVD 把构造场景扣到接近零的合规缺失，与契约设计一致：宁可缺失不可硬凑。恢复的余弦与无预算直接计算的数值逐一一致（检查 `restored_cosine_matches_direct`）。

## 4. 三网格 84 个余弦保持缺失（诚实性验收）

`attest_three_grid_direction.py` 读取现行权威归档 `2026-09-25_review_fixes/joint_roles_results.json`，核验 42 行 × layer/target 共 **84 个余弦全部为 null**（原原因 `direction_error_bound_not_declared`），另存细化原因 `fdtd_numerically_unresolved` 的证明文件到新目录；**归档本身未修改**。FDTD 模板在预算契约下不可声明预算，方向诊断继续关闭。

## 5. 回归与限制

- `python scripts/verify_workspace.py`：六组 157 → **八组 182 项**（新增 check_error_budget 20 项 + study_error_budget_restoration 5 项），全绿。仍纯数组，无 FDTD/实测/网络。
- 限制：窗预算携带全数组保守项（偏松不错误）；SVD 输入误差传播未认证（构造精确输入下不涉及）；增益曲线按存储值定义，曲线自身构造误差不入证书；预算不扩展为物理探测阈值；`training_enabled` 保持 false。
