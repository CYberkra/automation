# Line9 时窗收敛验证（4800 ns 延长窗对照）

日期：2026-10-12。状态：完成（3 站最差尾部样本）。

## 结论

1. 夜间 41 道批 2400 ns 窗口的尾部截断对本场景 501 音调复谱的影响很小：在最差 3 道（X170.0/170.5/171.0）上，截 2400 ns 与 4800 ns 全长的官方 direct 复谱相对 L2 差异为 0.09%–0.15%，单音调最大相对差异 0.26%–0.55%。
2. 4800 ns 窗口下官方尾部判据（末 5% 峰值相对电平）为 −71.7 至 −72.1 dB，已低于官方 −60 dB 建议；夜间批 2400 ns 时该指标为 −57.3 dB，是当时的未签认项。
3. 「夜间 2400 ns 原生求解」与「4800 ns 记录截 2400 ns」逐比特一致（相对 L2 = 0），截断等价性严格成立——差异只来自被截掉的晚到信号，不来自重跑数值噪声。

## 依据与边界

- 批次：ROG `E:\djh\line9_window_20261012_r1`，契约 `bcf1ca2b…`（3 组 4800 ns，除 `#time_window` 外与夜间批逐字节一致），3×约 717 s 完成，`completed_verification.json` PASS。
- 分析：`scripts/analyze_line9_window_convergence.py`；结果 `artifacts/research_checks/2026-10-12_line9_window_convergence_r1/analysis.json` 与 `window_convergence.png`。
- 夜间道哈希对照夜间交付 `analysis.json` 逐道核验。
- 边界：只覆盖尾部最差 3 道，不是全线收敛证书；但选题是全线最差尾部，其余 38 道夜间指标已优于这 3 道。
- 对地下解释的意义：截断造成的 ≤0.55% 单音调差异远小于后续任何处理/反演的误差预算，2400 ns 窗口对本场景（20 m 以浅 + 多次波族）可用；如需更保守可全线改用 4800 ns（成本约 2 倍）。

## 环境备注

本次运行前 ROG 的 conda Python 3.12.14（py312conda）被 Device Guard 策略封禁（2026-10-05 事故的同类，哈希封禁）。改用 python.org 官方签名 Python 3.12.10 可携版（`E:\djh\artifacts\local_checks\py312embed`）+ 既有两个 venv 的 site-packages（cp312 ABI 兼容，含重编 virtual_waveguide.pyd），冒烟验证 numpy/h5py/scipy/gprMax 4.0.1/PyCUDA/RTX4090 全部通过。旧路径未删除，策略恢复后仍可对照。
