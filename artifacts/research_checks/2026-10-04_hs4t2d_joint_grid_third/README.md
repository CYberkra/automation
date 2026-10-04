# X/Z 1/60m中心收缩检查（2道）

V4.0.0/CUDA double；同一36m域、15m航高、原台阶材料、源/接收器和物理PML厚度。实际材料数组独立核验等于5cm粗网格沿X/Z各重复3倍。两份原始H5与实际网格全部归档。

首次预检后、首次求解前，仅修正契约继承的用途说明与保守显存估计措辞；输入字节不变。原契约、原准备源码和原预检保留为 `execution_contract_preflight_v0_1.json`、`preparation_source_v0_1.py`、`preflight_verification_v0_1.json`。当前契约和预检用于实际执行。1.25cm显存描述是执行前保守外推，不是测得的硬件容量结论；本轮没有运行该档。

本组只认证中心相邻网格差异变化，不代表全13站位第三档收敛。分析见 `../2026-10-04_hs4t2d_joint_grid_third_results/`，包括有限的95MHz路径诊断；最早路径不等于最强回波。

复核：`python scripts/archive_hs4t2d_cause_controls.py --study artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_third`。

CPU重建：`python scripts/analyze_hs4t2d_joint_grid_third.py --third artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_third --scan-results artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_results --out <不存在的新目录>`。无求解、训练或实测。
