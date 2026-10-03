# X/Z 2.5cm中心先导（2道）

V4.0.0/CUDA double，起伏与全覆盖层各1道，复用已有5cm/36m域/左右2m PML中心基准。实际材料等于粗网格沿X/Z各重复2倍；收发位置和Y线源长度不变。输入/原始H5/全道网格验收与执行记录保留。

当前中心分析见 `../2026-10-04_hs4t2d_joint_grid_centre_results/`，全13站位分析见 `../2026-10-04_hs4t2d_joint_grid_results/`。中心数据与左右段组装，不再求解已消耗的中心attempt。两档网格不构成空间收敛证书。

复核：`python scripts/archive_hs4t2d_cause_controls.py --study artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre`。精确CPU重建用 `scripts/analyze_hs4t2d_joint_grid.py --centre <此目录> --out <不存在的新目录>`，官方固定V4环境，无求解。
