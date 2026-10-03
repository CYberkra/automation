# X/Z 2.5cm其余12站位（24道）

先导中心已独立验收后，另冻本契约并完成左/右各6站位的起伏和全覆盖层配对。中心61号不重跑，共用原有5cm/36m域/左右2m PML基准。24实际材料数组等于原粗网格沿X/Z各重复2倍。

全部24原始H5与每组代表网格随Git归档，20重复网格本机保留，原身份及执行时逐道检查见 `archive_selection.json` / `completed_verification.json`。克隆复核应使用归档检查，不把本机全道历史审计声称为克隆时重新检查遗漏网格。

复核：`python scripts/archive_hs4t2d_cause_controls.py --study artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_remaining`。

CPU分析：`python scripts/analyze_hs4t2d_joint_grid.py --centre artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre --remaining artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_remaining --out <不存在的新目录>`。保相位/符号及物理共用色标；无SVD/AGC/逐道归一化，非clean标签或物理验收。
