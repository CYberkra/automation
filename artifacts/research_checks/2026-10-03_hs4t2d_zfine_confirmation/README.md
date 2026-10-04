# HS4T2D 最后中心竖向网格配对（完成2道）

gprMax V4.0.0 / CUDA double；同一15m离地高度中心道，竖向1.25cm起伏/全覆盖层配对，X/Y仍5cm，PML物理厚度1m。两份原始Ey为float64且有限，两份实际材料网格独立核验通过。冻结契约、输入、日志、原始H5、网格及身份随Git归档。本阶段没有重复几何需省略。

三档中心比较见[总报告](../../../docs/research/2026-10-03_hs4t2d_cause_diagnostic.md)。相邻网格变化缩小，但只有中心竖向序列，不是全空间收敛。

使用本项目固定V4官方处理代码，在仓库根目录CPU复算到不存在的新目录，不调用求解器：

```powershell
python scripts/archive_hs4t2d_cause_controls.py --study artifacts/research_checks/2026-10-03_hs4t2d_zfine_confirmation
python scripts/analyze_hs4t2d_z_confirmation.py --contract artifacts/research_checks/2026-10-03_hs4t2d_zfine_confirmation/execution_contract.json --first-results artifacts/research_checks/2026-10-03_hs4t2d_cause_results --confirmation artifacts/research_checks/2026-10-03_hs4t2d_cause_confirmation --out artifacts/local_checks/hs4t2d_z_rebuild
```

归档分析为同级 `2026-10-03_hs4t2d_zfine_results/`。64倍补零仅时间插值；所有比较仍为20–170MHz/501点/Hann，无增益。表中的固定地下窗包络最大值不是经独立确认的界面真值到时。
