# HS4T2D 第二阶段源与竖向网格确认（完成32道）

gprMax V4.0.0 / CUDA double；6个Ricker95MHz中心道（12/24/36m域，起伏/全覆盖层配对）和26个竖向2.5cm的13站位配对。全部原始Ey为float64且有限，32份实际材料网格独立核验通过。细网格序列中心道与第一阶段细网格中心道逐元素相同。方案及结果见[总报告](../../../docs/research/2026-10-03_hs4t2d_cause_diagnostic.md)。

全部原始H5与每组一份实际网格随Git归档；重复几何本机保留，身份见 `archive_selection.json`，全道历史网格核验见 `completed_verification.json`。纯克隆当前可检查32份H5及8份代表材料网格。

使用固定V4官方处理源码在仓库根目录执行以下CPU复算；输出目录必须尚不存在。第一阶段归档数组会先验哈希。不要重跑已消耗求解契约。

```powershell
python scripts/archive_hs4t2d_cause_controls.py --study artifacts/research_checks/2026-10-03_hs4t2d_source_controls
python scripts/analyze_hs4t2d_source_controls.py --contract artifacts/research_checks/2026-10-03_hs4t2d_source_controls/execution_contract.json --first-results artifacts/research_checks/2026-10-03_hs4t2d_cause_results --first-controls artifacts/research_checks/2026-10-03_hs4t2d_cause_controls --out artifacts/local_checks/hs4t2d_source_rebuild
```

归档分析为同级 `2026-10-03_hs4t2d_cause_confirmation/`。源对照先除以实际保存的源谱，不能直接比较Ricker与impulse的原始幅度；地下固定窗和全记录复谱统计分开，不以窗最大峰跨度签認几何到时。
