# HS4T2D 第一阶段原因对照（完成67道）

gprMax V4.0.0 / CUDA double；13个匹配站位、全覆盖层参考、36m延拓域、2m离地高度及中心竖向2.5cm配对。输入、契约在执行前冻结，`execution.jsonl` 记录67道完成，无求解失败。方案与限制见[总报告](../../../docs/research/2026-10-03_hs4t2d_cause_diagnostic.md)。

101道中的本阶段67道实际材料网格在执行机器全部核验通过。所有原始H5及每组一份实际网格随Git归档；同组其他几何文件保留本机，其身份见 `archive_selection.json`，全道历史独立核验见 `completed_verification_v0_2.json`（较早checker结果保留，不覆盖）。纯克隆的当前检查覆盖67份原始H5和7份代表材料网格，不重新宣称其余几何已在该克隆读取。

在仓库根目录，用已安装本项目固定V4官方SFCW处理代码的Python执行以下CPU复算；选择一个尚不存在的输出目录。不要再次执行已经消耗的求解契约。

```powershell
python scripts/archive_hs4t2d_cause_controls.py --study artifacts/research_checks/2026-10-03_hs4t2d_cause_controls
python scripts/analyze_hs4t2d_cause_controls.py --contract artifacts/research_checks/2026-10-03_hs4t2d_cause_controls/execution_contract.json --out artifacts/local_checks/hs4t2d_cause_rebuild
```

归档分析为同级 `2026-10-03_hs4t2d_cause_results/`。复谱保留相位，20–170MHz/501点/Hann；先复谱差分后重建。低高度只补回固定空气时延，形态图只做每例一个矩阵标量归一化，不作幅度/物理验收。加宽改变外部层延续，不是仅改变PML参数。
