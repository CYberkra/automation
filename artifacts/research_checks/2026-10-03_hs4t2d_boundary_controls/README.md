# HS4T2D 固定几何的边界厚度对照（32道完成）

V4.0.0 / CUDA double，32份原始Ey为float64且有限，32份实际材料网格逐个与既有36m宽域基准相同。仅改变指定PML面的厚度；地层/站位/网格/源/时窗均未变。

- 左右20→40单元：起伏/全覆盖层各13道。
- 左右20→60单元：中心配对。
- 底部或顶部20→40单元：分别中心配对。

20单元基准复用 `../2026-10-03_hs4t2d_cause_controls/wide36_{rough,halfspace}/`，不重跑已消耗attempt。所有原始H5及每组代表网格归档，重复几何本机保留，其身份及历史32网格验收见 `archive_selection.json` / `completed_verification.json`。

第一次预检后、首次求解前，仅将间距说明11.6m更正为保守的11.1m，没有改动任何输入字节。原契约与原准备源码保留在 `execution_contract_preflight_v0_1.json` / `preparation_source_v0_1.py`，当前执行契约为 `execution_contract.json`。两份预检分别对应其时刻的契约，不能互换源码哈希。

当前分析为 `../2026-10-03_hs4t2d_boundary_results_v0_2/`。第一版分析在 `../2026-10-03_hs4t2d_boundary_results/` 保留，后续仅补绝对范数和40→60中心比较。复算到尚不存在的新目录，不调用求解器：

```powershell
python scripts/archive_hs4t2d_cause_controls.py --study artifacts/research_checks/2026-10-03_hs4t2d_boundary_controls
python scripts/analyze_hs4t2d_boundary_controls.py --contract artifacts/research_checks/2026-10-03_hs4t2d_boundary_controls/execution_contract.json --out artifacts/local_checks/hs4_boundary_rebuild
```

需要本项目固定V4官方SFCW源码。先做复谱匹配差分，保留符号/相位，固定20–170MHz/501点/Hann/8倍补零/200ns尾部渐消。宽带原始场与带限产品分开评价；有限厚度对照不提供绝对PML误差界，也不证明5cm网格充分。详见[报告](../../../docs/research/2026-10-03_hs4t2d_manual_and_boundary_findings.md)。
