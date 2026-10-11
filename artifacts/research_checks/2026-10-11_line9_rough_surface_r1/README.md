# 证据包：地表微起伏敏感批 r1（2026-10-11，x190，AGL 5 m）

报告：`docs/research/2026-10-11_rough_surface_results.md`

## 内容

- `execution_contract.json`（9286c39f…）/ `execution.jsonl` / `preflight_verification.json` / `completed_verification.json`：2 组契约与完成验证。
- `manifest.json`：准备清单（含分形参数与偏离披露）。
- `rough_stats.json`：分形微扰统计（RMS 5.4 cm，28191 体素）。
- `rough_analysis.json`：敏感度、差分稳健性、分带峰数值。
- `rough_traces.npz`：4 条迹线（rough/smooth × H0/H1）的 t/x/env/response。
- `figures/fig1`–`fig4`：中文标注图件。

## 复现

- ROG 侧分析：`python scripts/analyze_line9_rough_surface.py --rough execution --smooth E:\djh\line9_material_ablation_20261011_r1\combined_r1 --out analysis_r1`。
- 本地图件：`python scripts/render_line9_rough_figures.py --npz rough_traces.npz --analysis rough_analysis.json --out figures/`。

## 限制

数值收据一致性证据，不是物理/场地认证。分形参数为试验选择。原生 H5 留 ROG 批目录。
