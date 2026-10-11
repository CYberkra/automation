# 证据包：材料消融批 r1（2026-10-11，x190，AGL 5 m）

报告：`docs/research/2026-10-11_material_ablation_results.md`

## 内容

- `parent_execution_contract.json` / `parent_execution.jsonl` / `parent_preflight_verification.json`：原 8 组契约（3d174683…）与中断事件账（1 完成 + SESSION_LEASE_EXPIRED）。
- `continuation_execution_contract.json` / `continuation_execution.jsonl` / `continuation_preflight_verification.json` / `continuation_completed_verification.json` / `continuation_predecessor_partial_audit.json`：续跑 7 组契约（8c13bd10…）与完成验证。
- `execution_contract.json` / `completed_verification.json`：合并 8 组审计视图（分析直接读取）。
- `combined_certificate.json`：合并证书（cc064ed1…）。
- `ablation_analysis.json`：分带峰、RMS、C3 归因检验数值。
- `ablation_traces.npz`：8 算例的 t/x/env/response 数组（图件复现用）。
- `figures/fig1`–`fig5`：中文标注图件。

## 复现

- ROG 侧分析：`python scripts/analyze_line9_material_ablation.py --execution combined_r1 --out analysis_r1`（需 4.0.1 环境与 SFCW 工具箱）。
- 本地图件：`python scripts/render_line9_ablation_figures.py --npz ablation_traces.npz --analysis ablation_analysis.json --out figures/`。

## 限制

数值收据一致性证据，不是物理/场地认证。原生 H5 存本机 `artifacts/local_checks/2026-10-11_ablation_results_r1/combined/`（忽略目录）与 ROG `E:\djh\line9_material_ablation_20261011_r1\combined_r1\`。
