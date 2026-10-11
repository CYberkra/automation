# 归档：AGL 8 m 版模型（2026-10-11 被 5 m 版取代）

- 取代原因：2026-10-11 用户显性批准航高默认 8 m → 5 m（参考 Wang 2026 论文滑坡仿真航高 5 m 及其实测"11 m vs 16 m 差异很小"结论）。
- 内容：`cases/x190_H{0,1}_local/profile.in` 为 AGL 8 m 精细波场批（line9_fine_wavefield_20261011_r1）的诊断算例输入（含 31 探针×6 锚点和快照定义）；geometries/source/acquisition 与现行版逐字节相同（该次变更只动天线高度）。
- 对应 8 m 结果证据：`artifacts/research_checks/2026-10-11_line9_fine_wavefield_r2/` 等，报告 `docs/research/2026-10-11_fine_wavefield_results.md`。
- 注意：8 m 与 5 m 结果跨航高对比时须考虑空气路径差（单程约 10 ns）。
