# 新版 0.8 m 界面：一次中心道波场快照

用户明确要求后另冻契约；V4.0.0 / CUDA FP64 一次求解，九个时刻（请求 40–190 ns），t61 原始接收器 Ey 与未加快照输出逐元素相同。九份 H5 的 54 分量数组均 float64/有限值，时间交错独立核验通过。

`profile_snaps/` 是原始六分量场快照，轴为 [X,Y,Z]，E/H 有半步时间差；不能误当全部同时间场。原始宽带 impulse，不是 20–170 MHz 带限场。图/动画在同级 `2026-10-03_hs4t2d_relief08_snapshot_visual/`，共用带符号对数色标。

`execution_contract.json` 是执行前契约，`execution.json` 是完成与逐文件哈希，`independent_verification.json` 为独立复核。新增 `manifest.json` 追踪整个胶囊（不含清单本身）。[解释与限制](../../../docs/research/2026-10-03_hs4t2d_relief08_comparison.md)；当前未据该快照签认 B-scan 形态或识别唯一根因。
