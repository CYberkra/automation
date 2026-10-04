# 二维起伏 0.8 m：选择与完成状态

用户选择 relief2p0（全模型 0.8 m 高差）。原版、新版、深度 2.95 m 平界面各 121 道已完成 gprMax V4.0.0 / CUDA 双精度。所有原始 Ey 为 float64，原版 121 道与旧输出逐元素相同；363 实际材料网格与输入一致。

`proposal.json` / `manifest.json` 为执行前准备快照，0.6 m 推荐是当时建议，未执行。执行状态以 `execution_contract_cuda_cached_att4.json`、`completed_verification.json` 和 [完整报告](../../../docs/research/2026-10-03_hs4t2d_relief08_comparison.md) 为准。失败尝试和每次契约保留。新增文件由 `package_manifest.json` 单独追踪，旧准备清单不覆盖。

新版 `relief2p0/model.obj` 与 `hs4_scene.mtl` 为实际台阶剖面；可解压 `model_browser.zip` 后一起拖入 https://3dviewer.net/。显示不包含实际电磁场或三维有限目标，浏览器材质不是电性。

求解原始数据、处理结果、对比图分别在同级 `2026-10-03_hs4t2d_relief08_raw`、`...results`、`...visual_v0_2`。条带随起伏变化，但拱形对应仍不清楚；不能据输入网格正确宣布物理形态问题已解决。
