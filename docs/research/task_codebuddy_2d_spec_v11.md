# 任务书：2D 批量仿真规格书 v1.1 起草（固化 A0/ZFINE2/FINE2 收敛与 3D 校核证据）

> 委派 CodeBuddy 起草，kimi 逐条核验后入库。**你只允许新建一个文件**：`docs/research/2026-09-27_batch_2d_spec_v1_1_draft.md`（注意文件名带 `_draft`）。不得修改任何既有文件（v1.0 规格书是已执行批次的冻结引用，永不改写；v1.1 以新文件增补/取代未来批次口径）。不得运行命令、不得联网。

## 背景与必读材料（全部在仓库内，逐条引用须给出文件路径+节/锚点）

1. `docs/research/2026-09-26_batch_2d_spec_v1.md` —— v1.0 规格书（既有口径，v1.1 在其上增补）。
2. `docs/research/2026-09-27_a0_3d_results.md` —— A0 3D 校核结果（D10m 锚点：2D/3D 谱形状相关 0.9686/0.9729、远场变换 0.986/0.988 探索性、包络峰 3D 晚 +2.46/+18.13 ns、波形互相关弱 0.29/0.05、四类混淆未分离、硬限制六条）。
3. `docs/research/2026-09-26_dep3d_gold_results.md` —— D20m 类 3D 校核（0.935/0.958；3D 纯 dz 细化 0.9951）。
4. `docs/research/2026-09-26_batch2d_v1_results.md`（BASE 29/29，成本）与 `docs/research/2026-09-26_batch2d_v1_fine2_results.md`（FINE2 8/8，方向一致性：形状相关 0.99813–0.99955、符号一致率 88.0–93.8%、包络峰细档早 15.7–19.1 ns）。
5. `docs/research/2026-09-25_deep_zfine2_results.md` —— 垂向细化收敛链（全带 L2 66.6647%→17.2452%→4.1216%，相位 64.8431°→16.2741°→3.8847°，未做 Richardson 外推、收敛未认证）。
6. `configs/research/batch2d_v1/cases.json` 与 `configs/research/batch2d_v1/groups.json` —— 场景族/案例/分组事实。
7. `configs/research/project_context_v1.json` —— 设备与场景背景（SFCW 20–170 MHz、501 频点、约 20 m 浅层非显性滑坡）。

## v1.1 文档结构（逐节写，每节标注证据来源）

1. **版本关系**：v1.1 增补 v1.0、不 retroactive 改写；已执行批次仍引用 v1.0。
2. **网格依据**：BASE 25 mm 档的证据链（ZFINE2 收敛链数值 + FINE2 方向一致 + A0 3D 形状同形），逐条标注"收敛未认证、Richardson 未做"的限制原文。
3. **PML/时窗**：按 v1.0 与 cases.json 实际值照录。
4. **材料表**：按 cases.json/project_context 照录（含 Peplinski 频段限制警告——0.3–1.3 GHz 不直接作本项目低频真值）。
5. **单例成本**：实测值（BASE 29 例约 23.4 min 全批、FINE2 约 1320 s/例、3D 5CM 约 909–916 s/例），标注硬件（RTX 3060 Laptop, CUDA double）与"实测非 ETA 承诺"。
6. **场景族与种子分组**：按 groups.json 照录四族定义与开发/测试划分（{C1,C3}/{C5,C8}）。
7. **2D/3D 差异与适用范围（v1.1 新增核心节）**：综合 A0（D10m）与 dep3d（D20m）两锚点证据——谱形状类指标 2D 可代 3D 作批量主力（0.94–0.97 量级）；波形逐样本/绝对到时不可跨维直接比（波形互相关 0.05–0.74、到时方向两锚点相反不可迁移）；绝对幅值永不跨维比较；远场变换仅探索性。列出仍未覆盖的 3D 校核缺口（C5/C8 锚点、电性档）。
8. **不变限制**（照录）：reference_state=numerically_unresolved、网格收敛未认证、无物理/训练标签、阈值 null、结论限仿真域。
9. **批量跑启动条件清单**：未来新批次启动前必须满足的契约/gate 条件列表（冻结契约、gate 激活留痕、attempt 单次、监督器、预算回填）。

## 纪律

- 每个数值必须能从上述文件逐字溯源；不得发明任何数字；不确定处列"待 kimi 裁决"清单。
- 文档头标注：`codebuddy(glm-5.3-flash) 起草，kimi 验收前为草稿`。
- 中文正文、技术名词保留英文原词；硬限制/纪律句逐字照录不改写。
