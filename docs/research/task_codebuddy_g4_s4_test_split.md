# 任务书：G4 第 5 步 S4——`--split {dev,test}` 扩展（测试族 {C5,C8} 阶梯与能力导出数据生成）

> 委派 CodeBuddy 实现，kimi 验收。本任务书是唯一规格来源。**你只允许修改两个文件**：`scripts/run_damage_ladder.py`、`scripts/run_g4_mission_capability.py`。**不得**新建文件、不得改任何其他文件、不得运行任何命令（你没有 Bash）、不得触碰 configs/、docs/、artifacts/。

## 背景（只读上下文）

- G4 第 5 步锁定程序已冻结（`configs/research/g4_threshold_derivation_v1.0.json`），S1–S3 已完成（开发组阈值推导已执行并收口）。S4 = 在**测试族 {C5,C8}** 上生成损伤阶梯记录与能力导出（纯数据生成，**无阈值、无判定、无排名**）。
- 测试族数据事实：仅 MT 几何——C5/C8 各 3 个母模型（`B2D-C5m-BG` / `B2D-C5m-D10m-NC` / `B2D-C5m-D10m-W4m-T0.5m-E20-S0.02-MT33`，C8 同构）；CO 批只有 C3 族案例，测试族的 CO 组合自动为空（`expand_dev_rows` 现有逻辑已自然产生此结果，不要补造）。

## 改动一：`scripts/run_damage_ladder.py`

1. argparse 增加 `--split`，choices `['dev','test']`，**默认 `dev`**。
2. `split='dev'` 时行为必须与当前版本**逐字节一致**：同一 DEV_FAMILIES `('c1','c3')`、同一记录内容、同一 manifest 全部字段与取值（含 `dev_families`、`gates.dev_groups_only: True`、`hard_limits.test_groups_excluded`、`conclusions_scope` 原文）。kimi 将重跑 dev 全部 8 chunk 并比对合并剥离 SHA `6bc13f63c5985fa0b38b4724c2ffa0cf61716053751746520766ef6fa0e4d0d9`。
3. `split='test'` 时：族集合为 `('c5','c8')`；其余逻辑（损伤实例、矩形、指标、chunk 协议、确定性）完全不变。manifest 相应调整：
   - 新增 `'split': 'test'` 与 `'families_used': ['c5','c8']`；`'dev_families'` 字段保留但值为 `[]` 并在旁加 `'note'` 或改名字段需保持 schema 自洽（推荐：保留 `dev_families: []` + 新增字段，避免破坏既有读取方的键存在性假设）。
   - `gates.dev_groups_only` → `False`，新增 `gates.test_groups_only: True`。
   - `hard_limits.test_groups_excluded` 移除，替换为 `'dev_groups_excluded': ['B2D-C1m','B2D-C3m']`；`conclusions_scope` 字符串中的 "dev-group" 改为 "test-group"。
   - 其余 hard_limits 逐字保留（reference_state、training_eligible=false 等不变）。
4. docstring 顶部加一段：测试族模式属 S4 数据生成，阈值应用只在 S5 且只一次。

## 改动二：`scripts/run_g4_mission_capability.py`

1. argparse 增加 `--split`，choices `['dev','test']`，默认 `dev`。
2. `split='dev'`：逐字节一致（kimi 回归门：重跑结果与 `artifacts/research_checks/2026-09-26_g4_mission_capability_v3_r1/results.json` SHA-256 `2091f4171c16341731f57bb94c8fead0a72e47a1d58389caa74926286bf8e6c8` 相同）。注意当前文件的 CHUNK_GLOB 指向 r1 08 块、记录数断言 34086、DEV_FAMILIES 过滤——这些在 dev 模式全部不变。
3. `split='test'`：
   - 阶梯来源改为 `artifacts/research_checks/2026-09-27_damage_ladder_test_r1_c{:02d}`，N_CHUNKS 仍为 8；
   - 记录数 34086 断言仅对 dev 生效；test 模式断言 `len(recs) > 0` 并在结果中如实记录合并条数与来源 glob；
   - 族过滤改 `('c5','c8')`；
   - 结果 JSON 新增 `'split': 'test'`，`'ladder_source'` 字符串改为测试族 glob 并注明 r1（r2 字节一致由验收证明），`'dev_families'` 字段改为 `'families': ['c5','c8']`（dev 模式字段名与值保持原样不动）；hard_limits 全部逐字保留。
   - 聚合逻辑（mission_levels 映射、describe/p80、分层、nc 收集）一行不改。

## 纪律

- 不应用任何阈值；不触碰 ε；不排名不选优；不训练；不调求解器。
- 两脚本保持确定性（无随机、无时间戳进入结果 JSON 正文；manifest 的 total_wall_s 除外，它在既有 dev 产物中也存在）。
- 完成后报告：改动摘要、dev 路径为何逐字节不变（逐点说明）、test 模式 manifest/schema 差异清单、你无法验证的事项。
