# 任务书：方向核查扫描参考池候选清单（草稿）

你是本仓库的协作 agent。任务：**纯整理**，不新增任何外部检索、不做判断性研究、不改任何既有文件。

## 输入（只读）

- `artifacts/research_checks/2026-09-26_direction_check_litscan/q1_uav_gpr_landslide.csv`（10 条）
- `artifacts/research_checks/2026-09-26_direction_check_litscan/q2_gprmax_grid_convergence.csv`（10 条）
- `artifacts/research_checks/2026-09-26_direction_check_litscan/q3_clutter_suppression.csv`（10 条）
- `artifacts/research_checks/2026-09-26_direction_check_litscan/q4_detectability_metrics.csv`（10 条）
- `artifacts/research_checks/2026-09-26_direction_check_litscan/q5_gain_compensation.csv`（10 条）
- 背景（只读，供理解用途）：`docs/research/2026-09-26_direction_check_litscan.md` 与同名 `_ledger.json`

## 产出（唯一允许写的文件）

`docs/research/2026-09-26_direction_check_litscan_pool_draft.md`

## 格式要求

1. 文件头注明：**草稿，由 codebuddy（glm-5.3-flash）整理，待 kimi 验收**；日期 2026-09-26；性质：参考池候选清单（非结论）。
2. 按 Q1–Q5 分节，每节一张 Markdown 表，**逐条列出该 CSV 的全部 10 条记录，一条不漏**：标题 / 作者(前3人+et al.) / 年份 / 被引数 / 期刊或来源 / URL / 一句话相关性备注（用 CSV 的 abstract 片段归纳，不得脑补）。
3. 表末加"候选标记"列说明：哪些条目适合登记入后续参照池，按这两类标记：
   - `[G4-pool]`：步进频/SFCW 衰减补偿与处理增益直接先例（Q5 中的 Noon 1996、Liu 2018 必选此项）；
   - `[roadmap]`：UAV-GPR 滑坡成像终点先例（Q1 中的 Wang 2026 必选此项）；
   - 其余不标记。
4. 结尾加一节"覆盖核对"：声明 5 个 CSV × 10 条 = 50 条全部列入，无遗漏、无添加。

## 禁止事项

- 禁止修改/删除任何既有文件（含 CSV、ledger、报告、configs、scripts）。
- 禁止联网检索、禁止访问仓库外路径。
- 禁止写任何数值结论（能量比、相关系数等一律不写）。
- 禁止 git 操作。

完成后，最后一行输出：DONE 及产出文件路径。
