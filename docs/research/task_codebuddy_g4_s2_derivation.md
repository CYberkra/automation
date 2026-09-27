# 任务书：G4 第 5 步 S2 — 阈值推导脚本实现（草稿，kimi 验收后入库）

日期：2026-09-27。委派：CodeBuddy（hy4-preview）。性质：CPU 只读分析脚本实现；不跑求解器、不 git、不联网、不改任何 configs/ 下文件。

## 背景（你需要的全部上下文）

仓库 `E:\automation_djh\automation_repo`（Windows；Python 一律用 `artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe`，仓库根为工作目录）。G4 阈值校准已到第 5 步 S1：推导程序已冻结为 `configs/research/g4_threshold_derivation_v1.0.json`（SHA-256 `d114fd193259008a4882ba6a3e3af3ab9ed507fd608fe3cf72a2432d99283eca`）。**先读这个文件全文**——它是唯一规格来源。你的任务是实现 S2：`scripts/run_g4_threshold_derivation.py`，按冻结程序在开发组数据上执行推导。

## 输入（全部 sha256 门禁，脚本启动时逐一来核验，不符即退出）

| 输入 | 路径 | SHA-256 |
|---|---|---|
| 推导程序（冻结） | `configs/research/g4_threshold_derivation_v1.0.json` | `d114fd193259008a4882ba6a3e3af3ab9ed507fd608fe3cf72a2432d99283eca` |
| 任务容差 v0.2 | `configs/research/g4_mission_tolerance_v0.2.json` | `ee039fb1fab3f1147e8a5fe53809bbb8ae9d6f51aafc2feee0fb90071a6c57ea` |
| 稳定规则 | `configs/research/g4_fliprate_stability_rule_v0.1.json` | `8ee61d0b7777f188d84fd754e3302ade9fe678ed0eae1b7f3f0ade4c451b65f6` |
| 事件表 | `configs/research/batch2d_v1_event_table_v0.1.json` | `b0ad100334132cb6e6a706af4f5d8fbbff7cced26b6e07613c50be77a607f56c` |
| S4 划分 | `configs/research/s4_baseline_group_v0.1.json` | `8676f37d6aee603b7f2481779ffe361d2878735316863da0bed0db82e16b5593` |
| 能力导出 v3 | `artifacts/research_checks/2026-09-26_g4_mission_capability_v3_r1/results.json` | `2091f4171c16341731f57bb94c8fead0a72e47a1d58389caa74926286bf8e6c8` |
| 翻转率 v2 | `artifacts/research_checks/2026-09-26_g4_sensitivity_fliprate_v2_r1/results.json` | `f96d317a3e046dfcc9efb01da55d49f78ab09c6fdb562c37e4ce93c6480e558f` |
| 参考预算 | `artifacts/research_checks/2026-09-26_reference_uncertainty_budget_r1/results.json` | `38cf98bab723358dc8d75486ef6824057d42c32a9030b56f7101400fc87a8cc1` |
| 损伤阶梯合并 | `artifacts/research_checks/2026-09-26_damage_ladder_r1_c00..c07/records.json` | 合并剥离 resource 后 SHA-256 `6bc13f63c5985fa0b38b4724c2ffa0cf61716053751746520766ef6fa0e4d0d9`（合并算法照抄 `scripts/check_damage_ladder_acceptance.py` 的 `merged()`：按块序、逐记录 pop resource、sort_keys 规范 JSON） |

## 既有可复用代码

`scripts/run_g4_sensitivity_fliprate.py`：含 `pair_list / sign_compare / level_axis / window_axis / value_axis / build_strata / p80 / quantile / ordering` 等函数。**逐对重导必须复用这些函数或其逐字逻辑**（含 value_axis 的 `max(0.0, ...)` 截断、随机种子协议 `random.Random(f'{stratum_key}|{k}|{cid}')`、tie 精确相等排除口径），使逐对翻转计数与 v2 聚合完全一致。

关键数据结构（已核验）：
- 能力 v3 `results.json`：顶层 `capability`（2065 行：candidate_id/damage_type/damage_level/family/geometry/mission_class/D_p80/D_median/n 等）+ `negative_control`（324 行：candidate_id/family/geometry/nc_amplify_q/N_b_energy_ratio_p80/N_b_energy_ratio_median/n）。
- 翻转率 v2 `results.json`：`strata` 24 个，键 = (quantity, damage_type, mission_class, geometry, family)；每个含 `full_values`（候选→p80 值）、`full_ordering`、`level_axis`/`window_jackknife`/`value_jitter` 聚合（flips/comparisons/ties）、`value_jitter.epsilon`、`levels`。
- ε_stratum = 翻转率 v2 对应 stratum 的 `value_jitter.epsilon`。

## 要实现的功能（严格按冻结配置的 derivation_rules）

1. **D 约束**：对每个 (geometry, family) × 候选 × damage_type ∈ {amplitude_scale, polarity_flip, sample_shift, trace_deletion}：该类型任务档（mission_relevant）**每个冻结级别**（从任务容差 v0.2 `ladder_level_mapping.mission_relevant` 读，含符号对称移位网格）的 `D_p80 ≤ ε[geometry,family,damage_type,mission_relevant]`。候选在 (geometry,family) 通过 D 约束 ⟺ 4 个损伤类型全过。能力行缺失/不可用照录为 unavailable，不记零分、不记通过。
2. **N 约束**：每个 (geometry, family) × 候选：q=2 与 q=4 两条负控行的 `N_b_energy_ratio_p80` 都 ≤ 1.0 才通过（冻结配置的保守联合读法）。
3. **admissible** = D 约束 ∧ N 约束。某 stratum 无 admissible ⇒ 如实报 `no_admissible_candidate`，禁止放宽。
4. **稳定规则一次性应用**：仅对 admissible 候选、在翻转率 v2 的 21 个 D strata 内逐对重导三轴翻转（级别轴=冻结符号对称网格全级别对；window_jackknife=留一事件；value_jitter=K=200 确定性抖动，ε 用该 stratum 的 v2 值）。**内置验收门**：逐对重导的聚合必须逐 stratum 逐轴精确复现 v2 的 flips/comparisons/ties，任何不符即报错退出（exit≠0）。3 个 N_b strata 不做逐对稳定判定，只在报告中列 `cross_scale_caveat: true` 的说明行（N_b 约束判定本身照常）。
5. **stratum 判定**：通过约束的候选中，某候选对全部其他 admissible 候选 determined-better ⇒ `unique_stable_top`；否则 `set`/`partial_only`（对齐协议 v0.2 §7 语义）。只报结构，不做选型。
6. **能力表述（a80 式）**：每 (geometry,family,damage_type) × 候选输出任务档级别上 D_p80 的上界陈述 + **单列的系统偏差项**：该 stratum 事件窗在参考预算中的包络峰 BASE↔FINE2 偏移（ns）中位数（网格链方向量，只报数值不做深度换算、不混入统计项）。弱事件带级别永不进能力头条。
7. **hard_limits 逐字照录**冻结配置的 hard_limits 到输出；输出含 `g4_relieved: false`、`thresholds_are_simulation_domain_only: true`、构造参考声明。

## 输出

- `artifacts/research_checks/2026-09-27_g4_threshold_derivation_r1/results.json`（你跑两遍，r1/r2 目录，results.json 必须**逐字节一致**；manifest 只允许 run_tag/total_wall_s 差异）。
- 只读输入；不写 configs/；不动 gate；不调求解器。
- 脚本头部 docstring 写明口径来源（冻结配置 SHA）与"不排名不选优不训练"。

## 验收（kimi 执行，你不必做但须知道标准）

- r1/r2 逐字节一致；逐对重导精确复现 v2 聚合；identity 锚（B0_G1_BG）在所有 stratum 的 D 约束通过（D_p80=0）；抽样复算若干 stratum 的约束判定。
- 失败处理：交付前自查上述点；无法复现 v2 聚合时停下来报告差异，不要硬凑。
