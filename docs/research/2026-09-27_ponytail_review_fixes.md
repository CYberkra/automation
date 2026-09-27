# 最新审查问题修复

对应 `2026-09-27_ponytail_latest_review.md` 四项发现。保留其他agent的S4 dev/test扩展及规格草稿，只在其基础上修复；不改冻结物理/阈值配置、原始仿真或历史输出。

1. A0报告§6纠正指标等同推论，并同步v1.1草稿及起草任务书：幅度谱相关不证明实波形NRMSE、N_b或算子排序跨维稳定。原数字不变，草稿仍为草稿。
2. 能力导出在聚合前按冻结元数据复用阶梯的事件展开/损伤实例定义，核对完整record_id集合与族、几何、案例、候选、角色、损伤参数；重复/缺失报错，identity必须可运行且D=0。当前元数据推导开发34086条、测试12768条；不可用候选必须保留记录，不能静默消失。
3. 统一跨版本阶梯验收：源码身份始终真实保存；比较时仅排除resource与runner_script_sha256，manifest仅排除run_tag/total_wall_s。其他数值和输入溯源严格比较；同版本r1/r2仍比较runner身份。现有checker原来只打印记录相等与identity检查，现已改成失败即退出。支持`--split test`及左右8块路径模板，不创建新的验收框架。
4. START_HERE顶部只有一份有效当前状态及下一步；旧流水声明历史，底部接续步骤也指向当前S4/S5，不再导向已完成A0或ZFINE2。

## 核验

- `scripts/check_ponytail_review_fixes.py`：synthetic测试族完整集合通过；空输入、漏行、重复、漏族、漏identity、错族、identity不可用七类拒绝；跨版本只放行runner身份/资源差异，指标或输入hash变化仍失败；共享provenance字典也不会被原地修改。未读取测试族波形或计算其评价指标。
- `scripts/check_damage_ladder_acceptance.py`：旧开发r1/r2各34086条记录和manifest一致、覆盖及identity检查通过，合并SHA保持 `6bc13f63c5985fa0b38b4724c2ffa0cf61716053751746520766ef6fa0e4d0d9`。
- 开发组能力导出重跑SHA仍为 `2091f4171c16341731f57bb94c8fead0a72e47a1d58389caa74926286bf8e6c8`。
- 开发阶梯首个工作项：旧提交d95727d代码与新代码均执行chunk0/n_chunks598，57条结果及manifest在既定排除后相同。只是小范围烟雾检查，不冒充新版本全部8块数值回归；正式S4验收仍需按任务书执行该回归。
- 不新增FDTD、不重跑测试族阶梯、不应用测试阈值、不生成排名或训练标签。G4未解除。

## 接续命令

同版本新数据：`python scripts/check_damage_ladder_acceptance.py --split test --left-pattern 'artifacts/research_checks/2026-09-27_damage_ladder_test_{run}_c{k:02d}' --right-pattern 'artifacts/research_checks/2026-09-27_damage_ladder_test_{run}_c{k:02d}'`。

跨版本开发回归：同一checker加 `--cross-version`，左模板指向旧开发r1，右模板指向新开发输出；可以用固定`r1`字面量或`{run}`（左取r1、右取r2）。保留源输出，不改写旧哈希。只要缺少一块或记录不齐，验收失败；跨版本通过不能代替新版本r1/r2确定性证明。

最终验证：新增检查、原有182项纯数组回归、handoff检查通过。源码SHA与检查摘要归档于 `artifacts/research_checks/2026-09-27_ponytail_review_fixes/verification.json`。
