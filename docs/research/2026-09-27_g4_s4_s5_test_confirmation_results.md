# G4 S4 与 S5 测试族确认结果

日期：2026-09-27。S4 全量数据生成和 S5 唯一一次测试族确认评价已完成。所有步骤只读既有仿真并使用 CPU 后处理；未启动 FDTD、训练模型或解除 G4。

## S4：阶梯回归与测试族数据

- 新版开发阶梯 r1/r2 各 8 块、每遍 34,086 条记录。真实 runner SHA-256 为 `a178c983f2b67644f68ed3fceccf1c8633ca05ca78ae754abebfbb366864930d`。同版本验收保留 runner 身份，去除既定 resource 后数值记录规范化 SHA 为 `d84ae98f289dc3da4bb4127cb920deee1756252542dcd29be7e2c75ce32b2622`，覆盖、manifest 和 identity 检查通过。主控独立完成旧/新开发阶梯跨版本回归：34,086 条记录通过，SHA `f041d5788c2aaf416e2a1e37f38265773548a379990d3309b0b9253b4fb8b7c5`；只忽略 `resource`、`provenance.runner_script_sha256`，manifest 只忽略 `run_tag`、`total_wall_s`。
- 新版测试阶梯 r1/r2 各 8 块，每遍 12,768 条记录。保留真实源码身份后两遍同版本验收通过，SHA `cc9b29a5ef162f6865477721f7356f19f60069cda26917535e4cb93bfa4e8227`。完整性、manifest 和 identity 检查通过；族为 C5/C8，几何只有 MT，CO 无对应测试族案例。
- 开发能力导出重跑 2,065 条能力行和 324 条负控行，逐字节复现冻结 SHA `2091f4171c16341731f57bb94c8fead0a72e47a1d58389caa74926286bf8e6c8`。
- 测试能力导出实际分别读取测试阶梯 r1 与 r2。结果各含 1,482 条能力行、228 条负控行，57 个候选，只有 MT；r1/r2 结果逐字节一致，SHA `8c5bd82392cd0f457c678ec8bb8b9e66ded466c8af2f241f73b9d65ceeee3efe`。每遍 sidecar manifest 核实了实际 run 标签、8 个输入块的 records/manifest 路径与 SHA、导出器源码 SHA、命令和结果 SHA。

S4 产物位于 `artifacts/research_checks/2026-09-27_damage_ladder_{dev,test}_{r1,r2}_c00..07/` 与 `artifacts/research_checks/2026-09-27_g4_mission_capability_dev/`、`artifacts/research_checks/2026-09-27_g4_mission_capability_test_r1/`、`artifacts/research_checks/2026-09-27_g4_mission_capability_test_r2/`。测试能力导出验收脚本为 `scripts/check_g4_s4_test_split_acceptance.py`。

## S5：一次性确认评价

执行前冻结输入清单 SHA-256 为 `a265580988c2d982d922df6305f38a6d68540c7a3e8fd550a8a8419440e104a4`，覆盖 53 项源码、配置、开发参照、测试 r1/r2 能力结果与 16 个阶梯块。可变的 `decision_log.md` 在执行前逐字节快照至 `artifacts/research_checks/2026-09-27_g4_s5_test_confirmation_freeze/pre_execution_decision_log.md`，SHA 为 `e39d726c9630efeeee49e55a34482a4066483e9630e9c36e9af2886ff0702ffc`；清单记录了原始来源。执行器 SHA 为 `2d14cc8e67f19be19042d55acb6160df9fae4aeedcbd8df6a2770df0112f84aa`。

冻结执行命令（历史留痕，已消耗attempt，禁止重新执行；复核请用下方只读验收器）：

```powershell
artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe scripts/run_g4_test_confirmation.py --output-dir artifacts/research_checks/2026-09-27_g4_s5_test_confirmation_r1
```

命令 exit 0，且排他 attempt 标记只创建一次。只读验收命令 `artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe scripts/check_g4_test_confirmation_acceptance.py --output-dir artifacts/research_checks/2026-09-27_g4_s5_test_confirmation_r1` exit 0。

结果 SHA-256 为 `157c3d6b1e3ba5a61422e8eb17daae167aa9a50861f4b706c76a4d0d3b15bec3`；attempt 标记 SHA 为 `33ca389ac50305b29c32a6cb08649f472939d29ceb10b503b9f32aa1e2f29405`；运行 manifest SHA 为 `59bade2a3f3da097d7cb169e693e03a60143ac09baa257ba7811452a147da15b`。每个测试族组合都核算了 57 个候选的 9 个任务档 D 单元及 q=2/q=4 两个 N 单元，共 1,254 个候选×约束单元。任务档 ε 对 C5/C8 和四种损伤均取同 geometry、同 damage type 的开发 C1/C3 中较小值 `0.03003563001306807`（C1=`0.03003563001306807`，C3=`0.034122773235149714`）；这是开发标定值的保守沿用，不是测试族重标定。稳定性对每个完整 `(quantity, damage_type, mission_class, geometry, family)` 键应用冻结规则，共 14 个 D strata；弱事件档只进入稳定性检查，不进入 a80 能力头条。N_b 保留并列报告，因跨尺度 caveat 不应用稳定性规则。

| 测试组合 | D 通过候选数 | N 通过候选数 | 联合 admissible 数 |
|---|---:|---:|---:|
| MT / C5 | 1 / 57 | 43 / 57 | 1 / 57 |
| MT / C8 | 1 / 57 | 43 / 57 | 1 / 57 |

两组合唯一 admissible 都是 identity 锚 `B0_G1_BG`。14 个 D strata 均报告 `unique_stable_top`，且均只有一个 admissible 候选，故标记为 singleton 空真；这只记录冻结 pairwise 规则在空比较集上的结构，不证明某算法获胜或优越。8 个 a80 输出块覆盖 57 个候选，统计能力项与网格链系统偏差项分开记录；未将其换算成深度。

G4 仍未解除，冻结配置物理阈值仍为 null；没有训练资格或实测性能主张。测试结果是构造参考下的仿真域确认，不代表物理真值或现场泛化。

## 复现与留痕

逐块命令、日志、最终 S5 命令和验收结果位于 `artifacts/research_checks/2026-09-27_g4_s4_execution_logs/`。较长阶段会话为开发 r1 `26663`、开发 r2 `84675`、测试 r1 `34846`、测试 r2 `95536`；它们均以原会话完成。正式 S5 命令在工具返回前完成，没有分配长时会话 ID。首次非正式合成自检因导入常量路径错误而失败，错误记录与修正后通过的最终自检均保留在日志目录；该失败发生在 S5 冻结执行之前，不是正式测试阈值尝试。

主要文件：`scripts/run_g4_test_confirmation.py`、`scripts/check_g4_test_confirmation_acceptance.py`、`configs/research/g4_s5_test_confirmation_inputs_v1.0.json`、`artifacts/research_checks/2026-09-27_g4_s5_test_confirmation_r1/results.json`、`artifacts/research_checks/2026-09-27_g4_s5_test_confirmation_attempt.json`。
