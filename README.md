# UAV-GPR SFCW 自动处理研究

[当前状态与执行边界](START_HERE.md) · [跨机接续](docs/WORKFLOW.md) · [工作约定](AGENTS.md)

项目目标是为浅层非显性滑坡研究背景/相干杂波抑制和增益/衰减补偿。当前先核对 gprMax 仿真与 SFCW 采集处理的一致性；已有二维仿真与机制控制，不代表现场材料、真实天线或实测性能已验证。

设备频段为 20–170 MHz、步长 0.3 MHz，均匀含端点时共 501 点。当前已完成批次使用 gprMax 4.0.1；本轮按用户授权清理流程、修复 CPU 后处理，新方案或求解仍须显性批准。新输入激励已按用户要求改为官方内置 impulse，见[替换记录](docs/research/2026-10-10_official_impulse_source_preparation.md)；既有结果仍来自 Ricker。

当前仿真/SFCW 路径、官方依据和必要限制集中在[流程审查](docs/research/2026-10-10_simulation_sfcw_simplification.md)。原生资料保持只读；历史结果按原源码和哈希复现，不改写为新证据。

纯数组检查入口（不运行求解器、不读取实测或训练）：

```powershell
python -m pip install -r requirements.txt
python scripts/verify_workspace.py
```

它与 gprMax 4.0.1 环境分开；通过数以检查实际输出为准，不能当作电磁或实测验证。SFCW 修复的定向检查使用已安装官方库的环境运行 `scripts/check_line9_cover_sfcw.py`。

研究依据与历史：

- [评价与标签 v0.2](docs/research/2026-09-24_evaluation_and_labels_v0.2.md)
- [算子契约](docs/research/2026-09-24_cycle03_operator_contract.md)
- [首版范围](docs/research/2026-09-23_v1_scope.md)
- [持久决策](docs/research/decision_log.md)与[产物索引](docs/research/research_artifact_registry.json)
- [完整历史入口](START_HERE_history_20261010.md)

原始测线、地质资料、环境和大数组不随 Git 上传；本地历史聊天也不会随仓库迁移。
