# UAV-GPR SFCW 自动处理研究

仓库：<https://github.com/CYberkra/automation>。**换电脑或新会话先读 [START_HERE.md](START_HERE.md)**，再按[接续说明](docs/WORKFLOW.md)准备环境。当前已完成构造数据损伤试验与[V4 资料核对](docs/research/2026-09-24_gprmax_v4_review.md)；旧 cycle05 的损伤研究已完成首版，不要重做。

设备：**SFCW 20–170MHz，步长0.3MHz**，均匀且含端点时为501频点。地表覆盖层可以是粉质粘土，基岩一般为砂岩；具体电性与厚度待定。项目为浅层非显性滑坡，深度暂按地下约20m以内理解。任何gprMax仿真先与用户敲定，当前未批准。事实、假设与未知见[项目背景](configs/research/project_context_v1.json)，执行边界见[记录](configs/research/gprmax_v4_execution_gate.json)。

最小验证仅需 Python 3.10 与 `requirements.txt` 中的 NumPy：

```powershell
python -m pip install -r requirements.txt
python scripts/verify_workspace.py
```

当前预期通过 **136 项**构造数组检查（历史64＋损伤72）；不需要实测资料、gprMax、绘图库或训练环境。新结果写入 `artifacts/local_checks/`。原始资料、派生审计缓存和公开论文缓存不随 Git 上传，历史资料链接可能只在本机有效。V4 求解器需另用 Python3.11–3.13 环境，不能直接套用上述数组环境。

当前目标：依托现有非显性滑坡项目，先做场景特化的 B-scan 自动处理，自动选择背景抑制与增益的算法、参数和必要顺序，保留项目相关地下结构并减少人工调参。泛用性研究延期，首版不建设通用/特化双模式。

本项目已完成资料审计、首版范围约定及前四轮研究。2026-09-24 新增结构证据卡、材料与处理链复核、初步实验计划，以及 10 项代数/评价反例检查；尚未训练模型或运行正演。前轮综述完成 8 篇论文正文通读，本轮新增论文按各自台账记录阅读范围；SGC 等全文访问缺口仍保留。

2026-09-24 第二次续研进一步形成低频材料与源契约、M00–M04 案例提案及参考事件分级规则，并完成作者公开数据的 15 个 S2P 文件审计。公开谱尚未反演为材料参数，案例尚未执行。

第三次续研已将 27 个初始配置写为算子契约和小型数组实现，补充传统自适应基线依据，完成 23 项操作/参考反例检查。结果仅属数值机制证据；尚无 GPR 性能、配置学习或实测验证。

按用户最新要求，已进一步形成评价与标签 v0.2：分事件保真、幅度与干扰评价，增益参考条件，以及可行性/偏序/集合标签；26 项评价反例检查通过。指标和标签结构已明确，物理阈值及完整电磁参考仍待小规模校准，尚未生成正式训练集。

最新数据安排：**先用 gprMax 开展仿真开发与评价，现有测线留待后期使用。** 当前不使用测线训练、调参或提取杂波做混合训练。先定义参考真值、指标和传统基线，再决定选择器；仿真通过不等于实测适用性已验证。

用户已确定：**首版只研究“背景与相干杂波抑制”和“增益与衰减补偿”**，从简单方法逐步扩展。基线/仪器校正、频带滤波、航空/空间几何校正、去噪、子波整形/反卷积和成像延期，其余分类不纳入项目范围。

- [构造数据首轮实验：损伤、参考误差与候选能力](docs/research/2026-09-24_damage_pilot_findings.md) · [运行前计划](docs/research/2026-09-24_damage_pilot_plan.md) · [72 项检查与完整数据](artifacts/research_checks/2026-09-24_damage_pilot_r2/results.json)
- [项目工作约定](AGENTS.md)
- [当前评价与标签协议 v0.2](docs/research/2026-09-24_evaluation_and_labels_v0.2.md) · [机器契约](configs/research/evaluation_label_contract_v0.2.json) · [26 项检查结果](artifacts/research_checks/2026-09-24_cycle04_evaluation_labels_r2.json)
- [27 个配置的算子契约与参考反例（第三次续研）](docs/research/2026-09-24_cycle03_operator_contract.md) · [JSON 清单](configs/research/operator_catalogue_v0.1.json)
- [传统自适应基线与自动选秩的边界](docs/research/2026-09-24_cycle03_traditional_baselines.md) · [来源台账](docs/research/2026-09-24_cycle03_source_ledger.json)
- [23 项数组检查结果](artifacts/research_checks/2026-09-24_cycle03_operator_checks.json)
- [加权评分的局限与阶段状态（研究备忘，相关用户提议已撤回）](docs/research/2026-09-24_manual_weights_to_selector_plan.md)
- [低频材料、激励与观测量契约（第二次续研）](docs/research/2026-09-24_cycle02_materials_and_source.md)
- [P1 案例与参考事件契约](docs/research/2026-09-24_cycle02_pilot_reference_contract.md) · [JSON 提案](configs/research/p1_case_contract_v0.1.json)
- [第二次续研来源台账](docs/research/2026-09-24_cycle02_source_ledger.json)
- [结构证据卡、低频材料与处理链严谨性复核（2026-09-24）](docs/research/2026-09-24_structure_and_chain_rigor.md)
- [初步实验计划：阶段、候选配置、评价与预算](docs/research/2026-09-24_initial_experiment_plan.md)
- [本轮来源及阅读台账](docs/research/2026-09-24_rigor_source_ledger.json) · [代数检查结果](artifacts/research_checks/2026-09-24_processing_algebra.json)
- [持续研究进度与下一步](docs/research/continuous_research_status.md) · [任务书](docs/research/continuous_research_brief.md)
- [非显性滑坡与 gprMax 自动处理综述（Word）](docs/research/2026-09-23_landslide_literature_review.docx) · [PDF](docs/research/2026-09-23_landslide_literature_review.pdf) · [Markdown](docs/research/2026-09-23_landslide_literature_review.md)
- [文献阅读范围与来源校验记录](docs/research/2026-09-23_literature_reading_ledger.json)
- [首版范围决定与最小方案](docs/research/2026-09-23_v1_scope.md)
- [gprMax 仿真优先：相近文献、指标与测线留用安排](docs/research/2026-09-23_gprmax_simulation_first_evidence.md)
- [两类处理的跨领域 ML 架构与文献证据](docs/research/2026-09-23_cross_domain_ml_evidence.md)
- [背景抑制、增益与自动配置的评价协议草案](docs/research/2026-09-23_v1_evaluation_protocol.md)
- [第一轮研究与推荐路线](docs/research/2026-09-23_uav_sfcw_autoprocessing.md)
- [本地资料审计与待核实问题](docs/research/2026-09-23_material_audit.md)
- [26 份 CSV 审计汇总](artifacts/initial_audit/csv_summary.csv)
- [审计脚本](scripts/audit_materials.py)

原始资料位于 `探地雷达背景资料/`，保持只读。审计快照记录了工作期间曾存在、随后被移出的 8 个旧仿真/图片文件，当前是否存在以审计目录中的 `source_presence_at_completion.json` 和实际目录为准。
