# UAV-GPR SFCW 自动处理研究

仓库：<https://github.com/CYberkra/automation>。换电脑时先按[版本管理与跨电脑接续](docs/WORKFLOW.md)克隆、安装依赖并运行检查，再读[最新进度](docs/research/continuous_research_status.md)。当前研究接续点是[评价校准准备](docs/research/2026-09-24_cycle05_handoff.md)。

最小验证仅需 Python 3.10 与 `requirements.txt` 中的 NumPy：

```powershell
python -m pip install -r requirements.txt
python scripts/verify_workspace.py
```

当前预期通过 64 项构造数组检查；不需要实测资料、gprMax 或训练环境。新结果写入 `artifacts/local_checks/`。原始资料、派生审计缓存和公开论文缓存不随 Git 上传，因此下文部分历史资料链接仅在保留这些缓存的电脑可用。

当前目标：依托现有非显性滑坡项目，先做场景特化的 B-scan 自动处理，自动选择背景抑制与增益的算法、参数和必要顺序，保留项目相关地下结构并减少人工调参。泛用性研究延期，首版不建设通用/特化双模式。

本项目已完成资料审计、首版范围约定及前四轮研究。2026-09-24 新增结构证据卡、材料与处理链复核、初步实验计划，以及 10 项代数/评价反例检查；尚未训练模型或运行正演。前轮综述完成 8 篇论文正文通读，本轮新增论文按各自台账记录阅读范围；SGC 等全文访问缺口仍保留。

2026-09-24 第二次续研进一步形成低频材料与源契约、M00–M04 案例提案及参考事件分级规则，并完成作者公开数据的 15 个 S2P 文件审计。公开谱尚未反演为材料参数，案例尚未执行。

第三次续研已将 27 个初始配置写为算子契约和小型数组实现，补充传统自适应基线依据，完成 23 项操作/参考反例检查。结果仅属数值机制证据；尚无 GPR 性能、配置学习或实测验证。

按用户最新要求，已进一步形成评价与标签 v0.2：分事件保真、幅度与干扰评价，增益参考条件，以及可行性/偏序/集合标签；26 项评价反例检查通过。指标和标签结构已明确，物理阈值及完整电磁参考仍待小规模校准，尚未生成正式训练集。

最新数据安排：**先用 gprMax 开展仿真开发与评价，现有测线留待后期使用。** 当前不使用测线训练、调参或提取杂波做混合训练。先定义参考真值、指标和传统基线，再决定选择器；仿真通过不等于实测适用性已验证。

用户已确定：**首版只研究“背景与相干杂波抑制”和“增益与衰减补偿”**，从简单方法逐步扩展。基线/仪器校正、频带滤波、航空/空间几何校正、去噪、子波整形/反卷积和成像延期，其余分类不纳入项目范围。

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
