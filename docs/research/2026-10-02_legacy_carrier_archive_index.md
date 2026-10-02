# 旧载频（95 MHz）表示证据归档登记（2026-10-02）

用户指令："把以前错误的版本全部归档"。本表把**所有依赖旧带符号重建表示**（`env·cos(2π·95MHz·t + ∠env)`，载频误用带心、+75 MHz 频移、幅值为官方一半）的证据与脚本**原地归档**。

**归档纪律**：所有文件保留原位、一字节不动——冻结契约（容差 v0.1、参考窗 v0.1、权重 v0.2/v0.3）内嵌了其中多份证据的全文 SHA-256，载频修复的复现性核对也以这些归档值为基准。移动或改写会断历史复现链。**归档的含义是：自本表起，这些产物中的带符号波形数值、绝对 dB 刻度、实测峰登记、移位类指标不得再作为现行依据引用；作为历史记录与复现基准仍然有效。**

依据：[设计审查 §2](2026-10-02_model_design_review.md)、[载频修复 v0.2](2026-10-02_sfcw_carrier_fix_v0_2.md)、[独立复审](2026-10-02_sfcw_carrier_fix_rereview.md)、[复审修复与验收](2026-10-02_carrierfix_acceptance_repairs.md)。

## A 类：直接产出旧表示证据（归档，数值停用）

| 证据 | 来源脚本 | 失效范围 | 后继 |
|---|---|---|---|
| `2026-09-28_t3_damage_ladder_r1/r2.json` | `study_t3_damage_ladder.py` | 全部带符号指标（D/a/A/H/ρ/到时）；**按字段分列**：幅值/极性/删除档 D、a、H、ρ 经 v0.2 复算不变仍可引用；**加噪 Nb 作废**（v0.2 实测六项全变，如 S2TZX +6 dB 档 4.960311→5.056014，独立审查复核 S2X/S2TZX/S3X 两档各变）；**到时诊断作废**；移位档刻度作废 | `…_t3_ladder_carrierfix_v0_2`（开发）；v0_3 待全量重跑 |
| `2026-09-28_t3_operator_effects_r1/r2.json` | `study_t3_operator_effects.py` | 算子 D_e/Nb 绝对刻度；top-1 结构经 B2 间接复算稳定，完整排序以新版为准 | 待 v0_3 链重出 |
| `2026-09-29_t3_operator_effects_rpca_r1/r2.json` | `study_t3_operator_effects_rpca.py` | 同上（RPCA 算子） | 同上 |
| `2026-09-28_tolerance_sensitivity_r1/r2.json` | `study_tolerance_sensitivity.py` | 翻转率结论锚定在旧移位刻度上；τ_D 锚点需按官方表示重新登记（复审：2 样本 D=0.9557/0.9985/1.0106 仍全部挂，**无证据要求改阈值**） | 容差 v0.2 提案（待用户决策） |
| `2026-09-28_reward_weights_r1/r2.json`、`2026-09-29_reward_weights_v0_2_r1/r2.json` | `study_reward_weights_v0_1.py` 等 | 排名演示基于旧刻度；dB 恒等式构造本身不受影响 | 同上 |
| `2026-09-29_gain_ladder/gain_effects/gain_weights/clip_scale_r1/r2.json`（8 份） | `study_gain_*_v0_1.py`、`study_clip_scale_v0_1.py` | 增益/截幅全部带符号证据；**"2× 幅值在比值中抵消"不适用于载频变化，不可沿用**（复审 P2） | 待版本化重出 |
| `2026-09-29_reward_protocol_t3_first_run_r1/r2.json` | `run_reward_protocol_t3_v0_1.py` | 奖励首跑绝对分值 | 待 v0_3 链重出 |
| `2026-09-29_s1s3_ref_window_r1/r2.json` + `…_cache/` | `freeze_s1s3_reference_window_v0_1.py` | **t_measured 登记列、NC 与 floor 的绝对 RMS 均作废**（独立审查复算：S1X NC 均值 0.0004165575→0.0008200593，floor RMS 3.0384→6.0768，官方表示幅值×2 的直接后果）；费马几何窗不受影响；冻结断言（峰落窗内、NC 间隙为正）在官方表示下仍全过 | 契约 v0.2 更新（待用户决策） |
| `2026-10-01_reward_protocol_b2_pilot_r1/r2.json` | `run_reward_protocol_b2_pilot_v0_1.py` | 绝对 R 刻度；top-1 与可行性经 v0.2 复算零翻转仍成立；完整排序以新版为准 | `…_b2_reward_carrierfix_v0_2`（开发） |

## B 类：使用旧表示的增广与诊断产物（归档）

| 产物 | 说明 |
|---|---|
| `augment_sim2real_v0_1` … `v0_6` 全部版本输出 | 仿真侧输入经旧加载器；v0.6 振铃标定维持"经验候选、未验收"，其数值须按官方链重出后再议 |
| `fig_line9_vs_sim*.py` 产出图（含工作区根目录 `fig_line9_vs_sim_early.png`） | 仿真侧波形为旧表示；实测侧不受影响；图仅作历史展示 |
| `preview_b2_pilot_bscans.py`、`preview_bscan_operators.py`、`preview_ringing_bscan.py`、`diag_v04_display.py`、`diag_v04_t0rel.py` 产出图 | 展示性诊断，波形为旧表示 |
| `analyze_direct_fidelity.py`、`diag_rpca_catalogue_equivalence_cost.py` 产出 | 诊断性分析，带符号数值作废 |

## C 类：不受影响的资产（明确免责）

- **一切 `_attempt.json` 执行记录与归档 H5**：求解器原始输出，与后处理载频无关。
- **包络/复频响类证据**：`2026-09-27_c3_co_contamination_root_recompute.json`、`2026-09-27_g4_s6_root_review.json` 等基于复包络/频谱的分析，不经过带符号重建。
- `research_sfcw_unified_adapter.py`（09-27 统一适配器）：**一直使用官方 `real_bandpass`，本来就正确**，其产物不在归档范围。
- **冻结契约本体**（容差 v0.1、参考窗 v0.1、权重 v0.2/v0.3）：作为冻结历史契约保留；其数值锚点的版本升级是独立的用户决策流程，不因本归档自动失效。
- 3D 标杆批（benchmark3d_r2_co）全部冲激响应 H5 与几何、材料输入：原始仿真，无载频问题；其 SFCW 合成图应改用官方链重出（见修复报告 §5.3）。

## D 类：脚本处置

20 个引用旧加载器的脚本（3 个定义 + 17 个导入）**保留原样**（历史复现需要），自本表起视为**冻结历史脚本**：不得用于产出新证据；新研究统一走 `sfcw_official_loader_v0_2.py`（及通过 `sfcw_carrierfix_acceptance.py` 验收的 v0_3 入口）。清单见复审报告 §2 与本文 A/B 表来源列。
