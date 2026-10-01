# B1 任务定义冻结记录（2026-10-01）

**冻结产物**：`configs/research/task_definition_contract_v0.1.json`
**SHA-256**：`cd7bc8ed1280ebe0bed514fb9ced867501cf49902c0f600eea17fde286810e97`
**冻结脚本**：`scripts/freeze_task_definition_b1.py`（三遍执行字节一致：cd7bc8ed…）
**被冻结提案**：`docs/research/2026-10-01_contractual_adaptive_isp_proposal.md`（SHA-256 `bc3db0dbadbd03655c954f49f681bdb136c46d71ffce29fbc1442851722eed04`，已内嵌为漂移闸门 EXPECTED_PROPOSAL_SHA256）
**用户签认**：2026-10-01 13:23 "确认，先试试效果"（审阅页 `E:/automation_djh/B1_review_sheet.md` 过目后）

## 冻结语义

仅锁定**任务定义**：MDP 结构（T_max=4、γ=1.0、稀疏终端奖励）、动作空间（目录 v0.2 全 30 项枚举 + STOP、identity 锚 B0_G1_BG 恒在、不可行动作掩码、离散档纪律、无连续参数头）、奖励公式（R = R_contract − λc·cost）、开发/测试族隔离、5 条 Goodhart 护栏、§3.5 新族入场五步、验收标准草案、5 项训练前置。

**不授权任何仿真或训练**：契约自身携带 `training_enabled=false`、`fdtd_execution_enabled=false`、`g4_relieved=false`，并有断言守护。

## 逐键断言覆盖（冻结脚本内）

1. 五份输入证据哈希全部内嵌并核验存在；
2. 提案漂移闸门（EXPECTED_PROPOSAL_SHA256）；
3. 目录结构性事实：candidate_count=30、identity 锚在目录内、目录自身 training/fdtd 开关为 False；
4. 本契约三开关为 False；
5. 权重契约 supersede 关系核验（v0.3 声明 v0.2 背景类逐字不变）；
6. 容差契约 frozen 状态核验；
7. 奖励引用与实载文件一致；
8. 护栏 5 条、入场五步、两项不预留、验收与前置条目数核验；
9. 动作空间纪律旗标核验。

## 偏差说明（一处，已显性登记）

提案 §4 文字写"权重契约 v0.2"；冻结时绑定的是 `reward_weights_contract_v0.3.json`。依据：v0.3 的 supersedes 块自身声明"v0.2 background-class section unchanged; this contract adds the gain class"，故对背景类与提案文字等价，且额外覆盖增益类；该绑定理由已写入契约 `reward.weights_contract_note`。

## 纪律声明

本次冻结未新建/修改任何其他契约；未运行仿真或训练；G4 维持未解除；测试族未接触。训练启动仍需契约内 5 项前置全部满足（A1 G4 签认为总开关）。
