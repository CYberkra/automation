# B2 试点批奖励协议 baseline 判卷（P4-1）

- 日期：2026-10-01
- 前置任务台账条目：P4-1（docs/research/2026-10-01_pretraining_prerequisite_ledger.md）
- 执行脚本：`scripts/run_reward_protocol_b2_pilot_v0_1.py`
- 运行环境：`artifacts/local_checks/gprmax_v4_gpu_env`（与 09-29 首跑相同）
- 输出：`artifacts/research_checks/2026-10-01_reward_protocol_b2_pilot_r1.json`（r2 同字节，断言通过；r1 SHA256 前 16 位 `ffb57aa9c3a4c673`）

## 1. 契约基线（SHA 断言通过）

| 契约 | 文件 | SHA256 |
|---|---|---|
| 容差 v0.2 | reward_tolerance_contract_v0.2.json | 2b2b6aff6a6bb3710f6d361a61ca710946a52e6d1cfeb4dc8828a64d4c9a5fee |
| 权重 v0.3 | reward_weights_contract_v0.3.json | 282d083b34477bafff25473f853e3e8d47c6c44807c82e3c5f9d41170399e0d0 |
| 背景公式来源 v0.2 | reward_weights_contract_v0.2.json | 0427b6ed58328e5ae01d2e0bfbda90460ef1a68beea318808a3f2ae2ed69bfb0 |

容差数值：tau_A=0.2、tau_D=0.95、eps_Nb=1.5（arrival_drift 仅诊断）。背景类公式
`R_bg = contrast_db_delta - 20*log10(1/(1-D_e)) - max(10*log10(Nb_ratio), 0)` 来自 v0.2；
v0.3 的 supersedes 字段声明背景类一节不变，脚本对此做了断言。

## 2. 数据

| 族 | 母模型 | 归档日期 | 几何 | 界面 z |
|---|---|---|---|---|
| C1mX | B2D-C1mX-BG | 2026-10-01（B2 试点批） | 平地，覆土 1 m | 29.0 |
| S2X | B2D-C3mS2X-BG | 2026-09-28 | 坡地 | min(0.2y+22.625, 30) |
| S2TZX | B2D-C3mS2TZX-BG | 2026-09-28 | 坡地 | 同上 |
| C3mX | B2D-C3mX-BG | 2026-09-28 | 平地，覆土 3 m | 27.0 |

候选：目录 v0.2 背景类 8 算子（B0/B2/B3/B4/B5/B7/B8/B9_G1_BG）+ 目录外 local_mean_w11 参照。
B 层 4 个探针母模型（C1p5mS1X/C1p5mS3X/C2mS3X/C2mS3TZX）不判卷：按任务定义契约 scope_limits，
新场景须自建参考窗并重校准容差后方可使用冻结奖励。

## 3. 复现性核对

S2X / S2TZX / C3mX 三族的 selection、R_bg、ranking 及**每一行指标**与
2026-09-29 首跑（2026-09-29_reward_protocol_t3_first_run_r1.json）逐项相等。
契约升级（容差 v0.1→v0.2、权重 v0.2→v0.3）对背景类判卷行为零影响。

## 4. baseline 判卷结果

| 族 | 赢家 | R_bg (dB) | D_e | A | Δ对比度 (dB) | Nb_ratio | 被门排除 |
|---|---|---|---|---|---|---|---|
| **C1mX** | **B7_G1_BG** | **52.378** | 0.0015 | 0.0013 | +52.391 | ≈0 | B2(A)、B3(AD)、B4(AD)、B8(A)、B9(AD)、local_mean(AD)；B5 不可用 |
| S2X | B9_G1_BG | 27.017 | 0.0156 | 0.0002 | +27.154 | 0.0019 | 无（9 个全部可行） |
| S2TZX | B9_G1_BG | 25.311 | 0.0331 | 0.0006 | +25.603 | 0.0027 | 无 |
| C3mX | B7_G1_BG | 30.025 | 0.0141 | 0.0116 | +30.148 | 0.0010 | 同 C1mX 的排除模式 |

可行集完整排序（R_bg, dB）：

- C1mX：B7(52.38) > B0(0.00)；其余全部不可行
- S2X：B9(27.02) > local_mean(22.49) > B8(21.96) > B4≈B3(21.74) > B5(18.36) > B7(18.16) > B2(5.70) > B0(0)
- S2TZX：B9(25.31) > local_mean(22.03) > B4≈B3(21.69) > B8(21.63) > B5(18.62) > B7(17.12) > B2(5.66) > B0(0)
- C3mX：B7(30.02) > B0(0.00)；其余全部不可行

## 5. 观察

1. **赢家按几何分群**：平地族（C3mX、C1mX）选 B7，坡地族（S2X、S2TZX）选 B9。
   C1mX 作为全新平地族复现了这一模式——目录内不存在"通吃"算子，奖励协议能区分几何需求。
2. **C1mX 上 B9 直接毁掉事件**（D_e=1.0、A=1.0，对比度仅 +0.17 dB），被硬门排除；
   而 B9 在坡地族是赢家。平地上同相轴水平，B9 的道间参考机制把事件也消掉了。
3. **B8 在 C1mX 上对比度 +94.3 dB 但 A=0.44 超门**：增益类行为的幅度失真被保真门拦住，
   说明硬门确实在挡"只追对比度"的投机解。
4. C1mX 的 B7 对比度增益（+52.4 dB）明显高于 C3mX（+30.1 dB）：1 m 覆土下直达波/浅层
   杂波与事件的功率比更大，背景抑制收益更高。
5. B5 在两个平地族均报 `svd_cutoff_gap_unresolved`（SVD 截断间隙未解析），坡地族可用——
   平地数据奇异值谱更集中，该算子的自适应截断找不到间隙。维持不可用标记，不干预。

## 6. 注意事项

- C1mX 的分数是在 t3 冻结容差下的**诊断性 baseline**；容差标定于 t3 三族，
  C1mX 自身的容差校准检查在台账 P1-2（试点评审）中处理。
- 本运行无求解器执行、无测试族数据、无训练。
