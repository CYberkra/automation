# B2 试点批冻结记录（batch2d_b2_pilot_co，gate 未批准）

日期：2026-10-01
依据：B2 采样设计草案（§6 五点用户"都认"）+ 前置任务台账 P1-1（用户"可以，开始"）
冻结脚本：`scripts/freeze_batch2d_b2_pilot_co.py`（`--verify-only` 字节级复现已通过；全仓 verify_workspace 182/182 PASS）
Gate 状态：**`approved_to_simulate: false`**，396 个 run_id 全部在 `approved_run_ids_pending_agreement`，等待用户签认；执行前还须先重启机器修复 NVML（P0-1）。

## 1. 采样表（12 母模型 × 33 道 = 396 次 FDTD）

### A 层（奖励可求值：容差 v0.2 + 权重 v0.3 直接适用）

| # | 母模型 | 族/形态 | 场景 | 来源 | 备注 |
|---|---|---|---|---|---|
| 1 | B2D-C1mX-BG | C1 平地 | BG | v1 `B2D-C1m-BG` 色散替换 | 覆土 1 m（z 29–30） |
| 2 | B2D-C1mX-NC | C1 平地 | NC | v1 `B2D-C1m-D10m-NC` | 零衬度盒 = 岩性参数 |
| 3 | B2D-C1mX-D10m-W4m-T0.5m-E20-S0.02 | C1 平地 | TGT | v1 同名 | 锚点目标盒固定 z（顶 20.25 m） |
| 4 | B2D-C3mX-NC | C3 平地 | NC | v1 `B2D-C3m-D10m-NC` | BG 用 09-28 归档，不重跑 |
| 5 | B2D-C3mX-D10m-W4m-T0.5m-E20-S0.02 | C3 平地 | TGT | v1 同名 | 同上 |
| 6 | B2D-C3mS2X-NC | C3 T2坡·无TZ | NC | t2 `B2D-C3mS2-NC` 阶梯重建 | 阶梯与归档 BG 逐字节一致（已断言） |
| 7 | B2D-C3mS2X-D10m-W4m-T0.5m-E20-S0.02 | C3 T2坡·无TZ | TGT | t2 同名 | 同上 |
| 8 | B2D-C3mS2TZX-NC | C3 T2坡·有TZ | NC | t2 `B2D-C3mS2TZ-NC` | TGT 留主批 |

A 层复用存量（不重跑）：`B2D-C3mX-BG`、`B2D-C3mS2X-BG`、`B2D-C3mS2TZX-BG`（各 33 道，2026-09-28 t3_co 归档）。

### B 层（新族探针：机制/保真检查，冻结奖励暂不适用）

| # | 母模型 | 覆土（锚点） | 倾角 | 域宽 | 出露 / PML 余量 |
|---|---|---|---|---|---|
| 9 | B2D-C1p5mS1X-BG | 1.5 m | T1（θ_eff 5.7145°） | 50 m | 出露 y=30.75，余量 18.25 m |
| 10 | B2D-C1p5mS3X-BG | 1.5 m | T3（θ_eff 21.8014°） | 40 m | 出露 y=19.75，余量 19.25 m |
| 11 | B2D-C2mS3X-BG | 2.0 m | T3 | 40 m | 出露 y=21.0，余量 18.0 m |
| 12 | B2D-C2mS3TZX-BG | 2.0 m | T3·有TZ | 40 m | 同上（TZ 探针） |

**zc 规则披露**：A 层坡地母模型保持归档 t3 几何（zc=25.8，遗留 t2 规则，锚点覆土 4.2 m），保证与归档 BG 逐字节配对；B 层新档用 `zc = 30 − h`（锚点覆土=标称值）。两套规则不嵌套；B 层本来就是新族。覆土 2.5 m 档与 T1×h≥2.0 档在既有域宽下不满足出露/PML 余量 ≥2 m，留主批决策。

## 2. 几何与静态锚（static_check.json 全录）

- 阶梯重建与归档交叉验证：S2X/S2TZX 重建阶梯 box 行与归档 t3 BG 母模型**逐字节一致**（148 条带 ×2，已断言）。
- θ_eff 最小二乘拟合与名义档差 <0.05°；全部档出露位置、PML 余量、Fermat 静止点余量断言通过。
- 材料逐字对照 `dispersion_materials_v0.1.json`（色散 cover/tzone，基岩非色散），hash 锁定。
- 采集与 t3/t1t3 逐字一致：CO33，33 道 @0.25 m，Rx=Tx+1.30 m，z=45，仅 Ex，锚点 t17。
- NC−BG 配对差分逐位恒零为**跑后必验项**（新几何下重新验证，不假设继承）。

## 3. 预算与硬停止

| 域 | 母模型数 | 单道估 | 合计 |
|---|---|---|---|
| 32 m 平地 | 5 | 18–24 s | — |
| 40 m 坡地 | 6 | 24–29 s | — |
| 50 m S1 | 1 | 30–36 s | — |
| **总计** | **12（396 道）** | — | **估 2.4–3.6 h，保守上限 300 min** |

硬停止：单道 20 min、Job 6 GiB、输出 2 GiB、retries=0、失败即停、串行。执行序：平地 C1 三连 → C3 对 → S2X 对 → S2TZX NC → B 层 S3 探针 → S1 探针最后。

## 4. 纪律边界（照录草案 §7 并补充）

- 本冻结**不构成执行授权**：gate `approved_to_simulate:false`；执行须用户显式签认。
- GPU 备注（2026-10-01 重启后复验更新）：NVML/nvidia-smi 仍失败但 CUDA 计算路径实测完好（pycuda 设备在线、显存空闲 14.72/16 GiB ≥ 预检 4 GiB）；执行链无 NVML 依赖，**不构成阻塞**；NVML 修复降级为日后非阻塞项。
- 测试族 {C5,C8} 零接触；不用新数据回头调整已冻结结果；B 层数据在各自参考窗构造+容差重校准完成前不应用冻结奖励契约。
- 不生成物理/训练标签、无可探测性结论、不做算子选型；G4 维持未解除；`reference_state` 不升级。
- 每组按 group_id（母模型）划分，倾斜变体间不共享分组。

## 5. 产物清单

- 母模型：`configs/research/batch2d_b2_pilot_mothers/`（12 个 .in）
- 输入与契约：`configs/research/batch2d_b2_pilot_co/`（396 个 .in + cases.json + groups.json + budget.json + static_check.json）
- Gate：`configs/research/gprmax_v4_execution_gate.json`（batch2d_b2_pilot_co，未批准）
- 执行器（待命）：`scripts/run_approved_batch2d_b2_pilot_co.py`、`scripts/run_batch2d_b2_pilot_co_gpu.cmd`
- launcher_sha256：`27911e26696db3e928d8a685e7cedf2a32fbad9411623e7c9e84a901fefb8e20`
