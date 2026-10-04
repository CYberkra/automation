# HS4 跨机任务包本机执行记录（2026-10-04，五阶段全部完成）

> 状态：五阶段 20 道全部 COMPLETED 且独立验收 PASS；全部结果为**机制诊断证据**（`COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE`)，不是物理验收、绝对收敛或实测泛化结论。

## 背景与执行机

[跨机任务书](2026-10-04_hs4_cross_pc_task.md) 原定由另一台电脑执行；该机 RTX 3060/6 GiB 不满足三维 8 GiB 显存门槛，用户决定移交本机。本机：RTX 4090 Laptop(16376 MiB)、RAM 63.2 GB、venv `artifacts/local_checks/gprmax_v4_gpu_env`(Python 3.12.15,gprMax V4.0.0/CUDA float64;319/319 个冻结 .py 与旧契约逐字节一致，25 个 .pyd 为本机新编译，预检如实记录 `historical_binary_changes` + `CROSS_BUILD_REPLAY_REQUIRED`)。设计目录 `artifacts/research_checks/2026-10-04_hs4_cross_pc_design_r1/`(20 份输入 + 12 份已完成细网格控制）原样使用，未修改。

## 阶段1 replay(4 道，2D 中心 1.25 cm)

胶囊 `target_hs4_replay_r1`，结果 `target_hs4_replay_results_r1`。每道 73–80 s，峰值 RSS 1.84–1.86 GiB。

- **四道全部 `native_receiver_bitwise_equal=true`,501 复谱相对 L2 = 0.0**：本机不同编译 .pyd 的求解输出与旧机归档逐比特一致，跨构建复相位担忧闭环，`CROSS_BUILD_REPLAY_REQUIRED` 不转化为实际偏差。
- crest/slope early(160–180 ns)包络 L2 比 = 7.5561，与旧 47 道研究的 7.556 吻合，旧结论在本机运行时下成立。

## 阶段2 centre1cm(4 道，2D 中心 1 cm)

胶囊 `target_hs4_centre1cm_r1`，结果 `target_hs4_centre1cm_results_r1`。每道 156–163 s，峰值 RSS 2.73–2.74 GiB。分析参照为已归档 1/60 m(≈1.67 cm)中心配对 `2026-10-04_hs4_patch_finest`。

中心站位相邻网格差分观测量序列（带符号相对 L2,base,early 160–180 ns):

| 网格步 | 带符号相对 L2 |
|---|---|
| 5 cm → 2.5 cm | 61.291% |
| 2.5 cm → 1.67 cm | 16.596% |
| 1.67 cm → 1 cm（本轮） | **4.102%**(complex 4.292%) |

- 各情形一致：early crest 5.597%/5.563%(signed/complex)、slope 4.405%/6.058%;later(180–220 ns)最大项为 crest complex 7.634%;full(160–220 ns)各情形 4.33–5.60%。每档细化约 4 倍递减，**趋势继续下降但未走平：不宣布收敛，也不把相邻网格差当绝对误差界**。
- crest/slope early 比 = 7.3750(later 0.2424,full 2.1797)，与归档 7.556 一致：该不对称诊断在 1 cm 下依然成立，不是网格伪影。
- 全 501 复谱与归档相对 L2 = 7.29e-4–7.32e-4（四情形），与跨构建 replay 基座兼容。

## 阶段3 3d-centre(3 道，首个三维 inline 机制模型）

中断与续跑（按任务书纪律执行，现场全部保留）：宿主应用退出导致胶囊 `target_hs4_3d-centre_r1` 在 `centre_rough08` 完成（234.8 s,4.97 GiB,raw sha256 已入执行日志）后、`centre_flat` 求解中途被终止，无 FAILED 事件。监督器按设计拒绝续跑已消耗 attempt；依据任务书"另冻只需补做的具体任务，不重跑已完成组"，由新脚本 `scripts/freeze_hs4_cross_pc_continuation.py` 校验中断现场（恰好一道完成且哈希吻合、一道悬挂、一道未启动、无残留求解器进程）后冻结有界续跑契约 `target_hs4_3d-centre_cont_r1`(centre_flat 189.2 s + centre_halfspace 185.6 s，各 4.94 GiB)。随后由同一脚本 `assemble-analysis` 把三道字节核验副本合并为 `target_hs4_3d-centre_merged_r1` 并重新独立 audit(PASS 3 组：版本/网格/dt/float64/有限值/实际收发位置/全材料图），分析在合并胶囊上进行，结果 `target_hs4_3d-centre_results_r1`。

三视图人工核查（`comparison.png`)**通过**:

- 粗糙−平面：178–196 ns 主导单波 let，峰峰 ±6.5e-4，单一事件可分。
- 粗糙−全覆层（=粗糙界面下基覆界面响应本身）：±1.4e-4，弱约 5 倍——起伏界面散射把镜面反射能量打散、深层界面响应大幅削弱的机制在三维 inline 模型中直接可见；量级关系自洽（|rough−flat| ≤ |rough−halfspace| + |flat−halfspace|)。
- 原始：直达/耦合 ~146 ns 峰 5.5e-4(140–160 ns 窗 absmax)，深层区 ~2e-4，无不可解释形态。

三维为 inline 理想偶极机制模型：COL6 实际 X=1.60 m、中心 Tx/Rx=[1.6,5.6,27]/[1.6,6.9,27]、全三维面起伏 1.6 m（剖面 0.8 m 仅为 bin6 切片）、平面参照底界 z=9.05 m;**不认证设备 cross-track 行为或有限天线效应**。

## 阶段4 3d-ends(6 道，左右端站位）

胶囊 `target_hs4_3d-ends_r1`，结果 `target_hs4_3d-ends_results_r1`。每道 184.7–194.2 s，峰值 RSS 4.94 GiB。六视图核查通过：左右端与中心同模式——粗糙−平面主导单波 let 可分（左 ±7e-4 @176–200 ns、右 ±6e-4 @180–200 ns)，粗糙−全覆层弱约一个量级，原始道早期耦合正常。

## 阶段5 3d-xwide-centre(3 道，X 域 12→16 m，内部几何刚性 +2 m)

胶囊 `target_hs4_3d-xwide-centre_r1`，结果 `target_hs4_3d-xwide-centre_results_r1`（参照胶囊为 3d-centre **合并**胶囊；被中断的 r1 无 completed_verification，不能作参照，内容逐字节一致）。每道 250.2–253.4 s，峰值 RSS 6.43 GiB。

- 全 501 复谱：中心站位扩域前后相对 L2 = 8.6913e-4 / 8.6860e-4 / 8.6858e-4(rough08/flat/halfspace)，三模型几乎相同——边界贡献不依赖界面形态。
- 固定窗边界贡献（`window_boundary_contribution.json`,`scripts/compare_hs4_xwide_windows.py` 入库生成）:
  - **配对差分诊断量**（两侧共享边界几何，边界效应相消）:rough−flat early 0.248%/later 1.762%/full 1.639%（带符号）;rough−halfspace 0.712%/2.140%/2.004%。**配对差分对 12 m 域宽稳健**（≤2.2%，远小于研究的 5 倍级形态效应）。
  - **原始响应**:early 85.5%/later 24.0%/full 64.7%——原始深窗波形 (~1e-4 量级） 对域宽高度敏感，边界/PML 晚到贡献与深部信号同量级。**结论：深窗解释只能基于配对差分，不能基于原始波形绝对形态**；这也解释了此前在原始深窗上反复看到的不稳定形态。

## 结论与边界

1. 跨构建逐比特复现成立（replay L2=0.0)，后续全部阶段建立在该基座上。
2. 2D 中心网格收敛序列 61.291%→16.596%→4.102% 继续下降未走平；crest/slope 不对称稳定。
3. 三维 inline 机制：界面起伏使基覆界面响应削弱约 5 倍，端点与中心一致；该机制对域宽稳健（差分 ≤2.2%)。
4. 未认证：绝对收敛、唯一物理归因、cross-track/有限天线、实测泛化。20 道完成不自动等于问题完全解决；后续只能沿证据缺口另冻有界任务。

## 复现入口

```bash
# 环境变量(见任务书);设计目录 D=artifacts/research_checks/2026-10-04_hs4_cross_pc_design_r1
cmd //c "scripts\\run_hs4_cross_pc.cmd freeze --design $D --stage <stage> --replay-report artifacts\\research_checks\\target_hs4_replay_results_r1\\summary.json --out <capsule>"
cmd //c "scripts\\run_hs4_cross_pc.cmd run --out <capsule> --execute"
"$HS4_PYTHON" scripts/analyze_hs4_cross_pc.py --capsule <capsule> --out <results>
# 中断续跑与合并:
"$HS4_PYTHON" scripts/freeze_hs4_cross_pc_continuation.py freeze-continuation --prior <r1> --out <cont>
"$HS4_PYTHON" scripts/freeze_hs4_cross_pc_continuation.py assemble-analysis --prior <r1> --continuation <cont> --out <merged>
# 扩域窗口贡献:
"$HS4_PYTHON" scripts/compare_hs4_xwide_windows.py
```

原始 H5/输入/契约/验收/结果图入 Git;vtkhdf 几何导出（2D 约 30 MB、3D 约 146 MB 每份）与 cuda_cache.jsonl 留本机，SHA-256 台账在各胶囊 `completed_verification.json`(.gitignore 逐项记录）。108 m 大域 `ready_r2` 仍未消耗，可在本机另冻契约接续。G4/训练/保留实测/首版范围不变。
