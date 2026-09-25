# 垂向第二次细化（ZFINE2）结果：相邻网格变化降至约4%/3.9°，资格不升级

按[设计书](2026-09-25_deep_zfine2_design.md)执行：保持 dy12.5mm 不变，dz 由 6.25mm 降至 3.125mm，仅运行 DEP_BG_ZFINE2 与 DEP_20_ZFINE2 各一次。两次均 CUDA double、1 次 attempt、exit 0，无重试、无 CPU 回退。40 分钟墙钟与 20GiB Job 内存预算均未触及。

以下区分来源：墙钟与 exit 为监督器实测；"22m57s/22m13s"为求解器自报耗时；Job 峰值内存为监督记录；收敛链数值全部来自分析脚本对已归档 HDF5 输出的实测复算（`analysis_r1/results.json`），不是求解器声明。

## 执行记录

|运行|attempt|exit|墙钟（实测）|求解器自报|Job 峰值提交内存|
|---|---:|---:|---:|---:|---:|
|DEP_BG_ZFINE2|1|0|1438.6 s|22m57s|15.28 GB|
|DEP_20_ZFINE2|1|0|1334.8 s|22m13s|15.28 GB|

设计预算（单组约26.5分钟）与实测同量级。原始输出已核验：网格 1×2560×16000（40,960,000 单元，与 YFINE 总数相同）、dy12.5mm/dz3.125mm、Ex 全有限非零。场景、几何、PML、1200ns 时窗与官方单位 impulse/SFCW 20–170MHz 处理均与前轮一致，未增加其他参数扫描。

## 收敛链（目标减背景，全频带复数）

|相邻比较|全带相对 L2 变化|最大相位变化|最大幅度变化|中位复数误差|
|---|---:|---:|---:|---:|
|原2.5cm均匀→FINE（均匀1.25cm）|66.6647%|64.8431°|0.5940 dB|19.76%|
|FINE→ZFINE（dy12.5/dz6.25mm）|17.2452%|16.2741°|0.1425 dB|4.96%|
|ZFINE→ZFINE2（dy12.5/dz3.125mm）|4.1216%|3.8847°|0.0347 dB|1.18%|

趋势解读：

- 最近两步均为仅细化垂向（dy 保持 12.5mm），因此 17.2452/4.1216≈4.18、16.2741/3.8847≈4.19 的比值不再混入横向细化，是垂向二阶行为的更直接指示：垂向每减半，相邻变化约除以 4，符合二阶 FDTD 预期。此为实测相邻差的比值趋势，未做 Richardson 外推，相邻差也不等于相对精确解的误差。
- ZFINE2 相对最初 2.5cm 网格的复数 L2 差为 85.7514%、最大相位差 85.0018°。该距离增大不表示 ZFINE2 更差：粗格不是精确参考，判据应看相邻变化及其机制，而非靠近粗格。
- ZFINE2 的 200/400ns 尾窗相对 L2 为 0.0258%（前组 0.0259%），远小于网格间变化，尾窗拖尾不是当前主导差异。

## 结论与资格

- 明确结论：在保持 dy12.5mm 的前提下，垂向 6.25mm→3.125mm 引起的相邻变化已降至约 4%（L2）/3.9°（最大相位）量级，且相邻比值支持二阶收敛趋势。这不声明绝对收敛，也不做任意精度无限加密（设计书原则）；进一步加密必须由本次结果明确论证需求。
- 20m 全频带物理真值的 `reference_state` 是否由 `numerically_unresolved` 升级，由参考资格评估决定；本文保持 `numerically_unresolved`、physical_label_eligible=false、training_eligible=false，直到资格文档基于本证据更新。物理阈值仍为 null。

## 限制

- 不拟合源幅度、时间平移或相位补偿；各比较均在同一固定约定下直接求差。
- 尾窗诊断固定为 200/400ns，实波形诊断固定 440–600ns 窗；结论限于这些窗口设定。
- 保持全频带 20–170MHz，未靠裁掉高频获得数值改善；分带差异（120–170MHz 仍最大）如实保留。
- 全部为二维 YZ/TMx 域、沿飞行 X 无限延伸的研究假设模型；结论限定于该平层、4m 宽异常与既定几何，不构成真实 B-scan 或实测泛化证据。

## 归档与复现

原始运行归档：

- `artifacts/research_checks/2026-09-25_DEP_BG_ZFINE2/`（HDF5、`.in`、授权快照、日志、监督终态、record 哈希）
- `artifacts/research_checks/2026-09-25_DEP_20_ZFINE2/`（同上）
- attempt 记录：`2026-09-25_DEP_BG_ZFINE2_attempt.json`、`2026-09-25_DEP_20_ZFINE2_attempt.json`（含输入 SHA256：背景 `.in` 为 `1d21c51c…f259`，目标 `.in` 为 `06d528f2…e11`，完整值见文件）

分析与复算：

- `artifacts/research_checks/2026-09-25_deep_zfine2_analysis_r1/` 与 `_r2/`：两次连跑 `results.json` 逐字节一致、数组一致；输入哈希覆盖全部 8 个 HDF5 及 `2026-09-25_yz_depth_analysis/arrays.npz`（完整 SHA256 列表见 results.json `inputs` 字段；两个 ZFINE2 HDF5 哈希与 record.json 记录一致）。
- 复现入口：`scripts/analyze_deep_controls.py`（CPU 后处理，仅复算比较，不调用求解器；results.json 中 `solver_invoked=false`）。
- 输入与执行契约位于 `configs/research/deep_zfine2_v1/`；attempt 不可复用，旧会话不重启。

## 下一步

- 将 ZFINE2 证据并入批量仿真规格书的 P1–P10 占位项（各占位项在冻结具体输入与资源上限前不执行）。
- 着手 3D 基准链设计：[2026-09-26_dep3d_gold_design.md](2026-09-26_dep3d_gold_design.md)，用于校核有限目标、侧向散射与天线影响；2D 收敛证据为该设计提供网格参照，不替代三维验证。
