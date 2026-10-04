# 2026-10-04 迁移成像恢复大尺度界面形状 + 覆盖层介电常数扫描

> 对应用户指令："大尺度形状仍可成像、临界角随覆盖层介电常数变——能否从这些角度上尝试做好我们的模型？"
> 本单元分两步：A) 用已归档 15 m 十三站位数据做 Kirchhoff 迁移成像诊断（纯 CPU，无新求解）；B) 冻结并执行覆盖层 εr 扫描（18 组 78 道新求解），对每个 εr 重复同一迁移诊断，比较可恢复性。
> 范围声明：全部为二维 TM 不变 Y 仿真（HS4T2D 0.8 m 起伏界面、36 m 域、2.5 cm X/Z 网格、HORIPML），不是三维/实测结论；迁移只是诊断手段，不等于首版处理范围新增成像模块。

## A. 迁移成像诊断（无新求解）

- 脚本 `scripts/migrate_hs4_height15m_diff_v0_1.py`，输出 `artifacts/research_checks/2026-10-04_hs4_migration_diagnostic/`（summary.json + 两张图）。
- 数据：已归档 15 m 胶囊 13 站位（centre + remaining，全部 manifest 哈希核验），官方 SFCW 链重建（impulse 源、Hann、8x 零填充、200 ns 物理尾锥），锚点相对 L2 = 1.98e-11（门禁 1e-9）。
- 方法：Kirchhoff 叠加，13 道共偏移 1.3 m；旅行时为双层（空气/覆盖层 v=0.0706 m/ns）过平地表的费马最短时间（每腿对地表穿越点 600 点稠密扫描）。时间零点用各地面反射包络峰与空气双程预测的中位数校准（−0.3627 ns，13 站位恒定）。成像网格 x 13.5–22.5 m × z 8.2–12 m。
- 山脊提取：差分像 |rough−halfspace| 每列 50% 峰值以上质心（非 argmax，避免子波旁瓣跳变）。
- 结果（εr=18.017 色散基准）：
  - 山脊 vs 真实起伏：全频 corr **0.843**；>1.58 m 大尺度 corr **0.835**；<1.58 m 细节 corr **0.115**；RMSE 0.221 m（系统性偏深约 0.2–0.5 m，形状相关是重点）。
  - 模板匹配（高斯走廊归一相关）：真实 0.718 ≈ 低通 0.713 > 平移 1 m 0.655 > 平面 0.567 ≈ 平面+高频 0.566 > 平移 2 m 0.518 —— 成像像确实锁定真实大尺度形状，而非任意平面或移位版本。
  - 半空间对照界面带能量仅为差分的 0.99%（干净对照，差分像条带不是参考自身带入）。
- **结论 A：即使全部回波经过 13.6° 临界角滤波，大尺度（>1.58 m）起伏形状仍可从 15 m 机载数据经迁移成像恢复；亚分辨率细节不可恢复，与临界角分辨率上限 c/(2f)=1.58 m 的预测一致。**

## B. 覆盖层介电常数扫描（新求解）

- 契约/运行/验收脚本 `scripts/hs4_permittivity_scan_v0_1.py`，运行器 `scripts/run_hs4_permittivity_scan_v0_1.cmd`，胶囊 `artifacts/research_checks/2026-10-04_hs4_permittivity_scan/`。
- 分组：εr ∈ {6, 9, 18.017} × {rough, halfspace} × {centre(-n 1), left(-n 6), right(-n 6)} = 18 组 78 道。覆盖层改为**非色散**（删 Debye 行），σ=0.003 不变，rock（εr=9, σ=0.001）不变；εr=18.017 组（eps18nd）同时充当"去色散"对照，与已归档色散基准（eps18disp）分离 Debye 项的贡献。站位几何、网格、PML、impulse 源、时窗与 15 m 联合网格研究完全一致；319 个 V4 .py 与 8 m 冻结契约逐字节一致。
- 执行：V4.0.0 / CUDA double，GPU 锁+资源守卫，单次 attempt。18 组全部求解成功（约 10 分钟）。**已知过程瑕疵**：冻结的检查函数要求事件日志末尾为总 COMPLETED，而总事件在检查之后才写入，导致运行末尾抛 FAILED；求解器本身 18 组全部成功，原始输出完整保留。随后修复检查逻辑（改为逐组核对 COMPLETED 事件），以修复后脚本独立复验 **PASS**（输入不变量、覆盖层材料替换、站位格点、float64/有限值/迭代数、实际材料几何与基准 vtkhdf 逐元素一致），并把补验事件如实追加到 execution.jsonl。预检与复验的脚本哈希不同（检查逻辑修复），契约与全部输入未变。
- 迁移对比分析 `scripts/migrate_hs4_permittivity_scan_v0_1.py`，输出 `artifacts/research_checks/2026-10-04_hs4_permittivity_migration/`。对四种情形用**完全相同的**成像网格/度量/山脊提取，仅覆盖层速度 v=c/√εr 随情形变化；分辨率上限 1.58 m 按构造固定。色散基准重算与 A 单元逐数字一致（corr 0.8432 等，锚点 1.98e-11）。

### 关键数字

| 情形 | 临界角 | corr 全频 | corr >1.58m | corr <1.58m | RMSE (m) | 模板(真实) | 模板(平面) | 回波/直达波 | 回波峰时 |
|---|---|---|---|---|---|---|---|---|---|
| eps6 | 24.1° | 0.940 | 0.938 | 0.075 | 0.080 | 0.688 | 0.608 | **−51.9 dB** | 152 ns |
| eps9 | 19.5° | 0.962 | 0.949 | 0.155 | 0.068 | 0.707 | 0.575 | −70.7 dB | 163 ns |
| eps18nd（非色散） | 13.6° | **0.977** | **0.967** | 0.135 | **0.069** | **0.731** | 0.563 | **−50.3 dB** | 193 ns |
| eps18disp（色散基准） | 13.6° | 0.843 | 0.835 | 0.115 | 0.221 | 0.718 | 0.567 | −70.3 dB | 175 ns |

半空间/差分界面带能量比：0.0013 / 0.0501 / 0.0001 / 0.0099（均干净；eps9 最高因界面信号本身最弱）。

### 解读

1. **分辨率上限不随 εr 变（预测证实）**：四种情形 corr(<1.58 m) 全部 ≤0.155，亚分辨率细节在任何覆盖层介电常数下都不可恢复。c/(2f) 是硬底。
2. **大尺度形状在全部情形下可成像**（corr 大尺度 0.835–0.967），包括最差的原色散覆盖层——"机载+高介电覆盖层就看不到基覆界面形状"不成立，前提是做迁移成像而非看原始 B-scan。
3. **εr 的效果被界面阻抗反差混淆**：eps9 回波最弱（−70.7 dB）不是因为临界角，而是 cover εr=9 与 rock εr=9 几乎无介电反差（仅剩 σ 0.003 vs 0.001 的电导反差），界面反射系数趋零。即使如此，差分像仍锁定形状（corr 0.949）。**对真实场地，覆盖层 εr 不是自由参数；这个扫描说明的是物理机制，不是选材建议。**
4. **色散/损耗是比临界角更大的退化源**：eps18nd vs eps18disp，回波强 20 dB（−50.3 vs −70.3 dB）、corr 大尺度 0.967 vs 0.835、RMSE 0.069 vs 0.221 m。Debye 项的损耗把高角度分量进一步衰减。注意 eps18disp 山脊系统性偏深约 0.2 m——迁移用的非色散 v=c/√18.017 对色散覆盖层偏高估速度，速度模型误差直接转成深度偏差；真实场地处理必须实测/标定覆盖层速度。
5. **回波峰时随速度精确移动**（152/163/193 ns vs 预测 148/159/184 ns），色散情形因高频有效 εr 较低而提前（175 ns），与 Debye 谱一致。

### 对"能否做好模型"的回答（限本仿真范围）

- 可落地的改进：(a) 处理链以迁移成像（正确双层速度模型 + 时间零点校准）恢复大尺度界面形状，这一步对本项目全部已测覆盖层条件有效；(b) 速度模型必须来自标定而非假设——0.2 m 级深度偏差即由此产生；(c) 不要追求 <1.6 m 的界面细节，该信息在过地表时已物理丢失（任何 εr 下）。
- 不可行的方向：指望改变覆盖层 εr（场地属性不可选）；指望加大孔径恢复亚临界角细节（孔径只改善已透射锥内的聚焦，锥外分量仍困在覆盖层混响里）。
- 下一步候选（未冻结）：把迁移速度模型做成随深度/色散修正的版本并量化深度偏差收敛；或在保留实测数据介入前，用本扫描框架评估增益补偿对 −70 dB 级回波的可恢复余量。

## C. 可视化产物（2026-10-04 追加，用户要求"看各个 B-scan、波场快照"）

- **四情形 B-scan**（`artifacts/research_checks/2026-10-04_hs4_permittivity_bscans/`，脚本 `scripts/plot_hs4_permittivity_bscans.py`，纯 CPU 复用已验收 H5）：
  - `bscans_by_eps.png`：rough/halfspace/diff 三行共用固定标尺（±176/±176/±1.13）。rough/halfspace 行只见 ~100 ns 地表反射，地下全淹没；diff 行只有 eps6、eps18nd 条带可见且贴真实起伏双程时曲线。
  - `bscans_by_eps_selfnorm_dB.png`：各面板自归一 dB 包络（仅显示用，峰值已标注）。四种情形 diff 条带均跟随真实起伏曲线；峰值 eps18nd 1.13 > eps6 0.62 ≫ eps18disp 0.086 ≈ eps9 0.084。
- **各 εr 波场快照批**（胶囊 `2026-10-04_hs4_permittivity_wavefield`，脚本 `scripts/hs4_permittivity_wavefield_v0_1.py`）：15 m 航高、中心站位、ricker 95 MHz，快照 ROI/1.0024 ns 帧隔与 8m 批相同（599 帧/组），eps6/eps9/eps18nd × rough/halfspace 共 6 组，V4.0.0/CUDA double，探针/快照闭合 0.0、几何与已跟踪参考逐元素一致，验收 PASS。色散 15 m 属前机已消耗 attempt，不重跑；其波场参照用仓库已有 8m GIF。
- **GIF**（`artifacts/research_checks/2026-10-04_hs4_permittivity_wavefield_visual/`，脚本 `scripts/plot_hs4_permittivity_wavefield_gif.py`）：`wavefield_15m_eps6.gif`（2.9 MB）、`wavefield_15m_eps9.gif`（2.5 MB）、`wavefield_15m_eps18nd.gif`（4.4 MB）；每条 300 帧（2.0049 ns 步进），三面板（rough Ey / rough−halfspace 差分场 / 差分 log10 包络），**三情形共用全局固定标尺**（rough ±356，diff ±11.9）可直接对比。快照帧与 vtkhdf/cuda_cache 留本机（.gitignore），帧哈希在 completed_verification.json。

## 复现入口

- A：`artifacts/local_checks/hs4_budget_env/Scripts/python.exe scripts/migrate_hs4_height15m_diff_v0_1.py --out <新目录>`（输出目录必须不存在）。
- B 分析：同环境 `scripts/migrate_hs4_permittivity_scan_v0_1.py --out <新目录>`（要求扫描胶囊 completed_verification.json PASS）。
- C：同环境 `scripts/plot_hs4_permittivity_bscans.py --out <新目录>`；GIF 渲染 `scripts/plot_hs4_permittivity_wavefield_gif.py`（固定输出目录，需波场胶囊验收 PASS）。
- 182 项数组回归于本单元完成后复跑 PASS（`artifacts/local_checks/20261004T060039Z-83ec97b4`）。
- 求解环境纪律不变：.bat 内部设 env、完整继承环境启动（NVML 在剥离环境下不可初始化）、GPU 锁、单次 attempt、失败胶囊保留。
