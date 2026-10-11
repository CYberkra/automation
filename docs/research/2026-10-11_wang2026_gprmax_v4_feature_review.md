# Wang 2026 论文仿真方法与 gprMax v4 官方功能对照研究

- 日期：2026-10-11
- 性质：只读研究。不含任何模型/配置修改，不启动新求解；文末建议均需用户显性批准后执行。
- 论文来源：`探地雷达背景资料/Wang 等 - 2026 - Deep-Penetrating UAV-based GPR for Landslide Surveys.pdf`（IEEE GRSM 2026 接收版），阅读范围：全文 21 页（本地文本提取 `artifacts/local_checks/wang2026_fulltext.txt`）。
- 官方手册来源：docs.gprmax.com en/devel（v4）：`input_api.html`、`input_hash_cmds.html`、`features.html`、`accelerators.html`、`output.html`、`inc_SFCW.html`（逐页核读，检索 "OSG"/"huygens" 等关键词无命中处已注明）。

## 一、论文仿真部分事实（与我们项目的对应关系）

1. **论文不用 gprMax**。正文仿真为自研 OBFG-FDTD（Ryskin–Leal 贴体网格 + 曲线坐标 Maxwell 方程，式 1–9）；gprMax 仅作参考文献 [73] 出现。因此论文方法不能"直接用 gprMax 复现"，只能抽取思想。
2. **多区域策略（用户所称"计算空气传播时间以加速"）**：式 10–11。在天线到地表之间不建 FDTD 网格，在贴近地表处设虚拟入射面，入射场按 u_inc(t) = u(t − R/c) / R 解析给定（R = 源到入射面观测点的直线距离，τ = R/c 为空气传播延迟，1/R 为几何扩散），FDTD 只推进地下区域。
   - 缺口（论文未覆盖、对我们关键）：论文接收阵布在入射面附近（正演试验接收机在 y = 1 m），即不做"地下散射场回传到 8 m 空中的机载接收机"这一步。我们的 B-scan 接收机在空中，采用该思路还需补回传路径；且该方案只含直达几何扩散，不含空气层内天线–地表多次往返。
3. **OBFG 贴体网格**：消除起伏地表的阶梯伪影。论文定量证据：无噪声 RTM 下 OBFG 对矩形网格 SSIM 0.81 vs 0.77、LCC 0.53 vs 0.46、图像熵更低（图 10、11）。论文矩形网格阶梯步长 Δh ≈ 0.1–0.3 m，对 100 MHz、εr≈4 土壤（v ≈ 0.15 m/ns，λ ≈ 1.5 m）为 λ/5–λ/15，伪影显著。
4. **材料表（Table 1）**：滑坡模型 100×20 m：土层 εr=12/σ=0.003、泥岩 εr=18/σ=0.006、砂岩 εr=9/σ=0.001（即此前已提交用户的候选 1、2 的出处）。正演测试模型 20×15 m：土壤 εr=4/σ=0.001。
5. **采集与设备参数（Table 2 及实测节）**：20–170 MHz、501 频点、步进 300 kHz——与用户设备一致，确认用户系统即该论文系统。**新证据：时间采样 1.4 ns = 700 ns / (501−1)，支持 T/(N−1) 端点约定**（营山 CSV 头 700 ns/501 样本互洽；雅安 900 ns → 1.8 ns）。道间距约 0.1 m。仿真采集航高 5 m（滑坡模型），实测 11 m/16 m 两航高剖面高度一致、16 m 略高 SNR（图 17 及正文）。
6. **预处理链**（图 2，仅列名无公式）：直耦波去除、增益补偿、带通滤波。无参数细节，不能作为我们首版算子的直接依据。
7. **RTM/ERM 成像**：爆炸反射面模型单程化，v' = v/2、α' = α/4（式 12–16），较双程 RTM 省一半；速度模型 = 真模型 5×5 高斯平滑；速度不确定度 ±10% 计入 80% 概率不确定带（DP 界面提取，式 17–18）。实测 6 个钻孔深度对照误差 0.27%–9.33%（Table 4）。成像属本项目延期任务，仅存档。
8. **成像评价指标**：SSIM、LCC（局部互相关）、图像熵（能量聚焦）。可作后期成像评价参考，与我们评价 v0.2 分域逻辑不冲突。

## 二、gprMax v4 官方手册核对结果

1. **Python API 建模（用户点名项）**：v4 官方推荐方式为 `import gprMax` + `Scene` 对象式建模，`gprMax.run(scenes=[...])` 运行；**v3 的 `#python:`/`#end_python:` 输入文件内嵌块已弃用**（input_api.html "Introduction"）。2D 用 `DomainMode('TM')` + 某维 `inf`。逐站位循环天然由 Python 代码控制。我们现行"自写生成器输出 .in 文本"仍合法有效，但官方 API 是以下若干特性（HSG 子网格）的唯一入口。
2. **子网格（用户所称"HSG 其他形式"）**：v4 只有 **HSG（Huygens Subgridding）** 及其带人工损耗变体 **SHSG**；**OSG 在 v4 文档中不存在**。HSG 只能经 Python API 创建（`gprMax.SubGridHSG(p1, p2, ratio, id)`），支持 CPU/CUDA，ratio>1 强制双精度；ratio=1 模式不细化、用于局部色散材料存储。源/接收机/快照可放入子网格。子网格边界 PML 间距固定 ratio//2+2。（features.html "Subgridding"；input_api.html "Subgrid"；accelerators.html "Solver precision"）
3. **免除空气层的官方机制：不存在**。无 `#huygens_surface:` 入射注入命令；`#plane_wave_angles:` 是 TFSF 平面波注入（非偶极近场、波形受限）；NTFF 系列是场→远场的观测后处理，不是源注入。官方文档对机载天线大空气区无免除手段（input_hash_cmds.html 命令索引核实）。
4. **地形/非均质辅助命令**（官方、文本输入即可用）：`#fractal_box:`（分形属性体，可定种子）、`#add_surface_roughness:`（分形粗糙面）、`#add_surface_water:`、`#add_grass:`、`#soil_peplinski:`（**官方有效范围 0.3–1.3 GHz，不适用于我们 20–170 MHz**，与 AGENTS.md 既有判断一致）、`#material_range:`/`#material_crim:`。官方警告：分形材料贯穿 PML 时 PML 处理有限。
5. **性能/GPU**：CUDA 支持 2D 模型（官方示例即 `cylinder_Ascan_2D.in -gpu`）；FP64 用 `-gpu_precision double`（我们已 FP64）；多站位 B-scan 官方并行路径是 **MPI task farm**（`--taskfarm`，多 GPU 时每秩一卡）；逐站位输出合并官方工具 `toolboxes/Utilities/outputfiles_merge.py`。快照的 CUDA 支持状态文档无明文（我们的 GPU 快照批次已实际成功运行，属实证）。
6. **官方 SFCW 路径确认**：官方 SFCW 工具箱 = impulse 宽带正演 + 后处理卷积合成逐音调响应（inc_SFCW.html，GIA2023SFCW），与我们现行 impulse+direct 路线一致。
7. **色散材料语法**：v4 为 `#material:` + `#add_dispersion_debye:`（多极=同一材料多次 add）；v3 的 `#material_dispersive:` 在 v4 不存在。色散材料构建后有自动时间步稳定检查。

## 三、对现行模型（`model_current/`）的建议清单（全部待显性批准）

| # | 建议 | 依据 | 对我们的意义 | 预期代价/风险 |
|---|---|---|---|---|
| S1 | 维持现状：不迁移求解器、不引入自研 OBFG | 我们 2.5 cm 网格阶梯步长 ≈ λ/48（100 MHz、土壤 λ≈1.2 m），论文的阶梯问题（λ/5 起）对我们会小约一个量级；gprMax v4 也无贴体网格 | 无行动 | 无 |
| S2 | 空气层解析化（论文式 10–11）**不采用** | v4 官方无注入机制；论文方案不覆盖"回传到机载接收机"和空气层多次往返，而我们的关键发现（232/328 ns 覆盖层多次返回强带）正依赖完整空气–地表路径 | 自研成本极高且可能破坏已验证现象 | 建议永久搁置，除非官方未来支持 |
| S3 | 中期：把模型生成从自写 .in 生成器迁移到官方 Python API | `#python:` 块已弃用、Scene API 是官方方向；且是 HSG 的唯一入口 | 工程迁移，一次完成后续受益 | 需逐算例数值对照验证；单独批准 |
| S4 | 局部细化需求出现时（滑带/薄层）用官方 HSG ratio=2 | features.html/input_api.html | 替代全域 12.5 mm 加密（全域加密成本 ×8） | 依赖 S3；FP64 已满足 |
| S5 | 地表微起伏/覆盖层非均质敏感性试验：`#fractal_box:` + `#add_surface_roughness:` | 官方文本命令，无需 API | 检验"地表阶梯/起伏"对我们结论的实际影响，比 OBFG 便宜得多 | 一批新求解（约 12–13 min/对），单独批准 |
| S6 | 多站位整线扩展时用官方 `--taskfarm` + `outputfiles_merge.py` | accelerators.html "Task farm"；output.html | 替代我们自写调度/合并脚本 | 工程改动，单独批准 |
| S7 | CSV 端点约定暂采用 T/(N−1)（1.4 ns 互洽证据）并在文档中标注证据强度 | 论文 Table 2 + 营山头部 700 ns/501 | 影响时间/深度坐标换算 | 只是解释约定，不改数据；标注"单篇论文+头部互洽，非官方确认" |

此前已提交待批的四条候选（砂岩 εr9/σ0.001、覆盖层论文土 εr12/σ0.003、覆盖层 εr=18.4 证伪、航高 8→5 m）与本表独立，仍待你逐项决定。

## 四、证据边界

- 论文仿真结论出自其自研 OBFG-FDTD，不是 gprMax 行为；其 RTM/DP 指标为作者报告。
- v4 文档为 devel 版页面；HSG/Scene API 我们尚未实跑，S3/S4 属"官方支持但未验证于本项目"。
- S7 的 1.4 ns 是论文系统参数与营山 CSV 头部字段的互洽证据，不构成设备导出链的官方确认。
