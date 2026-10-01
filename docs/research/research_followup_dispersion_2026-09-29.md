# 面波色散实测数据方向研究进展跟进

跟进日期：2026-09-29（第二轮，修正调研范围）
调研范围：以**野外面波色散实测数据**为核心——数据采集（主动源/被动源/背景噪声）、色散曲线自动提取、反演成像方法、以及滑坡等近地表应用。第一轮误以 UAV-GPR 为范围（见 `research_followup_2026-09-29.md`），本轮纠正。

检索来源：Scholar 学术搜索引擎（2024–2026 年为主，按时间排序），五批检索：
- `scholar_dispersion_inversion.csv`：surface wave dispersion curve inversion field（5 条）
- `scholar_dispersion_landslide.csv`：Rayleigh wave dispersion landslide shallow shear wave velocity（3 条）
- `scholar_ml_dispersion.csv`：machine learning dispersion curve extraction near-surface（15 条）
- `scholar_ambient_dispersion.csv`：ambient noise surface wave dispersion tomography（1 条综述）
- `scholar_gpr_landslide.csv`（第一轮）：其中 Zhang et al. 2025 面波成果归入本报告

---

## 一、总体判断

色散实测数据方向近两年的主线非常清晰：**深度学习正在重塑"提取—反演"全流程**。传统流程（人工拾取频散曲线 → 单点 1D 反演 → 横向插值成图）正在让位于"数据选择/叠加自动化 → 端到端曲线提取 → 快速反演/直接成像"。同时滑坡场景出现了多条独立的实测验证案例，是与工程地质结合最紧密的应用出口。

## 二、实测数据采集进展

- **主被动源联合观测成为主流配置**：Zhang et al. (2025, Applied Geophysics) 在青藏高原深切割地貌用主动+被动面波联合探测做滑坡 2D/3D 成像（滑带划分、基岩界面刻画）；Akin & Sayil (2025, Pure Appl. Geophys.) 在滑坡易发区用主被动源面波 + Rayleigh 波椭圆率（RWE）获取更深 Vs 剖面。
- **线性台阵被动源观测**：Yin et al. (2026, EGU 预印本) 用轻量化 U-net 从线性台阵背景噪声中自动拾取频散曲线，追踪近地表内部结构演化。
- **背景噪声多模式层析方法成熟**：Nishida et al. (2024, Progress in Earth and Planetary Science) 的综述（被引 33）系统梳理了互相关计算到多模式频散层析的实用流程，是该方法目前的高被引基准。
- **新传感器形态出现**：Sun et al. (2026, IEEE) 提出 SA_ConvLSTM 从 DAS（分布式光纤）数据中提取面波频散曲线。

## 三、色散曲线自动提取（当前最热）

- **图像分割/检测类方法**：Gan et al. (2024, IEEE TGRS) 用 U-Net 从主动源面波记录自动提取频散曲线（被引 25）；Hu et al. (2025, J. Applied Geophysics) 图像分割自动拾取（被引 11）；Chamorro et al. (2024, Near Surface Geophysics) 从炮集直接深度学习提取并反演（被引 24）。
- **2026 年新趋势——半监督与注意力机制**：Nan et al. (2026, J. Applied Geophysics) 半监督迁移学习端到端提取；Chen et al. (2026, 地球物理学报) 注意力机制 + 多尺度特征融合，直接处理台站对背景噪声互相关。
- **值得注意的批评声音**：Naskar & Das (2026) 指出机器学习拾取"只在合成数据上表现好"，提出自适应波长依赖提取方法——做 ML 提取时这是一个必须回应的审稿人视角。
- **配套工具与数据质量控制**：DisperPy (2025) 监督+无监督混合 ML 自动拾取群速度曲线（地震与噪声互相关数据通用）；Wang et al. (2025) 用 ResNet-50 做被动源线性台阵数据筛选与叠加。

## 四、反演方法进展

- **从曲线反演走向全谱反演**：Zhang, Alkhalifah & Liu (2024, JGR) 全频散谱反演（被引 42），放弃"先拾曲线再反演"的两步范式。
- **多数据联合反演**：Wang Z. et al. (2024, Surveys in Geophysics) 导波 P + 面波频散联合反演同时估计 P、S 波速度（被引 16）；Sun et al. (2025) 对西山村滑坡用 Rayleigh 波 ZH 比 + Love 波群速度联合反演。
- **算法系统性比较**：Yang et al. (2024, Surveys in Geophysics) 对近地表频散曲线反演算法做全面比较（被引 32），是方法选型的实用参考。
- **稳健框架**：Pan et al. (2026, BSSA) 二次极值插值 + 随机分层多初值框架，用意大利 Mirandola 场地实测数据验证。
- **ML 快速反演**：Zhao et al. (2025)、Zhang & Song (2026) 深度学习反演；Mi et al. (2026, Big Data and Earth System) 综述深度学习在面波分析（采集优化—曲线提取—反演）全链条的应用与挑战。

## 五、滑坡应用案例（实测数据验证）

| 案例 | 年份/出处 | 方法要点 | 被引 |
|---|---|---|---|
| Majiagou 滑坡（三峡库区）| Cao et al. 2025, Bull. Eng. Geol. Env. | 地震面波 + 钻孔联合，2D Vs 结构 | 1 |
| 西山村滑坡 | Sun et al. 2025, J. Earth Sci. | ZH 比 + Love 群速度联合反演 | 0 |
| 青藏高原深切割地貌 | Zhang et al. 2025, Applied Geophysics | 主被动源联合 2D/3D 成像 | 3 |
| 滑坡易发区土体 | Akin & Sayil 2025, Pure Appl. Geophys. | 主被动源 + RWE 深部 Vs | 9 |
| 滑坡监测 ML 反演 | Anjom et al. 2025, EAGE | ML 频散反演用于滑坡监测 | 0 |

## 六、开放问题

1. **合成→实测的泛化鸿沟**：ML 提取方法多在高信噪比合成数据上训练，实测噪声/多模态干扰下鲁棒性存疑（Naskar & Das 2026 的批评即针对此点）；
2. **多模态分离**：高频段基阶与高阶模态混叠仍是实测数据反演不确定性的主要来源；
3. **二维效应**：1D 反演 + 横向插值的范式在横向强变速的滑坡体上有系统误差，直接 2D/全谱反演计算成本高；
4. **监测维度**：多数滑坡工作仍是单次结构探测，利用频散曲线时变特征做动态监测（如 Yin et al. 2026、Anjom et al. 2025）刚刚起步。

## 附：本轮关键文献清单

| 文献 | 年份/出处 | 主题 | 被引 |
|---|---|---|---|
| Nishida et al. | 2024, PEPS | 背景噪声多模式层析综述 | 33 |
| Yang et al. | 2024, Surv. Geophys. | 频散反演算法比较 | 32 |
| Gan et al. | 2024, IEEE TGRS | U-Net 自动提取频散曲线 | 25 |
| Chamorro et al. | 2024, NSG | 深度学习提取 + 反演 | 24 |
| Zhang et al. | 2024, JGR | 全频散谱反演 | 42 |
| Wang Z. et al. | 2024, Surv. Geophys. | 导波 P+面波联合反演 | 16 |
| Hu et al. | 2025, JAG | 图像分割自动拾取 | 11 |
| Akin & Sayil | 2025, PAG | 滑坡区主被动源+RWE | 9 |
| Zhang et al. | 2025, Appl. Geophys. | 高原滑坡主被动联合成像 | 3 |
| Cao et al. | 2025, BEEE | 马家沟滑坡 2D Vs | 1 |
| Anjom et al. | 2025, EAGE | ML 反演滑坡监测 | 0 |
| Pan et al. | 2026, BSSA | 稳健提取+反演框架 | 1 |
| Nan et al. | 2026, JAG | 半监督迁移学习提取 | 0 |
| Chen et al. | 2026, 地球物理学报 | 注意力机制提取 | 0 |
| Naskar & Das | 2026 | 自适应提取；批评 ML 泛化 | 1 |
| Sun et al. | 2026, IEEE | DAS 频散提取 | 0 |
| Yin et al. | 2026, EGU | 被动源自动拾取+结构演化 | 0 |
| Mi et al. | 2026, BDES | DL 面波分析综述 | 3 |

注：被引数为 Scholar 引擎单次检索快照，仅供参考；部分摘要为引擎截断片段，具体实验细节以原文为准。
