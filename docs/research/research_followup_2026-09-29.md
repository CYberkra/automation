# UAV 探地雷达（GPR）滑坡探测方向研究进展跟进

跟进日期：2026-09-29
基准文献：Wang et al. (2026), *Deep-Penetrating Uncrewed Aerial Vehicle-Based Ground-Penetrating Radar for Landslide Surveys*, IEEE Geoscience and Remote Sensing Magazine（工作目录本地 PDF）

检索来源：Scholar 学术搜索引擎（2024–2026 年，按时间排序），三批检索：
- `scholar_uav_gpr.csv`：UAV ground-penetrating radar aerial（15 条）
- `scholar_gpr_landslide.csv`：ground-penetrating radar landslide bedrock interface（15 条）
- `scholar_rtm_gpr.csv`：GPR reverse-time migration terrain topography（10 条，其中 2023 起）

---

## 一、总体判断

该方向 2024–2026 年处于明显升温期，文献集中在四个层面：系统硬件、成像算法、滑坡应用、综述与挑战。本地这篇 Wang et al. (2026, IEEE GRSM) 处于"低频深穿透 UAV-GPR + 曲线坐标 RTM + 基岩界面自动提取"这条主线的核心位置，且其团队（成都理工大学 Tianyang Li 组）在 2025–2026 年连续产出，是该方向最活跃的课题组之一。

## 二、系统与硬件进展

- **Guo et al. (2026, Journal of Applied Geophysics)**：面向危险区域的无人机载 GPR 系统设计与实现，强调危险地带非接触探测的工程化落地。
- **Alva Alarcon et al. (2025, Aerospace)**：面向民用基础设施巡检的无人机载穿透雷达系统设计、电磁建模与分析。
- **Wang Z. et al. (2025, IEEE Antennas and Wireless Propagation Letters)**：双 H 形超宽带天线，针对机载 GPR 载荷约束下的带宽/分辨率问题（被引 5）。
- **Qi et al. (2026, Transportation Geotechnics)**：无人机摄影测量 + GPR 联合评估青藏冻土公路路面–路基变形（被引 6），代表"UAV 光学 + GPR"多传感器融合的应用扩展。

## 三、成像算法进展（与基准文献最相关）

- **Wang W. et al. (2025, IEEE Geoscience and Remote Sensing Letters)**：同一团队的前置工作——RTM-TVD（逆时偏移 + 全变分去噪）并用多 GPU 加速，提升 UAV-GPR 偏移成像质量（被引 3）。
- **Mao et al. (2025, Remote Sensing)**：三维 FDTD(2,4) 子网格算法用于机载 GPR 滑坡模型探测，刻画基覆界面与滑坡体内含水层（被引 5）；另见 Mao (2025) 会议摘要 "Efficient 3D airborne GPR simulation and migration imaging for landslide investigation"。
- **Chi et al. (2024, Remote Sensing)**：复杂地形机载 GPR 成像技术，含起伏地形下基岩界面反射波成像（被引 10）。
- 本地 Wang et al. (2026) 在此基础上更进一步：曲线坐标系 RTM + 贴体正交网格（OBFG）消除阶梯误差 + 概率动态规划自动提取基岩界面及 80% 概率不确定带，四川两个地质场景验证，穿透深度达 20 m，钻孔验证绝对误差 <1 m、相对误差 <10%。

## 四、滑坡应用进展

- **Wang, Li & Yu (2026, 87th EAGE Annual Conference)**：上述 IEEE 论文的会议版延伸——"Deep-Penetrating UAV-GPR Imaging for Inapparent Landslide Investigation in Rugged Terrain"，面向隐蔽型（inapparent）滑坡，野外测线偏移处理 + 基岩界面提取。
- **Famiglietti et al. (2026, Drones)**：意大利南部 Melizzano 案例——低频 UAV-GPR + LiDAR 联合研究边坡变形过程，是独立团队在同一方向的并行验证。
- **Sperandio et al. (2026, Environmental Earth Sciences)**：城市热带区滑坡几何形态，GPR + 传统地球物理联合评估。
- **Zhang et al. (2025, Applied Geophysics)**：青藏高原深切割地层面波 2D/3D 成像（面波方法，GPR 的互补参照）。

## 五、综述、挑战与空白点

- **Catapano et al. (2026, IEEE Geoscience and Remote Sensing Magazine)**：UAV-GPR 三维地下成像的电磁建模、层析处理与开放挑战——目前最权威的方法学综述之一，与本地论文同刊同期，值得精读对标。
- **Giocoli et al. (2026)**：电阻率层析（ERT）+ GPR 在滑坡与桥梁研究中的文献计量综述。
- **Ebrahim et al. (2024, Remote Sensing)**：地下滑坡监测技术混合综述（被引 59），梳理多技术谱系。
- 明确的开放问题（综合各摘要）：
  1. 机载低频天线与载荷/续航的权衡（现有系统多偏高频浅目标）；
  2. 飞行高度、发射功率对数据质量的系统量化（本地论文已部分回答）；
  3. 三维成像与层析处理的计算成本（多 GPU 并行是现有解）；
  4. 结果不确定性的定量表达（本地论文的 80% 概率带是少见的尝试）。

## 六、与基准文献的关系定位

Wang et al. (2026, IEEE GRSM) 是当前少数同时具备「低频深穿透（20 m）+ 复杂山地验证 + 钻孔定量验证（误差 <10%）」的工作；其直接竞品/互补文献是 Mao et al. 2025（仿真层面）、Famiglietti et al. 2026（LiDAR 融合层面）和 Catapano et al. 2026（综述层面）。若做后续研究，自然的切入点是：三维层析而非二维剖面、多时相监测而非单次探测、以及与其他传感器（LiDAR/InSAR/ERT）的联合反演。

## 附：关键文献清单

| 文献 | 年份/出处 | 主题 | 被引 |
|---|---|---|---|
| Wang et al. | 2026, IEEE GRSM | 深穿透 UAV-GPR 滑坡探测（基准） | 0 |
| Wang, Li & Yu | 2026, EAGE 87th | 隐蔽滑坡深穿透成像（会议延伸） | 0 |
| Catapano et al. | 2026, IEEE GRSM | 3D UAV-GPR 建模与开放挑战综述 | 0 |
| Guo et al. | 2026, J. Appl. Geophys. | 无人机载 GPR 系统设计 | 5 |
| Famiglietti et al. | 2026, Drones | UAV-GPR+LiDAR 边坡变形（意大利） | 1 |
| Qi et al. | 2026, Transp. Geotech. | 冻土公路 UAV 摄影测量+GPR | 6 |
| Mao et al. | 2025, Remote Sensing | 3D FDTD 子网格滑坡 AGPR 仿真 | 5 |
| Wang W. et al. | 2025, IEEE GRSL | RTM-TVD 多 GPU 偏移 | 3 |
| Wang Z. et al. | 2025, IEEE AWPL | 双 H 形 UWB 天线 | 5 |
| Alva Alarcon et al. | 2025, Aerospace | 基础设施巡检机载雷达 | 4 |
| Chi et al. | 2024, Remote Sensing | 复杂地形机载 GPR 成像 | 10 |
| Tjoelker et al. | 2024, Remote Sensing | 埋藏冰无人机 GPR | 25 |
| Ebrahim et al. | 2024, Remote Sensing | 地下滑坡监测综述 | 59 |

注：被引数来自 Scholar 引擎单次检索快照，仅供参考；部分摘要为引擎截断片段，标注细节以原文为准。
