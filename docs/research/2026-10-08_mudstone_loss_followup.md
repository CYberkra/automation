# 泥岩损耗来源复核与非平层单因素对照

**本轮已完成3项非平v5 FP64控制并复用1项H0，独立审核通过。** 仅降低泥岩DC电导率，原非平站底砂差场增强7.34/7.40倍，H0/底砂由3.36/3.04降至0.296/0.257，原总场与底砂复相关从0.483/0.628升至0.968/0.979。另独立证明当前正Debye模型类中实部端点差3自带明显损耗下限。材料假设是经实际控制确认的重要因素，但较低σ不等于现场真值，正式参数不变，整线/三维/实测仍未解决。

用户在两机4.0.1配置与九项控制后要求“继续”。本单元保持正式材料，查明当前Debye损耗假设的约束，并在原非平v5站位验证敏感性；不用现场波形拟合参数。

## 执行前冻结设计

原v5站位198.6m/local x178.6、AGL8m，210×42.5m、2.5cm、1200ns、FP64、HORIPML80格=2m，40A/100MHz Ricker，收发及全部原非平体素保持。复用上一轮4.0.1已核验`nonflat_ricker_H0`，不重跑。新建三项：同材料非平Ricker H1、仅泥岩σDC .003→.0003的非平H0/H1。H1几何取已完成v5底砂配对的原H1；H0为仅底连通砂岩替换泥岩。低σ两项只改DC参数和相应诊断元数据，不改Debye/其他材料/源/域/时间或地形。

比较原冻结底砂窗359.756487–383.756487ns与完整20–170MHz/0.3MHz/501点复数SFCW，Hann和Blackman均报告H0/底砂差场、底砂增强、H0变化及总场/底砂复相关。无AGC、去背景、幅相拟合或按结果调窗；同站配置图明确不称空间测线。不设置“必须更清楚”的通过门，减弱或无改善也是有效结果。源/FP64/位置/时步/输入/几何/材料与复用native哈希均审核。每项一次attempt，失败保留，不重跑已完成配置。预计新增约11–12分钟ROG求解，冻结时复核容量与独占锁。

同时复核Wagner2013式33/Table4的分数阶谱与Römhild2019 Table2的DC定义；两份外地不同状态样品不能拼成现场已测配方。重点检查是否真实存在重复计入损耗，及正Debye谱在固定实部跨度下的最低极化损耗。材料预算不是实际探深或实测SNR，解析界限不等于FDTD收敛证书。结果、证据阅读范围和失败另追加。

## 来源复核：没有取得现场损耗真值

阅读台账见`docs/research/2026-10-08_mudstone_loss_sources.json`，新增公开论文缓存仍被忽略，不上传正文。只读相关章节及图表，不声称完整复现。

| 一手文献 | 本轮重新核对的位置 | 对当前配方的意义 |
|---|---|---|
| [Römhild2019](https://onlinelibrary.wiley.com/doi/10.1002/nsg.12073) | 材料/实验、Table2、普通泥岩LA/OL | 七个完整岩塞σ0为0.49–8.83 mS/m，当前3 mS/m在其量级内；样品NaCl饱和且各向异性。它是SIP/DC依据，不给本项目MHz介电谱。诊断0.3 mS/m低于这七个样本的最小值，不能借此论文称其为经实测支持的场地默认值。 |
| [Wagner2013](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/jgrb.50343) | 式33，Table4，§3.5/4；PDF4737–4738页 | 1MHz–10GHz COx黏土岩，分数阶GDR。中间状态表中的σ是微西门子/米量纲的“表观DC贡献”，作者明确其低频不确定性；很小的该项不等于完整MHz材料低损耗。不能只摘小σ、忽略较强弛豫过程，也不能把8ns当完整作者配方。 |
| [Josh与Clennell2015](https://doi.org/10.1190/GEO2013-0458.1) | D129–D132的损耗定义/仪器/制样，D138–D141的样本对比、Fig13 | 已取得原论文PDF。作者区分极化与欧姆贡献，完整岩样与重塑浆体不同，保存岩芯的10–100MHz谱有方向和样本差异。支持联合检查ε′/ε″、水状态、结构与方向；不提供可直接移植到营山的唯一Debye表。 |

上述实测研究没有证明营山必为低损耗，也没有证明目前3 mS/m错误。之前中心ε′95=12、端点差3是研究设计，8ns取另一黏土岩弛豫量级，3 mS/m取普通泥岩DC量级；它们**没有在同一岩样、同一状态下共同拟合**。单独选取每个量看似合理，不代表组成后的复介电谱已被实测支持。Fam1998出版方全文本轮仍不可访问，不从摘要补造数值曲线。Loewer2017相关式11与定义仅做辅助阅读，没有采用其土样参数。

4.0.1实际冻结材料JSON逐501点与官方`DispersiveMaterial.calculate_er`比较通过，相对差约1.6e-17：明确计算一次σDC项加一次Debye项，没有发现实现把同一个项重复加两遍。当前输入σDC来源也是DC而非用ε″反推的MHz有效电导，因此**未发现已发生的同一测量损耗重复记账**。更准确的问题是本构假设未联合校准。官方API一致只证明公式/参数使用一致，不认证所选材料。

## 新的解析约束：实部跨度3本身携带最低损耗

对任意正参数普通Debye极点，令ω1/ωc/ω2分别对应20/95/170MHz，

\[
d(\tau)=\frac{1}{1+\omega_1^2\tau^2}-\frac{1}{1+\omega_2^2\tau^2}>0,\qquad
\ell(\tau)=\frac{\omega_c\tau}{1+\omega_c^2\tau^2}.
\]

单极对实部端点差贡献是Δε·d，对95MHz极化虚部贡献是Δε·ℓ。若正极点总实部差固定为3，

\[
\epsilon''_{\rm pol}(95)=\sum_k \Delta\epsilon_k\,d(\tau_k)\frac{\ell(\tau_k)}{d(\tau_k)}
\ge3\min_{\tau>0}\frac{\ell(\tau)}{d(\tau)}=1.2411779574.
\]

该推导适用于非负极点强度、正时间和非负DC的普通Debye和（含正弛豫时间分布），不是任意因果材料的一般下界。常量ε∞不影响实部差或极化损耗；DC只增加ε″。频点固定后比值在τ→0和∞均发散，导数分子化为τ²的三次多项式。独立精确有理数Sturm序列认证它仅有一个正驻点；独立标量黄金分割得到同一全局最低值。最低点τ=7.45563ns，可用ε∞=11.72110、Δε=5.80249实现ε′95=12，因此此约束在当前物理参数类可取到。

当前8ns的ε″极化=1.243896，距离最低值仅约0.22%。保持这两个实部约束而只换τ，或增加更多正Debye极点，不能把当前95MHz极化吸收大幅消掉；不应把更多极点数当成自动改善探深的方法。这里只讨论连续体物性，不保证离散FDTD精度。

| 当前τ8ns，95MHz | σDC .003 | σDC .0003 | σDC=0解析负控 |
|---|---:|---:|---:|
| 极化ε″ | 1.243896 | 1.243896 | 1.243896 |
| DC ε″ | 0.567635 | 0.056763 | 0 |
| 总有效σ，mS/m | 9.5741 | 6.8741 | 6.5741 |
| 穿过7m泥岩往返的场幅吸收 | −63.128dB | −45.387dB | −43.412dB |

即使采用数学最低点且DC=0，此7m控制层往返场幅吸收也约−43.317dB。7m仅为本控制站；不能把它当全线厚度。预算不含扩散、反射、透射、多次/非局部散射和仪器噪声，不是现场SNR。σ=0这一列只做解析，未新增该FDTD配置。

入口`scripts/review_mudstone_loss_bound.py`；独立`scripts/audit_mudstone_loss_bound.py`与官方API核查`scripts/audit_mudstone_material_v401.py`。图及全部数值/审核：`artifacts/research_checks/2026-10-08_mudstone_loss_bound_r1/mudstone_positive_debye_loss_bound.png`及同目录JSON。没有改正式材料，没有为追求好看结果调整实部跨度。

## 实际非平模型：材料敏感性并非平层特例

三项新增均一次完成，baseline_H1 225.843秒、low_H0 219.734秒、low_H1 221.406秒，合计666.983秒约11.1分钟；RSS峰约4.12GB。原3项输入/几何/材料在4.0.1、CUDA13.3、实际FP64下完成；同域同站位，未缩域、减时窗或重跑已完成H0。没有新增快照观察器。

| 原冻结底砂窗359.756487–383.756487ns | Hann | Blackman |
|---|---:|---:|
| 原H0/底砂差场范数 | 3.3554 | 3.0357 |
| 低σ对照H0/底砂差场范数 | 0.2955 | 0.2566 |
| 对照底砂差场/原底砂差场 | 7.3448 | 7.3997 |
| 对照H0/原H0 | 0.6468 | 0.6255 |
| 原总场/原底砂复相关 | 0.4826 | 0.6279 |
| 对照总场/对照底砂复相关 | 0.9680 | 0.9794 |

较宽的原300–450ns窗也成立：底砂增强7.3449/7.4149倍，H0/底砂从3.8917/3.2574降至0.3897/0.2863，总场/底砂复相关从0.4157/0.6024升至0.9424/0.9730。没有按峰位置重选窗口或拟合幅相。两个窗一致支持材料损耗敏感性，仍不是整条测线轮廓认证。

H0包括真实地质响应，不称纯噪声；σ改变复介电虚部，既影响传播吸收，也影响界面复反射和多层交互，不能把7.34倍全部归为纯吸收收益。底砂H1−H0只对应本次底部替换定义；H0低σ还会改变替代底砂区域的泥岩及其PML延续，因此不是全部实测clean，也不据此删除H0。

执行合同SHA256 `09f1ad819bc1858d995f8d71794a091257eb076739bb6de46049a16ad6fdbbb3`；传回ZIP `9ffcb4e8ff3df6d30fe79c0a245b9006ea852f7aa1ab0a698adffee373966a3e`。复用H0原生SHA256 `e8e90acc5d39440a19af71bf48006e7a436777604c54b7dd4bce1206bacae5eb`。新增原生哈希分别为`d7b6f2104717a73467612daf36501dcdeb95ab689239f23ce1a6bfc8a3e84cef`、`5ab2a44865f6002707cca4994bc9ec8cc1dbf9107dad0175bdef92fe3a197109`、`adf5c7bbd51a9d678c673af000cc6fb114e4855da1c1e1201f1fa84a5901ba7d`。

独立审核重新核对原生版本、FP64、20352样本、时步/整数格位/物理位置、源样本逐位相同；输入除标题外完全相同；高低σ的H1几何逐字节相同、H0几何逐字节匹配复用源；逐列独立验证H0只替换自底边向上的5,606,675个砂岩体素。恢复唯一DC修改和诊断元数据后全部材料相等。独立分块直接DFT原场相对差约1.05e-13，较弱底砂相减后约4.66e-11/9.06e-12；直接逆求和复算全部指标通过。这是处理/来源核对，不把相减后的误差当作FDTD整体误差证书。

准备入口`scripts/line9_v401_nonflat_material.py`；分析`scripts/analyze_line9_v401_nonflat_material.py`；独立审核`scripts/audit_line9_v401_nonflat_material.py`。图件：

- `artifacts/research_checks/2026-10-08_v401_nonflat_material_r1/v401_nonflat_conductivity_hann_gray.png`
- `artifacts/research_checks/2026-10-08_v401_nonflat_material_r1/v401_nonflat_conductivity_blackman_gray.png`
- 同目录`analysis.json`、`independent_audit.json`、冻结合同、原生完成审核、事件日志和准备清单。

完整几何/native/复谱保持在私有`artifacts/local_checks/2026-10-08_v401_nonflat_material_prepared_r1`、`2026-10-08_v401_nonflat_material_results_r1`；ROG根`E:\line9_basal_pair_20261008_r1\v401_nonflat_material_r1`。没有新GIF，因为本批没有快照。不上传原测线、PDF、私有几何或原生H5。

## 当前判断与下一单元

已经有受控证据表明，当前非平v5站位的底砂在仿真中存在，所假设的泥岩复介电损耗会显著削弱它；较强真实地质响应与之叠加，造成总场不随底砂形态直观显示。网格相位误差也独立存在；4.0.1升级没有消除两者。不同v3站位还存在远处夹层侧向响应，不合并成同一个根因。

不能据此宣布“泥岩参数错了”或“UAV高航高必然不可行”。本轮缺少营山泥岩同状态复介电谱，文献亦有强损耗页岩，无法认证0.0003 S/m更接近现场。实测图像可见某条纹理也不独立证明它就是本模型14.3m底砂界面。数学损耗下限只约束当前材料类和假设的端点差，不给现场普适探深极限。

下一步先做**正Debye实部跨度与DC的分开敏感性设计**：保持95MHz中心ε′12、几何/航高不变，分别改变有标签的极化跨度或DC，防止再次只定实部却忽略其必带损耗；范围是机制研究而非场地置信区间。保留现3/.003配方和低σ诊断作为对照，不按图像选择“真实”配方。随后用少量不同层厚/坡度站位检查底砂是否呈现正确相对到时和脉络，再决定是否值得扩大测线。完整1.25cm整线和有限3D仍需独立资源/误差预算，当前不启动。新单元须另冻具体参数/站位合同，不能沿用本批3项额度扩跑。

## 只读复算入口

数学界限可在有NumPy/SciPy/Matplotlib的环境从Git参数重新计算，输出必须新建；此操作不启动求解：

```powershell
.\scripts\run_gprmax_v401.cmd scripts/review_mudstone_loss_bound.py --out artifacts/local_checks/mudstone_bound_recompute
.\scripts\run_gprmax_v401.cmd scripts/audit_mudstone_loss_bound.py --analysis artifacts/local_checks/mudstone_bound_recompute/analysis.json --out artifacts/local_checks/mudstone_bound_recompute/independent_audit.json
```

非平结果独立复核依赖上面列出的私有native和准备包；另一台电脑仅拉Git没有这些原始文件，须使用ROG既有完成目录或按哈希另行复制，禁止为了重现图而重跑已消费attempt：

```powershell
.\scripts\run_gprmax_v401.cmd scripts/audit_line9_v401_nonflat_material.py --source artifacts/local_checks/2026-10-08_v401_nonflat_material_results_r1 --package artifacts/local_checks/2026-10-08_v401_nonflat_material_prepared_r1 --public artifacts/research_checks/2026-10-08_v401_nonflat_material_r1 --out artifacts/local_checks/nonflat_material_independent_recompute.json
```

本轮独立审核通过不代表回归训练/实测验证通过；未修改通用处理算子，不重复运行其182项数组回归。ROG科学批完成后求解树退出，不留后台续扫。
