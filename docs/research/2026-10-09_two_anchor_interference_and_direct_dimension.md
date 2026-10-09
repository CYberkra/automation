# 两配对锚点、二维/三维直达场与本批中断记录

本单元继续回答仿真为何看不出实测中类似的地层脉络。新证据支持材料损耗和相干叠加共同影响深部可见性；尚不能把全部差异归于PML、二维源或后处理。没有更换正式材料，没有实测标定、三维地质求解或训练。

## 新配对锚点告诉我们什么

已独立审核的24新增+2复用快照包含190→187.5m的11站双损耗H1，以及190/187.5m的双损耗H0配对。公开保存于`artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r5/`。精确20–170MHz/0.3MHz/501点、4.0.1原生FP64、源与Yee空间偏移归一化保持；无AGC、逐道归一化或幅相拟合。

H0仅把底部连通砂岩替换为泥岩，仍有真实覆盖层和全部传播交互。δ=H1−H0表示底砂材料对比的完整场差，不能叫纯一次波，也不能把H0全称噪声。

187.5m的主评价窗为两配方事前模板窗并集340.628077–369.618097ns，不沿用190m的345.618097–374.608117ns。固定300–450ns宽窗作为敏感性检查。

|187.5m主窗量|Hann|Blackman|
|---|---:|---:|
|低损耗δ范数 / 高损耗δ范数|750.857|774.749|
|低损耗H0范数 / 高损耗H0范数|0.99073|0.91816|
|高损耗H0范数 / δ范数|1.38439|0.71686|
|低损耗H0范数 / δ范数|0.001827|0.000850|

低损耗总场与δ在约356.786ns共同达到包络峰；高损耗宽窗总场主峰约318.530/317.698ns，而δ峰为351.796ns。190m也有强总场约330/329ns与底砂差场约357ns分离。因而宽窗中最显眼的波包不能直接标成深部砂泥界面。低损耗使底砂明显增强，并不表示H0绝对幅度消失，也不证明低损耗配方代表现场。

## 相消是局部且随站位变化的

新脚本逐一复算复数组恒等式：

`||H1||² = ||δ||² + ||H0||² + 2 Re<δ,H0>`。

187.5m高损耗主窗的交叉项除以两分量平方范数之和为−0.20925（Hann）/−0.39970（Blackman）；相位约−170.68°/−175.66°，相干幅度约0.223/0.423，支持该分解中有局部相消。190m同窗交叉项约+0.04124，两站不能统一称作反相抵消。187.5m整个300–450ns宽窗相同比值只有−0.01057/−0.01282，不能把局部20–40%的平方范数减少推广到整段波形。

这些是已存复数组上的描述性分解，不是物理能量标定或认证的方向余弦；FDTD误差预算仍未闭合。相位接近180°本身不足以说明强相消，必须同时看相干幅度、范数和实际交叉项。

独立审核以原生直接DFT/直接逆求和复算，另用补零IFFT核对。初次审核发现低幅值内积的相位对逆变换舍入敏感，没有简单放宽一个相位容差：改用两种逆变换的实际数组差、Cauchy–Schwarz界与浮点求和界给出CPU方法一致性范围。192项标量检查通过，最大缩放误差1.67e−10；相位差均在计算的一致性范围内。此范围不替代FDTD来源误差证书。

结果与两张中文图：`artifacts/research_checks/2026-10-09_v401_dense_anchor_interference_r1/`。

## 改为三维点源会自动解决直达波尾部吗

官方[#hertzian_dipole说明](https://docs.gprmax.com/en/latest/input_hash_cmds.html)明确二维是线源；本轮阅读了该段、[Hertzian解析比较章节](https://docs.gprmax.com/en/latest/comparisons_analytical.html)及[官方解析实现中的相关函数](https://raw.githubusercontent.com/gprMax/gprMax/master/testing/analytical_solutions.py)。在线latest/master不冒充当前安装4.0.1的相同源码。读取范围、缓存SHA见`source_reading.json`；下载文件未执行。

CPU解析比较使用真实190m收发间距1.302162m，空气传播时间4.343544ns。二维Hankel线源与包含辐射、感应和静电项的三维点偶极分别合成精确501点SFCW。二维按既有0.025m参考源尺度比较；三维与二维源分布、量纲不同，另图各自按0–20ns复数L2范数归一，只比较形状，不能比较发射功率、场强或地质反射幅度。没有从生产数据扣除解析波形。

原生H0早段与二维自由空间直达场未拟合误差约0.0321%（Hann）/0.0222%（Blackman），支持当前二维早期直达实现内部一致；不能由此认证弱晚时场精度。三维解析直达场仍有有限频带尾部；在190m底砂窗，其自身早段归一后的尾部范数约为二维的1.39–1.40倍。故“三维点源一定消除这些尾波”不成立。实际高损耗330ns强波仍明显超过解析直达尾部，结合既有原生分段证据，不能全部归为早段旁瓣。

独立验证包括：完整三维点偶极场沿不变方向积分复现二维Hankel场（全部501点相对误差约3e−14/9.8e−14），替代横向dyadic表达式一致至4.43e−16，20/95/170MHz自适应标量Green积分至5.42e−16，两份原生DFT至1.05e−13，68项直接逆变换指标检查通过。它们仅认证公式/数组一致性，未包含三维地表、真实天线、端口或仪器。

结果与两张中文图：`artifacts/research_checks/2026-10-09_v401_direct_dimensionality_r1/`。

## 46项批次实际中断，后续不能沿旧状态直接启动

本次实时核查发现远端28/46新增项完成，随后`low_x18675_H1`已STARTED但未COMPLETED，最终事件为`FAILED / SESSION_LEASE_EXPIRED`；owned Python进程为0，GPU已空闲。会话未及时续上600s心跳，属于调度中断，不是已证实的数值发散或物理模型失败。已完成28项与复用2项取回，ZIP SHA256为`4243b0deee56206a67a97e03ea7c7ee9e9eaefca9d34aa169dd573ef79014211`。这不意味着46项完整完成。

28项包含13站双H1（190→187m）及190/187.5m双H0。终止快照及其独立审核/六张灰度图另存`artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r6/`；r5保留新相消分析依赖的原24项字节证据，不用新哈希改写旧研究结果。

13站低损耗峰相对局部模板提前约0.832ns，仅188m站提前约1.663ns，Hann与Blackman一致；不是每站误差都完全相同。沿用图中的“三个锚点”指计划位置190/187.5/185m，当前实际仅前两站有配对，185m与其他未计算列均为蓝灰色空白。终止记录、29项已消耗attempt及17项未开始清单见r6的`continuation_inventory.json`。

已完成28项不重跑。中断项的attempt已消耗，不能原卡自动重试；17项从未STARTED的剩余组应另冻结有界接续契约，复用既有结果并明确缺失站位。不得把更新heartbeat或旧run重新执行当作合法恢复。新的四项波场控制仍未冻结/运行，原门禁要求完整46项证书，当前必须拒绝；后续须审查接续及组合审核规则后再冻结新契约，不直接把门禁改成true。

研究尚未完成：强H0晚到路径没有唯一归因，正式材料没有现场联合复介电谱，有限三维/真实天线与实测导出链尚未匹配。当前证据已否定“只要降航高/换三维/改一幅图就自然解决”的充分性，下一项仍是保留attempt审计的有界接续与覆盖层波场控制。

## 图件与复现

本次新增10张图（r5六张、相消两张、解析两张），终止快照r6另六张。所有未计算列留空、不插值。公开目录均有逐文件SHA交付清单；原生H5、完整复谱与原始资料仅留忽略目录。

- `artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r5/basal_zoom/dense_loss_hann_basal_zoom.png`
- `artifacts/research_checks/2026-10-09_v401_dense_loss_partial_r6/basal_zoom/dense_loss_hann_basal_zoom.png`
- `artifacts/research_checks/2026-10-09_v401_dense_anchor_interference_r1/anchor_interference_hann.png`
- `artifacts/research_checks/2026-10-09_v401_direct_dimensionality_r1/direct_dimensionality_hann.png`

生成/独立审核分别为`diagnose_line9_dense_anchor_interference.py` / `audit_line9_dense_anchor_interference.py`和`diagnose_line9_direct_dimensionality.py` / `audit_line9_direct_dimensionality.py`，各脚本`--help`给出输入路径参数。依赖输入SHA与审核SHA在结果中绑定。密采样仍使用既有`analyze_line9_dense_loss_line.py`、`audit_line9_dense_loss_results.py`、`compare_line9_dense_common_gates.py`及`audit_line9_dense_common_gates.py`。
