# 夹层敏感区域的密集平滑源波场：2/2完成

2026-10-08，接续用户自主研究和波场快照授权。源/时窗及两站控制9/9审核完成后，实际冻结并执行本有界波场单元。本文保留求解前设计，并追加完成证据。

原冲激快照覆盖全域但时间抽样约2ns，激励含强高频网格场；不能只凭影片的环纹判定20–170MHz内反射来源。本对照使用已核查与冲激在335ns目标窗频响一致的Ricker100MHz，每17个原生时间步记录一次，约1.002ns；仅缩小观察器视野，不改变FDTD网格或整个110×40m求解域。

两个模型都是c0008、AGL10.275m、原H0、Ricker100MHz/40A、800ns、2.5cm、HORIPML80格、V4原生FP64；只比较H0与x45–60m上覆砂岩→泥岩。两组相同观察器：local x40–85m、y10–37m、z0–0.025m，每0.1×0.1×0.025m记录六分量，共799帧，每帧450×270×1。单模型六场观察器历史4.33975GiB，加6GiB求解储备，必须在目标机冻结时仍能容纳；完整两组快照约9.3GB留ROG，回传审核、固定色标GIF与少量静帧。

准备入口`scripts/line9_interbed_wavefield.py`，已通过同源材料/输入哈希、局部替换不变量、两组相同输入卡及V4原生解析检查。私有包`artifacts/local_checks/2026-10-08_line9_interbed_wavefield_prepared_r1`；公开准备快照`artifacts/research_checks/2026-10-08_line9_interbed_wavefield_prepare_r1/preparation.json`。运行时另冻实际设备/构建/代码/输入及容量、GPU锁、USER_STOP、会话心跳，不排队绕过当前批锁、不重跑已完成的源/时窗批。

执行审核包括认证H0原生接收序列与无观察器Ricker800ns一致，审核全部快照的次数/迭代/原点/间隔/六分量float64及哈希；并比较Ez原场与ROI材料差场的传播。GIF已标注实际Ricker源延迟、观察窗和独立差场色标，不能把脉冲加权影片称完整501频点SFCW影片，不能把局部替换产生的散射自动认证为某个尖灭点的一次绕射。需另附本批同站灰度/源归一化复谱对照，原始场与差场分行标清。

此v3夹层机制仍不能替代第一包v5/8m的验证：v5采用另一套重排层序、泥岩路径更厚且源距顶PML更近，需分别核查FP64弱场、顶边距离及内部多次，而非把v3的远处夹层结论直接套过去。全线/有限三维/真实天线与实测适用性继续保持未完成。

## 执行与独立交付审核

独立运行目录`E:/line9_basal_pair_20261008_r1/interbed_wavefield_r1`；部署ZIP SHA256 `7cfbf4cfbb717f014d85b41d51620a8d10fa74c186dd634414775e7c878e9493`。完整执行合同与原生审核分别为：

- contract SHA256 `d725771b445e2c20f78982d0a8e3568dc9f4881f11c994473e3bd5439ef7a606`
- completed verification SHA256 `4ff5f24a7aab61dd923a162d2c8593f010e64e656a1789446d6e1a70ffbf7153`
- 回传结果ZIP SHA256 `38f29ac4156b74396e1160162f4a176df44a5f79e232e064de9e4e4a164c9632`。

|模型|耗时 / s|峰值所属进程RSS / GiB|
|---|---:|---:|
|H0|85.375|10.385|
|far_removed|85.922|9.937|

2次启动/2次完成，无失败或重试；GPU使用约7177MiB，冻结时保守预算为观察器历史4.33975GiB+原生储备6GiB。799×2=1598帧均在ROG核查全部六分量float64、形状、原点/间距、迭代/时间、有限值及文件SHA256。H0原生接收H5与无观察器的Ricker800ns **整文件字节相同**，不只是波形看似相同。

本地独立审核复算原生源样本/收发位置/网格/事件/代码哈希/媒体哈希和两频窗逆变换；最大逆变换误差2.22e-14。本地核对所有快照收据身份、动画选帧哈希与原生审核一致，不声称把约9.3GB完整六场又下载复算一遍。完整快照留ROG上述`prepared/<case>/cases/full2d_c0008/profile_snaps/`。

## 结果与传播解释

|窗|替换后 / 原H0目标窗范数|新Ricker替换场相对历史冲激替换场变化|材料差场主峰 / ns|
|---|---:|---:|---:|
|hann|0.064980|0.165178%|336.826347|
|blackman|0.081192|0.054626%|336.826347|

平滑源重新正演仍给出6.5%/8.1%的目标窗剩余，源归一化差场主峰336.826ns；上轮冲激的局部归因不依赖原冲激高频图样。这里的差场是“原H0−ROI替换”，不同于此前“原H1−去底砂H0”的底砂差场，不混用峰时或把它当clean。

静帧与动画显示：约201.49ns时材料差场主要在夹层区域内及其附近；随后向观察窗上方和右侧传播，348.84ns时进入接收点附近空气区。这与“远处夹层相关响应返回并在底砂窗叠加”相容，结合既有顶/侧/底边距离控制使PML主导解释不获支持。该观察只描述场的空间传播，不认证唯一发射点/波型/传播阶次。

重要限制：ROI的x60m新截断也会产生散射，差场既含原夹层响应变化，也含替换几何引入的新响应和交互。动画里的圆弧中心不能直接当原地质尖灭点；两模型的实际场、差场必须同时看。图中绿线统一为原H0层界指引（替换图亦如此），紫框为ROI；观察窗以外的场没有显示。没有将地下某条环纹标为已认证绕射/多次。

## 中文灰度、动画与复现

- [本批同站SFCW灰度与包络](../../artifacts/research_checks/2026-10-08_line9_interbed_wavefield_results_r2/wavefield_pair_sfcw_gray.png)：精确501点20–170MHz/0.3MHz；Hann/Blackman并排，无AGC/移时/拟合。三个配置列不冒称连续测线。
- [波场动画167帧](../../artifacts/research_checks/2026-10-08_line9_interbed_wavefield_results_r2/interbed_wavefield.gif)：固定对称对数灰度，原场/替换场共享范围，差场另有固定范围；采集约1.002ns一帧，动画每三帧显示一帧（约3.007ns，至499.209ns）。16.06MB，可入Git。
- [201.49ns静帧](../../artifacts/research_checks/2026-10-08_line9_interbed_wavefield_results_r2/wavefield_03417.png)、[348.84ns静帧](../../artifacts/research_checks/2026-10-08_line9_interbed_wavefield_results_r2/wavefield_05916.png)、[451.09ns静帧](../../artifacts/research_checks/2026-10-08_line9_interbed_wavefield_results_r2/wavefield_07650.png)。

原场固定绝对限幅7091.870974V/m，差场102.980966V/m，对称对数线性阈值各自范围的1e-5；用于传播可视化，不能从颜色直接比较两个色标的强弱。下方曲线是原生Ricker接收序列，含14.142ns源延迟；上面的SFCW灰度已按实际源复谱相除。

渲染脚本`scripts/render_line9_interbed_wavefield.py`只读已审核原始快照，在ROG生成动画/收据；本地入口`scripts/analyze_line9_interbed_wavefield_delivery.py`与独立审核`scripts/audit_line9_wavefield_delivery.py`。源代码、审核及动画均随本批入库；无需重跑FDTD即可复现数值与媒体审核。

后续优先分别核查v5/8m第一包的弱场、顶边距离及层序多次；对本v3模型继续区分实际夹层散射与ROI新边缘，并检查159.25m底砂差场偏离局部模板的原因。不将这两次快照认作195m全线或有限三维/真实天线/实测验证。

图件复核纠正：r1草稿的包络纵轴受早时直耦峰牵引，目标窗被压扁；r2按300–385ns可见区间的共同绝对幅度设纵轴，未归一化或改变数据。r1草稿保存在私有`artifacts/local_checks/2026-10-08_line9_interbed_wavefield_results_r1/public_draft_r1/`，r2数值指标和原生哈希逐项不变，独立审核重做通过；没有重跑FDTD。
