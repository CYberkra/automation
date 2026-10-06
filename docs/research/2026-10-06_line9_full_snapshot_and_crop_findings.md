# 九号线大域快照及缩域对照：实际执行结果

本单元只做X220完整域波场观察及200/160m两个横向域宽控制；没有启动99/391道、三维或训练。SSH恢复后沿用同一ROG、V4.0.0与既有FP64/CUDA环境；旧关机中断五站只复用已完成X220，不重试原attempt。

## 完整域被动快照已完成

400×75×0.025m、2.5cm求解网格、15m地形跟随、四材料、HORIPML80格、1200ns、原z向Ricker100MHz源与1.3m沿线二维代理全保留，仅追加快照。599帧，空间20cm、时间2.004856ns，六分量原始float64、形状2000×375×1全部核查。接收Ez与源波形逐点相同；原始接收文件SHA-256也与原无快照结果完全一致：`2b9af6e225437d9516d6bdf8b358bd09a224729cead0e39b3560442ded2e7689`。

实测658.391s（10.97min），峰值自有RSS30.048GiB；原参考631.171s，含快照本批增加约4.31%时间。原生转存记录确认599帧/21,564,000,000bytes历史，逐帧GPU→内存。快照不是网格加粗，原数值接收没有被观察器改变。

[波场GIF](../../artifacts/research_checks/2026-10-06_line9_large_snapshot_rog_r3/line9_large_domain_wavefield.gif)（17.86MB、201帧）、[402.98ns静态总场](../../artifacts/research_checks/2026-10-06_line9_large_snapshot_rog_r3/large_domain_wavefield_example.png)、[可读单站灰度图](../../artifacts/research_checks/2026-10-06_line9_large_snapshot_rog_r3/single_station_details/single_station_gray.png)、[原数值分窗A-scan](../../artifacts/research_checks/2026-10-06_line9_large_snapshot_rog_r3/single_station_details/single_station_windowed_ascans.png)。原整线轴单站图只有左端0.5m宽，难看清，因此额外给出明确标注的一站横向展开图；没有补道或插值。

GIF全帧共同SymLog尺度，动画隔3帧约6ns，只用于浏览波包及路径；不能从条纹运动认证100MHz相位传播方向。准确时序须查599帧原始H5（约20.08GiB），留在ROG `E:\line9_rog_snapshot_accounted_20261006\artifacts\local_checks\line9_large_snapshot_rog_r3\cases\full2d_pilot_0\profile_snaps`。原2ns保存仍不是完整源高频谱的认证。动画是Ricker瞬时总场，非SFCW波场。

## 两次失败及修正据实保留

r1为不支持的旧式CLI参数，argparse退出，没有FDTD；r2这份安装V4源码的延迟初始化导致容量检查时快照nbytes为零，未触发原生转存，资源保护停止，无完整输出。旧监督器未记录r2触发瞬间RSS，不能回填实际峰值。r3单独worktree/新契约，在进程入口只于原生容量检查期间补入实际尺寸/dtype/分量推导的nbytes，再恢复供原生分配；没有修改安装目录或场更新/快照内核，也没有放宽40GiB RSS门槛。运行后的原始接收哈希一致与转存日志核验通过。失败契约、日志与worktree保留，见[失败记录与检查](../../artifacts/research_checks/2026-10-06_line9_domain_acceleration_r1/)。

## SFCW与观感的处理因素

仍为20–170MHz、0.3MHz步长、501个真实复数频点，按实际源半步时间轴除源谱与0.025m电流矩尺度，不称端口S21；矩形及均值归一化Hann、8倍零填充、无额外时间平移/增益/背景扣除。7个频点独立连续DFT相对L2为1.55e−13，独立逆求和矩形/Hann为7.24e−14/2.83e−14，仅认证实现一致性。

原始300–600ns峰值0.6064V/m，仅为整道峰值3738.21V/m的1.62e−4（约−75.8dB）。整体L2/全幅灰度容易由直耦主导，不能代替弱事件保真核查。相同复数谱在600–1200ns的矩形/Hann带符号峰值23.394/0.09737，相差240.25倍（47.61dB）；因此晚时段观感受频窗明显影响，不能把所有振荡都归于侧边PML。Hann改变谱权和分辨率，这个峰值差不是SNR改善或现场导出链验证；设备实际频窗仍未知。

全域顶侧PML空气镜像路径估计约138ns/854ns；200m/160m最近侧PML约654ns/520ns，不含源延迟。这些只是路径估计，不能由条纹或到时单独证明边界反射的幅度/主导性。瞬时图显示空气波、地表反射与地下传播同时存在；没有物理天线本体，不能把理想Hertzian源的图认作天线本体反复散射证据。

## 缩域实测结果

两道均完成，原始FP64、源样本/参数、空间位置、时间网格和完整接收长度审核通过；相同保留区域再次逐体素复验。**不能签认200m或160m为完整1200ns等价替代，保持400m基准，不放行整线缩域。**

| 域宽 | 单道实测 | 相对400m原无快照 | RSS峰值 | 300–600ns原始L2差 | 300–600ns Hann差 | 600–1200ns Hann差 |
|---|---:|---:|---:|---:|---:|---:|
| 400m参考 | 631.171s | 1.00× | 10.00GiB | — | — | — |
| 200m | 326.781s | 1.93× | 5.17GiB | 0（逐点一致） | 6.54e−6（0.000654%） | 67.82% |
| 160m | 275.922s | 2.29× | 4.17GiB | 8.78% | 9.33% | 56.84% |

全501点复谱全局L2差只有9.00e−6/1.66e−5，弱晚窗却相差约57–68%，直接证明全谱/直耦主导的总体误差不能代替弱事件评价。200m在X220的0–600ns原始波形逐点一致，但有限带宽重建把晚时段变化以旁瓣形式带入早窗，故Hann差非零；这不授权删掉600ns之后或改变1200ns任务。160m已在目标相关弱窗产生变化，不推荐用作当前正式批。

[六联单站灰度](../../artifacts/research_checks/2026-10-06_line9_x220_crop_controls_rog_r1/crop_single_station_bscan.png)、[300–600/600–1200ns实幅及差值](../../artifacts/research_checks/2026-10-06_line9_x220_crop_controls_rog_r1/independent_local_checks/weak_windows_and_residuals.png)、[数值报告](../../artifacts/research_checks/2026-10-06_line9_x220_crop_controls_rog_r1/comparison_report.json)。本地独立取回三份原始H5，SHA复验，全部501点按H5真实源/接收半步时钟直接DFT，复数L2误差约1.4e−13；76个逆求和时刻及每窗fsum独立L2复算通过。原始H5留忽略目录，图/数值/复现脚本入库。

两个缩域同时移动PML且移除域外地质散射，差异只能归因于横向域截断敏感性，不能独立判为PML回弹。下一步先定位599帧中晚路径与域外地表散射，再用保留散射区、加大外侧缓冲的更保守候选及其他站位核查；无条件“不影响质量”现在没有证据。99站同硬件线性外推400/200/160m约17.36/8.99/7.59h，仅为速度预算，后两者没有完整质量放行。
