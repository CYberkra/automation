# 第8道底砂配对与波场快照：SSH执行记录

2026-10-08。用户明确授权“通过ssh连接我们那个机器然后开跑吧”。本批为两个单道V4.0.0/CUDA原生FP64模型，串行；实际两份模型完成、一次H0被心跳中断并保留，另建一次H0恢复attempt。未扩道、未改频带/网格，不混用三包作为纯航高控制。

## 冻结方案

第二组`full2d_c0008`，剖面横坐标179.25m，native局部Tx=(78.6,35,0)、Rx=(79.9,35,0)，名义z=0.0125m折叠为二维唯一平面；收发中点离地10.275m。域110×40m、网格2.5cm、TMz、800ns、HORIPML80格，沿用实际40A impulse源、原四材料及色散/界面平均设置。

- H1：原始几何，完整保留底部连续砂岩。
- H0：逐列仅将与模型底边连续相连的ID3砂岩改为ID2泥岩，保留全部上覆夹层；包括该区域在侧/底PML内的材料延续，不是单独PML对照。
- 两份输入卡逐字节相同（原卡仅追加快照），各自导入不同几何。实际替换2,341,672体素；未替换体素逐点不变。原输入、原native与原几何均只读，新私有包另存。
- 两次各256帧；空间0.1×0.1×0.025m、1100×400×1。前480ns约2.005ns间隔，之后约20.049ns；六分量历史每次5.035GiB。原生求解网格/时间步不变，动画不作完整脉冲源频谱的相位传播认证。
- 原始接收与源必须float64，实际时间步/样本数/源样本/源收位置与既有FP32参考相同。与FP32的比较仅报告精度/构建敏感性，不假设FP32一定不可信。
- 源归一化采用实际源采样与Yee偏移，精确20–170MHz/0.3MHz/501点，无尾渐消、AGC或逐道归一化；H1−H0必须先在复数域相减。差分包含材料改变引起的传播及多次相互作用，不冒充纯一次回波或实测clean。

## 目标机与保护

SSH连接ROG@172.27.22.61，主机LAPTOP-83GRKUNE，公钥指纹与用户确认的`SHA256:P5pIrBpnIXSicC5Zs5kw1WUa5/j9vqQX1OlgScJwA+w`一致；凭据不入库。首读4090 Laptop16GiB、约48GiB可用RAM、GPU约1%且无Python进程。实际使用既有V4环境`E:/automation_djh/artifacts/local_checks/gprmax_v4_gpu_env2`；旧E盘MSVC junction在SSH中不能遍历，使用已核验的C盘真实BuildTools/CUDA13.3路径，不改用户junction或安装目录。

新运行在独立目录部署，保留目标机`E:/automation_djh`的未提交资料与分支。共享GPU文件锁仍沿用其`artifacts/local_checks/hs4_gpu_exclusive.lock`；首读及每道启动前要求空闲VRAM≥12GiB、可用RAM≥20GiB，并显式要求快照历史+6GiB原生数组预留可容纳。RSS≤24GiB、系统剩余RAM≥2GiB、单道≤30min、两道≤60min。会话心跳120s过期或USER_STOP只终止本任务树；失败attempt保留、不原地重试。

原生快照历史本批可容纳于显存，使用既有缓存CUDA入口；不套用大域强制streaming补丁，不修改场更新。目标freeze绑定全部代码、实际V4源码/扩展、编译器/Python和输入哈希。准备不等于启动，启动不等于完成。

## 准备与复现

[准备数值](../../artifacts/research_checks/2026-10-08_line9_basal_pair_prepare_r1/preparation.json)、[实际两份地质图](../../artifacts/research_checks/2026-10-08_line9_basal_pair_prepare_r1/geometry_pair.png)。本批脚本`scripts/line9_basal_pair.py`复用既有Windows监督器，独立audit适配本批域/impulse/快照。底连通区域/保留夹层/异常拒绝的数组检查通过，V4输入解析通过。

```powershell
D:\gprmax_v4_gpu_env\Scripts\python.exe scripts/line9_basal_pair.py prepare --package artifacts/local_checks/2026-10-07_line9_three_packages_review_r1/line9_pkg2_x160-110_agl35_142st --review artifacts/research_checks/2026-10-07_line9_three_packages_review_r3/line9_pkg2_x160-110_agl35_142st_audit.json --out artifacts/local_checks/line9_basal_prepare_reproduce --evidence artifacts/local_checks/line9_basal_prepare_check_reproduce
```

私有输入包及FP32原始参考需单独带到目标机；Git只含代码、聚合证据及图。目标机通过`run_line9_2d_v4.cmd`设置实际GPRMAX_PYTHON/GPRMAX_VCVARS/GPRMAX_CUDA_BIN后调用`freeze`、`run`、`verify`。输出目录要求新建，不能重用消耗过的attempt。

## 完成与恢复记录

实际独立目录为ROG上的`E:/line9_basal_pair_20261008_r1`。H1完成91.5s；首轮H0在会话上下文整理期间因120s心跳过期被监督器停止，属于本次监控配置误触发，未认证为完成，也未覆盖原attempt。`resume_line9_basal_pair.py`先核对H1已完成事件及native哈希，逐字节复制H0科学输入到新目录，仅运行新的H0一次，完成82.813s；H1不重跑。恢复心跳上限改为600s，仍保留取消文件、30min单道上限、资源门禁和独占锁，没有改求解/材料/快照配置。

累计3次求解启动，2份不同模型完成，1次中断保留。两次完成的进程树峰值RSS约6.96GiB；实际求解GPU利用率曾达100%、显存约7.7GiB。原生版本4.0.0、Ez/源样本float64、13569样本、dt=5.896635841874211e−11s，与历史FP32参考的网格、位置、时钟及实际源逐点一致。两份各256帧、六分量float64/有限值/形状/时间/逐文件哈希核验通过，最终恢复验证同时包含旧H1和新H0。完整快照约10GiB留ROG，不入Git。

原始合同SHA256为`8ea9118f43c60bc8ba5c78a4934b78446cc8174b51af9eaa1cc9c74f317abfc0`；恢复合同`ba2857169e7b1af735786a0a6226d1a48fbb11d46ef3db25b165475bfb6fed18`；两模型联合验证合同`f30be51909bd066f29fa458b854764b26bad5c0193cb3ad527d6d93959e12ede`。独立部署源码基于973eed5，远端原仓库未提交资料与分支保留。

## 此次明确缩小的原因范围

在预先按完整Debye一次反射近似定义的底砂332.31–356.31ns窗内：

| 量 | Hann | Blackman |
|---|---:|---:|
| H1总场与底砂主反射模板的复相关 | 0.275 | 0.453 |
| H1−H0与模板的复相关 | 0.962 | 0.976 |
| 模板峰 / ns | 344.31 | 344.31 |
| 差场包络峰 / ns | 342.65 | 342.65 |
| 差场L2 / H1总场L2 | 0.319 | 0.340 |
| H0剩余场L2 / 差场L2 | 3.70 | 3.63 |
| H0与差场复内积相位 / ° | 172.08 | 173.04 |
| FP64−原FP32的L2 / 材料差场L2 | 0.00557 | 0.00592 |

**该站底砂相关响应确实存在，总场的强纹理却主要由仍留在H0中的其他回波决定，并与底砂差场相消。** 去掉底砂后，约335ns的较大包络峰仍在；差场则在约343ns呈现更接近底砂预测的波形。复内积相位是整个窗内的相对关系，不是每个时间样本/频点都固定172°。这些L2比是指定窗内幅度范数，不是SNR、杂波能量占比或整线探测率。

覆盖层底与首砂窗的材料差场/H1仅约8.8e−6和2.2e−5（Hann），符合此次主要改动深部区域的定位；它们小于该处FP64−FP32差异，不能把微小浅时差分当作可靠目标事件。底砂窗的FP64−FP32配置差异仅为底砂差场约0.56%，本次不支持“换成FP64就能修好主要形态”。该比较同时含构建/舍入与被动快照配置差异，没有隔离出通用FP32误差底限。

本次支持**底砂响应被其他相干回波覆盖并部分相消**，尚未把H0的较强回波独立归因给层间多次、地形散射或PML。材料替换包含侧/底PML内的相应材料延续，差场还包含传播和相互作用；不能称为纯底界一次反射、实测clean或独立PML认证。站位来自前期诊断选择，预测窗也利用已知模型；没有盲测/整线/实测验证。下一步应围绕H0在335ns附近的剩余事件做路径及边界控制，不能直接据本图去多次或删除生产信号。

## 图件、数值与复现

聚合目录：`artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4/`。

- [现有142道B-scan及本次站位](../../artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4/existing_bscan_station_marked.png)：仍是历史FP32整线，红线为本次第8道，未假造新整线。
- [中文单站灰度对照](../../artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4/single_station_grayscale.png)：H1/H0/复数差场/FP64−FP32；全窗、地下段与底砂段分行，各行四列共用色标。
- [包络对照](../../artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4/paired_complex_envelopes.png)：先作复数差分再取模；Hann地下段及Hann/Blackman底砂段，共同绝对幅度。
- [波场GIF](../../artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4/wavefield/wavefield_pair.gif)：原场、对照与差场；96帧取自256帧，原存0.1m、显示0.2m，时间并非均匀；两总场全程共同固定对称对数色标，差场独立固定色标。它展示全带impulse场，不能把稀疏快照当20–170MHz窄带相位证书。
- 同目录`wavefield/wavefield_120ns.png`至`wavefield_790ns.png`为七个固定时刻，附H1地质线、收发点及PML接口。对照列仍叠H1地质线供定位，不表示H0仍有该底界。
- [处理数值](../../artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4/analysis.json)、[执行及独立复算](../../artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4/execution_summary.json)、`native_verification.json`及两份execution日志给出源/原生/快照/图件哈希、完成与失败事件。私有501点复频谱/两窗复数时域存于本机`artifacts/local_checks/2026-10-08_line9_basal_pair_results_r1/paired_sfcw_r2.h5`。

精确7个非FFT频点独立DFT相对误差约1.2e−12；独立复数直接求和逆变换误差约5.5e−14，H1=H0+(H1−H0)复数组重组误差小于2e−21。仅属处理与执行身份验证；没有新增物理误差预算。底连通mask/夹层保留/异常拒绝检查及三个执行脚本编译检查通过；原核心代码未改，不重复已完成的182项回归。

复现入口为`scripts/analyze_line9_basal_pair.py traces --help`及`wavefield --help`；用本机/ROG所列私有文件作参数，输出目录要求全新。`scripts/audit_line9_basal_delivery.py --private artifacts/local_checks/2026-10-08_line9_basal_pair_results_r1 --numerical artifacts/local_checks/2026-10-08_line9_basal_pair_results_r1/paired_sfcw_r2.h5 --analysis artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4 --wavefield artifacts/research_checks/2026-10-08_line9_basal_pair_results_r4/wavefield`可复核汇总，不启动求解器。早期可视化草稿r1–r3保存在本机忽略目录，最终图为r4；没有追加求解来重画图。
