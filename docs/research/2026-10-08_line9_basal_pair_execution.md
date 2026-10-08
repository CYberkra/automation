# 第8道底砂配对与波场快照：SSH执行记录

2026-10-08。用户明确授权“通过ssh连接我们那个机器然后开跑吧”。本批限两次单道V4.0.0/CUDA原生FP64求解，串行；不扩道、不改频带/网格，不混用三包作为纯航高控制。

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

## 当前执行状态

已完成本机准备和SSH只读资源/身份核查，尚未将启动记为事实。启动后追加实际execution.jsonl/原生输出审计；完成后交付中文单站对照、现有B-scan标位和波场GIF，缺道仍缺，不生成假整线。
