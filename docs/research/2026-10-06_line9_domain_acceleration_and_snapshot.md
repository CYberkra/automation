# 九号线域宽加速与大域被动快照

用户提出“太大了，我们能调整模型在不影响质量的前提下加速吗”，随后明确要求“大域的等会再跑一个带波场快照的版本”。本单元保留已有完整大域结果，准备只改变水平计算域的候选；优先追加一条完整大域被动快照。**候选已制作不等于缩域质量已验证；只检查总场全图误差不能保护被直耦掩盖的弱地下回波。**

## 执行中断与可复用基准

原五站契约`line9_2d_pilots_rog_r1`已消耗attempt。原生执行日志确认`full2d_pilot_0`（X220）完成，631.171s、进程树RSS峰值10,739,511,296字节，原生H5 SHA-256 `2b9af6e225437d9516d6bdf8b358bd09a224729cead0e39b3560442ded2e7689`。第二站已启动后用户关机，随后重启且IP改变；重新连接核对保留文件，只有第一站原生输出，未完成第二站/其余三站，也未产生整批completed_verification。

不覆盖旧契约/执行日志，不在旧attempt原地重跑。原工作目录保持旧代码身份；新快照在独立worktree接续、复用既有V4环境与共享GPU锁。快照改选已完成X220（不再使用未完成X196.75），将主接收Ez和原生源samples与已有结果逐元素核对。只读completed-group审核视图明确标记`READONLY_COMPLETED_GROUP_AUDIT`，不能作为可执行契约，更不能补造“五站完成”记录。

## 优先加速因素：水平移动域

已制作200×75m、160×75m两种宽度，在五个既定试点各有一份输入，共10份；私有包`artifacts/local_checks/2026-10-06_line9_domain_acceleration_r1`。由原始完整H5直接按整数格切片，所有保留体素逐格复验一致；域窗口限制在原X−50–350m以内，超界时整体平移窗口，不发明延续地质。收发器坐标作相同整数格平移，其绝对位置和Yee网格相位保持；材料、75m高度、2.5cm网格、2m HORIPML、15m航高、1200ns、原生FP64、源电流元与501个20–170MHz频点均保留。

| 水平域宽 | 单元数 | 相对原域工作量 | 理想工作量倍率 | 以首道外推的单道 / 99道 |
|---|---:|---:|---:|---:|
| 原400m | 4800万 | 100% | 1 | 10.52min / 17.36h |
| 200m | 2400万 | 50% | 2 | 5.26min / 8.68h |
| 160m | 1920万 | 40% | 2.5 | 4.21min / 6.94h |

外推由[计算脚本证据](../../artifacts/research_checks/2026-10-06_line9_domain_acceleration_r1/workload_estimates.json)记录，**不是缩域实测耗时或保证**。PML/导入/编译和GPU利用率会改变倍率。与网格放粗相比，这个因素避免同时改变空间分辨率、时间步长和地层栅格，但会移近人工边界并删除远处地形/地层，不能未经对照称为无损。

建议先在有完整400m基准的站位对照160m；边界差异不可接受时退到200m或400m。保持源/接收时间原点和绝对幅相，比较501点复频响、两种窗的复轮廓及带符号波形，不允许逐道归一化、后移相位、拟合增益再宣称一致。全记录和预声明的地表后弱事件/晚时段分别比较，不能让直耦主导的全图相对L2掩盖地下差异；低于可信误差底的弱窗标为不可判定。具体可接受误差和保护窗须在新对照契约执行前冻结，不能用事后ROI或未声明eps调整结论。

不在本轮顺带修改网格、精度、材料损耗或记录时窗。5cm网格可能降低二维单元数和步数，但不具备现成的无损依据；缩短时窗还会改变源归一化复谱，不能因为B-scan只展示前半段就认为安全。官方依据：[FDTD网格、CFL及边界说明](https://docs.gprmax.com/en/latest/gprmodelling.html)；当前安装仍固定4.0.0，未升级到文档站显示的4.0.1。

## 用户授权的大域快照

新私有准备包`artifacts/local_checks/2026-10-06_line9_large_snapshot_prepared_r3`，参考`full2d_pilot_0`。原始输入按字节保留，**仅追加**599条快照命令，不改变求解网格或物理因素。保存整域400×75m、间隔34个原生时间步（约2.005ns），从0到约1198ns；保存间距XYZ=0.2/0.2/0.025m，每帧2000×375×1。原生H5仍输出完整FP64 Ez与源samples。

安装V4的hash快照语法只接受11参数，默认六分量；其CUDA `dtoh_snapshot_array`即使最终只查看Ez也保留六个host数组。因此预算按**六分量**约20.08GiB历史计算，不误按一分量估计。V4的`utilities.host_info.mem_check_device_snaps`在不含快照的模型能装入显存时，自动设置内部`snapsgpu2cpu=True`，使GPU只持有一帧约36MB附加缓存；无需命令行选项或修改求解器。最低可用RAM约40GiB、空闲VRAM约11.7GiB、额外磁盘约26.1GiB，最终以生成snapshot_plan为准；监督器自有RSS上限40GiB、剩余系统RAM下限1.5GiB、单道/全批60min、no_retry。重启后只读预检为可用RAM53.52GiB、空闲显存15104MiB、E盘空闲1362.96GiB，实际freeze/run再次核查。共享GPU锁沿用原`E:\automation_djh\artifacts\local_checks\hs4_gpu_exclusive.lock`。

首次快照胶囊`line9_large_snapshot_rog_r1`使用了V4不存在的`-snapsgpu2cpu`命令行选项，argparse退出，未加载模型、未执行FDTD、没有原始输出。失败契约、源码worktree与日志保留；新启动必须使用独立worktree和新r2胶囊，不能原地重试。该错误不代表模型或容量失败。

快照20cm/约2ns仅为保存采样，2.5cm FDTD求解不变。图用于传播路径与大尺度波前诊断，不能认证源全部高频分量/精细地下相位；GIF再隔3帧播放，必须标注动画采样。灰度单站B-scan保留未计算区域；GIF标原始Ricker瞬时总场、共享SymLog色标、真实地表/岩性接触线、收发位置和时间，不能把瞬时绝对值称解析包络，不能把脉冲波场称SFCW波场。实际回波身份仍需路径与对照证据，不能仅据动画颜色判断“PML反弹”。

交付：原生被动核对、599帧E/H时标/原始dtype/形状/有限值/哈希、官方复源归一化矩形/Hann SFCW单站图、GIF及聚合报告。大体积原生快照留目标机；GIF超出合理Git体量时留本机并说明路径，配套静态帧与B-scan图入库。

## 复现入口与已做检查

`scripts/prepare_line9_domain_acceleration.py`只准备精确切片；`scripts/line9_large_domain_snapshot.py`分prepare/freeze/run/verify；`scripts/analyze_line9_large_wavefield.py`核查后生成实际GIF和SFCW报告。监督器拒绝附加求解参数，快照转存使用冻结版本V4的原生自动机制；SFCW分析器支持快照审核及正确的一站留白图。**不要把新监督器/分析代码覆盖到旧冻结worktree，旧源码哈希必须保留。**

检查见[准备和回归证据](../../artifacts/research_checks/2026-10-06_line9_domain_acceleration_r1/)：10个切片完整体素、材料键、空气收发点及坐标平移核对通过，V4原生命令均解析；599帧命令解析与内存/磁盘计算完成；原生/SFCW数组回归与5项被动观测、错误极性、FP32快照、E/H半步、违规输入修改的反例通过。均为准备/数组证据，不是缩域回波质量签认或已完成大域快照。

目标机新worktree中，设置与既有任务书一致的真实Python/MSVC/CUDA路径后：

```powershell
# 生成独立准备包；旧大域模型包保持只读。
scripts/run_line9_2d_v4.cmd scripts/line9_large_domain_snapshot.py prepare --package E:/line9_rog_first_20261006/artifacts/local_checks/2026-10-06_line9_material_model_r1 --out artifacts/local_checks/line9_large_snapshot_prepared_r3 --evidence artifacts/research_checks/line9_large_snapshot_prepare_target --reference-id full2d_pilot_0
scripts/run_line9_2d_v4.cmd scripts/line9_large_domain_snapshot.py freeze --package artifacts/local_checks/line9_large_snapshot_prepared_r3 --out artifacts/local_checks/line9_large_snapshot_rog_r1 --pilots E:/line9_rog_first_20261006/artifacts/local_checks/line9_2d_pilots_rog_r1
scripts/run_line9_2d_v4.cmd scripts/line9_large_domain_snapshot.py run --out artifacts/local_checks/line9_large_snapshot_rog_r1
scripts/run_line9_2d_v4.cmd scripts/analyze_line9_large_wavefield.py --study artifacts/local_checks/line9_large_snapshot_rog_r1 --out artifacts/research_checks/line9_large_snapshot_rog_r1_analysis
```

先完成本次大域快照，再决定缩域对照；99站/391站和3D没有追加启动。未经新的缩域执行契约，不用完整400m运行器直接跑这些不同形状输入。
