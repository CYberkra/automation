# PEC镜像参考与提前响应排查

## Material Passport
2026-09-24；AI辅助确定性数值验证；用户授权自主设计与GPU仿真。分析证据已归档，不涉及统计推断、实测调参或现场材料认定。状态：单次求解完成，解析比较完成；网格收敛及介质提前响应因果验证尚未完成。

## 结果

[预声明设计](2026-09-24_PEC_design.md)：仅把原介质半空间改为内置PEC，保持15m高度、横向1.3m基线、X源/Ex、5cm网格、24m立方域、800ns记录和20格HORIPML。复用空气T800，不再求解空气。CUDA double 8310步正常完成，447.078秒；主机作业内存峰值22431076352 bytes，低于32GiB。一次GPU采样99%/10108MiB，不是整程峰值。原始10项检查通过，目前无求解任务。

独立参考是无限PEC平面的反向水平镜像电流元，镜像距离sqrt(30²+1.3²)m，包含近场项，无拟合。该规则见[Purdue ECE604 Lecture32](https://engineering.purdue.edu/wcchew/ece604s20/Lecture%20Notes/Lect32.pdf)印刷第316页；场表达式沿用已核对的[Cornell Lecture28](https://courses.cit.cornell.edu/ece303/Lectures/lecture28.pdf)。

| 800ns记录尾渐消 | 最大复数相对误差 | 中位复数相对误差 | 最大幅度差 | 最大相位差 |
|---|---:|---:|---:|---:|
| 无 | 107.68% | 22.80% | 27.83dB | 89.04° |
| 200ns | 9.3954% | 1.6325% | 0.003552dB | 5.3862° |
| 400ns | 9.3954% | 1.6325% | 0.003552dB | 5.3862° |

比较的是PEC总场减空气场后的反射频谱，不是被强直达波主导的总场。501点20–170MHz由官方direct转换；没有Ricker替代或拟合时移。幅度很接近，但最大复数误差仍约9.4%，不能说整体精度高或已经收敛。

## 相位偏差的独立预算

在读取PEC结果前保存了轴向Yee平面波预算：sin(ωdt/2)=(c dt/dx)sin(k_num dx/2)，对近垂直30m空气路径估算传播相位。公式来源[Schneider第7章](https://eecs.wsu.edu/~schneidj/ufdtd/chap7.pdf)7.4节。

170MHz预计相位滞后5.41755°，实际5.38619°；全频带观测相位与该无拟合预算最大残差0.03135°。这支持网格色散是PEC相位偏差主要来源，但轴向平面波预算不是完整三维点源离散解析解，也没有代替网格细化验证。交付频谱未作相位修正。

## 提前响应的范围已经缩小

- PEC减空气在80ns前逐样本为零；εr9介质减空气在0–80ns峰值0.04820V/m，约77ns。
- 两者完整差分去均值/Hann周期图主峰均约5.16756GHz。共同带外振荡并不能解释为何只有介质对照有提前分量。
- 官方Hann频率窗、8倍补零重构后，PEC和介质包络峰都在100.632ns，解析PEC峰99.800ns。重构网格间隔约0.832ns，峰值取样有量化；补零不提高物理分辨率，带限前旁瓣不能作为违反因果性的证据。这项波形比较不证明介质半空间幅相准确。

进一步只读本地V4实现并核对与已安装Python源hash一致：`grid/fdtd_grid.py`中`_calculate_average_pml_material_properties`对整个侧面材料截面取平均，`pml.py`中`calculate_sigmamax`含1/sqrt(er*mr)。小规模调用已安装官方`pml_average_er_mr`（未创建仿真网格、未调用求解器）得到：

| 下方4m材料 | 侧截面平均εr | 自动sigma_max相对空气 |
|---|---:|---:|
| 空气 | 1 | 1 |
| PEC（内置名义εr=1） | 1 | 1 |
| εr9介质 | 7/3 | 0.654654 |

所以空气/介质两算例并非只改变物理地层，还隐式改变了整个侧面PML自动参数。源接收位置到侧PML内外界面的几何反射路径约73.05–80.39ns，与77ns现象相符。**这是源码和时间对应支持的候选原因，尚无固定PML参数对照，不能写成已证实因果或V4程序缺陷。** 不修改gprMax源码，也不把原始提前分量裁掉以制造好结果。

## 接续与复现

当前推荐下一项：先核对官方固定PML参数接口，构建空气/介质共享同一明确PML参数的有限对照，以隔离自动参数变化。这样比直接扩大地层扫描更能回答当前问题。之后再决定横域或网格细化；5cm网格的长路径相位误差已明确需要预算。当前不生成训练集或宣称20m可探测。

完整原始归档：`artifacts/research_checks/2026-09-24_PEC_T800/`；分析：`2026-09-24_PEC_analysis/`；静态预测与源码探针：`2026-09-24_PEC_sources/`。各目录record.json含文件/脚本hash。当前attempt已消耗，不重跑。

不求解复现（输出目录/文件须不存在）：
```powershell
& artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe scripts/analyze_pec_reference.py --pec artifacts/research_checks/2026-09-24_PEC_T800/PEC_T800.h5 --air artifacts/research_checks/2026-09-24_T800/M00_x_3d.h5 --dielectric artifacts/research_checks/2026-09-24_M01_T800/M01_T800.h5 --output artifacts/local_checks/PEC_replay
& artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe scripts/budget_yee_phase.py --output artifacts/local_checks/phase_replay.json
& artifacts/local_checks/gprmax_v4_gpu_env/Scripts/python.exe scripts/audit_pml_material_average.py --output artifacts/local_checks/pml_average_replay.json
```

phase_budget_comparison.json复算方法：CSV取taper_ns=400，angle(PEC_reflected/theory,deg=True)减去phase_budget_before_result.json的predicted_response_phase_error_deg，取绝对值最大；不拟合也不修正输出。
