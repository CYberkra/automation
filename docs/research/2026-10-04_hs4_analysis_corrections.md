# HS4 角谱与形态指标修正、CPU复算和接续

用户在[独立复审](2026-10-04_hs4_independent_height_attribution_review.md)后同意推进。本单元修正两项分析、复算已归档资料并检查后续执行资源；没有新增求解。历史结果与冻结胶囊原样保留，以下新证据限定旧文档的“61.5%功率”“1.58m硬底”“细节全部丢失”表述。

## 实现

- `analyze_hs4_height8m_surface_filter.py`逐频使用`|kx| <= 2πf/c_air`；结果命名为电场谱平方统计，不当成上下行能流或全反射比例。移除固定PASS及固定13.6°精确归因，明示有限时空窗、泄漏、不同高度波包和无方向分解的限制。
- 角谱入口要求`--out`新目录，可通过`--capsule`指向保留原始帧的胶囊；核对完成状态、契约哈希、双侧完整帧清单和相同迭代序列，校验所选原始帧哈希/时间/形状/原生float64/有限值及接收器身份。不从GIF重建数值场，不覆盖旧输出。
- 两份迁移脚本共用`hs4_analysis_metrics.py`。提取山脊与真值使用相同空间FFT掩码，低/高频分量分别给相关、RMSE、相对L2误差、幅度范数比和带符号投影增益；保留整曲线的平均深度偏差。1.577855m为诊断分频尺度，取消“测得的硬分辨率”命名与物理PASS。
- 迁移输出新增完整数组与源码/指标源码哈希，供独立计算和有限区间敏感性分析；每项小分量的数值误差预算仍为null，`recoverability_certified=false`。

## 复算结果

同一批归档H5、同一官方20–170MHz/501点复谱链、同一迁移图及山脊提取算法，仅修正评价。新[四情形结果](../../artifacts/research_checks/2026-10-04_hs4_permittivity_matched_metrics_r1/summary.json)如下：

| 覆盖层 | 整曲线相关 | 匹配低频相关 | 匹配高频相关 | 高频相对L2误差 | 高频幅度范数比 |
|---|---:|---:|---:|---:|---:|
| εr6 | 0.940 | 0.950 | 0.479 | 1.053 | 1.060 |
| εr9 | 0.962 | 0.970 | 0.751 | 1.048 | 1.565 |
| εr18.017，无Debye | 0.977 | 0.983 | 0.753 | 0.952 | 1.441 |
| εr18.017，原Debye | 0.843 | 0.850 | 0.614 | 0.838 | 0.895 |

原模型旧高频相关约0.115改为0.614，说明旧指标低值不能证明细节全部消失。但新高频相对误差约84%，低频相对误差约56%，整曲线平均偏深0.179m、RMSE0.221m，仍不能认证准确恢复；相关不能代替幅度、误差或物理资格。无Debye组高频相关0.753伴随44%范数放大和95%相对误差，也不能把较高相关单独称为保真通过。

整曲线相关与旧四情形结果差至多约1.53e−13，RMSE差至多约9.51e−14；差异属于当前CPU环境重建的数值变化，没有改迁移行为。两份修正脚本的原Debye基准山脊及三幅迁移数组逐元素完全一致，官方复谱对既有局部扰动基准的锚点相对L2为0。

[空间分频/区间敏感性](../../artifacts/research_checks/2026-10-04_hs4_profile_sensitivity_r1/summary.json)固定迁移数组，比较分频尺度1.0/1.577855/2.0m及原14.75–21m、内裁15.25–20.5m两区间。原Debye高频相关0.077–0.614、相对误差0.838–1.379，表明小分量受有限区间周期延拓影响明显。该分频不是独立测得的物理分辨率，也没有消除网格/速度/成像算法影响。后续需受控细节扰动和独立误差预算，不能用此表直接选物理“可恢复”阈值。

## 验证与本机限制

- 新12项数组反例通过：完美细节、删除细节、极性翻转、过度放大、深度偏移、常量方向缺失、逐频传播界限/正负kx对称、非法频率/尺度、非有限值及形状不一致。常数增益保持相关但不保持幅度的反例已覆盖。这些是数学/实现验证。
- 182项现有纯数组回归全部通过（`artifacts/local_checks/20261004T062909Z-bcf71737`）；不代表正演/实测验证。
- 实际角谱运行因本机rough/halfspace各599帧缺失被拒绝，未创建结果目录；以历史结果目录为输出同样被拒绝，旧角谱summary SHA-256未变。**角谱代码已修正，真实角谱百分比尚未复算**。
- 当前V4.0.0包与既有细网格包的源码/编译库身份相同，预检`SAME_ARCHIVED_PACKAGE_FILES`。CUDA11.8/nvcc可见，普通PowerShell未激活cl；如运行求解须使用既有MSVC环境。
- CPU分析中预检可用RAM约1.04GiB；结束后约3.09GiB，空闲显存约4.56GiB。已有1cm中心门槛3.75GiB RAM，三维中心8GiB RAM；均未满足，没有冻结或启动新求解，不停用户程序或降低物理/精度条件绕过预算。
- [本单元验证与接续状态](../../artifacts/research_checks/2026-10-04_hs4_analysis_correction_checks/summary.json)。G4、训练、保留实测与首版处理范围未改变；迁移仍仅为诊断。

## 接续命令

下列输出目录必须不存在。在保留8m原始快照的机器同步Git后，先做无需求解的角谱重算：

```powershell
python scripts/analyze_hs4_height8m_surface_filter.py --capsule artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c --out artifacts/research_checks/target_hs4_aircone_corrected_r1
```

这一步只需NumPy/h5py/Matplotlib；不依赖GPU或导入求解器。先复核帧完整性和源码哈希，再看逐频电场统计；要进一步认定全反射与多次波，仍需要E/H对齐、方向分解和法向通量，而非只读新的百分比。

本机CPU复算可在同样官方V4 SFCW源码身份的Python环境执行：

```powershell
python scripts/check_hs4_analysis_metrics.py
python scripts/migrate_hs4_permittivity_scan_v0_1.py --out artifacts/research_checks/target_hs4_matched_metrics_r1
python scripts/analyze_hs4_profile_metric_sensitivity.py --arrays artifacts/research_checks/target_hs4_matched_metrics_r1/migration_arrays.npz --out artifacts/research_checks/target_hs4_metric_sensitivity_r1
```

随后按[已有跨电脑任务](2026-10-04_hs4_cross_pc_task.md)在足够资源的机器接续网格及三维控制，先预检、按阶段另冻契约。不同构建须先replay；历史attempt不重跑，已完成对照不因本次改指标重新求解。
