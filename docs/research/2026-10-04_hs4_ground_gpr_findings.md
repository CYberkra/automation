# 地面雷达对照批次：同 36m 模型天线落地（z=12.025m）13 站稠密 B-scan（2026-10-04）

> 机制诊断证据，不是物理验收、网格收敛认证或实测泛化结论。状态：
> `COMPLETED_DIAGNOSTIC_NOT_PHYSICAL_ACCEPTANCE`。

## 任务与授权

用户 2026-10-04："再跑个地面雷达的吧"。在已完成的 15m 航高 13 站稠密 B-scan
（`2026-10-04_hs4_large_domain_dense_r2`）基础上做**地面耦合对照**：同一 HS4T2D 36m
域模型、同 13 站位（源 x=14.6+0.5k，rx 偏移 +1.3m）、同 2.5cm X/Z 网格、同 Debye 材料、
同 PML、同 600ns impulse、同原生 float64、同源归一化 501 点 20–170MHz 处理链；
**唯一改动是天线高度 z 27m → 12.025m**（地表 z=12m 上方一个网格步）。

## 执行与验收

- 冻结契约 `artifacts/research_checks/2026-10-04_hs4_ground_gpr_r1/execution_contract.json`
  （26 组=13 站×起伏/全覆盖层），preflight PASS：逐组输入重生成比对、raster 材料图与
  36m 归档原图逐元素相同、间距 [0.025,0.05,0.025]、形状 [1440,1,1320]。
- 求解：V4.0.0/CUDA 原生 float64，26/26 COMPLETED，墙钟合计 327.485s（约 12.6s/道），
  最大自有 RSS 0.666 GiB（预算上限 2.35 GiB）。
- 完成验收 PASS：逐道 H5 版本 4.0.0/float64/有限值/dt/Iterations 与 15m 参照一致、
  收发格点位置与新 tx/rx 一致、源激励样本与归一化属性与参照逐字节一致、
  vtkhdf 全材料图与 raster 一致。
- 独立复核：全部 52 条响应（26 新 + 26 归档参照）实际源直接 DFT 七点复算
  相对 L2 ≤8.4e-13（门禁 1e-9）。

## 结果（`2026-10-04_hs4_ground_gpr_results_r1/summary.json`）

| 指标 | 地面（z=12.025m） | 15m 航高 |
|---|---|---|
| 界面对比包络峰时与真实起伏双程时 Pearson 相关 | **0.9358** | 0.6652 |
| 峰时跨度 ns（真实起伏跨度 18.41ns） | 14.14 | 21.62（含边缘离群） |
| 峰时−双程时中位偏差 ns | −3.70 | −1.05 |
| 界面对比最大幅度（固定窗） | 22.24 dB 增益（地面/15m） | 基准 |
| 动态范围 total/contrast dB | **15.25** | **48.77** |

图 `ground_vs_15m_bscan.png`：2 行（total/contrast）× 2 列（地面/15m），每行跨两列共享
symlog 标尺，叠真实起伏双程时曲线（忽略色散的直射线估计，覆盖层 εr=18.017）。

## 解读（机制层面，限定条件内）

- 落地后界面对比条带**显著增强（+22.2dB）且逐站峰时紧贴真实起伏双程时曲线**
  （相关 0.936）；15m 时同一条带弱约 22dB、峰时与地形相关仅 0.665，边缘站位峰跳到
  ~190ns 晚事件。
- 动态范围从 48.8dB 降至 15.2dB：地面直达/地表耦合波仍主导 total，但界面对比只低
  约 15dB，共享色标下直接可见——与 15m 航高下"界面对比淹没在直达波下约 49dB"
  的既有结论定量一致，把"抬高天线损失约 22dB 界面回波相对幅度"落实为单因素测量。
- 这与此前结论链一致：地下传输无问题；高航高的损失来自空气路径衰减与宽照射叠加，
  不是模型或边界错误。

## 边界与限制

- 二维线源、36m 域、13 站段；地面天线是"距地表一个网格"的赫兹偶极，**不是真实地面
  耦合天线模型**；无损耗/色散变化（与 15m 共用材料）。
- 峰时窗为预声明固定窗（地面 40–160ns、15m 160–220ns），包络峰≠到时 onset，
  中位 −3.7/−1.05ns 偏差属显示尺度，不作物理到时结论。
- 不认证网格收敛、唯一物理归因、三维/实测泛化。vtkhdf/cuda_cache 留本机
  （.gitignore 通配 `*hs4_ground_gpr*`），哈希台账在胶囊验收文件。
- G4/训练/保留实测/首版范围不变。

## 复现入口

```bash
# 冻结（已完成，勿重复冻结同目录）
python scripts/hs4_ground_gpr_dense.py freeze --out artifacts/research_checks/2026-10-04_hs4_ground_gpr_r1
# 求解（需完整 vcvars64+CUDA 环境，见 scripts/run_hs4_ground_gpr_dense.cmd）
scripts\run_hs4_ground_gpr_dense.cmd run --execute --out artifacts\research_checks\2026-10-04_hs4_ground_gpr_r1
# 分析与图
python scripts/analyze_hs4_ground_gpr_dense.py --capsule artifacts/research_checks/2026-10-04_hs4_ground_gpr_r1 --out <new_results_dir>
```
