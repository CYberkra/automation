# 2026-10-03 SFCW 仿真与后处理链独立审查及修复记录

> 审查人：kimi（应用户要求独立审查最近仿真模型及后处理链路，对照 gprMax v4 官方手册）。
> 范围：benchmark3d_r2_co 39 道、verify_runs 近期批（halfspace_standard / hs4t2d / iface_sensitivity / bf2d_scan）、官方 SFCW 工具箱后处理链、sim2real 增广链。
> 依据：本地官方树 `E:\gprMax-v.4.0.0\gprMax-v.4.0.0`（`docs/source`、`gprMax/toolboxes/SFCW/`），官方 `processing.py` 当前 SHA-256 仍为钉住值 `adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b`（全树 24 条登记哈希逐一重算一致）。未运行求解器、未消耗 attempt、未触碰冻结族 {C5, C8}。

## 1. 审查结论总览

仿真模型与官方 SFCW 处理链本体**无数值错误**，可作证据基础。修复了 2 个高优先级口径问题与 2 个低优先级脚本缺陷；2 个 legacy 批次建模缺陷不可追溯修复（attempt 一次性、批次已冻结），按 §4 声明登记。

### 核验通过项（摘要，均有文件/数值证据）

- 全部受审 `.in` 为 gprMax v4.0.0 合法语法；`#waveform: impulse` 源谱在 20–170 MHz×501 全点高于 −100 dB 底限，符合官方 SFCW 工具箱"内置 impulse + 源启动 t=0"推荐（`waveforms.py:157-162`、SFCW README）。
- 几何网格对齐无半格错位；39 个 `.in` 道位置与 `sample_table.json` 全量一致；界面场按生成规则重算 0 处不符。
- 频网格 `linspace(20e6,170e6,501)` 500 个差值严格相等（官方 `rtol=1e-10` 必过）；taper 公式 `(round(200ns/dt)−0.25)/N` 两链等价（实际锥化 199.99 ns）；1200 ns 截断对 600 ns 记录为 no-op。
- 载波=`frequencies[0]`=20 MHz、`real_bandpass=2·Re(…)`、IFFT×zpf 保幅、Hann 均值归一——与官方实现逐字一致；平坦 H=1 重建峰值精确为 1（四组合 atol 1e-15 实测）。
- 文档对官方代码/手册的 24 条哈希与全部函数签名/默认值/公式声明核验成立；`sfcw_r1`/`sfcw_r2` manifest+npz 字节一致；carrierfix 回归 17/17 通过（修复后复跑仍 17/17）。

## 2. 已修复问题（本记录同日执行）

| # | 问题 | 修复 |
|---|---|---|
| P1（高） | sim2real 增广主链输入为 legacy 95 MHz 载波（`run_reward_protocol_b2_pilot_v0_1.load_bscan`），90–125 MHz 标定带对仿真/实测落在不同频带 | 新建 `scripts/augment_sim2real_v0_7.py`：v0.6 算子逐行不动，输入改经 `sfcw_official_loader_v0_2.load_bscan_both`（官方 20 MHz 载波），legacy 路径同种子重跑仅作连续性对照；重标定产出 `augment_v07_metrics.json`（指标见 §3）。v0.6 指标保留为历史记录，registry 标记 superseded |
| P2（高） | `plot_sfcw_bscan.py` 实为 house FFT 带通合成（无源谱归一/非 501 网格/Hann 未均值归一），图题却写"官方 SFCW 链路产物" | 图题/docstring 改为"house FFT 带通合成（非官方 SFCW 工具箱）"，注明幅度口径不用于定量，官方链产物指向 `artifacts_check/sfcw_r1` |
| P3（中） | `render_sfcw_compare.py` 色标参考（0–50 ns 直达窗）与打印指标参考（0–400 ns 全局）不一致 | 统一为逐道 0–50 ns 窗 max\|env\|，docstring 钉死 "dB re direct" 定义；本数据下数值不变（直达峰即全局峰，−81.1/−70.1/−68.5 dB 复跑一致） |
| P5（低） | `render_v6_ringing.py` 直达窗相关用 13 道峰位并集掩码按道计算；`t0_of` 定义未用 | 改为逐道掩码（直达窗相关复跑 1.0000/1.0000，符合起始门设计预期）；删除未用函数 |
| —（低） | `analyze_direct_fidelity.py` 未声明 legacy 载波口径 | docstring 补登记：相对保真指标同口径成立，跨链比较带载波偏差，官方口径见 v0.7 |

## 3. v0.7 重标定结果（官方载波，种子 20261001 不变）

- Line9 实侧标定不变：q p25–p75 = 0.053–0.091（中位 0.0713）。
- 直达保真门通过：corr_med 0.9950 / corr_min 0.9907 / 峰值变化 0.00% / 到时漂移 0；β 触顶 0/33。
- 横向纹理 diff_rms：v0.7 = 0.163（实测 0.179；v0.6 为 0.190）。
- 包络剖面 v7/实测：1.19× / 0.97× / 0.84× / 0.88× / **1.82×**（t0+140~200 ns 尾窗超出 0.8–1.25 目标带；该窗绝对量级接近显示下限，与 v0.6 尾窗 1.73× 同性质已知项，未修饰）。
- 产物：`augment_v07_metrics.json`、`fig_augment_v07_bscan.png`、`fig_augment_v07_profile.png`。

## 4. 不可追溯修复项（声明登记，跨批引用时必读）

- **M1**：benchmark3d_r2_co 的 2D/3D 配对界面不完全相同——3D bin 为 5×5 单元平均、2D 取单列，40 个 y-bin 中 10 个相差 1 格（0.05 m），与 `sample_table.json` "differ ONLY by dimension" 表述不符。HS4 起流程已修复（2D 从最终 binned table 提取，`build_hs4_rough_interface_v0_1.py`）。引用该批 2D/3D 配对差分结论时须注明 ±0.05 m 台阶差异混入。
- **M2**：benchmark3d_r2_co 与 BF2D 全部 51 个 `.in` 材料不贯通 PML（x/y 侧与底部 1 m 空气环恰等于 PML 厚度），晚期响应可能含 PML 内界面物理反射。批内对比可用（缺陷系统性抵消）；与 HS 修复版跨批比较晚期响应时须声明。HS1–HS4 已修复。
- 已知取舍重申：5 cm 组 170 MHz 处 cover 内 ≈8.3 单元/λ（低于官方 ≥10 经验，2.5 cm 组为收敛对照）；`pml_cfs` σmax 按 5 cm 空气调优跨网格复用（次优非错误）；2D 线源与 3D 偶极幅度不可直接比（配对限运动学/归一化特征）。

## 5. 待用户确认项（未动）

1. cover 单极 Debye 参数（ε∞=18.017、Δε=7.878、τ=6.4567 ns，σ=0.003 S/m）出处与频段适用性；BF2D clutter 的 Δε 按 εr 缩放而 σ 恒定的依据。
2. 色散材料界面未开 `#dispersive_averaging: y`（逐格硬跳变）是否为预期。
3. `render_sfcw_compare.py` 硬编码 t_surf=100.2 / t_iface=185.1 ns 出处未注明；官方要求记录末 5% 尾部 < −60 dB，手头"窗尾能量 ≤0.04%"口径建议按官方幅值口径复核一次（200 ns taper 与无 taper 敏感性档均已保留，风险可控）。
4. 未来做仿真-实测**绝对幅值**比对前，须统一官方因子 2 与 legacy 1× 幅度口径（当前全为比值/归一化比较）。
5. v0.7 是否取代 v0.6 作为后续增广基线（建议：是；v0.6 保留为 legacy 载波历史记录）。

## 6. 验证记录

- `augment_sim2real_v0_7.py` 于 venv（gprMax 4.0.0）运行完成，官方运行时 SHA 校验通过。
- `plot_sfcw_bscan.py` / `render_sfcw_compare.py` / `render_v6_ringing.py` 修复后复跑，图件按新口径重生成。
- `check_sfcw_carrierfix_acceptance.py` 修复后复跑 17/17 OK（`artifacts/local_checks/carrierfix_acceptance_checks.json`）。
