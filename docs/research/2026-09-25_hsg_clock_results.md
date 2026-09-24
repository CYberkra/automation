# HSG 三维时钟对照结果

本轮是完整三维 6×6×6 m 全空气数值验证，不是二维剖面，也不是已经建立的滑坡模型。用户要求的 15 m 飞行高度、0–20 m 地层目标属于后续地质模型；当前小对照不包含这些尺度。

## 新证据

已按[结果前设计](2026-09-25_hsg_clock_design.md)完成两个 CUDA double 对照：ratio1 用时 11.547 秒、ratio5 用时 36.125 秒，均 exit 0，无重试。Job 峰值提交内存分别 2,493,296,640 和 2,961,784,832 bytes；不是显存峰值。本轮没有记录可用 GPU 峰值采样，以启动检查和官方求解日志确认 CUDA 路径。

|细化倍率|170 MHz 内部相位差|源码时钟预测|全频带最大预测残差|
|---|---:|---:|---:|
|1|约 0°|0°|2.49×10⁻¹²°|
|3（复用旧结果）|−1.93232°|−1.96435°|0.05586°|
|5|−2.32587°|−2.35722°|0.05742°|

ratio1 内外接收的最大复数相对差均低于 1.2×10⁻¹³。ratio5 外接收最大变化 0.03870 dB / 0.03829°；内接收最大变化 0.06235 dB / 2.32587°。时钟假说的倍率变化得到支持，比笼统称为“加密网格误差”更具体。

原版 hsg_1/hsg_2 的符号执行显示输出计数连续，但耦合 precursor 时间相对标称细网格时间落后 (r−1)/(2r) 个主时间步。170 MHz 测得等效滞后：ratio3 31.57380 ps、ratio5 38.00446 ps；预测分别 32.09722、38.51666 ps。预测不含拟合参数，原始频谱未平移或相位校正。

## 限制与下一步

这支持时序机制，不单凭小型全空气实验判定官方实现有缺陷。ratio1 自动采用不同的过滤/PML/插值设置；细化还改变空间误差。普通参考相对无限空气仍有有限域/有限窗差异，三组比较不等于绝对精度认证。官方高层 direct 函数不接收混合 dt，因此内部使用官方底层 DFT 的项目诊断组合；相同 dt 路径已逐点核对 direct。

下一步核对 HSG 原设计对“滞后细网格”的约定，尤其内部电流源、外部耦合、回馈主网格是否共用一致时钟，再决定是否需要最小内部源对照。不能仅重写接收元数据或移动时间轴就宣称修复；暂未运行局部有损地质候选。当前无求解任务。

## 复现与归档

原始输出、输入、授权快照、日志、终态监督和 SHA256 分别保存在 `artifacts/research_checks/2026-09-25_HSG_ratio1/`、`2026-09-25_HSG_ratio5/`。本轮设计门禁和 attempt 已消耗，不重复运行。`live_status.json` 是历史快照，终态看 `supervision.json`。

分析保存在 `artifacts/research_checks/2026-09-25_hsg_clock/`，含符号追踪、两组分析、三倍率比较和 replay。两组分析复算 JSON/CSV/PNG 字节一致、NPZ 数组逐值一致。分析器仅增加输入路径/倍率选项，旧默认行为不变。输入准备时一次 GBK 解码失败及门禁拒绝日志已留档，该阶段未运行 FDTD。

使用 GPU 环境 Python（仅后处理使用 CPU）：

```powershell
python scripts/analyze_hsg_smoke.py --ratio 1 --hsg-path artifacts/research_checks/2026-09-25_HSG_ratio1/HSG_ratio1.h5 --output <new-ratio1-directory>
python scripts/analyze_hsg_smoke.py --ratio 5 --hsg-path artifacts/research_checks/2026-09-25_HSG_ratio5/HSG_ratio5.h5 --output <new-ratio5-directory>
python scripts/compare_hsg_clock.py
```

比较脚本读取归档路径，写入 comparison.json；不运行求解器。安装版求解器与源代码未修改。
