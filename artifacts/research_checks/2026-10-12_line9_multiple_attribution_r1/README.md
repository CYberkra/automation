# 2026-10-12 Line9 强斜带归因（A1+A2，纯分析）

- 报告：`docs/research/2026-10-12_line9_multiple_attribution.md`
- 脚本：`scripts/analyze_line9_multiple_attribution.py`（预声明常量在文件头）
- `attribution.json`：41 道逐站几何/拾取/残差/路径比 + 斜率行 + 输入指针
- 图：`bscan_multiple_overlay.png`（双面板叠加）、`multiple_residual_ratio.png`、`slope_comparison.png`
- 结论：强斜带全线属覆盖层多次波族（斜率比 0.999±0.056），第二程侧偏缩短，校准后到时提前 10.47±0.54 ns。
- 复现命令见报告；不调用求解器。
