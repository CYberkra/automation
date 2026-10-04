# -*- coding: utf-8 -*-
"""Generate 3 figures for the cross-domain research report."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "Noto Sans CJK SC", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# Cohesive palette
C = {
    "blue":   "#3B6EA5",
    "teal":   "#2E8B8B",
    "orange": "#D9822B",
    "red":    "#B5443B",
    "purple": "#7A5C9E",
    "green":  "#5C8A4D",
    "gray":   "#8C8C8C",
    "light":  "#F2F4F7",
}

# ---------------------------------------------------------------- Figure 1
# Framework taxonomy: six paradigms observed across fields
fig, ax = plt.subplots(figsize=(13.5, 7.2), dpi=150)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

rows = [
    ("① 端到端直接反演", "数据 → 网络 → 参数场", "GPRInvNet/3DInvNet(GPR) · InversionNet(地震) · CNN-3D-ERT(ERT)", C["blue"]),
    ("② 提取—反演两步走", "数据 → 自动拾取中间量 → 物理反演", "Surf-Net/DisperNet(面波) · 双曲线检测(GPR) · PhaseNet(到时)", C["teal"]),
    ("③ 物理引导学习", "损失函数内嵌正演/PDE 约束", "物理引导DC电阻率反演 · PINN弹性FWI · FF-PINN波场模拟", C["orange"]),
    ("④ 混合式精化", "DL 给初值/低波数 → 传统反演精修", "DNN初值+FWI(仅15%迭代) · 低频FWI+井约束+DL补高波数", C["purple"]),
    ("⑤ 生成式反演 + UQ", "潜空间重参数化 → 多解采样 → 不确定性", "VAE-ERT生成反演 · MC-Dropout · 贝叶斯DREAM(GPR)", C["red"]),
    ("⑥ 绕过中间量", "全频散谱/全波形直接反演", "全频散谱反演(JGR'24) · 全波形FWI(GPR/地震)", C["green"]),
]

ax.text(50, 97, "“采集 — 提取 — 反演”全流程的六类范式框架（跨领域归纳）",
        ha="center", va="top", fontsize=15, fontweight="bold", color="#222")

y0 = 88; dy = 14.2
for i, (name, flow, ex, col) in enumerate(rows):
    y = y0 - i * dy
    ax.add_patch(FancyBboxPatch((2, y-5.2), 20, 9.5, boxstyle="round,pad=0.4",
                                fc=col, ec="none", alpha=0.92))
    ax.text(12, y-0.45, name, ha="center", va="center", fontsize=11.5,
            color="white", fontweight="bold")
    ax.add_patch(FancyBboxPatch((25, y-5.2), 30, 9.5, boxstyle="round,pad=0.4",
                                fc=C["light"], ec=col, lw=1.5))
    ax.text(40, y-0.45, flow, ha="center", va="center", fontsize=10.5, color="#333")
    ax.annotate("", xy=(25, y-0.45), xytext=(22.6, y-0.45),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=2))
    ax.add_patch(FancyBboxPatch((58, y-5.2), 40, 9.5, boxstyle="round,pad=0.4",
                                fc="white", ec="#CCCCCC", lw=1))
    ax.text(59.5, y-0.45, ex, ha="left", va="center", fontsize=8.6, color="#444")
    ax.annotate("", xy=(58, y-0.45), xytext=(55.6, y-0.45),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=2))

ax.text(3, 6.5, "范式名称", fontsize=9, color=C["gray"])
ax.text(27, 6.5, "数据流结构", fontsize=9, color=C["gray"])
ax.text(60, 6.5, "代表工作（领域）", fontsize=9, color=C["gray"])
ax.text(50, 1.5, "我们的面波项目当前以 ② 为主线；④⑤⑥ 是文献中已验证的升级路径",
        ha="center", fontsize=10.5, color=C["red"], fontweight="bold")

plt.tight_layout()
plt.savefig("fig1_framework_taxonomy.png", bbox_inches="tight", facecolor="white")
plt.close()

# ---------------------------------------------------------------- Figure 2
# Metric-family × pipeline-stage applicability matrix
stages = ["L0 数据质量\n(采集/预处理)", "L1 中间量提取\n(频散曲线/双曲线)",
          "L2 反演重建\n(Vs/介电常数)", "L3 工程验证\n(钻孔/独立数据)", "L4 不确定性\n(UQ/校准)"]
families = ["检测分类类\nP/R/F1·mAP·AUC", "图像相似度类\nSSIM·PSNR·MSSIM",
            "数值误差类\nRMSE·MAE·R²", "物理一致类\nmisfit·数据残差",
            "校准类\nPICP·ECE·CRPS", "效率类\n推理时间·加速比"]

# applicability strength 0-3 (source: literature synthesis)
M = np.array([
    [0, 2, 1, 0, 0, 1],   # L0: SNR/RIHSNR(图像类), 噪声测试
    [3, 2, 2, 0, 1, 2],   # L1
    [1, 3, 3, 2, 1, 3],   # L2
    [2, 1, 3, 3, 1, 0],   # L3
    [1, 0, 2, 2, 3, 1],   # L4
])
notes = [
    ["", "RIHSNR", "扰动噪声实验", "", "", ""],
    ["Surf-Net阈值法", "频散谱IoU", "速度差μ/σ", "", "拾取置信度", "80%人工↓"],
    ["缺陷检出F1", "介电常数场", "OpenFWI三件套", "misfit<1", "", "0.59s vs 40min"],
    ["钻孔界面检出", "", "深度误差<10%", "RMS misfit↓50-70%", "80%概率带", ""],
    ["", "", "Vs30 CoV", "覆盖率vs真值", "PICP/MIW/ECE", ""],
]

fig, ax = plt.subplots(figsize=(13.5, 6.4), dpi=150)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("m", ["#FFFFFF", "#C7DCEC", C["blue"]])
im = ax.imshow(M, cmap=cmap, vmin=0, vmax=3, aspect="auto")
ax.set_xticks(range(len(families))); ax.set_xticklabels(families, fontsize=9.5)
ax.set_yticks(range(len(stages))); ax.set_yticklabels(stages, fontsize=10)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if notes[i][j]:
            ax.text(j, i, notes[i][j], ha="center", va="center", fontsize=8,
                    color="#1a1a1a" if M[i][j] < 3 else "white")
ax.set_xticks(np.arange(-.5, len(families), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(stages), 1), minor=True)
ax.grid(which="minor", color="#DDDDDD", lw=1.2)
ax.tick_params(which="minor", bottom=False, left=False)
ax.set_title("GPR 指标体系 → 面波项目流水线：五层级 × 六族指标适配矩阵\n（颜色越深 = 文献中应用越成熟；格内为该处的典型用法）",
             fontsize=13, fontweight="bold", pad=14)
cbar = plt.colorbar(im, ax=ax, shrink=0.75, ticks=[0, 1, 2, 3])
cbar.ax.set_yticklabels(["未见使用", "零星", "常见", "成熟主流"], fontsize=9)
plt.tight_layout()
plt.savefig("fig2_metric_matrix.png", bbox_inches="tight", facecolor="white")
plt.close()

# ---------------------------------------------------------------- Figure 3
# Relative-error benchmark comparison across communities
items = [
    ("GPR 滑坡界面深度\n(钻孔验证, UAV-GPR)", 10, C["blue"], "<10% 相对误差, 绝对<1 m"),
    ("GPR 全波形反演振幅\n(2D近似引入误差)", 30, C["blue"], "警示值: 维数简化代价"),
    ("面波频散拾取 RMSE\n(MASW 堤防实测)", 16.2, C["teal"], "13.6–21% 站点均值"),
    ("面波 Vs30 变异系数\n(MASW 主动源)", 5.5, C["teal"], "侵入式仅 1–3%"),
    ("ERT 单元电阻率\n(DL 反演, 合成)", 10, C["orange"], "R²>0.75 同时达成"),
    ("医学超声声速图\n(DL, 物理体模)", 1.0, C["purple"], "RMSE 15.2 m/s ≈1%"),
    ("导波相速度重建\n(增强相位谱法)", 3, C["green"], "经典法 ±11% 作对照"),
    ("PINN 材料参数\n(合成, 分层模型)", 1.7, C["gray"], "合成数据下限参考"),
]
fig, ax = plt.subplots(figsize=(12.5, 7.0), dpi=150)
y = np.arange(len(items))[::-1]
vals = [it[1] for it in items]
cols = [it[2] for it in items]
bars = ax.barh(y, vals, color=cols, alpha=0.88, height=0.62)
for yi, it, v in zip(y, items, vals):
    ax.text(v + 0.4, yi + 0.12, f"{v}%", va="center", fontsize=10, fontweight="bold", color="#333")
    ax.text(v + 0.4, yi - 0.22, it[3], va="center", fontsize=8.3, color="#777")
ax.set_yticks(list(y))
ax.set_yticklabels([it[0] for it in items], fontsize=9.5)
ax.set_xlabel("典型误差水平（%，越低越好）", fontsize=11)
ax.set_xlim(0, 42)
ax.axvline(10, color=C["red"], ls="--", lw=1.4, alpha=0.7)
ax.text(10.4, -0.62, "10%：工程可接受的经验分界线", fontsize=9, color=C["red"])
ax.set_title("各领域可比的相对误差参考值（文献实测/验证值）", fontsize=13.5, fontweight="bold", pad=12)
handles = [mpatches.Patch(color=C["blue"], label="GPR"),
           mpatches.Patch(color=C["teal"], label="面波/地震"),
           mpatches.Patch(color=C["orange"], label="ERT"),
           mpatches.Patch(color=C["purple"], label="医学超声"),
           mpatches.Patch(color=C["green"], label="导波NDT"),
           mpatches.Patch(color=C["gray"], label="PINN(合成)")]
ax.legend(handles=handles, loc="lower right", fontsize=9, framealpha=0.9)
plt.tight_layout()
plt.savefig("fig3_error_benchmarks.png", bbox_inches="tight", facecolor="white")
plt.close()

print("done")
