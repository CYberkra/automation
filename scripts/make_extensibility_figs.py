# -*- coding: utf-8 -*-
"""Figures for pipeline-extensibility research report (2026-10-01)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "Noto Sans CJK SC", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

OK   = "#5C8A4D"   # 已吸收
TODO = "#D9822B"   # 建议吸收
NO   = "#B5443B"   # 明确不采用
NEU  = "#3B6EA5"

# ============ Fig 1: 机制借鉴矩阵 ============
lineages = ["ISP 族\n(ReconfigISP/\nAdaptiveISP)",
            "AutoML 管线合成\n(AlphaD3M/\nML-ReinBo)",
            "可微 DSP/音频\n(DeepAFx/灰盒建模)",
            "认知雷达\n(RL 波形选择)",
            "语音增强\n(门控专家/Quality-Net)",
            "MoE/自适应计算\n(动态路由/PonderNet)"]
mechs = ["模块池+注册表", "语法/掩码约束\n动作空间", "效率惩罚\nλc·cost", "防复选/防坍缩\n机制", "STOP/动态停止",
         "动作空间扩展\n(权重续用)", "金牌参照\nGoodhart防护", "可微代理\n(黑盒算子入场)"]
# 0=不适用/空白 1=建议吸收(TODO) 2=已吸收(OK) 3=不采用(NO)
M = np.array([
    [2, 2, 2, 2, 2, 0, 2, 3],   # ISP 族
    [2, 2, 0, 0, 0, 1, 0, 0],   # AutoML
    [0, 0, 0, 0, 0, 0, 0, 3],   # 可微DSP (代理不采用)
    [2, 1, 1, 0, 0, 0, 0, 0],   # 认知雷达
    [1, 0, 0, 1, 0, 0, 1, 0],   # 语音
    [0, 1, 0, 2, 2, 1, 2, 0],   # MoE/自适应计算
])
cmap = {0: "#EDEFF2", 1: TODO, 2: OK, 3: NO}

fig, ax = plt.subplots(figsize=(13.5, 7.2), dpi=150)
for i in range(len(lineages)):
    for j in range(len(mechs)):
        v = M[i, j]
        ax.add_patch(plt.Rectangle((j, len(lineages)-1-i), 0.94, 0.9, fc=cmap[v], ec="white", lw=2))
        if v == 2: ax.text(j+0.47, len(lineages)-1-i+0.45, "已吸收", ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        elif v == 1: ax.text(j+0.47, len(lineages)-1-i+0.45, "建议吸收", ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        elif v == 3: ax.text(j+0.47, len(lineages)-1-i+0.45, "不采用", ha="center", va="center", fontsize=9, color="white", fontweight="bold")
ax.set_xlim(0, len(mechs)); ax.set_ylim(0, len(lineages))
ax.set_xticks([j+0.47 for j in range(len(mechs))])
ax.set_xticklabels(mechs, fontsize=9.5)
ax.set_yticks([len(lineages)-1-i+0.45 for i in range(len(lineages))])
ax.set_yticklabels(lineages, fontsize=9.5)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_title("六大谱系 × 可借鉴机制矩阵（对本项目契约化自适应管线的映射，2026-10-01）",
             fontsize=13.5, fontweight="bold", pad=14)
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color=OK, label="已吸收进提案/既有设计"),
                   mpatches.Patch(color=TODO, label="建议吸收（登记为后续参照）"),
                   mpatches.Patch(color=NO, label="明确不采用（附理由）"),
                   mpatches.Patch(color="#EDEFF2", label="该谱系无此机制")],
          loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=4, fontsize=10, frameon=False)
plt.tight_layout()
plt.savefig("fig_borrow_matrix.png", bbox_inches="tight", facecolor="white")
plt.close()

# ============ Fig 2: 决策粒度 × 摊销成本 ============
fig, ax = plt.subplots(figsize=(11.5, 7), dpi=150)
pts = [
    # (名称, 粒度x 1=每任务全局 2=每数据集 3=逐输入 4=逐token/逐阶段, 单次决策摊销成本log10(ms), 颜色, 注释)
    ("ReconfigISP\n(NAS/DARTS)", 1.0, 6.0, NEU, "每任务离线搜索一次\n代理+剪枝"),
    ("AlphaD3M\n(MCTS+自对弈)", 2.0, 5.0, NEU, "每数据集合成一条管线\n快传统AutoML约10×"),
    ("贝叶斯优化\n(Auto-sklearn/BOHB)", 1.7, 4.3, NEU, "每任务串行试配\n适合昂贵小预算"),
    ("认知雷达RL\n(波形目录选型)", 2.85, 2.1, OK, "冻结策略LUT在线选型·毫秒级"),
    ("AdaptiveISP\n(逐阶段贪心RL)", 3.8, 0.08, OK, "1.2 ms/阶段·逐输入生成管线"),
    ("MoE 动态路由", 4.25, -0.6, NEU, "逐token门控·亚毫秒"),
    ("本项目目标\n(契约化逐输入选型)", 3.25, 0.75, TODO, "≤数ms/阶段·目录离散档"),
]
for name, x, y, c, note in pts:
    ax.scatter(x, y, s=560, color=c, alpha=0.85, edgecolor="white", lw=2, zorder=3)
    if name.startswith("认知雷达"):
        off1, off2, ha2 = (-16, 26), (-130, -6), "center"
        ax.annotate(name, (x, y), textcoords="offset points", xytext=off1,
                    ha="center", fontsize=10, fontweight="bold", color="#1a1a1a")
        ax.annotate(note, (x, y), textcoords="offset points", xytext=off2,
                    ha=ha2, fontsize=8.2, color="#444")
    elif name.startswith("MoE"):
        ax.annotate(name, (x, y), textcoords="offset points", xytext=(-40, 26),
                    ha="center", fontsize=10, fontweight="bold", color="#1a1a1a")
        ax.annotate(note, (x, y), textcoords="offset points", xytext=(-55, -40),
                    ha="center", fontsize=8.2, color="#444")
    else:
        ax.annotate(name, (x, y), textcoords="offset points", xytext=(0, 26),
                    ha="center", fontsize=10, fontweight="bold", color="#1a1a1a")
        ax.annotate(note, (x, y), textcoords="offset points", xytext=(0, -40),
                    ha="center", fontsize=8.2, color="#444")
ax.axvspan(2.6, 4.4, color=OK, alpha=0.07)
ax.text(3.5, 6.6, "逐输入/逐阶段自适应区\n（本项目目标区间）", ha="center", fontsize=10.5,
        color=OK, fontweight="bold")
ax.set_xlim(0.4, 4.6); ax.set_ylim(-1.4, 7.2)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(["每任务全局\n(离线一次)", "每数据集\n(元学习)", "逐输入\n(在线)", "逐token/阶段\n(在线)"], fontsize=10)
ax.set_ylabel("单次决策摊销成本  log10(ms)", fontsize=11)
ax.set_xlabel("决策粒度（越靠右越细）", fontsize=11)
ax.set_title("管线配置方法的决策粒度 × 摊销成本全景（定性定位，依据各文献自报数据）",
             fontsize=13, fontweight="bold", pad=12)
ax.grid(alpha=0.25, zorder=0)
plt.tight_layout()
plt.savefig("fig_granularity_map.png", bbox_inches="tight", facecolor="white")
plt.close()
print("done")
