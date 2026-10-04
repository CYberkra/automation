# -*- coding: utf-8 -*-
"""Flowchart: contractual Adaptive-ISP training pipeline (current status + target architecture)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "Noto Sans CJK SC", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

OK   = "#5C8A4D"   # 已建成/冻结
TODO = "#D9822B"   # 待建
GATE = "#B5443B"   # 闸门
NEU  = "#3B6EA5"   # 中性/文献继承
LIGHT= "#F2F4F7"

fig, ax = plt.subplots(figsize=(16, 11.5), dpi=150)
ax.set_xlim(0, 160); ax.set_ylim(0, 115); ax.axis("off")

def box(x, y, w, h, text, fc, ec=None, fs=9.5, tc="white", bold=True, lw=1.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5",
                                fc=fc, ec=ec or fc, lw=lw))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", linespacing=1.45)

def arrow(x1, y1, x2, y2, color="#555", lw=1.8, style="-|>", curve=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 color=color, lw=lw,
                                 connectionstyle=f"arc3,rad={curve}",
                                 shrinkA=2, shrinkB=2))

def band(y, h, label, color):
    ax.add_patch(FancyBboxPatch((1, y), 13, h, boxstyle="round,pad=0.3",
                                fc=color, ec="none", alpha=0.25))
    ax.text(7.5, y + h/2, label, ha="center", va="center", fontsize=10.5,
            color="#333", fontweight="bold", rotation=90)

ax.text(80, 112.5, "契约化自适应处理管线：训练方法与网络架构全景（2026-10-01）",
        ha="center", fontsize=16, fontweight="bold", color="#1a1a1a")

# ---------------- Band 1: 数据生成层 ----------------
band(88, 20, "① 数据生成", OK)
box(16, 95, 26, 9, "场景族参数采样设计\n【待冻结 B2】覆盖层厚度/Debye色散/\n界面倾角/目标档", TODO, fs=8.5)
box(46, 95, 26, 9, "gprMax V4 FDTD 批量仿真\n【已建成】GPU·t3基线\n99例/批·约44 min", OK, fs=8.5)
box(76, 95, 26, 9, "BG/NC/TGT 配对道集\n【已建成】NC负控零差分\n2D/3D一致性已校核", OK, fs=8.5)
box(106, 95, 24, 9, "官方 SFCW 链\n【已建成】20–170 MHz\nHann窗·带限B-scan", OK, fs=8.5)
box(134, 95, 24, 9, "训练输入：B-scan\n(dB刻度·固定尺寸)", NEU, fs=9)
arrow(42, 99.5, 46, 99.5); arrow(72, 99.5, 76, 99.5)
arrow(102, 99.5, 106, 99.5); arrow(130, 99.5, 134, 99.5)
arrow(146, 95, 146, 84)

# ---------------- Band 2: 策略网络层 ----------------
band(52, 34, "② 策略网络\n（待建）", TODO)
box(16, 74, 26, 10, "状态 s_t\n当前B-scan + EC通道\n(阶段计数t·算子历史)", NEU, fs=8.5)
box(46, 74, 30, 10, "特征骨干【继承AdaptiveISP】\n4×[Conv4×4/s2→BN→LReLU0.2]\n通道32→256 + FC-128\n参数量约1–2M", TODO, fs=8.5)
box(80, 80, 24, 6.5, "算子头 softmax\n30项目录 + STOP", TODO, fs=8.5)
box(80, 72, 24, 6.5, "参数档头 softmax\n(离散档·沿用目录标定)", TODO, fs=8.5)
box(80, 64, 24, 6.5, "价值头 v(s)\nPPO优势估计", TODO, fs=8.5)
box(110, 74, 22, 10, "算子作用于道集\n确定性数组运算\n(mean/SVD/RPCA/增益)", OK, fs=8.5)
box(136, 74, 22, 10, "逐阶段循环\nt ≤ T_max=4\n多数输入2–3阶段", NEU, fs=8.5)
arrow(42, 79, 46, 79)
arrow(76, 79, 80, 82.5)
arrow(76, 76, 80, 75.2, curve=0.15)
arrow(76, 73, 80, 67.2, curve=0.25)
arrow(104, 83, 108, 80.5); arrow(104, 75, 108, 77.5, curve=0.1)
arrow(132, 79, 136, 79)
# 循环回流
arrow(147, 74, 147, 62, color=NEU); arrow(147, 62, 29, 62, color=NEU)
arrow(29, 62, 29, 74, color=NEU)
ax.text(40, 63.5, "输出 STOP 或达到 T_max → 结算终端奖励", fontsize=9,
        color=NEU, ha="center", fontweight="bold")
# 扩展回路：目录升版 → 动作空间同步
box(78, 52.5, 56, 7.5, "扩展回路【预留】算子目录升版 vX.Y → 动作空间同步扩列\n新族五步入场（提案§3.5）· 旧权重续用不重构", TODO, fs=8)
ax.add_patch(FancyArrowPatch((92, 60.2), (92, 79.8), arrowstyle="-|>",
                             color=TODO, lw=1.6, linestyle=(0, (4, 2)),
                             shrinkA=1, shrinkB=2))
arrow(120, 74, 120, 47)

# ---------------- Band 3: 奖励层 ----------------
band(30, 20, "③ 契约奖励\n（已冻结）", OK)
box(16, 34, 30, 12, "冻结契约求值【已建成】\n容差契约v0.2\nτ_A=0.20/τ_D=0.95/ε=1.5\n权重契约v0.2(dB恒等式)", OK, fs=8.5)
box(52, 34, 26, 12, "R_contract\nD保真 − N_b惩罚 − 截幅项\n硬约束违规→钳制到\nidentity锚之下", OK, fs=8.5)
box(84, 34, 24, 12, "结构锚在线校准\nidentity恒读−due\noracle恰读+due", OK, fs=8.5)
box(114, 34, 22, 12, "成本惩罚\n−λc·cost(算子时延)\n推理期可调λc", NEU, fs=8.5)
box(140, 34, 18, 12, "终端奖励 R\n（稀疏·回合末）", GATE, fs=9)
arrow(120, 47, 65, 34.5, curve=0.0)
arrow(46, 40, 52, 40); arrow(78, 40, 84, 40)
arrow(108, 40, 114, 40); arrow(136, 40, 140, 40)

# ---------------- Band 4: 训练层 ----------------
band(12, 18, "④ 训练\n（闸门后）", GATE)
box(16, 15, 28, 11, "PPO 更新\n策略/价值网络\n奖励无需可微", TODO, fs=9)
box(50, 15, 28, 11, "数据纪律\n仅用开发族 {C1,C3}\n测试族 {C5,C8} 零接触", OK, fs=8.5)
box(84, 15, 30, 11, "监控与护栏\n奖励曲线·选择频率熵(防坍缩)\n多次训练报std·固定种子\nr1/r2双遍", TODO, fs=8.5)
box(120, 15, 38, 11, "【启动闸门】G4解除(用户签认)\n+ S1/S3容差重校准完成\n+ 采样设计冻结 + 奖励链验收 + PyTorch环境锁定", GATE, fs=8.5)
arrow(149, 34, 40, 26.5, color=GATE, lw=2, curve=0.15)

# ---------------- Band 5: 验证层 ----------------
band(0, 10, "⑤ 验证", OK)
box(16, 2.5, 34, 6, "测试族一次性评价（S5式纪律）\n基线：identity / 协议首跑全局赢家B9·B7", OK, fs=8)
box(56, 2.5, 34, 6, "物理验收（前置中）\n20–170MHz全带真实几何3D\n网格/PML/记录长度", TODO, fs=8)
box(96, 2.5, 34, 6, "实测桥接与钻孔验证\n域差桥接(pretrain-to-alignment)\n对齐 <1m/<10% 表述", TODO, fs=8)
arrow(33, 15, 33, 8.5); arrow(73, 15, 73, 8.5); arrow(90, 8.5, 96, 5.5)

# legend
handles = [mpatches.Patch(color=OK, label="已建成/已冻结"),
           mpatches.Patch(color=TODO, label="待建/待冻结"),
           mpatches.Patch(color=GATE, label="闸门/硬约束"),
           mpatches.Patch(color=NEU, label="中性设计（继承文献）")]
ax.legend(handles=handles, loc="upper right", fontsize=10, framealpha=0.95,
          bbox_to_anchor=(0.995, 0.985))

plt.tight_layout()
plt.savefig("fig_training_pipeline.png", bbox_inches="tight", facecolor="white")
plt.close()
print("done")
