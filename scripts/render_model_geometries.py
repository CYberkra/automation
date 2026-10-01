# -*- coding: utf-8 -*-
"""
把每个标准模型的 gprMax .in 输入逐箱解析并栅格化成真实材料分布图。
输出: fig_model_<short>.png (E:\\automation_djh\\) + model_gallery.html
"""
import os, re, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = r"E:\automation_djh"
CFG_PILOT = os.path.join(ROOT, r"automation_repo\configs\research\batch2d_b2_pilot_co")
CFG_T3 = os.path.join(ROOT, r"automation_repo\configs\research\batch2d_slope_t3_co")

MODELS = [
    # (short, 目录, 文件名, 中文说明) —— 全部为无目标的纯背景（BG）母模型
    ("C1mX_BG", CFG_PILOT, "B2D-C1mX-BG-CO33-t17.in",
     "平地 · 覆土 1 m · 无目标背景"),
    ("C3mX_BG", CFG_T3, "B2D-C3mX-BG-CO33-t17.in",
     "平地 · 覆土 3 m · 无目标背景"),
    ("C3mS2X_BG", CFG_T3, "B2D-C3mS2X-BG-CO33-t17.in",
     "坡地（阶梯 0.25 m/格）· 覆土 3 m · 无目标背景"),
    ("C3mS2TZX_BG", CFG_T3, "B2D-C3mS2TZX-BG-CO33-t17.in",
     "坡地 + 过渡带（TZ，εr=12.5 半强度色散）· 覆土 3 m · 无目标背景"),
    ("C1p5mS1X_BG", CFG_PILOT, "B2D-C1p5mS1X-BG-CO33-t17.in",
     "B 层缓坡（S1X）· 覆土 1.5 m（天线处）· 无目标背景"),
    ("C1p5mS3X_BG", CFG_PILOT, "B2D-C1p5mS3X-BG-CO33-t17.in",
     "B 层陡坡（S3X）· 覆土 1.5 m（天线处）· 无目标背景"),
    ("C2mS3X_BG", CFG_PILOT, "B2D-C2mS3X-BG-CO33-t17.in",
     "B 层陡坡（S3X）· 覆土 2 m（天线处）· 无目标背景"),
    ("C2mS3TZX_BG", CFG_PILOT, "B2D-C2mS3TZX-BG-CO33-t17.in",
     "B 层陡坡 + 过渡带（TZ）· 覆土 2 m · 无目标背景"),
]

# 材料类别配色（kimi 风格、全套图一致）
CAT_STYLE = {
    "air":    ("#F7F9FB", "空气"),
    "rock":   ("#9AA5B1", "基岩"),
    "cover":  ("#6BAED6", "覆土"),
    "target": ("#E4572E", "目标体"),
    "tz":     ("#74C69D", "过渡带"),
    "ncpos":  ("#F2B134", "目标对照位（与基岩同质）"),
}
OTHER_COLOR = "#F2B134"

def categorize(name):
    n = name.lower()
    if "nullcontrast" in n:
        return "ncpos"
    for k in ("cover", "rock", "target", "tz"):
        if k in n:
            return k
    return "other"

def parse_in(path):
    d = {"Y": None, "Z": None, "materials": {}, "dispersive": set(),
         "boxes": [], "tx": None, "rx": None, "dx": 0.025, "pml": None}
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#domain:"):
                p = line.split(":")[1].split()
                d["Y"], d["Z"] = float(p[1]), float(p[2])
            elif line.startswith("#dx_dy_dz:"):
                d["dx"] = float(line.split(":")[1].split()[0])
            elif line.startswith("#pml_cells:"):
                d["pml"] = [int(v) for v in line.split(":")[1].split()]
            elif line.startswith("#material:"):
                p = line.split(":")[1].split()
                d["materials"][p[4]] = {"er": float(p[0]), "sigma": float(p[1])}
            elif line.startswith("#add_dispersion_debye:"):
                d["dispersive"].add(line.split()[-1])
            elif line.startswith("#box:"):
                p = line.split(":")[1].split()
                d["boxes"].append((float(p[1]), float(p[2]), float(p[4]), float(p[5]), p[6]))
            elif line.startswith("#hertzian_dipole:"):
                p = line.split(":")[1].split()
                d["tx"] = (float(p[2]), float(p[3]))
            elif line.startswith("#rx:"):
                p = line.split(":")[1].split()
                d["rx"] = (float(p[1]), float(p[2]))
    return d

def rasterize(d, res=0.02):
    Y, Z = d["Y"], d["Z"]
    ny, nz = int(round(Y / res)), int(round(Z / res))
    grid = np.zeros((nz, ny), dtype=np.int16)  # 0 = air
    cats = ["air"]
    for y1, z1, y2, z2, mat in d["boxes"]:
        cat = categorize(mat)
        if cat not in cats:
            cats.append(cat)
        iy1 = max(0, int(round(y1 / res))); iy2 = min(ny, int(round(y2 / res)))
        iz1 = max(0, int(round(z1 / res))); iz2 = min(nz, int(round(z2 / res)))
        grid[iz1:iz2, iy1:iy2] = cats.index(cat)
    return grid, cats

def render(short, path, desc):
    d = parse_in(path)
    grid, cats = rasterize(d)
    Y, Z = d["Y"], d["Z"]
    res = Y / grid.shape[1]

    # 颜色表
    rgba = np.ones((grid.shape[0], grid.shape[1], 3))
    color_of = {"air": "#F7F9FB", "other": OTHER_COLOR}
    for i, cat in enumerate(cats):
        hexc = CAT_STYLE.get(cat, (OTHER_COLOR, ""))[0]
        rgb = tuple(int(hexc[j:j+2], 16) / 255 for j in (1, 3, 5))
        rgba[grid == i] = rgb

    fig, ax = plt.subplots(figsize=(11, 6.2), dpi=150)
    ax.imshow(rgba, extent=[0, Y, 0, Z], origin="lower",
              interpolation="nearest", aspect="auto")
    ax.invert_yaxis()  # z=0 在下(深)，z=Z 在上(空气/天线)

    # PML 边界带（1 m = 40 cells × 0.025）
    if d["pml"]:
        pml_m = d["pml"][1] * d["dx"]
        for spine_pos in ("left", "right", "top", "bottom"):
            pass
        from matplotlib.patches import Rectangle
        for (x0, y0, w, h) in [(0, 0, pml_m, Z), (Y - pml_m, 0, pml_m, Z),
                               (0, 0, Y, pml_m), (0, Z - pml_m, Y, pml_m)]:
            ax.add_patch(Rectangle((x0, y0), w, h, fill=False,
                                   edgecolor="#888888", linestyle=":", linewidth=1.0, alpha=0.9))

    # 天线
    if d["tx"]:
        ax.plot(*d["tx"], marker="v", color="#C0392B", markersize=11,
                markeredgecolor="white", markeredgewidth=1.2, zorder=5)
        ax.annotate("Tx", d["tx"], textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=10, color="#C0392B", fontweight="bold")
    if d["rx"]:
        ax.plot(*d["rx"], marker="^", color="#1F6FEB", markersize=11,
                markeredgecolor="white", markeredgewidth=1.2, zorder=5)
        ax.annotate("Rx", d["rx"], textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=10, color="#1F6FEB", fontweight="bold")

    # 图例：材料（带 εr/σ/色散）+ 天线 + PML
    handles = []
    for cat in cats:
        if cat == "air":
            continue
        hexc, zh = CAT_STYLE.get(cat, (OTHER_COLOR, cat))
        # 找到对应材料的参数
        prm = ""
        for mname, mp in d["materials"].items():
            if categorize(mname) == cat:
                prm = f"  εr={mp['er']:g}, σ={mp['sigma']:g} S/m"
                if mname in d["dispersive"]:
                    prm += "（Debye 色散）"
                break
        handles.append(Patch(facecolor=hexc, edgecolor="none", label=f"{zh}{prm}"))
    handles.append(Line2D([0], [0], marker="v", color="w", markerfacecolor="#C0392B",
                          markersize=9, label=f"Tx 发射 ({d['tx'][0]:.2f}, {d['tx'][1]:.2f}) m"))
    handles.append(Line2D([0], [0], marker="^", color="w", markerfacecolor="#1F6FEB",
                          markersize=9, label=f"Rx 接收 ({d['rx'][0]:.2f}, {d['rx'][1]:.2f}) m"))
    handles.append(Line2D([0], [0], color="#888888", linestyle=":", label="PML 吸收边界（1 m）"))
    ax.legend(handles=handles, loc="lower right", fontsize=9, framealpha=0.92)

    ax.set_xlabel("y（测线方向，m）")
    ax.set_ylabel("z（m）")
    ax.set_title(f"标准模型几何 · {short}    域 {Y:g} m × {Z:g} m，网格 {d['dx']*1000:.0f} mm",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = os.path.join(ROOT, f"fig_model_{short}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out, desc, d

def main():
    entries = []
    for short, cfg, fname, desc in MODELS:
        path = os.path.join(cfg, fname)
        if not os.path.exists(path):
            print("MISSING:", path)
            continue
        out, desc, d = render(short, path, desc)
        entries.append({"short": short, "png": os.path.basename(out), "desc": desc,
                        "Y": d["Y"], "Z": d["Z"],
                        "mats": list(d["materials"].keys()),
                        "nboxes": len(d["boxes"])})
        print("OK", short, "boxes=", len(d["boxes"]), "mats=", list(d["materials"].keys()))

    # HTML 画廊
    cards = []
    for e in entries:
        cards.append(f"""
      <div class="card">
        <h3>{e['short']}</h3>
        <p class="meta">域 {e['Y']:g} m × {e['Z']:g} m · 几何箱 {e['nboxes']} 个 · 材料：{', '.join(e['mats'])}</p>
        <p>{e['desc']}</p>
        <img src="{e['png']}" alt="{e['short']}" />
      </div>""")
    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8" />
<title>标准仿真模型几何画廊</title>
<style>
  body {{ font-family: "Microsoft YaHei", sans-serif; margin: 0; background: #f5f7fa; color: #1f2933; }}
  header {{ padding: 24px 32px 8px; }}
  header h1 {{ margin: 0 0 6px; font-size: 22px; }}
  header p {{ margin: 0; color: #5b6470; font-size: 13px; }}
  main {{ display: grid; grid-template-columns: 1fr; gap: 20px; padding: 20px 32px 40px; }}
  .card {{ background: #fff; border: 1px solid #e3e8ee; border-radius: 10px; padding: 16px 20px; box-shadow: 0 1px 3px rgba(16,24,40,.06); }}
  .card h3 {{ margin: 0 0 4px; font-size: 16px; color: #0f62fe; }}
  .card .meta {{ margin: 0 0 6px; font-size: 12px; color: #8a94a6; }}
  .card p {{ margin: 0 0 10px; font-size: 13px; }}
  .card img {{ width: 100%; border-radius: 6px; border: 1px solid #eef1f4; }}
</style>
</head>
<body>
<header>
  <h1>标准仿真模型几何画廊（试点批 8 族 · 按 .in 逐箱真实还原）</h1>
  <p>生成于 2026-10-01 · 数据源：automation_repo/configs/research/ 各 t17 道输入文件 · 阶梯坡地、目标体、过渡带均为输入文件原样栅格化</p>
</header>
<main>{''.join(cards)}
</main>
</body>
</html>"""
    with open(os.path.join(ROOT, "model_gallery.html"), "w", encoding="utf-8") as f:
        f.write(html)
    print("HTML written")

if __name__ == "__main__":
    main()
