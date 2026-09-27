#!/usr/bin/env python3
"""tilt_target_v0 倾斜目标阶梯切片生成器（设计草案，未冻结）。

把水平平板目标（4 m 沿测线 × 0.5 m 厚）绕板心在 y-z 平面倾斜，
用轴对齐 #box 条带做阶梯逼近。全部计算在整数格点上进行：
条带 y 边界与 z 偏移都是格距整数倍，从构造上保证网格对齐。

对称整数偏移规则（n 偶数条带、步长 s 偶数格）：
    dz_offset_i = -(2i + 1 - n) * (s // 2)   [格]
即条带偏移关于板心反对称（+y 方向下倾），成对 ±(1,3,5,...) 格。

角度档位由格点量化决定，如实报告有效角度：
    theta_eff = atan(s * dz / w)
    S1: s=2, w=0.5(2D)/1.0(3D)  -> atan(0.1)   = 5.7106 deg
    S2: s=2, w=0.25(2D)/0.5(3D) -> atan(0.2)   = 11.3099 deg
（名义 5°/10° 不可格点整除，按纪律调整到可整除档并记录。）

几何来源（逐值取自真实配置，SHA 随产物登记）：
  2D 锚点 configs/research/batch2d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02.in
    目标 y 14-18, z 19.75-20.25, dy=dz=0.025, 材料名 target
  3D 锚点 configs/research/a0_3d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM.in
    目标 x 3-5, y 3-7, z 19.75-20.25, dy=dz=0.05, 材料名 wet

纪律：本脚本只生成几何提案，不冻结契约、不生成 gate、不产生物理/训练标签。
"""
import hashlib
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

SRC_2D = "configs/research/batch2d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02.in"
SRC_3D = "configs/research/a0_3d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM.in"
OUT = REPO / "configs/research/tilt_target_v0/slices.json"

# (tier_key, spacing, y_center, z_center, x_range_or_None, material_name)
GEOMS = {
    "2d": dict(d=0.025, yc=16.0, zc=20.0, xr=None, mat="target",
               ylim=(0.0, 32.0), pml_cells_y=40, pml_cells_z=40),
    "3d": dict(d=0.05, yc=5.0, zc=20.0, xr=(3.0, 5.0), mat="wet",
               ylim=(0.0, 10.0), pml_cells_y=20, pml_cells_z=20),
}
L_M = 4.0       # 板水平投影长
T_M = 0.5       # 板厚（竖直方向，垂直厚度近似见 JSON notes）
COVER_BOTTOM_Z = 27.0

# (variant_suffix, strip_width_m, step_cells)
VARIANTS = {
    "S1": dict(w=0.5, s=2),   # theta_eff = atan(2d / 0.5)
    "S2": dict(w=0.25, s=2),  # theta_eff = atan(2d / 0.25)（3D 档 w 见下）
}
# 3D 的条带宽与 2D 不同（格距粗，为保持同一有效角度）：S1 w=1.0, S2 w=0.5
VARIANT_W_3D = {"S1": 1.0, "S2": 0.5}


def sha256_of(rel):
    return hashlib.sha256((REPO / rel).read_bytes()).hexdigest()


def build(tier, suffix):
    g = GEOMS[tier]
    d = g["d"]
    w = VARIANT_W_3D[suffix] if tier == "3d" else VARIANTS[suffix]["w"]
    s = VARIANTS[suffix]["s"]
    w_c = round(w / d)
    assert abs(w_c * d - w) < 1e-12, "strip width not grid aligned"
    n = round(L_M / w)
    assert n * w_c == round(L_M / d), "strip count does not tile plate length"
    assert n % 2 == 0 and s % 2 == 0, "symmetric integer offsets need even n and s"
    y0_c = round((g["yc"] - L_M / 2) / d)
    zc_c = round(g["zc"] / d)
    t2_c = round((T_M / 2) / d)
    assert abs(t2_c * d - T_M / 2) < 1e-12

    slices = []
    for i in range(n):
        off_c = -(2 * i + 1 - n) * (s // 2)
        ya_c = y0_c + i * w_c
        yb_c = ya_c + w_c
        za_c = zc_c - t2_c + off_c
        zb_c = zc_c + t2_c + off_c
        slices.append(dict(i=i, off_cells=off_c,
                           y0=round(ya_c * d, 6), y1=round(yb_c * d, 6),
                           z0=round(za_c * d, 6), z1=round(zb_c * d, 6)))
    # 反对称检查：strip i 与 strip n-1-i 偏移相反
    for i in range(n):
        assert slices[i]["off_cells"] == -slices[n - 1 - i]["off_cells"]
    # 连续性：y 向无缝无重叠
    for i in range(1, n):
        assert slices[i]["y0"] == slices[i - 1]["y1"]
    # 包络与约束
    z_top_max = max(s["z1"] for s in slices)
    z_bot_min = min(s["z0"] for s in slices)
    assert z_top_max < COVER_BOTTOM_Z, "slice intrudes into cover layer"
    pml_y_m = g["pml_cells_y"] * d
    pml_z_m = g["pml_cells_z"] * d
    assert slices[0]["y0"] >= pml_y_m and slices[-1]["y1"] <= g["ylim"][1] - pml_y_m
    assert z_bot_min >= pml_z_m
    if g["xr"] is not None:
        assert g["xr"][0] >= pml_y_m and g["xr"][1] <= 8.0 - pml_y_m
    theta_eff = math.degrees(math.atan2(s * d, w))
    max_dev = (n - 1) * (s // 2) * d

    if tier == "2d":
        box_lines = ["#box: inf %.6g %.6g inf %.6g %.6g %s" % (
            s_["y0"], s_["z0"], s_["y1"], s_["z1"], g["mat"]) for s_ in slices]
    else:
        box_lines = ["#box: %.6g %.6g %.6g %.6g %.6g %.6g %s" % (
            g["xr"][0], s_["y0"], s_["z0"], g["xr"][1], s_["y1"], s_["z1"],
            g["mat"]) for s_ in slices]

    return dict(
        tier=tier, variant=suffix, grid_m=d,
        strip_width_m=w, step_cells=s, n_strips=n,
        theta_eff_deg=round(theta_eff, 6),
        theta_nominal_deg=5.0 if suffix == "S1" else 10.0,
        center_y=g["yc"], center_z=g["zc"],
        max_z_deviation_m=round(max_dev, 6),
        z_range=[z_bot_min, z_top_max],
        x_range=g["xr"], material=g["mat"],
        direction="z decreases along +y (沿测线下倾)；符号对称档未做，留待后续决定",
        slices=slices, box_lines=box_lines,
    )


def main():
    out = {
        "contract_id": "tilt_target_v0",
        "status": "design_draft_unfrozen",
        "created_by": "scripts/gen_tilt_target_slices.py",
        "source_inputs": {
            "2d_anchor": {"file": SRC_2D, "sha256": sha256_of(SRC_2D)},
            "3d_anchor": {"file": SRC_3D, "sha256": sha256_of(SRC_3D)},
        },
        "base_plate": {"length_y_m": L_M, "thickness_z_m": T_M,
                       "note": "倾斜后 y 向水平投影保持 4 m；条带竖直厚度保持 0.5 m，"
                               "垂直厚度=0.5*cos(theta)，5.71°/11.31° 下差异 0.5%/1.9%，"
                               "为已知近似并在此声明"},
        "variants": {},
        "hard_limits": [
            "设计草案未冻结，不构成执行契约；执行前须另行冻结 gate 与用户签认",
            "阶梯逼近是格点离散的标准做法，角度档位为格点量化后的有效角度，名义 5/10 度不可整除",
            "倾斜档与水平档共用同一 BG/NC 母模型，背景配对关系不变",
            "不生成物理/训练标签；倾角档位为工程假设非场地实测",
            "测试族 {C5,C8} 已在 S5 消费，本扩展只面向开发族 {C1,C3}，不据此调参",
        ],
    }
    for tier in ("2d", "3d"):
        for suffix in ("S1", "S2"):
            v = build(tier, suffix)
            out["variants"]["%s_%s" % (tier, suffix)] = v
    OUT.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(out, ensure_ascii=False, indent=1)
    OUT.write_text(text, encoding="utf-8", newline="\n")
    for key, v in out["variants"].items():
        print("%s: n=%d w=%.3g m step=%d cells theta_eff=%.4f deg z_range=%s" % (
            key, v["n_strips"], v["strip_width_m"], v["step_cells"],
            v["theta_eff_deg"], v["z_range"]))
    print("written:", OUT.relative_to(REPO))


if __name__ == "__main__":
    sys.exit(main())
