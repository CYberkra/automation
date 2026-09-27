#!/usr/bin/env python3
"""slope_family_v0 坡地族（倾斜基覆界面）几何生成器（设计草案，未冻结）。

方案 A（本脚本）：地表保持水平 z=30、采集几何不变（天线 z=45），
基覆界面（覆盖层底）沿测线 y 阶梯化倾斜；方案 B（地表同步倾斜+天线跟随）
因天线净空/域高约束留待设计书分析，本脚本不生成。

做法：y 向等宽条带，每条带内界面水平，条带间界面按格点取整递变：
    z_if_i = round((zc + (y_strip_center - y_mid) * tan_theta) / d)   [格]
实现斜率 theta_eff 由实现界面的最小二乘拟合如实报告（量化抖动 <= 0.5 格/条带）。

中心深度约束（最低覆盖层厚度 >= 1 m，即界面最高点 <= 29 m）：
    zc = min(27, 29 - (YSPAN/2) * tan_theta)
T2/T3 档为满足该约束中心界面下移（覆盖层中心变厚），逐档如实记录。

过渡带（TZ 变体）：界面下方 1 m 厚强风化过渡带，材料 tzone eps_r=12 sigma=0.005
——插值假设（覆盖层 16/0.01 与砂岩 9/0.001 之间），非实测。

几何来源（SHA 随产物登记）：
  2D 锚点 configs/research/batch2d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02.in
  3D 锚点 configs/research/a0_3d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM.in

纪律：只生成几何提案；不冻结契约、不生成 gate、不产生物理/训练标签。
"""
import hashlib
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC_2D = "configs/research/batch2d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02.in"
SRC_3D = "configs/research/a0_3d_v1/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-3D5CM.in"
OUT = REPO / "configs/research/slope_family_v0/geometry.json"

SURF_Z = 30.0
MIN_COVER_M = 1.0          # 最薄覆盖层约束（文献强风化带下限量级，工程取整）
TZ_THICK_M = 1.0           # 过渡带厚度
TZ_MATERIAL = ("tzone", 12.0, 0.005)   # 插值假设，非实测

TIERS = {"T1": 0.1, "T2": 0.2, "T3": 0.4}   # tan(theta)

GEOMS = {
    "2d": dict(d=0.025, y0=0.0, y1=32.0, strip_w=0.25, pml_m=1.0, xr=None,
               domain_z=50.0),
    "3d": dict(d=0.05, y0=0.0, y1=10.0, strip_w=0.5, pml_m=1.0, xr=(0.0, 8.0),
               domain_z=50.0),
}


def sha256_of(rel):
    return hashlib.sha256((REPO / rel).read_bytes()).hexdigest()


def fit_slope(ys, zs):
    n = len(ys)
    my = sum(ys) / n
    mz = sum(zs) / n
    sxy = sum((y - my) * (z - mz) for y, z in zip(ys, zs))
    sxx = sum((y - my) ** 2 for y in ys)
    return sxy / sxx


def build(tier_key, tan_theta, with_tz, g):
    d = g["d"]
    w = g["strip_w"]
    w_c = round(w / d)
    assert abs(w_c * d - w) < 1e-12
    yspan = g["y1"] - g["y0"]
    n = round(yspan / w)
    assert n * w_c == round(yspan / d)
    ymid = (g["y0"] + g["y1"]) / 2
    half = yspan / 2 * tan_theta
    zc = min(27.0, SURF_Z - MIN_COVER_M - half)
    zc_c = round(zc / d)
    assert abs(zc_c * d - zc) < 1e-12, "interface center not grid aligned"

    strips = []
    for i in range(n):
        ya_c = round(g["y0"] / d) + i * w_c
        yb_c = ya_c + w_c
        ycen = (ya_c + yb_c) / 2 * d
        zif_c = zc_c + round((ycen - ymid) * tan_theta / d)
        strips.append(dict(i=i, y0=round(ya_c * d, 6), y1=round(yb_c * d, 6),
                           z_if=round(zif_c * d, 6)))
    # 连续性 + 单调性
    for i in range(1, n):
        assert strips[i]["y0"] == strips[i - 1]["y1"]
        assert strips[i]["z_if"] >= strips[i - 1]["z_if"]
    z_min = strips[0]["z_if"]
    z_max = strips[-1]["z_if"]
    cover_min = SURF_Z - z_max
    cover_max = SURF_Z - z_min
    assert cover_min >= MIN_COVER_M - 1e-9, "min cover violated"
    assert z_min >= g["pml_m"] + (TZ_THICK_M if with_tz else 0.0), "PML margin violated"
    # 目标区（锚点 y 14-18(2D) / 3-7(3D)，z 19.75-20.25）须在界面之下（砂岩侧）
    if g["xr"] is None:
        tgt_y0, tgt_y1 = 14.0, 18.0
    else:
        tgt_y0, tgt_y1 = 3.0, 7.0
    z_if_at_target = [s["z_if"] for s in strips
                      if s["y1"] > tgt_y0 and s["y0"] < tgt_y1]
    assert min(z_if_at_target) > 20.25, "target would intrude into cover"

    slope = fit_slope([(s["y0"] + s["y1"]) / 2 for s in strips],
                      [s["z_if"] for s in strips])
    theta_eff = math.degrees(math.atan(slope))
    max_dev = max(abs(s["z_if"] - (zc + ((s["y0"] + s["y1"]) / 2 - ymid) * tan_theta))
                  for s in strips)

    if g["xr"] is None:
        cover_lines = ["#box: inf %.6g %.6g inf %.6g %.6g cover" % (
            s["y0"], s["z_if"], s["y1"], SURF_Z) for s in strips]
        tz_lines = ["#box: inf %.6g %.6g inf %.6g %.6g %s" % (
            s["y0"], round(s["z_if"] - TZ_THICK_M, 6), s["y1"], s["z_if"],
            TZ_MATERIAL[0]) for s in strips] if with_tz else []
    else:
        cover_lines = ["#box: %.6g %.6g %.6g %.6g %.6g %.6g cover" % (
            g["xr"][0], s["y0"], s["z_if"], g["xr"][1], s["y1"], SURF_Z)
            for s in strips]
        tz_lines = ["#box: %.6g %.6g %.6g %.6g %.6g %.6g %s" % (
            g["xr"][0], s["y0"], round(s["z_if"] - TZ_THICK_M, 6),
            g["xr"][1], s["y1"], s["z_if"], TZ_MATERIAL[0])
            for s in strips] if with_tz else []

    return dict(
        tier=tier_key, tan_theta_nominal=tan_theta, transition_zone=with_tz,
        grid_m=d, strip_width_m=w, n_strips=n,
        zc_m=round(zc_c * d, 6),
        z_interface_range=[z_min, z_max],
        cover_thickness_range=[round(cover_min, 6), round(cover_max, 6)],
        theta_eff_deg=round(theta_eff, 6),
        max_quantization_deviation_m=round(max_dev, 6),
        direction="interface rises along +y（基覆界面沿 +y 抬升）；符号对称档未做",
        note_zc=("zc=27 锚点原值" if zc == 27.0 else
                 "zc 自 27 下移至 %.6g 以满足最薄覆盖层 >= %.6g m" % (zc, MIN_COVER_M)),
        strips=strips, cover_box_lines=cover_lines, tz_box_lines=tz_lines,
    )


def main():
    out = {
        "contract_id": "slope_family_v0",
        "status": "design_draft_unfrozen",
        "created_by": "scripts/gen_slope_family_geometry.py",
        "scheme": "A（地表水平 z=30、采集几何不变、仅基覆界面阶梯倾斜）；方案B（地表同步倾斜+天线跟随）未生成",
        "source_inputs": {
            "2d_anchor": {"file": SRC_2D, "sha256": sha256_of(SRC_2D)},
            "3d_anchor": {"file": SRC_3D, "sha256": sha256_of(SRC_3D)},
        },
        "tier_definition": {k: {"tan_theta": v} for k, v in TIERS.items()},
        "constraints": {
            "min_cover_m": MIN_COVER_M, "surface_z_m": SURF_Z,
            "tz_thickness_m": TZ_THICK_M,
            "tz_material_assumed": {"name": TZ_MATERIAL[0], "eps_r": TZ_MATERIAL[1],
                                    "sigma": TZ_MATERIAL[2],
                                    "basis": "覆盖层(16,0.01)与砂岩(9,0.001)之间的插值假设，非实测"},
        },
        "required_material_line": "#material: 12 0.005 1 0 tzone",
        "variants": {},
        "hard_limits": [
            "设计草案未冻结，不构成执行契约；执行前须另行冻结 gate 与用户签认",
            "方案 A 地表水平、天线 z=45 不变；仅基覆界面倾斜，不与真实地形的地表同步倾斜混淆",
            "T2/T3 档界面中心下移（覆盖层中心厚于 3 m 锚点），跨档比较须标注该几何差异",
            "倾角档位参考 2026-09-27_cover_bedrock_interface_research.md 的片段级文献值，非场地实测",
            "过渡带材料为插值假设；倾斜档的 BG/NC 须随界面几何各自重建（背景配对不再与水平族共用）",
            "不生成物理/训练标签；本扩展面向开发族 {C1,C3}，测试族 {C5,C8} 已消费不得据此调参",
        ],
    }
    for gname, g in GEOMS.items():
        for tk, tt in TIERS.items():
            for tz in (False, True):
                key = "%s_%s%s" % (gname, tk, "TZ" if tz else "")
                out["variants"][key] = build(tk, tt, tz, g)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1),
                   encoding="utf-8", newline="\n")
    for key, v in out["variants"].items():
        print("%s: n=%d zc=%.2f if_range=%s cover=%s theta_eff=%.4f maxdev=%.4f" % (
            key, v["n_strips"], v["zc_m"], v["z_interface_range"],
            v["cover_thickness_range"], v["theta_eff_deg"],
            v["max_quantization_deviation_m"]))
    print("written:", OUT.relative_to(REPO))


if __name__ == "__main__":
    sys.exit(main())
