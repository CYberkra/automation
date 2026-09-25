"""dep3d_gold_v1 静态预算与网格对齐核验（纯自算，不调用 gprMax 求解器）。

两种模式：
  python scripts/check_dep3d_gold_static.py                 # 只打印核验结果
  python scripts/check_dep3d_gold_static.py --write-budget  # 额外写入 budget.json

预算方法见 docs/research/2026-09-26_dep3d_gold_design.md 第 7 节；
dt 公式与 gprMax V4 grid/fdtd_grid.py::FDTDGrid.calculate_dt 一致
（2D TMx 只用 dy、dz；之后按 round_float(prec-1) 向下舍入），
并用 histories/round 的 12 位 str repr 做"%s" 输出对齐性检查。
"""

from __future__ import annotations

import argparse
import decimal as d
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
CASE_DIR = REPO / "configs" / "research" / "dep3d_gold_v1"

C = 299792458.0
TIME_WINDOW = 1200e-9

# 已敲定的人为常数字段（设计书 §7）
BYTES_PER_CELL_USER = 131            # 用户给定
BYTES_PER_CELL_DIM23 = 122           # dim23 3D 实测 3.8 GB / 3120 万
BYTES_PER_CELL_2D = 183              # DEP 2D 实测 7.5 GB / 4096 万（悲观上界，维度不可比）
CELL_STEPS_PER_S = 792e6             # 保守下界（dim23 3D 反算）
JOB_RATIO = (1.9, 2.1)               # Job 提交 / 设备内存 实测比
CELL_CAP_3D = 64_000_000             # §11-3 批准
BUILD_OVERHEAD_S = (158.0, 173.0)    # dim23 3D 实测几何构建/编译/写盘开销

# dt/步数口径的三条交叉验证基准（均已跑过求解器，日志在 artifacts/research_checks/）
KNOWN_STEPS = [
    # 名称, 维度, 参与 CFL 的格距, 时窗(s), 求解器 stdout 报的 iterations
    ("dim23 3D 5cm 800ns", "3D", (0.05, 0.05, 0.05), 800e-9, 8310),
    ("dim23 2D 5cm 800ns", "2D", (0.05, 0.05), 800e-9, 6785),
    ("DEP 2D ZFINE2 1200ns", "2D", (0.0125, 0.003125), 1200e-9, 118665),
]

NUMBER = re.compile(r"^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$")


def round_value(value: float, decimalplaces: int) -> float:
    """与 gprMax utilities.round_float 一致：向下舍入到指定小数位。"""
    precision = f"1.{'0' * decimalplaces}"
    return float(d.Decimal(value).quantize(d.Decimal(precision), rounding=d.ROUND_FLOOR))


def parse_in(path: Path) -> dict:
    """提取内核 YAML 片段（忽略 #python:/#end_python: 块）。"""
    lines, keep = [], []
    in_py = False
    for ln, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        s = raw.strip()
        if s.startswith("#python:"):
            in_py = True
            continue
        if in_py and s.startswith("#end_python:"):
            in_py = False
            continue
        if not in_py:
            lines.append(raw)
            keep.append(ln)
    cmds: list[tuple[int, list[str]]] = []
    for ln, raw in zip(keep, lines):
        s = raw.strip()
        if not s or s.startswith("##"):
            continue
        if not s.startswith("#"):
            raise ValueError(f"{path.name}:{ln} 非 '# ' 前缀条目: {s!r}")
        body = s[1:].strip()
        head, _, param = body.partition(":")
        cmds.append((ln, [head.strip()] + param.strip().split()))
    return {"name": path.name, "cmds": cmds, "text": "\n".join(lines)}


def get_cmd(doc: dict, key: str, n: int) -> list[str]:
    hits = [p for _, c in doc["cmds"] if c[0] == key for p in [c[1:]]]
    if len(hits) != 1:
        raise ValueError(f"{doc['name']}: # {key} 出现 {len(hits)} 次（期望 1）")
    if len(hits[0]) != n:
        raise ValueError(f"{doc['name']}: # {key} 参数个数 {len(hits[0])} != {n}")
    return hits[0]


def as_float(tok: str, doc: dict, key: str) -> float:
    if not NUMBER.match(tok):
        raise ValueError(f"{doc['name']}: # {key} 非法数值 {tok!r}")
    return float(tok)


def extract(path: Path) -> dict:
    doc = parse_in(path)
    domain = get_cmd(doc, "domain", 3)
    dstep = get_cmd(doc, "dx_dy_dz", 3)
    mode_cmds = [c[1:] for _, c in doc["cmds"] if c[0] == "domain_mode"]
    mode = mode_cmds[0][0] if mode_cmds else ""

    dxdydz = [as_float(t, doc, "dx_dy_dz") for t in dstep]
    dims: list[float | None] = []
    for tok in domain:
        if tok == "inf":
            dims.append(None)
        else:
            dims.append(as_float(tok, doc, "domain"))

    # PML 含在域内（gprMax _validate_pml_thickness 按 nx/ny/nz 校验）
    pml = [int(t) for t in get_cmd(doc, "pml_cells", 6)]

    src = get_cmd(doc, "hertzian_dipole", 5)     # pol x y z waveform
    rx = get_cmd(doc, "rx", 5)                   # x y z measurement Ex
    boxes = [c[1:-1] for _, c in doc["cmds"] if c[0] == "box"]      # x0 y0 z0 x1 y1 z1 <id>
    mats = [c[1:] for _, c in doc["cmds"] if c[0] == "material"]

    n_pml = [pml[0] + pml[3], pml[1] + pml[4], pml[2] + pml[5]]

    return {
        "name": path.stem,
        "title": " ".join(" ".join(c[1:]) for _, c in doc["cmds"] if c[0] == "title"),
        "mode": mode,
        "domain": domain,
        "dims": dims,
        "steps_tuple": tuple(dxdydz),
        "pml": pml,
        "n_pml_axis": n_pml,
        "src_pol": src[0],
        "src": src[1:4],
        "rx": rx[0:3],
        "rx_field": rx[4],
        "boxes": boxes,
        "materials": mats,
        "time_window": get_cmd(doc, "time_window", 1)[0],
        "pml_formulation": get_cmd(doc, "pml_formulation", 1)[0],
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def dt_from_spacing(spacings: tuple[float, ...]) -> tuple[float, float]:
    """返回 (原始 dt, gprMax round_float 之后的 dt)。2D TMx 只用 (dy, dz)。"""
    ss = sum((1.0 / s) ** 2 for s in spacings)
    raw = 1.0 / (C * math.sqrt(ss))
    return raw, round_value(raw, d.getcontext().prec - 1)


def analytic(cells_axis: list[int], n_pml_axis: list[int]) -> tuple[int, int]:
    """返回 (主网格单元数, PML 格数)；PML 含在域内，故用补集计数。"""
    cells = 1
    interior = 1
    for n, p in zip(cells_axis, n_pml_axis):
        cells *= n
        interior *= max(0, n - p)
    return cells, cells - interior


def check_step_rule() -> int:
    """用三条已跑过的求解器日志校验 dt 口径与 ceil+1 步数规则。"""
    print("=" * 78)
    print("A0. dt / 步数口径回归（gprMax cmds_singleuse.py TimeWindow.build: "
          "iterations = ceil(time/dt) + 1）")
    print("=" * 78)
    bad = 0
    for name, dim, sp, tw, expected in KNOWN_STEPS:
        _, dt = dt_from_spacing(sp)
        got = int(np.ceil(tw / dt)) + 1
        ok = got == expected
        bad += 0 if ok else 1
        print(f"[{'PASS' if ok else 'FAIL'}] {name:22s} dt={dt*1e12:>9.4f} ps  "
              f"复算={got:<7d} 求解器={expected}")
    return bad


def grid_alignment(e: dict) -> list[str]:
    """逐坐标核验：坐标 / 对应格距 必须为整数格（相对域原点 0）。"""
    errs: list[str] = []
    d = e["steps_tuple"]

    def chkx(label: str, coord: str, axis: int):
        if coord == "inf":
            return
        v = float(coord)
        if e["dims"][axis] is None:
            errs.append(f"{e['name']}: {label} 给了有限坐标 {v}，但该轴域为 inf")
            return
        q = v / d[axis]
        if abs(q - round(q)) > 1e-9:
            errs.append(
                f"{e['name']}: {label} = {v} 不在格距 {d[axis]} 上（商 {q:.6f}）"
            )
        gi = int(round(q))
        if gi < 0 or gi > int(round(e["dims"][axis] / d[axis])):
            errs.append(f"{e['name']}: {label} = {v} 超出域范围")

    for tok, ax in zip(e["src"], (0, 1, 2)):
        chkx(f"源 {tok}", tok, ax)
    for tok, ax in zip(e["rx"], (0, 1, 2)):
        chkx(f"接收 {tok}", tok, ax)
    for b in e["boxes"]:
        toks = b[0:6]
        name = b[6] if len(b) > 6 else "box"
        for tok, ax in zip(toks, (0, 1, 2, 0, 1, 2)):
            chkx(f"{name} 坐标 {tok}", tok, ax)

    # PML 厚度：有限方向必须各 1 m（格数 x 格距），且 2*thickness < 该轴格数
    # （gprMax FDTDGrid._validate_pml_thickness 会据此抛 ValueError）
    names = ("x0", "y0", "z0", "xmax", "ymax", "zmax")
    for i, (a, b) in enumerate(((0, 3), (1, 4), (2, 5))):
        if e["dims"][i] is None:
            continue
        ncells = int(round(e["dims"][i] / d[i]))
        for j in (a, b):
            t = e["pml"][j]
            thick = t * d[i]
            if abs(thick - 1.0) > 1e-9:
                errs.append(
                    f"{e['name']}: PML {names[j]} = {t} 格 x {d[i]} m = {thick:.4f} m != 1 m"
                )
            if 2 * t >= ncells:
                errs.append(
                    f"{e['name']}: PML {names[j]} = {t} 格，2*t={2*t} >= 该轴格数 {ncells}"
                )
    return errs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-budget", action="store_true")
    args = ap.parse_args()

    files = sorted(CASE_DIR.glob("*.in"))
    if len(files) != 8:
        print(f"期望 8 份 .in，实际 {len(files)}", file=sys.stderr)
        return 2

    total_bad = check_step_rule()

    print()
    print("=" * 78)
    print("A. 坐标网格对齐核验（所有有限坐标必须整除对应格距，并落在域内）")
    print("=" * 78)
    all_errs: list[str] = []
    ex = {}
    for f in files:
        e = extract(f)
        ex[f.stem] = e
        errs = grid_alignment(e)
        all_errs += errs
        ok = "PASS" if not errs else "FAIL"
        print(f"[{ok}] {f.name:20s} domain={e['domain']}  step={e['steps_tuple']}")
        for m in errs:
            print(f"        {m}")
    print(f"\n合计: {len(files)} 份输入，{len(all_errs)} 项对齐问题")

    # ---- 预算（设计书 §7 方法） ----
    print()
    print("=" * 78)
    print("B. 预算自算（设计书 §7 方法；estimate，非承诺）")
    print("=" * 78)
    budget_cases, rows = [], []
    for name in [
        "DEP3D_5CM_BG", "DEP3D_5CM_TGT",
        "DEP3D_ANISO_BG", "DEP3D_ANISO_TGT",
        "B2D5CM_BG", "B2D5CM_TGT",
        "B2DANISO_BG", "B2DANISO_TGT",
    ]:
        e = ex[name]
        dx, dy, dz = e["steps_tuple"]
        is3d = not e["mode"]
        if not is3d and e["domain"][0] != "inf":
            raise ValueError(f"{name}: 2D 情形 x 轴应为 inf")
        dt_raw, dt = dt_from_spacing((dx, dy, dz) if is3d else (dy, dz))

        n_cells_axis = []
        for i, dimv in enumerate(e["dims"]):
            if dimv is None:
                n_cells_axis.append(1)          # inf 方向 1 格
            else:
                n_cells_axis.append(int(round(dimv / e["steps_tuple"][i])))
        nx, ny, nz = n_cells_axis
        cells, pml_cells = analytic(n_cells_axis, e["n_pml_axis"])

        tw = float(e["time_window"])
        nsteps = int(np.ceil(tw / dt)) + 1

        cell_steps = cells * nsteps
        solve_s = cell_steps / CELL_STEPS_PER_S
        solve_s_fast = solve_s / 2.0

        rec = {
            "run_id": name,
            "file": name + ".in",
            "sha256": e["sha256"],
            "dimension": "3D" if is3d else "2D TMx",
            "domain_m": [None if v is None else v for v in e["dims"]],
            "grid_spacing_m": [dx, dy, dz],
            "n_cells_axis": [nx, ny, nz],
            "cells": cells,
            "pml_cells": pml_cells,
            "pml_fraction_of_main_grid": pml_cells / cells,
            "dt_s": dt,
            "dt_ps": dt * 1e12,
            "time_window_s": tw,
            "time_steps_analytic_raw": tw / dt,
            "time_steps": nsteps,
            "cell_steps": cell_steps,
            "device_memory_GiB_at_131B_per_cell": cells * BYTES_PER_CELL_USER / 2**30,
            "device_memory_GiB_at_122B_per_cell": cells * BYTES_PER_CELL_DIM23 / 2**30,
            "device_memory_GiB_at_183B_per_cell": cells * BYTES_PER_CELL_2D / 2**30,
            "job_commit_GiB_estimate_range": [
                cells * BYTES_PER_CELL_USER / 2**30 * JOB_RATIO[0],
                cells * BYTES_PER_CELL_USER / 2**30 * JOB_RATIO[1],
            ],
            "job_commit_GiB_pessimistic": cells * BYTES_PER_CELL_2D / 2**30 * JOB_RATIO[1],
            "solve_s_at_792M_cellsteps_per_s": solve_s,
            "solve_s_at_2x_rate": solve_s_fast,
            "wall_clock_s_estimate_range": [
                solve_s + BUILD_OVERHEAD_S[0], solve_s + BUILD_OVERHEAD_S[1]
            ],
            "within_3d_cell_cap": cells <= CELL_CAP_3D if is3d else None,
        }
        budget_cases.append(rec)
        rows.append(rec)

    hdr = (f"{'run_id':16s} {'cells':>12s} {'nx x ny x nz':>16s} {'dt_ps':>10s} "
           f"{'steps':>7s} {'s@792M':>9s} {'GiB@131B':>9s} {'Job GiB':>9s}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(
            f"{r['run_id']:16s} {r['cells']:>12,d} "
            f"{r['n_cells_axis'][0]}x{r['n_cells_axis'][1]}x{r['n_cells_axis'][2]}"
            f"{'  ' + r['dimension'] if r['dimension'] != '3D' else '':<8s}"
            f"{r['dt_ps']:>10.4f} "
            f"{r['time_steps']:>7,d} {r['solve_s_at_792M_cellsteps_per_s']:>9.1f} "
            f"{r['device_memory_GiB_at_131B_per_cell']:>9.2f} "
            f"{r['job_commit_GiB_estimate_range'][0]:>5.1f}-{r['job_commit_GiB_estimate_range'][1]:<5.1f}"
        )

    total = sum(r["cells"] for r in rows)
    total_s = sum(r["wall_clock_s_estimate_range"][1] for r in rows)
    print("-" * len(hdr))
    print(f"{'内核总量':16s} {total:>12,d}   墙钟上界合计 {total_s/60:>7.1f} min（8 例串行）")

    cap_violations = [r["run_id"] for r in rows if r["within_3d_cell_cap"] is False]
    print(f"3D 单元上限 {CELL_CAP_3D:,}：{'全部满足' if not cap_violations else cap_violations}")

    if args.write_budget:
        out = CASE_DIR / "budget.json"
        payload = {
            "batch": "dep3d_gold_v1",
            "generated_by": "scripts/check_dep3d_gold_static.py",
            "method": "docs/research/2026-09-26_dep3d_gold_design.md §7",
            "notes": [
                "cells = prod(domain/step)，PML 含在域内（gprMax _validate_pml_thickness 依据）。",
                "dt = 1/(c*sqrt(sum 1/d_i^2))；2D TMx 只用 dy、dz；随后按 gprMax round_float(prec-1) 向下舍入。",
                "time_steps = ceil(time_window/dt) + 1。该式取自 gprMax 源码 "
                "gprMax/user_objects/cmds_singleuse.py TimeWindow.build；本脚本以三条已跑过的求解器日志回归通过"
                "（dim23 3D 800ns->8310、dim23 2D 800ns->6785、DEP 2D ZFINE2 1200ns->118665）。"
                "设计书 2026-09-26 早期草稿所写的 12463/17625 未计入 round_float 后的 dt 差异，"
                "正确值为 12464/17626，已据源码与上述三条日志更正。",
                "cell-steps/s = 792e6 为保守下界（dim23 3D 反算，RTX 3060 Laptop 实测；本机 RTX 4090 Laptop 按 1-2 倍估计）。",
                "设备内存三档：131 B/单元（用户给定）、122 B/单元（dim23 3D 实测）、183 B/单元（DEP 2D 实测，维度不可比，悲观上界）。",
                "Job 提交内存 = 设备内存 x 1.9-2.1（dim23 3D 7.25 GiB、DEP 2D 15.27 GiB 实测比）。",
                "墙钟估计含 dim23 3D 实测 +158~173 s 的几何构建/编译/写盘开销；2D 对照未单独校准开销。",
                "全部为自算估计，非承诺；实际值以执行时监督器记录为准。",
            ],
            "constants": {
                "c_m_s": C,
                "time_window_s": TIME_WINDOW,
                "bytes_per_cell_user": BYTES_PER_CELL_USER,
                "bytes_per_cell_dim23_3d_measured": BYTES_PER_CELL_DIM23,
                "bytes_per_cell_dep2d_measured_pessimistic": BYTES_PER_CELL_2D,
                "cell_steps_per_s_conservative_lower_bound": CELL_STEPS_PER_S,
                "job_commit_over_device_ratio_range": list(JOB_RATIO),
                "cell_cap_3d": CELL_CAP_3D,
                "build_overhead_s_range": list(BUILD_OVERHEAD_S),
                "resource_caps_per_case": {
                    "wall_minutes": 40,
                    "job_commit_GiB": 20,
                    "free_ram_GiB_preflight": 24,
                    "free_vram_GiB_preflight": 8,
                    "output_GiB": 1,
                    "retries": 0,
                    "fdtd_runs": 1,
                    "threads": 8,
                    "backend": "CUDA",
                    "precision": "double",
                    "device": 0,
                    "cpu_fallback": False,
                },
            },
            "totals": {
                "cases": len(rows),
                "total_cells": total,
                "total_wall_clock_upper_s": total_s,
                "execution_order": [r["run_id"] for r in rows],
            },
            "cases": rows,
        }
        out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"\n已写入 {out.relative_to(REPO)}")

    return 1 if (all_errs or total_bad) else 0


if __name__ == "__main__":
    sys.exit(main())
