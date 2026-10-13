"""Independently audit and compare the authorized V3.1 diagnostic study.

Consumes unchanged official SFCW outputs; no FDTD, filtering or gain is run.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from study_line9_v31_smoothing_domain import load, save, sha


def relative(candidate, reference):
    den = float(np.linalg.norm(reference))
    return float(np.linalg.norm(candidate - reference) / den) if den > 0 else None


def read_sf(path, index=None, frequency_spec=None):
    with h5py.File(path) as h:
        def take(name):
            x = h[name][:]
            return x if index is None else x[:, index]
        out = dict(f=h["frequency"][:], H=take("response"), t=h["time_response/time"][:] * 1e9,
                   env=take("time_response/complex_envelope"), real=take("time_response/real_bandpass"),
                   weights=h["time_response/weights"][:], tail_db=float(h.attrs["ReceiverTailRelativeDB"]))
        assert h.attrs["Method"] == "direct"
        assert h["time_response"].attrs["Window"] == "gaussian"
        assert h["time_response"].attrs["ZeroPadFactor"] == 4
        assert h.attrs["TailTaperFraction"] == 0 and h["time_response"].attrs["TimeShift"] == 0
        assert np.all(h["source_valid"][:])
    start, stop, count = frequency_spec or (20e6, 170e6, 501)
    assert np.array_equal(out['f'], np.linspace(start, stop, count))
    assert out['env'].shape == (4 * count,) and out['H'].shape == (count,)
    assert np.isfinite(out["env"]).all() and np.isfinite(out["H"]).all()
    return out


def independent_check(raw, sf):
    with h5py.File(raw) as h:
        ds = h["rxs/rx1/Ez"]
        y = np.asarray(ds[:], dtype=float)
        excitation = h["srcs/src1/excitation"]
        x = np.asarray(excitation["samples"][:], dtype=float)
        dt = float(ds.attrs["SampleInterval"])
        ry = np.arange(len(y)) * dt + float(ds.attrs["TimeSampleOffset"])
        sx = np.arange(len(x)) * float(excitation.attrs["SampleInterval"]) + float(excitation.attrs["TimeSampleOffset"])
    count = len(sf['f'])
    idx = np.linspace(0, count - 1, 11, dtype=int)
    numerical = np.array([np.sum(y * np.exp(-2j * np.pi * f * ry)) /
                          np.sum(x * np.exp(-2j * np.pi * f * sx)) for f in sf["f"][idx]])
    dft_error = relative(numerical, sf["H"][idx])
    tidx = (np.unique(np.r_[np.arange(0, 2004, 67), 2003]) if count == 501
            else np.linspace(0, len(sf['t']) - 1, 31, dtype=int))
    recon = np.array([np.sum(sf["H"] * sf["weights"] * np.exp(2j * np.pi *
                     (sf["f"] - sf["f"][0]) * sf["t"][i] * 1e-9)) / count for i in tidx])
    inverse_error = relative(recon, sf["env"][tidx])
    assert dft_error < 1e-10 and inverse_error < 1e-10
    return dict(independent_11_tone_DFT_relative_l2=dft_error,
                independent_31_time_inverse_relative_l2=inverse_error)


def audit(root, stage):
    contract = load(root / "contract.json")
    started = load(root / "execution_start.json")
    assert sha(root / "contract.json") == started["contract_sha256"]
    assert started["gprMax"] == "4.0.1" and "4090" in started["gpu"]
    stages = ["A"] if stage == "A" else ["A", "B"] if stage == "B" else ["A", "B", "C"]
    results, receipts, checks = {}, {}, {}
    for s in stages:
        complete = load(root / f"stage_{s}_complete.json")
        assert complete["status"] == "COMPLETE"
        for r in complete["receipts"]:
            assert r["status"] == "PASS"
            for f in r["files"]:
                p = root / f["path"]
                assert p.stat().st_size == f["bytes"] and sha(p) == f["sha256"], p
            receipts[r["case"]] = r
            case = next(c for c in contract["cases"] if c["id"] == r["case"])
            if case["geometry_only"]:
                continue
            result = read_sf(root / "cases" / r["case"] / "sfcw.h5")
            checks[r["case"]] = independent_check(root / "cases" / r["case"] / "impulse.h5", result)
            results[r["case"]] = result
    recovery = load(root / "audit_recovery.json")
    assert sha(root / "cases/full_x150.0_y/impulse.h5") == recovery["raw_sha256"]
    assert sha(root / "cases/full_x150.0_y/solver.log") == recovery["log_sha256"]
    return contract, results, receipts, checks


def compare(candidate, reference, case, screen):
    assert np.array_equal(candidate["t"], reference["t"])
    assert np.array_equal(candidate["f"], reference["f"])
    strong = np.abs(reference["H"]) >= np.abs(reference["H"]).max() * .001
    phases = np.angle(candidate["H"][strong] / reference["H"][strong], deg=True)
    spec_error = relative(candidate["H"], reference["H"])
    windows = [{"label": "common", "lo_ns": lo, "hi_ns": hi} for lo, hi in screen["common_windows_ns"]]
    for ev in case["vertical_reference_events"]:
        center = ev["time_ns"]
        if 0 < center < 1000:
            windows.append(dict(label=ev["label"], lo_ns=center-screen["event_halfwidth_ns"], hi_ns=center+screen["event_halfwidth_ns"]))
    stats = []
    for win in windows:
        mask = (reference["t"] >= win["lo_ns"]) & (reference["t"] < win["hi_ns"])
        t = reference["t"][mask]
        r, q = reference["env"][mask], candidate["env"][mask]
        rp, qp = float(np.abs(r).max()), float(np.abs(q).max())
        stats.append({**win, "complex_relative_l2": relative(q, r),
                      "envelope_peak_relative_error": abs(qp/rp - 1) if rp else None,
                      "envelope_peak_shift_ns": float(t[np.argmax(np.abs(q))] - t[np.argmax(np.abs(r))])})
    phase_rms = float(np.sqrt(np.mean(phases**2)))
    passed = (spec_error is not None and spec_error <= screen["complex_band_relative_l2_max"] and
              phase_rms <= screen["strong_bin_phase_rms_degrees_max"] and all(
                  s["complex_relative_l2"] is not None and s["complex_relative_l2"] <= screen["time_window_relative_l2_max"] and
                  s["envelope_peak_relative_error"] is not None and s["envelope_peak_relative_error"] <= screen["event_peak_relative_error_max"] and
                  abs(s["envelope_peak_shift_ns"]) <= screen["event_peak_shift_ns_max"] for s in stats))
    return dict(complex_band_relative_l2=spec_error, strong_bin_phase_rms_degrees=phase_rms,
                windows=stats, passes_frozen_domain_screen=bool(passed))


def coefficients(root, cid):
    path = root / "cases" / cid
    with h5py.File(path / "surface_yee.h5") as h:
        keys = [k.decode() for k in h["material_keys"][:]]
        ids = h["ID"][2, :, :, 0]
        voxels = h["data"][:, :, 0]
        spacing = float(h.attrs["dx_dy_dz"][0])
    db = load(path / "surface_yee_materials.json")["materials"]
    er = np.array([db[k]["base"]["relative_permittivity"] for k in keys], dtype=float)
    sigma = np.array([db[k]["base"]["electric_conductivity_s_per_m"] for k in keys], dtype=float)
    return dict(epsilon=er[ids], sigma=sigma[ids], voxels=voxels, dl=spacing,
                mixed_ez_edge_count=int(np.sum((er[ids] > 1) & (er[ids] < 12))),
                unique_ez_epsilon=np.unique(er[ids]).tolist())


def sparse(ax, values, stations, t, title, lim):
    dt = float(t[1]-t[0])
    for i, station in enumerate(stations):
        ax.imshow(values[:, i, None], origin="upper", aspect="auto", cmap="gray", vmin=-lim, vmax=lim,
                  interpolation="nearest", extent=[station-.3, station+.3, t[-1]+dt/2, t[0]-dt/2])
    ax.set_xlim(149.5, 198.5)
    ax.set_ylim(t[-1], t[0])
    ax.set(title=title, ylabel="时间（ns）", xlabel="测线里程（m）；空白＝未计算")


def plot_sparse_compare(out, name, stations, reference, candidate, time, title, reference_events=None):
    mask = (time >= 20) & (time < 450)
    time = time[mask]
    reference, candidate = reference[mask], candidate[mask]
    lim = float(np.quantile(np.abs(reference[(time >= 100)]), .995))
    fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True, constrained_layout=True)
    sparse(axes[0], reference, stations, time, "基准总场｜各站共同色标", lim)
    sparse(axes[1], candidate, stations, time, "候选总场｜与基准同一色标", lim)
    sparse(axes[2], candidate-reference, stations, time, "候选－基准｜同一色标；此差量不称噪声", lim)
    if reference_events is not None:
        points = [(s, ev['time_ns']) for s, events in reference_events.items() for ev in events
                  if time[0] <= ev['time_ns'] <= time[-1]]
        for ax in axes:
            ax.scatter([x for x,y in points], [y for x,y in points], marker='_', color='#d88900',
                       s=70, linewidths=1, label='基准体素垂向双程时参考；非唯一传播路径')
            ax.legend(loc='upper right', fontsize=8)
    fig.suptitle(title+f"\n官方 direct SFCW real bandpass；±{lim:.4g} 共同显示色标；没有增益、归一化或背景抑制", fontsize=12)
    fig.savefig(out / name, dpi=170); plt.close(fig)


def analyze_A(a, contract, results, receipts):
    anchors = contract["anchors_m"]
    old_complete = load(a.baseline / "completion.json")
    oldsf = a.baseline / "results/line9_v31_bscan_sfcw.h5"
    assert sha(oldsf) == old_complete["official_sfcw"]["sha256"]
    rows, base, on = [], [], []
    for station in anchors:
        case = next(c for c in contract["cases"] if c["id"] == f"full_x{station:.1f}_y")
        ref = read_sf(oldsf, case["baseline_trace_index"]-1)
        cand = results[case["id"]]
        old_receipt = old_complete["receipts"][case["baseline_trace_index"]-1]
        oldraw = a.baseline / "results" / f"line9_v31_bscan{case['baseline_trace_index']}.h5"
        assert sha(oldraw) == old_receipt["output_sha256"]
        with h5py.File(oldraw) as old, h5py.File(a.root / "cases" / case["id"] / "impulse.h5") as new:
            assert np.array_equal(old["srcs/src1/excitation/samples"][:], new["srcs/src1/excitation/samples"][:])
            assert np.array_equal(old["rxs/rx1"].attrs["Position"], new["rxs/rx1"].attrs["Position"])
            values = old["rxs/rx1/Ez"][:]
            n_tail = int(np.ceil(.05 * len(values)))
            baseline_tail = float(20*np.log10(float(np.abs(values[-n_tail:]).max())/float(np.abs(values).max())))
        row = compare(cand, ref, case, contract["domain_screen"])
        rows.append(dict(station_m=station, **row, baseline_tail_db=baseline_tail,
                         smoothed_tail_db=receipts[case["id"]]["tail_db"],
                         tail_change_db=receipts[case["id"]]["tail_db"]-baseline_tail))
        base.append(ref["real"]); on.append(cand["real"])
    plot_sparse_compare(a.out, "smoothing_sparse_bscan.png", anchors, np.array(base).T, np.array(on).T,
                        ref["t"], "V3.1 7 个同站位｜原生平均 n → y｜代表站位对照，非完整连续扫描")
    save(a.out / "full_model_smoothing.json", rows)

    normalized = {}
    for dl in [.025, .0125]:
        air = results[f"air_d{dl:g}"]["H"]
        assert np.min(np.abs(air)) > np.max(np.abs(air)) * 1e-8
        for avg in ["n", "y"]:
            normalized[f"d{dl:g}_{avg}"] = (results[f"surface_d{dl:g}_{avg}"]["H"] - air) / air
    fine_y = normalized["d0.0125_y"]
    errors = {k: relative(v, fine_y) for k, v in normalized.items()}
    fine_spread = relative(normalized["d0.0125_n"], fine_y)
    fine_n = normalized["d0.0125_n"]
    result = dict(diagnostic_station_m=contract["surface_diagnostic_station_m"],
                  definition="R=(H_surface-H_air)/H_air; dimensionless surface scattering diagnostic",
                  relative_l2_to_fine_smoothed_proxy=errors,
                  fine_n_vs_y_relative_l2=fine_spread,
                  coarse_n_to_fine_n_relative_l2=relative(normalized["d0.025_n"], fine_n),
                  coarse_y_to_fine_n_relative_l2=relative(normalized["d0.025_y"], fine_n),
                  coarse_y_closer_to_fine_y_than_coarse_n=errors["d0.025_y"] < errors["d0.025_n"],
                  error_reduction_factor=errors["d0.025_n"] / errors["d0.025_y"] if errors["d0.025_y"] else None,
                  warning="Fine proxy shares numerical assumptions and is not an analytic truth. The n/y difference is not itself isolated nonphysical scattering.")
    save(a.out / "surface_only_grid_reference.json", result)
    labels = {"d0.025_n": "2.5 cm，关闭平滑", "d0.025_y": "2.5 cm，原生平滑",
              "d0.0125_n": "1.25 cm，关闭平滑", "d0.0125_y": "1.25 cm，原生平滑参考"}
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), constrained_layout=True)
    for k, v in normalized.items():
        axes[0].plot(ref["f"] / 1e6, 20*np.log10(np.maximum(np.abs(v),1e-15)), label=labels[k])
        if k != "d0.0125_y":
            axes[1].plot(ref["f"] / 1e6, np.abs(v-fine_y)/max(np.linalg.norm(fine_y)/np.sqrt(501),1e-15)*100, label=labels[k])
    axes[0].set(title="同一地表、地下均为原覆盖层｜空气参照消除二维源网格缩放", ylabel="|R|（dB）")
    axes[1].set(title="相对 1.25 cm、平滑 y 代理的差异｜不称解析误差", xlabel="频率（MHz）", ylabel="逐频差／代理 RMS（%）")
    for ax in axes: ax.legend(); ax.grid(alpha=.2)
    fig.savefig(a.out / "surface_only_proxy.png", dpi=170); plt.close(fig)

    n = coefficients(a.root, "full_x174.0_n_geometry_audit")
    y = coefficients(a.root, "full_x174.0_y")
    assert np.array_equal(n["voxels"], y["voxels"])
    audit = dict(same_voxel_geometry=True, n_unique_ez_epsilon=n["unique_ez_epsilon"],
                 y_unique_ez_epsilon=y["unique_ez_epsilon"], n_mixed_edges=n["mixed_ez_edge_count"],
                 y_mixed_edges=y["mixed_ez_edge_count"], coefficient_changes=int(np.sum(n["epsilon"] != y["epsilon"])))
    assert audit["coefficient_changes"] > 0
    save(a.out / "native_coefficient_audit.json", audit)
    audit_case = next(c for c in contract["cases"] if c["id"] == "full_x174.0_y")
    lower, upper = audit_case["surface_coefficient_export_bounds"]
    ox, oy = audit_case["origin_global_x_elevation_m"]
    half = n["dl"] / 2
    extent = [ox + lower[0] - half, ox + upper[0] + half,
              oy + lower[1] - half, oy + upper[1] + half]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, field, title in zip(axes, [n["epsilon"], y["epsilon"]], ["n：未做介电平均", "y：原生介电平均"]):
        im=ax.imshow(field.T, origin="lower", aspect="auto", extent=extent, vmin=1, vmax=12, cmap="viridis", interpolation="nearest")
        ax.set(title=title, xlabel="地质里程（m）", ylabel="高程（m）")
    fig.colorbar(im, ax=axes, label="Ez 边上的有效相对介电常数")
    fig.suptitle("174 m 地表小块：体素不变，Yee 系数改变")
    fig.savefig(a.out / "native_smoothing_coefficients.png", dpi=170); plt.close(fig)
    return dict(full_model=rows, surface_only=result, coefficients=audit)


def analyze_domains(a, contract, results, receipts):
    screen = contract["domain_screen"]
    stage = a.stage
    controls, rows = [], []
    for case in contract["cases"]:
        if case["family"] not in ["expanded_domain", "reduced_domain"] or case["id"] not in results:
            continue
        station=case["station_m"]
        baseid=f"full_x{station:.1f}_y"
        baseline=results[baseid]
        cand=results[case["id"]]
        # Use the larger domain as reference for its guard comparison.
        measurement=compare(baseline,cand,case,screen) if case["family"]=="expanded_domain" else compare(cand,baseline,case,screen)
        coefficient_ref=coefficients(a.root,baseid)
        coefficient_cand=coefficients(a.root,case["id"])
        same_coeff=np.array_equal(coefficient_ref["epsilon"],coefficient_cand["epsilon"]) and np.array_equal(coefficient_ref["sigma"],coefficient_cand["sigma"])
        assert same_coeff,"retained surface coefficient model changed with crop"
        row=dict(case=case["id"],family=case["family"],station_m=station,domain_m=case["domain_m"],
                 same_retained_surface_coefficients=same_coeff,
                 fdtd_speedup_vs_full=receipts[baseid]["official_fdtd_seconds"]/receipts[case["id"]]["official_fdtd_seconds"],
                 whole_run_speedup_vs_full=receipts[baseid]["official_build_and_solver_seconds"]/receipts[case["id"]]["official_build_and_solver_seconds"],**measurement)
        (controls if case["family"]=="expanded_domain" else rows).append(row)
    guard_pass=all(r["passes_frozen_domain_screen"] for r in controls) and len(controls)==3
    summary=[]
    for dimensions in sorted({tuple(r["domain_m"]) for r in rows}):
        group=[r for r in rows if tuple(r["domain_m"])==dimensions]
        summary.append(dict(domain_m=list(dimensions),stations=len(group),
                            all_station_screens_pass=len(group)==7 and all(r["passes_frozen_domain_screen"] for r in group),
                            expanded_domain_guard_pass=guard_pass,
                            accepted_candidate=len(group)==7 and all(r["passes_frozen_domain_screen"] for r in group) and guard_pass,
                            max_complex_band_error=max(r["complex_band_relative_l2"] for r in group),
                            max_window_error=max(w["complex_relative_l2"] for r in group for w in r["windows"]),
                            median_fdtd_speedup=float(np.median([r["fdtd_speedup_vs_full"] for r in group])),
                            median_whole_run_speedup=float(np.median([r["whole_run_speedup_vs_full"] for r in group]))))
    combined=next(s for s in summary if s["domain_m"]==[40.0,45.0])
    save(a.out/"stage_C_screen.json",dict(proceed=bool(combined["accepted_candidate"]),condition=contract["stage_C_condition"],
                                         guard_pass=guard_pass,combined=combined,contract_sha256=sha(a.root/"contract.json")))
    save(a.out/"domain_comparisons.json",dict(expanded_controls=controls,reduced_domains=rows,summary=summary,
                                            interpretation="equivalence to a finite-domain numerical reference, not absolute physics validation"))
    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    for item in summary:
        label=f"{item['domain_m'][0]:g}×{item['domain_m'][1]:g} m"
        group=[r for r in rows if r['domain_m']==item['domain_m']]
        axes[0].plot([r['station_m'] for r in group],[100*max(w['complex_relative_l2'] for w in r['windows']) for r in group],'-o',label=label)
        axes[1].plot([r['station_m'] for r in group],[r['whole_run_speedup_vs_full'] for r in group],'-o',label=label)
        ids=[next(c for c in contract['cases'] if c['id']==r['case']) for r in group]
        stations=[c['station_m'] for c in ids]
        refmatrix=np.array([results[f"full_x{s:.1f}_y"]['real'] for s in stations]).T
        newmatrix=np.array([results[c['id']]['real'] for c in ids]).T
        plot_sparse_compare(a.out,f"domain_{item['domain_m'][0]:g}x{item['domain_m'][1]:g}_sparse_bscan.png",
                            stations,refmatrix,newmatrix,baseline['t'],f"原生平滑 y｜50×50 m → {label}｜只显示已计算的7站")
    axes[0].axhline(1,color='black',ls='--',label='本轮1%窗口容差')
    axes[0].set(yscale='log',xlabel='站位（m）',ylabel='最差参考窗口复响应差（%）',title='保留弱深部与晚期响应的误差检查')
    axes[1].axhline(1,color='black',ls='--')
    axes[1].set(xlabel='站位（m）',ylabel='同站完整求解加速比',title='官方构建＋求解日志计时；大于1为更快')
    for ax in axes: ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.savefig(a.out/'domain_errors_and_speed.png',dpi=170);plt.close(fig)
    return dict(expanded_controls=controls,summary=summary)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--stage',choices=['A','B','C'],required=True)
    a=p.parse_args();a.root=a.root.resolve();a.baseline=a.baseline.resolve();a.out=a.out.resolve()
    a.out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    contract,results,receipts,checks=audit(a.root,a.stage)
    result=analyze_A(a,contract,results,receipts) if a.stage=='A' else analyze_domains(a,contract,results,receipts)
    save(a.out/'audit.json',dict(stage=a.stage,status='FILES_HASHES_AND_INDEPENDENT_SUMS_PASS',
                               solved_outputs=len(results),checks=checks,contract_sha256=sha(a.root/'contract.json'),
                               script_sha256=sha(Path(__file__)),source_recovery=load(a.root/'audit_recovery.json')))
    if a.stage == 'A':
        brief = dict(surface_only=result['surface_only'], coefficients=result['coefficients'],
                     full_model_spectrum_error_range_pct=[100*min(r['complex_band_relative_l2'] for r in result['full_model']),
                                                          100*max(r['complex_band_relative_l2'] for r in result['full_model'])])
    else:
        brief = dict(summary=result['summary'], expanded_guard=[dict(station_m=r['station_m'],passed=r['passes_frozen_domain_screen']) for r in result['expanded_controls']])
    print(json.dumps(brief,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
