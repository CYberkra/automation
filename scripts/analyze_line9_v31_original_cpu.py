"""Audit the single original V3.1 CPU result and plot existing quantities.

Run with the matching gprMax 4.0.1 toolboxes installed. No solver is launched.
The requested original Gaussian/zero-pad-four product remains unchanged.
"""
from pathlib import Path
import argparse
import hashlib
import json

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(root, out):
    out.mkdir(parents=True, exist_ok=False)
    completion = json.loads((root / "completion.json").read_text(encoding="utf-8"))
    for row in completion["files"]:
        p = root / row["path"]
        assert p.stat().st_size == row["bytes"] and sha(p) == row["sha256"], p
    contract = json.loads((root / "contract.json").read_text(encoding="utf-8"))
    package = root / "package"
    for name, expected in contract["original_file_sha256"].items():
        assert sha(package / name) == expected
    raw = package / "results/line9_pilot_impulse.h5"
    processed = package / "results/line9_pilot_sfcw.h5"
    src = load_source(raw, "/srcs/src1")
    rec = load_receiver(raw, "/rxs/rx1", "Ez")
    with h5py.File(raw, "r") as h:
        native = h["rxs/rx1/Ez"][:]
        dt = float(h.attrs["dt"])
        tx = h["srcs/src1"].attrs["Position"].tolist()
        rx = h["rxs/rx1"].attrs["Position"].tolist()
        assert str(native.dtype) == "float32"
        assert int(h.attrs["Iterations"]) == contract["expected_native_samples"] == len(native)
        assert np.isfinite(native).all()
    with h5py.File(processed, "r") as h:
        freq = h["frequency"][:]
        response = h["response"][:]
        t = h["time_response/time"][:]
        env = h["time_response/complex_envelope"][:]
        signed = h["time_response/real_bandpass"][:]
        weights = h["time_response/weights"][:]
        tail_db = float(h.attrs["ReceiverTailRelativeDB"])
        assert h.attrs["Method"] == "direct"
        assert h["time_response"].attrs["Window"] == "gaussian"
        assert int(h["time_response"].attrs["ZeroPadFactor"]) == 4
        assert h.attrs["TailTaperFraction"] == 0
        assert h["time_response"].attrs["TimeShift"] == 0
        assert np.all(h["source_valid"][:])
    assert np.array_equal(freq, np.linspace(20e6, 170e6, 501))
    assert len(t) == 2004 and np.isfinite(response).all()
    # Independent defining sums with official loader time offsets retained.
    def dft(signal):
        axis = np.arange(len(signal.samples)) * signal.dt + signal.time_offset
        return np.array([np.sum(np.asarray(signal.samples, dtype=np.float64) *
                                np.exp(-2j * np.pi * f * axis)) for f in freq])
    independent = dft(rec) / dft(src)
    relative_spectrum = float(np.linalg.norm(independent - response) / np.linalg.norm(response))
    assert relative_spectrum < 1e-10
    # Independent inverse sum at all published time samples, same weights/scaling.
    inverse = np.array([np.sum(response * weights * np.exp(2j * np.pi *
                       (freq - freq[0]) * delay)) / len(freq) for delay in t])
    relative_inverse = float(np.linalg.norm(inverse - env) / np.linalg.norm(env))
    assert relative_inverse < 1e-10
    expected_time = np.arange(2004) / (2004 * 300000.)
    assert np.array_equal(expected_time, t)
    assert np.allclose(signed, 2 * (env * np.exp(2j * np.pi * freq[0] * t)).real,
                       rtol=1e-12, atol=np.max(np.abs(signed)) * 1e-14)
    contacts = json.loads((package / "model/contacts_m.json").read_text(encoding="utf-8"))
    keys = ["surface", "clay_base", "upper_mudstone_base", "lower_mudstone_top", "lower_mudstone_base"]
    labels = ["地表", "覆盖层底", "上泥岩底", "下泥岩顶", "下泥岩底"]
    elev = []
    for key in keys:
        a = np.asarray(contacts["contact_points_x_elevation_m"][key])
        elev.append(float(np.interp(190., a[:, 0], a[:, 1])))
    c = 299792458.
    times = [2 * 5. / c * 1e9]
    for upper, lower, er in zip(elev[:-1], elev[1:], [12, 18, 9, 18]):
        assert upper >= lower
        times.append(times[-1] + 2 * (upper - lower) * np.sqrt(er) / c * 1e9)
    references = [dict(label=l, elevation_m=z, vertical_two_way_ns=float(delay))
                  for l, z, delay in zip(labels, elev, times)]
    reference_amplitudes = []
    for reference in references:
        amplitude = float(np.interp(reference["vertical_two_way_ns"], t * 1e9, np.abs(env)))
        reference_amplitudes.append({"reference": reference["label"],
                                     "vertical_approx_ns": reference["vertical_two_way_ns"],
                                     "envelope_at_reference": amplitude,
                                     "db_to_global_envelope_peak": float(20 * np.log10(amplitude / np.max(np.abs(env))))})
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
                         "axes.unicode_minus": False, "font.size": 11})
    colors = ["#edf4fa", "#dfbc72", "#8e827d", "#d98077"]
    cmap = ListedColormap(colors)
    with h5py.File(package / "model/line9_pilot_geometry.h5", "r") as h:
        geology = h["data"][:, :, 0]
    fig, ax = plt.subplots(figsize=(12, 7))
    im = ax.imshow(geology.T, origin="lower", extent=[165, 215, 410, 460], aspect="equal",
                   interpolation="nearest", cmap=cmap, norm=BoundaryNorm(np.arange(-.5, 4.5), 4))
    for label, point in [("Tx", tx), ("Rx", rx)]:
        ax.scatter(point[0] + 165, point[1] + 410, marker="v", color="black", s=40)
        ax.annotate(label, (point[0] + 165, point[1] + 410), xytext=(0, 14 if label == "Tx" else -18),
                    textcoords="offset points", ha="center")
    cb = fig.colorbar(im, ax=ax, ticks=range(4), pad=.02)
    cb.ax.set_yticklabels(["空气", "覆盖层 εr12 / σ0.003", "泥岩 εr18 / σ0.006", "砂岩 εr9 / σ0.001"])
    ax.set(xlabel="地质里程（m）", ylabel="高程（m）", title="V3.1独立纵横比例修正版｜190m单站二维模型\n5m离地、2.5cm网格；颜色表示材料，不表示雷达幅度")
    fig.tight_layout(); fig.savefig(out / "geometry_and_station.png", dpi=160); plt.close(fig)
    fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    native_t = (np.arange(len(native)) * dt + rec.time_offset) * 1e9
    axes[0].plot(native_t, native, lw=.6, color="#3c5c7b")
    axes[0].set(title="① 原生CPU单精度总场 A-scan｜未去背景、未加增益", ylabel="Ez（V/m）")
    axes[1].plot(t * 1e9, signed, lw=.8, color="#3c5c7b")
    axes[1].set(title="② 官方direct SFCW带符号时域｜20–170MHz、501点、Gaussian、补零4", ylabel="官方重建幅值\n非端口S21")
    env_db = 20 * np.log10(np.maximum(np.abs(env) / np.max(np.abs(env)), 1e-12))
    axes[2].plot(t * 1e9, env_db, lw=1, color="#3c5c7b")
    axes[2].set(title="③ 同一重建复包络｜以全记录最大值为0dB；仅对数显示", ylabel="相对包络（dB）",
                xlabel="物理时间（ns）；虚线为垂直近似到时，不是回波归因", ylim=(-90, 3), xlim=(0, 800))
    for ax in axes:
        ax.grid(alpha=.18)
        for delay in times:
            ax.axvline(delay, color="#9e6041", ls=":", lw=.8)
    axes[2].legend([f"原生末5%尾峰 {tail_db:.1f} dB；−60dB建议{'通过初筛' if tail_db <= -60 else '未满足'}"], loc="lower right")
    fig.text(.5, .01, "垂直到时（ns）：" + "；".join(f"{l} {d:.1f}" for l, d in zip(labels, times)), ha="center", fontsize=10)
    fig.tight_layout(rect=[0, .035, 1, 1]); fig.savefig(out / "ascan_sfcw_review.png", dpi=160); plt.close(fig)
    fig, ax = plt.subplots(figsize=(3.4, 8))
    keep = (t >= 0) & (t <= 800e-9)
    limit = float(np.max(np.abs(signed)))
    ax.imshow(signed[keep, None], cmap="gray", vmin=-limit, vmax=limit, interpolation="nearest",
              aspect="auto", extent=[0, 1, t[keep][-1] * 1e9, 0])
    ax.set(xticks=[.5], xticklabels=["唯一站位190m"], ylabel="时间（ns）",
           title="单道灰度显示\n不是移动B-scan")
    fig.tight_layout(); fig.savefig(out / "single_station_gray.png", dpi=160); plt.close(fig)
    audit = {"status": "PASS_ORIGINAL_IDENTITY_AND_OFFICIAL_SFCW_NOT_PHYSICAL_CONVERGENCE",
             "native_dtype": str(native.dtype), "native_samples": len(native), "dt_s": dt,
             "tx_local_m": tx, "rx_local_m": rx, "references": references,
             "reference_assumption": "vertical two-way, low-loss speed c/sqrt(er), 5m AGL; ignores bistatic offset/slope and loss correction",
             "source_time_offset_s": src.time_offset, "receiver_time_offset_s": rec.time_offset,
             "independent_DFT_relative_L2": relative_spectrum,
             "independent_inverse_relative_L2": relative_inverse, "tones": len(freq),
             "inverse_samples": len(t), "envelope_period_ns": 1 / 300000. * 1e9,
             "native_last5percent_tail_db": tail_db, "official_minus60db_screen_met": tail_db <= -60,
             "raw_sha256": sha(raw), "processed_sha256": sha(processed),
             "remote_completion_sha256": sha(root / "completion.json"),
             "archive_sha256": contract["archive_sha256"], "retrieved_verified_files": len(completion["files"]),
             "matched_background_available": False, "Bscan_available": False,
             "window_boundary_precision_convergence_certified": False,
             "timing_peaks_are_not_geological_labels": True,
             "reference_time_amplitudes_not_detection": reference_amplitudes,
             "processed_frequency_response_is_not_port_S21": True,
             "independent_sum_tolerance": 1e-10,
             "solver_wall_seconds": next(r["wall_seconds"] for r in completion["receipts"] if r["stage"] == "03_solve_cpu_single")}
    (out / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.root.resolve(), args.out.resolve())
