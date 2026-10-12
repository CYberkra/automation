"""Read-only comparison of the approved V3.1 800/1100 ns CPU runs.

Uses the matching official gprMax 4.0.1 SFCW loaders. No simulation, taper,
gain, alignment, normalization change or replacement SFCW product is applied.
Output diagnostics use shared plotting scales; time intervals are explicitly
reported and are not geological labels or physical acceptance thresholds.
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
from gprMax.toolboxes.SFCW.processing import load_receiver, load_source


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(a, ref):
    return float(np.linalg.norm(a) / np.linalg.norm(ref))


def read(root, sample_count_erratum=None):
    contract = json.loads((root / "contract.json").read_text(encoding="utf-8"))
    completion = json.loads((root / "completion.json").read_text(encoding="utf-8"))
    assert completion["version"] == "4.0.1"
    for row in completion["files"]:
        p = root / row["path"]
        assert p.stat().st_size == row["bytes"] and sha(p) == row["sha256"], p
    for row in completion["receipts"]:
        assert row["exit_code"] == 0
        assert sha(root / "logs" / (row["stage"] + ".txt")) == row["log_sha256"]
    assert sha(root / "original.zip") == contract["archive_sha256"]
    expected = dict(contract["original_file_sha256"])
    if contract["record_ns"] in (1100, 1600):
        expected["config/pilot_ascan.json"] = contract["modified_pilot_sha256"]
    for name, digest in expected.items():
        assert sha(root / "package" / name) == digest, name
    raw = root / "package/results/line9_pilot_impulse.h5"
    processed = root / "package/results/line9_pilot_sfcw.h5"
    src = load_source(raw, "/srcs/src1")
    rec = load_receiver(raw, "/rxs/rx1", "Ez")
    with h5py.File(raw, "r") as h:
        native = h["rxs/rx1/Ez"][:]
        assert native.dtype == np.dtype("float32") and np.isfinite(native).all()
        dt = float(h.attrs["dt"])
        expected_samples = contract["expected_native_samples"]
        if sample_count_erratum is not None:
            erratum = json.loads(sample_count_erratum.read_bytes())
            assert sha(root / "contract.json") == erratum["frozen_contract_sha256"]
            assert expected_samples == erratum["original_expected_native_samples"]
            # Official TimeWindow(time=...) uses ceil(time/dt) + 1.
            derived = int(np.ceil(contract["record_ns"] * 1e-9 / dt)) + 1
            assert derived == erratum["corrected_expected_native_samples"]
            expected_samples = derived
        assert len(native) == int(h.attrs["Iterations"]) == expected_samples
        position = {name: h[name].attrs["Position"].tolist() for name in ("srcs/src1", "rxs/rx1")}
    with h5py.File(processed, "r") as h:
        assert h.attrs["Method"] == "direct" and h.attrs["TailTaperFraction"] == 0
        g = h["time_response"]
        assert g.attrs["Window"] == "gaussian" and g.attrs["ZeroPadFactor"] == 4
        assert g.attrs["TimeShift"] == 0 and np.all(h["source_valid"][:])
        freq, response = h["frequency"][:], h["response"][:]
        t, envelope, signed, weights = [g[name][:] for name in ("time", "complex_envelope", "real_bandpass", "weights")]
        tail_db = float(h.attrs["ReceiverTailRelativeDB"])
    assert np.array_equal(freq, np.linspace(20e6, 170e6, 501))
    assert len(t) == 2004 and np.isfinite(response).all()
    assert np.array_equal(t, np.arange(2004) / (2004 * 300000.))

    def dft(signal, start=0):
        axis = np.arange(start, len(signal.samples)) * signal.dt + signal.time_offset
        samples = np.asarray(signal.samples[start:], dtype=np.float64)
        return np.array([np.sum(samples * np.exp(-2j * np.pi * f * axis)) for f in freq])

    source_dft, receiver_dft = dft(src), dft(rec)
    defining_response = receiver_dft / source_dft
    dft_error = relative(defining_response - response, response)
    independent_envelope = np.array([np.sum(response * weights * np.exp(2j * np.pi *
                                    (freq - freq[0]) * delay)) / len(freq) for delay in t])
    inverse_error = relative(independent_envelope - envelope, envelope)
    assert dft_error < 1e-10 and inverse_error < 1e-10
    assert np.allclose(signed, 2 * (envelope * np.exp(2j * np.pi * freq[0] * t)).real,
                       rtol=1e-12, atol=np.max(np.abs(signed)) * 1e-14)
    tail = native[-int(np.ceil(len(native) * .05)):]
    computed_tail_db = float(20 * np.log10(np.max(np.abs(tail)) / np.max(np.abs(native))))
    assert abs(computed_tail_db - tail_db) < 1e-4
    return dict(root=root, contract=contract, completion=completion, raw=raw, processed=processed,
                native=native, dt=dt, position=position, src=src, rec=rec, freq=freq,
                response=response, t=t, envelope=envelope, signed=signed, weights=weights,
                tail_db=tail_db, dft_error=dft_error, inverse_error=inverse_error,
                source_dft=source_dft, receiver_dft=receiver_dft, dft=dft)


def compare(baseline, extended, out):
    a, b = read(baseline), read(extended)
    assert a["contract"]["record_ns"] == 800 and b["contract"]["record_ns"] == 1100
    assert sha(baseline / "completion.json") == b["contract"]["baseline_completion_sha256"]
    assert sha(a["raw"]) == b["contract"]["baseline_raw_sha256"]
    ca, cb = [json.loads((root / "package/config/pilot_ascan.json").read_text(encoding="utf-8"))
              for root in (baseline, extended)]
    cb["time_window_ns"] = ca["time_window_ns"]
    assert ca == cb, "Unexpected physical configuration change"
    assert a["dt"] == b["dt"] and a["position"] == b["position"]
    assert a["completion"]["gprMax_init_sha256"] == b["completion"]["gprMax_init_sha256"]
    for name in ("model/contacts_m.json", "model/line9_wang2026_constant_v2.json"):
        assert sha(baseline / "package" / name) == sha(extended / "package" / name)
    assert a["src"].time_offset == b["src"].time_offset
    assert a["rec"].time_offset == b["rec"].time_offset
    for key in ("freq", "t", "weights"):
        assert np.array_equal(a[key], b[key]), key
    geometry_checks = {}
    for filename in ("line9_350m_voxels.h5", "line9_pilot_geometry.h5"):
        pa, pb = [root / "package/model" / filename for root in (baseline, extended)]
        with h5py.File(pa, "r") as ha, h5py.File(pb, "r") as hb:
            assert set(ha.attrs) == set(hb.attrs)
            for key in ha.attrs:
                assert np.array_equal(ha.attrs[key], hb.attrs[key]), key
            assert np.array_equal(ha["data"][:], hb["data"][:])
        geometry_checks[filename] = dict(material_array_and_attributes_equal=True,
                                        bytes_equal=sha(pa) == sha(pb))
    n = len(a["native"])
    prefix_difference = b["native"][:n].astype(float) - a["native"].astype(float)
    source_prefix_equal = np.array_equal(a["src"].samples, b["src"].samples[:len(a["src"].samples)])
    prefix_equal = np.array_equal(a["native"], b["native"][:n])
    assert prefix_equal and source_prefix_equal, "Uncontrolled native prefix difference"
    assert not np.any(b["src"].samples[len(a["src"].samples):])
    delta = b["response"] - a["response"]
    # Prefix is identical and the added source samples are zero. Thus the
    # additional receiver samples alone must account for the spectrum change.
    added_dft = b["dft"](b["rec"], n)
    predicted_delta = added_dft / b["source_dft"]
    delta_closure = relative(predicted_delta - delta, delta)
    assert delta_closure < 1e-7
    magnitude_db = 20 * np.log10(np.abs(b["response"]) / np.abs(a["response"]))
    phase_deg = np.angle(b["response"] * np.conj(a["response"]), deg=True)
    common_peak = float(np.max(np.abs(b["envelope"])))
    envelope_delta = b["envelope"] - a["envelope"]
    intervals = []
    for lo, hi, label in [(0, 100, "直达与地表参考区"), (150, 400, "地下界面参考区"),
                          (400, 800, "旧窗晚期区"), (800, 1100, "延长记录区")]:
        keep = (b["t"] * 1e9 >= lo) & (b["t"] * 1e9 < hi)
        intervals.append(dict(label=label, start_ns=lo, stop_ns=hi,
                              complex_envelope_relative_L2=relative(envelope_delta[keep], b["envelope"][keep]),
                              peak_difference_db_to_common_global_peak=float(20 * np.log10(
                                  np.max(np.abs(envelope_delta[keep])) / common_peak))))
    endpoint_checks = []
    for center in (800, 1100):
        keep = np.abs(b["t"] * 1e9 - center) <= 15
        endpoint_checks.append(dict(center_ns=center, half_width_ns=15,
            peak_db_to_common_global_peak=[float(20 * np.log10(np.max(np.abs(x["envelope"][keep])) / common_peak)) for x in (a, b)]))
    references = json.loads(Path("artifacts/research_checks/2026-10-12_line9_v31_original_cpu_r1/audit.json").read_text(encoding="utf-8"))["references"]
    out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
                         "axes.unicode_minus": False, "font.size": 10})
    colors = ["#54799c", "#c15f34"]
    labels = ["800 ns（复用）", "1100 ns（本轮）"]

    def db(values, peak):
        return 20 * np.log10(np.maximum(np.abs(values) / peak, 1e-12))

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    native_peak = float(np.max(np.abs(b["native"])))
    for result, color, label in zip((a, b), colors, labels):
        axis = np.arange(len(result["native"])) * result["dt"] * 1e9
        axes[0].plot(axis, result["native"], color=color, lw=.55, label=label)
        axes[1].plot(axis, db(result["native"], native_peak), color=color, lw=.55,
                     label=f"{label}末5%尾峰：{result['tail_db']:.2f} dB")
    axes[0].set(title="原生总场｜CPU float32；前800 ns逐比特相同", ylabel="Ez（V/m）")
    axes[1].set(title="原生幅度的共同对数显示｜未滤波、未加尾窗", ylabel="相对原生峰（dB）",
                xlabel="时间（ns）", ylim=(-100, 3), xlim=(0, 1100))
    axes[1].axhline(-60, color="gray", ls="--", label="官方尾部初筛建议 −60 dB")
    for ax in axes:
        ax.axvline(800, color="black", ls=":", lw=1)
        ax.legend(loc="upper right"); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(out / "native_window_comparison.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    spectrum_peak = np.max(np.abs(b["response"]))
    for result, color, label in zip((a, b), colors, labels):
        axes[0].plot(result["freq"] / 1e6, db(result["response"], spectrum_peak), color=color, label=label)
    axes[0].set(title="官方direct复数响应幅度｜共同参考；不是端口S21", ylabel="相对幅度（dB）")
    axes[0].legend()
    axes[1].plot(b["freq"] / 1e6, magnitude_db, color=colors[1])
    axes[1].set(title="逐频幅度差＝20 log10（|H1100| / |H800|）", ylabel="幅度差（dB）")
    axes[2].plot(b["freq"] / 1e6, phase_deg, color=colors[1])
    axes[2].set(title="逐频包裹相位差＝arg（H1100 × conj H800）｜弱谱点相位需谨慎", ylabel="相位差（°）", xlabel="频率（MHz）")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(out / "complex_spectrum_comparison.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    for result, color, label in zip((a, b), colors, labels):
        axes[0].plot(result["t"] * 1e9, result["signed"], lw=.8, color=color, label=label)
        axes[1].plot(result["t"] * 1e9, db(result["envelope"], common_peak), color=color, lw=1, label=label)
    axes[0].set(title="190 m单道｜官方direct / Gaussian / 补零4；共同幅度标尺", ylabel="带符号重建幅值")
    axes[0].legend()
    axes[1].set(title="复包络｜相同全记录参考峰；未逐道归一化", ylabel="相对包络（dB）", ylim=(-100, 3))
    axes[2].plot(b["t"] * 1e9, db(envelope_delta, common_peak), color="#65518d", label="复包络差（1100−800）")
    axes[2].set(title="时窗改变引起的差异｜差分仅作数值诊断，不是地下净回波", ylabel="差异 / 共同峰（dB）",
                xlabel="时间（ns）；虚线为几何垂直近似到时，不是已检出界面", ylim=(-120, -20), xlim=(0, 1100))
    for ax in axes:
        ax.axvspan(150, 400, color="#ddbd72", alpha=.12)
        ax.axvline(800, color="black", ls=":")
        for ref in references:
            ax.axvline(ref["vertical_two_way_ns"], color="#888888", ls=":", lw=.7)
        ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(out / "sfcw_window_comparison.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for result, color, label in zip((a, b), colors, labels):
        axes[0].plot(result["t"] * 1e9, db(result["envelope"], common_peak), color=color, label=label)
    axes[0].set(title="地下参考时段细看｜两条复包络基本重合；共同全记录参考峰", ylabel="包络（dB）", ylim=(-100, -40))
    axes[0].legend()
    axes[1].plot(b["t"] * 1e9, 20 * np.log10(np.maximum(np.abs(envelope_delta) / common_peak, 1e-18)), color="#65518d")
    axes[1].set(title=f"复包络差值｜150–400 ns相对L2差{100 * intervals[1]['complex_envelope_relative_L2']:.4f}%；差值不代表物理可信底限", ylabel="差值 / 共同峰（dB）", xlabel="时间（ns）", ylim=(-150, -120), xlim=(150, 400))
    for ax in axes:
        for ref in references[1:]:
            ax.axvline(ref["vertical_two_way_ns"], color="#888888", ls=":", lw=.8)
        ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(out / "underground_reference_zoom.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(5, 8), sharey=True)
    keep = b["t"] <= 1100e-9
    limit = np.max(np.abs(b["signed"]))
    for result, ax, label in zip((a, b), axes, labels):
        ax.imshow(result["signed"][keep, None], cmap="gray", vmin=-limit, vmax=limit,
                  interpolation="nearest", aspect="auto", extent=[0, 1, b["t"][keep][-1] * 1e9, 0])
        ax.axhline(800, color="#bd704a", ls=":")
        ax.set(xticks=[.5], xticklabels=["190 m唯一站位"], title=label)
    axes[0].set_ylabel("时间（ns）")
    fig.suptitle("单道共同灰度显示｜不是B-scan\n800 ns曲线之后为有限频点周期重建，不是新增实测响应", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, .94]); fig.savefig(out / "single_station_gray_pair.png", dpi=150); plt.close(fig)

    audit = dict(status="CONTROLLED_PAIR_VERIFIED_NOT_ASYMPTOTIC_OR_PML_CONVERGENCE",
                 authorization=b["contract"]["authorization"], new_solver_runs=1, dt_s=b["dt"],
                 native_samples=[len(a["native"]), len(b["native"])], native_dtype="float32",
                 receiver_prefix_bitwise_equal=prefix_equal, source_prefix_bitwise_equal=source_prefix_equal,
                 prefix_max_absolute_difference=float(np.max(np.abs(prefix_difference))), geometry=geometry_checks,
                 only_physical_change="800ns to1100ns FDTD record length", native_tail_last5percent_db=[a["tail_db"], b["tail_db"]],
                 official_minus60db_screen_met=[a["tail_db"] <= -60, b["tail_db"] <= -60],
                 source_spectrum_relative_difference=relative(b["source_dft"] - a["source_dft"], b["source_dft"]),
                 complex_spectrum_relative_L2=relative(delta, b["response"]),
                 frequency_intervals=[dict(start_MHz=lo, stop_MHz=hi,
                    complex_response_relative_L2=relative(delta[keep], b["response"][keep]))
                    for lo, hi in [(20, 40), (40, 95), (95, 170)]
                    for keep in [(b["freq"] >= lo * 1e6) & (b["freq"] <= hi * 1e6)]],
                 extra_recorded_tail_explains_spectrum_delta_relative_residual=delta_closure,
                 max_abs_magnitude_change_db=float(np.max(np.abs(magnitude_db))),
                 max_abs_wrapped_phase_change_deg=float(np.max(np.abs(phase_deg))),
                 max_abs_magnitude_change_frequency_MHz=float(b["freq"][np.argmax(np.abs(magnitude_db))] / 1e6),
                 max_abs_phase_change_frequency_MHz=float(b["freq"][np.argmax(np.abs(phase_deg))] / 1e6),
                 phase_change_abs_percentile_deg=np.percentile(np.abs(phase_deg), [50, 95, 100]).tolist(),
                 intervals=intervals, endpoint_checks=endpoint_checks, references=references,
                 reference_limitation="vertical low-loss approximate timing; not isolated interface echoes",
                 independent_DFT_relative_L2=[a["dft_error"], b["dft_error"]],
                 independent_inverse_relative_L2=[a["inverse_error"], b["inverse_error"]],
                 solver_wall_seconds=[next(r["wall_seconds"] for r in x["completion"]["receipts"] if r["stage"] == "03_solve_cpu_single") for x in (a, b)],
                 raw_sha256=[sha(a["raw"]), sha(b["raw"])], processed_sha256=[sha(a["processed"]), sha(b["processed"])],
                 completion_sha256=[sha(baseline / "completion.json"), sha(extended / "completion.json")],
                 analysis_script_sha256=sha(Path(__file__)),
                 executed_supervisor_sha256=sha(extended / "run_line9_v31_original_cpu.py"),
                 verified_package_file_counts=[len(a["completion"]["files"]), len(b["completion"]["files"])],
                 inverse_period_ns=1e9 / 300000., Bscan_available=False, matched_background_available=False,
                 asymptotic_window_convergence_certified=False, PML_convergence_certified=False,
                 precision_convergence_certified=False, physical_detection_certified=False)
    (out / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    np.savez(out / "comparison_vectors.npz", frequency_hz=b["freq"], response800=a["response"], response1100=b["response"],
             time_s=b["t"], envelope800=a["envelope"], envelope1100=b["envelope"], magnitude_change_db=magnitude_db, phase_change_deg=phase_deg)
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--extended", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    compare(args.baseline.resolve(), args.extended.resolve(), args.out.resolve())
