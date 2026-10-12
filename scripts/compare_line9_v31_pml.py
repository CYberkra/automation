"""Audit the approved 1100 ns V3.1 PML80/PML120 CPU-single pair.

Reuses the independently checked official SFCW loader/sum audit. No solver,
filter, taper, gain, realignment or substitute SFCW product is applied.
"""
from pathlib import Path
import argparse
import json

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_line9_v31_windows import read, sha, relative


def compare(baseline, candidate, out):
    a, b = read(baseline), read(candidate)
    assert b["contract"]["task_id"] == "v31_pml80_vs120_at1100ns"
    assert a["contract"]["record_ns"] == b["contract"]["record_ns"] == 1100
    assert sha(baseline / "completion.json") == b["contract"]["baseline_completion_sha256"]
    assert sha(a["raw"]) == b["contract"]["baseline_raw_sha256"]
    ca, cb = [json.loads((root / "package/config/pilot_ascan.json").read_bytes()) for root in (baseline, candidate)]
    assert ca["pml"]["cells_by_face"] == [80, 80, 0, 80, 80, 0]
    assert cb["pml"]["cells_by_face"] == [120, 120, 0, 120, 120, 0]
    cb["pml"]["cells_by_face"] = ca["pml"]["cells_by_face"]
    assert ca == cb, "Unexpected configuration change beyond PML thickness"
    assert a["dt"] == b["dt"] and a["position"] == b["position"]
    assert a["completion"]["gprMax_init_sha256"] == b["completion"]["gprMax_init_sha256"]
    for key in ("freq", "t", "weights"):
        assert np.array_equal(a[key], b[key]), key
    assert np.array_equal(a["src"].samples, b["src"].samples)
    assert a["src"].time_offset == b["src"].time_offset
    assert a["rec"].time_offset == b["rec"].time_offset
    assert a["native"].shape == b["native"].shape
    assert np.array_equal(a["source_dft"], b["source_dft"])
    geometry = {}
    for name in ("line9_350m_voxels.h5", "line9_pilot_geometry.h5"):
        paths = [root / "package/model" / name for root in (baseline, candidate)]
        with h5py.File(paths[0]) as ha, h5py.File(paths[1]) as hb:
            assert set(ha.attrs) == set(hb.attrs)
            for key in ha.attrs:
                assert np.array_equal(ha.attrs[key], hb.attrs[key]), key
            for key in ("data", "material_keys"):
                assert np.array_equal(ha[key][:], hb[key][:]), key
        geometry[name] = dict(arrays_attributes_equal=True, bytes_equal=sha(paths[0]) == sha(paths[1]))
    for name in ("contacts_m.json", "line9_wang2026_constant_v2.json"):
        assert sha(baseline / "package/model" / name) == sha(candidate / "package/model" / name)

    delta = b["response"] - a["response"]
    native_delta = b["native"].astype(float) - a["native"].astype(float)
    envelope_delta = b["envelope"] - a["envelope"]
    source = b["source_dft"]
    # Equal source history, sampling and offsets: the transform of the native
    # receiver difference divided by that source is the frequency-response change.
    native_time_s = np.arange(len(native_delta)) * b["dt"] + b["rec"].time_offset
    difference_dft = np.array([np.sum(native_delta * np.exp(-2j * np.pi * f * native_time_s)) for f in b["freq"]])
    closure = relative(difference_dft / source - delta, delta)
    assert closure < 1e-7
    phase = np.angle(b["response"] * np.conj(a["response"]), deg=True)
    magnitude = 20 * np.log10(np.abs(b["response"]) / np.abs(a["response"]))
    t_ns = b["t"] * 1e9
    native_t_ns = np.arange(len(b["native"])) * b["dt"] * 1e9
    peak = float(np.max(np.abs(b["envelope"])))
    raw_peak = float(np.max(np.abs(b["native"])))
    intervals = []
    for lo, hi, label in [(0, 100, "直达与地表参考区"), (150, 400, "地下界面参考区"),
                          (400, 800, "晚期区"), (800, 1100, "记录末段")]:
        keep = (t_ns >= lo) & (t_ns < hi)
        raw_keep = (native_t_ns >= lo) & (native_t_ns < hi)
        intervals.append(dict(label=label, start_ns=lo, stop_ns=hi,
            complex_envelope_relative_L2=relative(envelope_delta[keep], b["envelope"][keep]),
            native_trace_relative_L2=relative(native_delta[raw_keep], b["native"][raw_keep]),
            peak_envelope_difference_db_to_common_peak=float(20 * np.log10(np.max(np.abs(envelope_delta[keep])) / peak))))
    nz = np.flatnonzero(native_delta)
    references = json.loads(Path("artifacts/research_checks/2026-10-12_line9_v31_original_cpu_r1/audit.json").read_bytes())["references"]
    reference_windows = []
    for ref in references[1:]:
        center = ref["vertical_two_way_ns"]
        keep = np.abs(t_ns - center) <= 15
        at = [np.interp(center, t_ns, x.real) + 1j * np.interp(center, t_ns, x.imag)
              for x in (a["envelope"], b["envelope"])]
        reference_windows.append(dict(label=ref["label"], center_ns=center, half_width_ns=15,
            complex_envelope_relative_L2=relative(envelope_delta[keep], b["envelope"][keep]),
            complex_change_at_reference_relative=float(abs(at[1] - at[0]) / abs(at[1])),
            window_peak_db_to_common_peak=[float(20 * np.log10(np.max(np.abs(x["envelope"][keep])) / peak)) for x in (a, b)]))
    # Diagnostic FFT of an unchanged subset of the existing raw records.
    # This is not the source-normalised official SFCW response. Its finite
    # rectangular support introduces spectral leakage and limits resolution.
    tail_keep = native_t_ns >= 800
    tail_spectra = [np.fft.rfft(x["native"][tail_keep].astype(float)) for x in (a, b)]
    tail_frequency = np.fft.rfftfreq(int(np.sum(tail_keep)), b["dt"])
    tail_bin_weights = np.full(len(tail_frequency), 2.)
    tail_bin_weights[0] = 1
    if np.sum(tail_keep) % 2 == 0:
        tail_bin_weights[-1] = 1
    tail_in_band = (tail_frequency >= 20e6) & (tail_frequency <= 170e6)
    tail_diagnostics = []
    for spectrum in tail_spectra:
        power = abs(spectrum) ** 2 * tail_bin_weights
        dominant = float(tail_frequency[np.argmax(abs(spectrum))])
        tail_diagnostics.append(dict(FFT_peak_Hz=dominant, FFT_bin_spacing_Hz=float(tail_frequency[1]),
            finite_record_spectral_energy_fraction_20_170MHz=float(np.sum(power[tail_in_band]) / np.sum(power)),
            air_cells_per_wavelength_at_peak=299792458. / dominant / .025))
    out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
                         "axes.unicode_minus": False, "font.size": 10})
    colors, labels = ["#54799c", "#c15f34"], ["PML 80格（复用）", "PML 120格（本轮）"]

    def db(values, reference):
        return 20 * np.log10(np.maximum(np.abs(values) / reference, 1e-15))

    def reference_lines(ax):
        for ref in references:
            ax.axvline(ref["vertical_two_way_ns"], color="gray", ls=":", lw=.7)
        ax.grid(alpha=.2)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for result, color, label in zip((a, b), colors, labels):
        axes[0].plot(native_t_ns, db(result["native"], raw_peak), color=color, lw=.55,
                     label=f"{label}；末5%峰{result['tail_db']:.2f} dB")
    axes[0].axhline(-60, color="gray", ls="--", label="官方尾部初筛建议 −60 dB")
    axes[0].set(title="1100 ns原生CPU float32总场｜未滤波、未加尾窗；共同标尺", ylabel="原生幅度 / 共同峰（dB）", ylim=(-100, 3))
    axes[0].legend()
    axes[1].plot(native_t_ns, db(native_delta, raw_peak), color="#65518d", lw=.6)
    axes[1].set(title="原生接收差值（120格−80格）｜不是地下净回波", ylabel="差值 / 共同峰（dB）", xlabel="时间（ns）", ylim=(-140, 0), xlim=(0, 1100))
    for ax in axes:
        reference_lines(ax)
    fig.tight_layout(); fig.savefig(out / "native_pml_comparison.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    frequency = b["freq"] / 1e6
    frequency_peak = np.max(np.abs(b["response"]))
    for result, color, label in zip((a, b), colors, labels):
        axes[0].plot(frequency, db(result["response"], frequency_peak), color=color, label=label)
    axes[0].set(title="官方direct复数响应｜20–170 MHz、501点；不是端口S21", ylabel="相对幅度（dB）")
    axes[0].legend()
    axes[1].plot(frequency, magnitude, color=colors[1])
    axes[1].set(title="逐频幅度差＝20 log10（|H120| / |H80|）", ylabel="幅度差（dB）")
    axes[2].plot(frequency, phase, color=colors[1])
    axes[2].set(title="包裹相位差＝arg（H120 × conj H80）｜弱谱点相位需谨慎", ylabel="相位差（°）", xlabel="频率（MHz）")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(out / "complex_spectrum_pml_comparison.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(12, 11))
    for result, color, label in zip((a, b), colors, labels):
        axes[0].plot(t_ns, db(result["envelope"], peak), color=color, label=label)
        axes[1].plot(t_ns, db(result["envelope"], peak), color=color, label=label)
    axes[0].set(title="190 m同站总场｜1100 ns、官方direct / Gaussian / 补零4", ylabel="复包络 / 共同峰（dB）", ylim=(-100, 3), xlim=(0, 1100))
    axes[0].legend()
    axes[1].set(title="地下参考时段细看｜共同全记录参考峰；虚线不是检出标签", ylabel="复包络 / 共同峰（dB）", ylim=(-100, -40), xlim=(150, 400))
    axes[2].plot(t_ns, db(envelope_delta, peak), color="#65518d")
    axes[2].set(title=f"复包络差值｜150–400 ns相对L2差{100 * intervals[1]['complex_envelope_relative_L2']:.4f}%", ylabel="差值 / 共同峰（dB）", xlabel="时间（ns）", ylim=(-130, 0), xlim=(0, 1100))
    for ax in axes:
        reference_lines(ax)
    fig.tight_layout(); fig.savefig(out / "sfcw_pml_comparison.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(5, 8), sharey=True)
    keep = t_ns <= 1100
    limit = float(np.max(np.abs(b["signed"])))
    for result, ax, label in zip((a, b), axes, labels):
        ax.imshow(result["signed"][keep, None], cmap="gray", vmin=-limit, vmax=limit,
                  interpolation="nearest", aspect="auto", extent=[0, 1, t_ns[keep][-1], 0])
        ax.set(xticks=[.5], xticklabels=["190 m唯一站位"], title=label)
    axes[0].set_ylabel("时间（ns）")
    fig.suptitle("1100 ns单道共同灰度显示｜不是移动B-scan", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, .96]); fig.savefig(out / "single_station_gray_pair.png", dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 5))
    tail_peak = np.max(abs(tail_spectra[0]))
    for spectrum, color, label in zip(tail_spectra, colors, labels):
        ax.plot(tail_frequency / 1e9, db(spectrum, tail_peak), color=color, lw=.8, label=label)
    ax.axvspan(.020, .170, color="#ddbd72", alpha=.4, label="正式SFCW工作带20–170 MHz")
    ax.set(title="已有原生800–1100 ns尾段的FFT｜只读诊断；未改变正式SFCW\n有限矩形记录有谱泄漏；曲线不是传递响应或边界反射率", xlabel="频率（GHz）", ylabel="FFT幅度 / 80格尾谱峰（dB）", xlim=(0, tail_frequency[-1] / 1e9), ylim=(-100, 3))
    ax.legend(); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(out / "tail_frequency_diagnostic.png", dpi=150); plt.close(fig)

    audit = dict(status="CONTROLLED_PML_CONFIGURATION_SENSITIVITY_VERIFIED_NOT_FULL_CONVERGENCE",
        authorization=b["contract"]["authorization"], authorization_context=b["contract"]["authorization_context"],
        new_solver_runs=1, record_ns=1100, native_dtype="float32", native_samples=len(b["native"]), dt_s=b["dt"],
        geometry=geometry, full_source_history_bitwise_equal=True, unchanged_physical_configuration_except_PML=True,
        PML_cells_by_face=[[80, 80, 0, 80, 80, 0], [120, 120, 0, 120, 120, 0]],
        domain_m=[50, 50], PML_entrance_shift_m=1,
        native_tail_db=[a["tail_db"], b["tail_db"]], official_minus60db_screen_met=[a["tail_db"] <= -60, b["tail_db"] <= -60],
        first_observed_bitwise_native_difference_ns=float(native_t_ns[nz[0]]) if len(nz) else None,
        first_bitwise_difference_is_not_a_physical_arrival=True,
        native_trace_relative_L2=relative(native_delta, b["native"]),
        complex_spectrum_relative_L2=relative(delta, b["response"]),
        frequency_intervals=[dict(start_MHz=lo, stop_MHz=hi,
            complex_response_relative_L2=relative(delta[keep], b["response"][keep]))
            for lo, hi in [(20, 40), (40, 95), (95, 170)]
            for keep in [(b["freq"] >= lo * 1e6) & (b["freq"] <= hi * 1e6)]],
        max_abs_magnitude_change_db=float(np.max(np.abs(magnitude))),
        max_abs_wrapped_phase_change_deg=float(np.max(np.abs(phase))),
        phase_change_abs_percentile_deg=np.percentile(np.abs(phase), [50, 95, 100]).tolist(),
        native_difference_spectrum_closure_relative_L2=closure,
        intervals=intervals, references=references, reference_windows=reference_windows,
        reference_windows_are_not_isolated_interface_echoes=True,
        tail_FFT_diagnostic=tail_diagnostics, tail_FFT_range_ns=[800, 1100],
        tail_FFT_limitation="finite rectangular segment; leakage and3.3MHz bins; not official source-normalised SFCW, reflection coefficient or physical validation",
        independent_DFT_relative_L2=[a["dft_error"], b["dft_error"]],
        independent_inverse_relative_L2=[a["inverse_error"], b["inverse_error"]],
        solver_wall_seconds=[next(r["wall_seconds"] for r in x["completion"]["receipts"] if r["stage"] == "03_solve_cpu_single") for x in (a, b)],
        raw_sha256=[sha(a["raw"]), sha(b["raw"])], processed_sha256=[sha(a["processed"]), sha(b["processed"])],
        completion_sha256=[sha(baseline / "completion.json"), sha(candidate / "completion.json")],
        analysis_script_sha256=sha(Path(__file__)), shared_audit_script_sha256=sha(Path(__file__).with_name("compare_line9_v31_windows.py")),
        executed_supervisor_sha256=sha(candidate / "run_line9_v31_original_cpu.py"),
        verified_package_file_counts=[len(a["completion"]["files"]), len(b["completion"]["files"])],
        inverse_period_ns=1e9 / 300000., Bscan_available=False, matched_background_available=False,
        asymptotic_PML_convergence_certified=False, asymptotic_window_convergence_certified=False,
        individual_boundary_reflection_coefficient_identified=False, physical_detection_certified=False)
    (out / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    np.savez(out / "comparison_vectors.npz", frequency_hz=b["freq"], response80=a["response"], response120=b["response"],
             time_s=b["t"], envelope80=a["envelope"], envelope120=b["envelope"], magnitude_change_db=magnitude, phase_change_deg=phase,
             tail_fft_frequency_hz=tail_frequency, tail_fft80=tail_spectra[0], tail_fft120=tail_spectra[1])
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    compare(args.baseline.resolve(), args.candidate.resolve(), args.out.resolve())
