"""Apply the installed official gprMax SFCW API to existing C3 CO/MT records only."""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

import h5py
import numpy as np
import scipy
from gprMax.toolboxes.SFCW import processing
from gprMax.toolboxes.SFCW.processing import (
    SampledSignal,
    direct_frequency_response,
    load_receiver,
    load_source,
    reconstruct_time_response,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "artifacts/research_checks/2026-09-27_sfcw_unified_dev_adapter_mt33_rerun"
FREQUENCIES = np.linspace(20e6, 170e6, 501)
TAPER_NS = (None, 200.0)
ZERO_PAD_FACTOR = 4


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inputs() -> list[tuple[str, str, Path]]:
    base = ROOT / "artifacts/research_checks"
    rows = [
        ("MT", "BG", base / "2026-09-26_B2D-C3m-BG-MT33/B2D-C3m-BG-MT33.h5"),
        ("MT", "TGT", base / "2026-09-26_B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-MT33/B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-MT33.h5"),
    ]
    for role, mother in (
        ("BG", "B2D-C3m-BG"),
        ("TGT", "B2D-C3m-D10m-W4m-T0.5m-E20-S0.02"),
    ):
        rows.extend(
            ("CO", role, base / f"2026-09-26_{mother}-CO11-t{k:02d}/{mother}-CO11-t{k:02d}.h5")
            for k in range(1, 12)
        )
    return rows


def synthetic_check() -> dict:
    """Exercise official phase, delay, reconstruction and multitrace linearity."""
    df = 300e3
    frequencies = FREQUENCIES
    dt = 1e-9
    pulse = np.zeros(64)
    pulse[0] = 1.0
    src = SampledSignal("synthetic-source", pulse, dt, 0.5 * dt)

    # A half-step source origin must appear with the correct positive phase in receiver/source.
    same = SampledSignal("synthetic-receiver", pulse, dt, 0.0)
    half = direct_frequency_response(src, same, frequencies)
    expected_half = np.exp(2j * np.pi * frequencies * (0.5 * dt))
    assert np.all(half.source_valid)
    assert np.max(np.abs(half.response - expected_half)) < 1e-12

    # A delay aligned to the zero-padded official reconstruction grid lands at the expected sample.
    delay_steps = 92
    delay = delay_steps / (len(frequencies) * df * ZERO_PAD_FACTOR)
    delayed = SampledSignal("synthetic-delayed", pulse, dt, 0.5 * dt + delay)
    delayed_fr = direct_frequency_response(src, delayed, frequencies)
    recon = reconstruct_time_response(
        delayed_fr, window="rectangular", zero_pad_factor=ZERO_PAD_FACTOR, time_shift=0.0
    )
    assert int(np.argmax(np.abs(recon.complex_envelope))) == delay_steps
    assert abs(recon.time[delay_steps] - delay) < 1e-15

    # Multitrace response preserves the independent trace axis and is linear in receiver data.
    r1 = pulse.copy()
    r1[3] = 0.5
    r2 = np.roll(pulse, 5)
    pair = SampledSignal("synthetic-pair", np.column_stack((r1, r2)), dt, 0.5 * dt)
    summed = SampledSignal("synthetic-sum", r1 + r2, dt, 0.5 * dt)
    pair_fr = direct_frequency_response(src, pair, frequencies)
    sum_fr = direct_frequency_response(src, summed, frequencies)
    assert pair_fr.response.shape == (len(frequencies), 2)
    assert np.max(np.abs(pair_fr.response.sum(axis=1) - sum_fr.response)) < 1e-12
    return {
        "status": "passed",
        "frequency_points": int(len(frequencies)),
        "half_step_source_phase_max_abs_error": float(np.max(np.abs(half.response - expected_half))),
        "delay_s": float(delay),
        "delay_peak_index": int(np.argmax(np.abs(recon.complex_envelope))),
        "delay_expected_index": delay_steps,
        "multitrace_shape": list(pair_fr.response.shape),
        "linearity_max_abs_error": float(np.max(np.abs(pair_fr.response.sum(axis=1) - sum_fr.response))),
    }


def metadata(path: Path, family: str, role: str, src, rx, receiver_attrs, receiver_names, source_attrs) -> dict:
    return {
        "input_path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "input_sha256": sha256(path),
        "family": family,
        "role": role,
        "receiver_path": "name:mt01..mt33" if family == "MT" else "name:measurement",
        "component": "Ex",
        "source_path": "src1",
        "receiver_names": receiver_names,
        "receiver_positions_xyz_m": [np.asarray(a["Position"], dtype=float).tolist() for a in receiver_attrs],
        "source_position_xyz_m": np.asarray(source_attrs["Position"], dtype=float).tolist(),
        "coordinate_note": "2D x-invariant input; x coordinate is an H5 placeholder, not a physical cross-track/flight position",
        "source_type": src.source_type,
        "source_quantity": src.quantity,
        "source_units": src.units,
        "source_spatial_scale": float(src.spatial_scale),
        "source_dt_s": float(src.dt),
        "source_time_offset_s": float(src.time_offset),
        "source_n_samples": int(src.samples.shape[0]),
        "receiver_quantity": rx.quantity,
        "receiver_units": rx.units,
        "receiver_dt_s": float(rx.dt),
        "receiver_time_offset_s": float(rx.time_offset),
        "receiver_n_samples": int(rx.samples.shape[0]),
        "source_time_s_array": "source physical sample axis = source_time_offset_s + n*source_dt_s",
        "receiver_time_s_array": "receiver physical sample axis = receiver_time_offset_s + n*receiver_dt_s",
        "response_units": "receiver sample units / source sample units; no extra spatial-scale normalization",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--self-test-only", action="store_true")
    args = parser.parse_args()

    check = synthetic_check()
    if args.self_test_only:
        print(json.dumps(check, indent=2))
        return
    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite existing output directory: {args.output_dir}")
    paths = inputs()
    missing = [str(p) for _, _, p in paths if not p.is_file()]
    if missing:
        raise SystemExit("missing required existing C3 H5 inputs:\n" + "\n".join(missing))

    records: list[dict] = []
    arrays: dict[str, np.ndarray] = {}
    for family, role, path in paths:
        src = load_source(path)
        with h5py.File(path, "r") as h5:
            rx_keys = sorted(h5["rxs"].keys(), key=lambda key: int(key[2:]))
            receiver_attrs = [dict(h5["rxs"][key].attrs) for key in rx_keys]
            source_attrs = dict(h5["srcs/src1"].attrs)
        if family == "MT":
            assert len(rx_keys) == 33
            receiver_names = [str(a["Name"]) for a in receiver_attrs]
            receivers = [
                load_receiver(path, receiver_path=f"name:mt{k:02d}", component="Ex")
                for k in range(1, 34)
            ]
            assert all(r.dt == receivers[0].dt and r.time_offset == receivers[0].time_offset for r in receivers)
            rx = replace(receivers[0], samples=np.column_stack([r.samples for r in receivers]))
        else:
            receiver_names = [str(receiver_attrs[0]["Name"])]
            rx = load_receiver(path, receiver_path="name:measurement", component="Ex")
        assert src.dt == rx.dt, f"sample interval mismatch: {path}"
        base = metadata(path, family, role, src, rx, receiver_attrs, receiver_names, source_attrs)
        base["source_sample_count"] = int(len(src.samples))
        base["receiver_sample_count"] = int(len(rx.samples))
        key_base = f"{family}_{role}_{path.stem}"
        arrays[f"{key_base}_raw_source_time_s"] = src.times
        arrays[f"{key_base}_raw_receiver_time_s"] = rx.times
        for taper_ns in TAPER_NS:
            taper_fraction = 0.0 if taper_ns is None else (round(taper_ns * 1e-9 / rx.dt) - 0.25) / len(rx.samples)
            fr = direct_frequency_response(src, rx, FREQUENCIES, tail_taper_fraction=taper_fraction)
            assert bool(np.all(fr.source_valid)) and bool(np.all(np.isfinite(fr.response))), str(path)
            tr = reconstruct_time_response(
                fr, window="rectangular", zero_pad_factor=ZERO_PAD_FACTOR, time_shift=0.0
            )
            taper_tag = "no_taper" if taper_ns is None else f"tail_{int(taper_ns)}ns"
            prefix = f"{family}_{role}_{path.stem}_{taper_tag}"
            arrays[f"{prefix}_frequency_hz"] = fr.frequency
            arrays[f"{prefix}_response"] = fr.response
            arrays[f"{prefix}_source_valid"] = fr.source_valid
            arrays[f"{prefix}_source_spectrum"] = fr.source_spectrum
            arrays[f"{prefix}_receiver_spectrum"] = fr.receiver_spectrum
            arrays[f"{prefix}_envelope_time_s"] = tr.time
            arrays[f"{prefix}_complex_envelope"] = tr.complex_envelope
            arrays[f"{prefix}_complex_bandpass"] = tr.complex_bandpass
            arrays[f"{prefix}_real_bandpass"] = tr.real_bandpass
            arrays[f"{prefix}_window_weights"] = tr.weights
            rec = dict(base)
            rec.update(
                taper=taper_tag,
                tail_taper_fraction=float(taper_fraction),
                tail_relative_db_before_taper=float(fr.receiver_tail_relative_db),
                frequency_response_array=f"{prefix}_response",
                official_reconstruction=dict(
                    function="gprMax.toolboxes.SFCW.processing.reconstruct_time_response",
                    window="rectangular",
                    zero_pad_factor=ZERO_PAD_FACTOR,
                    time_shift_s=0.0,
                    time_period_s=float(1.0 / np.diff(FREQUENCIES)[0]),
                    time_sample_interval_s=float(tr.time[1] - tr.time[0]),
                    representation="complex_envelope preferred; real_bandpass is twice real(carrier-modulated envelope)",
                    periodic_boundary="complex_envelope repeats every 1/df; real_bandpass endpoint continuity is not assumed",
                    output_keys={
                        "frequency_hz": f"{prefix}_frequency_hz",
                        "response": f"{prefix}_response",
                        "source_spectrum": f"{prefix}_source_spectrum",
                        "receiver_spectrum": f"{prefix}_receiver_spectrum",
                        "source_valid": f"{prefix}_source_valid",
                        "envelope_time_s": f"{prefix}_envelope_time_s",
                        "complex_envelope": f"{prefix}_complex_envelope",
                        "complex_bandpass": f"{prefix}_complex_bandpass",
                        "real_bandpass": f"{prefix}_real_bandpass",
                        "window_weights": f"{prefix}_window_weights",
                    },
                ),
            )
            records.append(rec)

    args.output_dir.mkdir(parents=True)
    npz_path = args.output_dir / "responses_and_reconstructions.npz"
    np.savez_compressed(npz_path, **arrays)
    py_mod = Path(processing.__file__).resolve()
    script_path = Path(__file__).resolve()
    manifest = {
        "schema": "research-sfcw-unified-dev-adapter/1",
        "scope": "C3 2D mechanism-model development proxy only; not clean field data or physical acceptance",
        "status": "completed",
        "record_decay_certified": False,
        "record_decay_note": "full raw receiver tails are retained; no taper is treated as a convergence test, and archived tail levels do not certify decay",
        "acquisition_contract": {
            "path": "docs/research/2026-09-27_acquisition_and_sfcw_data_contract.md",
            "sha256": sha256(ROOT / "docs/research/2026-09-27_acquisition_and_sfcw_data_contract.md"),
        },
        "inputs_n": len(paths),
        "records_n": len(records),
        "tapers": ["no_taper", "tail_200ns"],
        "frequency_hz": FREQUENCIES.tolist(),
        "frequency_step_hz": float(FREQUENCIES[1] - FREQUENCIES[0]),
        "record_window": "full archived receiver history; no truncation",
        "reconstruction": {
            "window": "rectangular",
            "zero_pad_factor": ZERO_PAD_FACTOR,
            "time_shift_s": 0.0,
            "purpose_of_zero_padding": "interpolates time sampling to 601.2 MHz so real_bandpass can represent the 170 MHz upper carrier; does not increase bandwidth or resolution",
            "complex_envelope_period_s": float(1.0 / np.diff(FREQUENCIES)[0]),
            "complex_envelope_preferred": True,
            "real_bandpass_periodic_continuity_assumed": False,
        },
        "api": {
            "module": "gprMax.toolboxes.SFCW.processing",
            "calls": ["load_source", "load_receiver", "direct_frequency_response", "reconstruct_time_response"],
            "processing_source_path": str(py_mod),
            "processing_source_sha256": sha256(py_mod),
            "adapter_script_path": str(script_path.relative_to(ROOT)).replace("\\", "/"),
            "adapter_script_sha256": sha256(script_path),
        },
        "runtime_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "h5py": h5py.__version__,
            "gprMax": importlib.metadata.version("gprMax"),
        },
        "synthetic_check": check,
        "array_archive": {
            "path": npz_path.name,
            "sha256": sha256(npz_path),
            "keys_n": len(arrays),
        },
        "records": records,
    }
    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir), "records_n": len(records), "npz_sha256": manifest["array_archive"]["sha256"], "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
