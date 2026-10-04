"""Apply the official gprMax SFCW API to benchmark3d_r2_co records (39 runs).

House chain: 20-170 MHz x 501 points, tail taper {none, 200 ns}, rectangular +
Hann reconstructions (zero-pad 4). Assembles per-group B-scans of the complex
envelope. Deterministic; run twice (r1/r2) and byte-compare the manifest.

Usage: python sfcw_adapter_benchmark3d.py <out_dir>
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path
import sys

import h5py
import numpy as np
from gprMax.toolboxes.SFCW.processing import (
    direct_frequency_response,
    load_receiver,
    load_source,
    reconstruct_time_response,
)

ROOT = Path(r'E:\automation_djh\automation_repo')
SIMS = ROOT / 'artifacts' / 'simulations'
FREQUENCIES = np.linspace(20e6, 170e6, 501)
TAPER_NS = (None, 200.0)
ZERO_PAD_FACTOR = 4
GROUPS = {
    '3D_5cm': 'B3D5CM-C3mR2-BG-CO13',
    '2D_5cm': 'B2D5CM-C3mR2-BG-CO13',
    '2D_2p5cm': 'B2D-C3mR2-BG-CO13',
}
TRACES = [f't{k:02d}' for k in range(1, 14)]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    records: list[dict] = []
    for gname, prefix in GROUPS.items():
        for t in TRACES:
            run_id = f'{prefix}-{t}'
            path = SIMS / f'2026-10-01_{run_id}' / f'{run_id}.h5'
            src = load_source(path)
            with h5py.File(path, 'r') as h5:
                rx_key = sorted(h5['rxs'].keys())[0]
                rx_name = str(h5['rxs'][rx_key].attrs['Name'])
            rx = load_receiver(path, receiver_path=f'name:{rx_name}', component='Ex')
            assert src.dt == rx.dt
            tag = f'{gname}_{t}'
            for taper_ns in TAPER_NS:
                taper_fraction = 0.0 if taper_ns is None else (round(taper_ns * 1e-9 / rx.dt) - 0.25) / len(rx.samples)
                fr = direct_frequency_response(src, rx, FREQUENCIES, tail_taper_fraction=taper_fraction)
                assert bool(np.all(fr.source_valid)) and bool(np.all(np.isfinite(fr.response)))
                taper_tag = 'no_taper' if taper_ns is None else f'tail_{int(taper_ns)}ns'
                for window in ('rectangular', 'hann'):
                    tr = reconstruct_time_response(fr, window=window, zero_pad_factor=ZERO_PAD_FACTOR, time_shift=0.0)
                    p = f'{tag}_{taper_tag}_{window}'
                    arrays[f'{p}_envelope_time_s'] = tr.time
                    arrays[f'{p}_complex_envelope'] = tr.complex_envelope
                    arrays[f'{p}_real_bandpass'] = tr.real_bandpass
                    if taper_ns is None and window == 'rectangular':
                        arrays[f'{p}_frequency_hz'] = fr.frequency
                        arrays[f'{p}_response'] = fr.response
            records.append(dict(group=gname, trace=t, input_sha256=sha256(path),
                                rx_name=rx_name, receiver_n=int(len(rx.samples)),
                                source_n=int(len(src.samples))))
        print(f'{gname}: 13 traces processed')

    npz_path = out_dir / 'sfcw_responses.npz'
    np.savez_compressed(npz_path, **arrays)
    manifest = dict(schema='sfcw-adapter-benchmark3d/1', status='completed',
                    groups=list(GROUPS), n_traces=len(records),
                    frequencies_hz=[float(FREQUENCIES[0]), float(FREQUENCIES[-1]), int(len(FREQUENCIES))],
                    tapers_ns=[None, 200.0], windows=['rectangular', 'hann'],
                    zero_pad_factor=ZERO_PAD_FACTOR, records=records,
                    npz='sfcw_responses.npz')
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding='utf-8')
    print('wrote', npz_path)


if __name__ == '__main__':
    main()
