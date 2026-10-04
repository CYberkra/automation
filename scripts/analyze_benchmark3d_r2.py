"""Assemble benchmark3d_r2_co B-scans and compute draft-§6 acceptance metrics.

Deterministic: sorted keys, fixed formatting, no timestamps. Run twice into
separate output dirs and byte-compare the JSON (house r1/r2 convention).

Usage: python analyze_benchmark3d_r2.py <out_dir>
"""
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(r'E:\automation_djh\automation_repo')
SIMS = ROOT / 'artifacts' / 'simulations'
CONFIGS = ROOT / 'configs' / 'research' / 'benchmark3d_r2_co'

GROUPS = {
    '3D-5cm': 'B3D5CM-C3mR2-BG-CO13',
    '2D-5cm': 'B2D5CM-C3mR2-BG-CO13',
    '2D-2.5cm': 'B2D-C3mR2-BG-CO13',
}
TRACES = [f't{i:02d}' for i in range(1, 14)]


def parse_in(path):
    """Extract rx y-position and run tag from a frozen .in file."""
    txt = path.read_text(encoding='utf-8', errors='replace')
    m = re.search(r'^#rx:\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s*$', txt, re.M)
    if not m:
        raise ValueError(f'no #rx line in {path}')
    y = float(m.group(1))
    tag = m.group(2)
    field = m.group(3)
    ms = re.search(r'^#hertzian_dipole:\s+\S+\s+\S+\s+(\S+)\s+\S+', txt, re.M)
    src_y = float(ms.group(1)) if ms else None
    return dict(rx_y=y, src_y=src_y, tag=tag, field=field)


def load_trace(run_dir, run_id):
    h5_path = run_dir / f'{run_id}.h5'
    with h5py.File(h5_path, 'r') as f:
        arr = np.asarray(f[f'rxs/rx1/{next(iter(f["rxs/rx1"].keys()))}']).ravel()
        dt = float(f.attrs['dt'])
    return arr, dt


def window_idx(t_ns, lo, hi):
    return (t_ns >= lo) & (t_ns < hi)


def main():
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    report = {'groups': {}, 'traces_missing': []}

    for gname, prefix in GROUPS.items():
        traces = {}
        for t in TRACES:
            run_id = f'{prefix}-{t}'
            run_dir = SIMS / f'2026-10-01_{run_id}'
            if not (run_dir / f'{run_id}.h5').exists():
                report['traces_missing'].append(run_id)
                continue
            cfg = parse_in(CONFIGS / f'{run_id}.in')
            arr, dt = load_trace(run_dir, run_id)
            traces[t] = dict(cfg=cfg, arr=arr, dt=dt)
        if not traces:
            continue
        dts = {v['dt'] for v in traces.values()}
        assert len(dts) == 1, f'{gname}: dt mismatch {dts}'
        nsteps = {len(v['arr']) for v in traces.values()}
        assert len(nsteps) == 1, f'{gname}: nsteps mismatch {nsteps}'
        dt = dts.pop()
        n = nsteps.pop()
        t_ns = np.arange(n) * dt * 1e9

        order = sorted(traces, key=lambda t: traces[t]['cfg']['rx_y'])
        B = np.stack([traces[t]['arr'] for t in order])  # 13 x nsteps

        per_trace = {}
        for t in order:
            a = traces[t]['arr']
            fin = bool(np.isfinite(a).all())
            direct_w = window_idx(t_ns, 0, 50)
            surf_w = window_idx(t_ns, 60, 140)
            iface_w = window_idx(t_ns, 140, 350)
            tail_w = t_ns >= t_ns[-1] * 0.9
            e_total = float(np.sum(a ** 2))
            i_idx = np.argmax(np.abs(a[iface_w]))
            iface_t = float(t_ns[iface_w][i_idx])
            per_trace[t] = dict(
                rx_y=traces[t]['cfg']['rx_y'],
                finite=fin,
                peak=float(np.max(np.abs(a))),
                e_direct=float(np.sum(a[direct_w] ** 2)),
                e_surface=float(np.sum(a[surf_w] ** 2)),
                e_interface=float(np.sum(a[iface_w] ** 2)),
                e_total=e_total,
                tail_ratio=float(np.sum(a[tail_w] ** 2) / e_total),
                iface_peak_t_ns=iface_t,
                iface_peak=float(np.abs(a[iface_w][i_idx])),
            )

        # cross-trace interface-window correlation
        iw = B[:, iface_w]
        iw = iw - iw.mean(axis=1, keepdims=True)
        norms = np.linalg.norm(iw, axis=1)
        valid = norms > 0
        corrs = []
        if valid.sum() > 1:
            nw = iw[valid] / norms[valid][:, None]
            C = nw @ nw.T
            iu = np.triu_indices(C.shape[0], 1)
            corrs = C[iu].tolist()
        iface_ts = [per_trace[t]['iface_peak_t_ns'] for t in order]
        g = dict(
            n_traces=len(order),
            dt_ps=round(dt * 1e12, 4),
            nsteps=n,
            window_ns=round(float(t_ns[-1]), 4),
            rx_y=[per_trace[t]['rx_y'] for t in order],
            per_trace={t: per_trace[t] for t in order},
            iface_corr_mean=float(np.mean(corrs)) if corrs else None,
            iface_corr_min=float(np.min(corrs)) if corrs else None,
            iface_t_mean_ns=float(np.mean(iface_ts)),
            iface_t_std_ns=float(np.std(iface_ts)),
            direct_e_mean=float(np.mean([per_trace[t]['e_direct'] for t in order])),
            surface_e_mean=float(np.mean([per_trace[t]['e_surface'] for t in order])),
            interface_e_mean=float(np.mean([per_trace[t]['e_interface'] for t in order])),
            tail_ratio_max=float(np.max([per_trace[t]['tail_ratio'] for t in order])),
        )
        report['groups'][gname] = g
        tag = gname.replace('-', '_').replace('.', 'p')
        np.save(out_dir / f'bscan_{tag}.npy', B)
        np.save(out_dir / f'rx_y_{tag}.npy', np.array([per_trace[t]['rx_y'] for t in order]))
        np.save(out_dir / f't_ns_{tag}.npy', t_ns)

    # ---- cross-group paired metrics (3D vs 2D-5cm vs 2D-2.5cm) ----
    paired = {}
    gs = report['groups']
    if '3D-5cm' in gs and '2D-5cm' in gs:
        d3 = np.array([gs['3D-5cm']['per_trace'][t]['e_direct'] for t in TRACES if t in gs['3D-5cm']['per_trace']])
        d2 = np.array([gs['2D-5cm']['per_trace'][t]['e_direct'] for t in TRACES if t in gs['2D-5cm']['per_trace']])
        if len(d3) == len(d2) and len(d3):
            paired['direct_energy_ratio_3d_over_2d5cm'] = (d3 / d2).tolist()
    if '2D-5cm' in gs and '2D-2.5cm' in gs:
        common = [t for t in TRACES if t in gs['2D-5cm']['per_trace'] and t in gs['2D-2.5cm']['per_trace']]
        if common:
            i5 = np.array([gs['2D-5cm']['per_trace'][t]['iface_peak_t_ns'] for t in common])
            i25 = np.array([gs['2D-2.5cm']['per_trace'][t]['iface_peak_t_ns'] for t in common])
            paired['grid_refinement_interface_dt_ns'] = (i25 - i5).tolist()
    report['paired'] = paired

    (out_dir / 'acceptance_metrics.json').write_text(
        json.dumps(report, indent=2, sort_keys=True), encoding='utf-8')
    print(f'wrote {out_dir / "acceptance_metrics.json"}')
    print('missing:', report['traces_missing'])


if __name__ == '__main__':
    main()
