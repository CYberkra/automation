"""Independent acceptance re-computation for benchmark3d_r2_co (2026-10-02).

Verifier-side (house rule: never trust the runner's narrative, recompute from raw h5).
Uses theoretical windows instead of argmax over a wide window:
  direct 0-50 ns | surface 85-125 | coda 130-160 | interface(theory) 165-215
  late rock-bottom check 320-370 | tail = last 10% of trace.
Also correlates interface-echo peak time with cover thickness parsed from each .in
(box column under the Tx-Rx midpoint), and checks the 2D group with and without
the unverified t11 (h5 exists but supervision.json missing => burned attempt).
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
OUT = Path(r'E:\automation_djh\artifacts_check\verify')

GROUPS = {
    '3D-5cm': ('B3D5CM-C3mR2-BG-CO13', [f't{i:02d}' for i in range(1, 14)]),
    '2D-5cm': ('B2D5CM-C3mR2-BG-CO13', [f't{i:02d}' for i in range(1, 12)]),  # t12/t13 not run
}
BURNED = {('2D-5cm', 't11')}  # no supervision.json -> treat as unverified

W = dict(direct=(0, 50), surface=(85, 125), coda=(130, 160),
         iface=(165, 215), late=(320, 370))


def hilbert_env(a):
    n = len(a)
    X = np.fft.fft(a)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1
        h[1:n // 2] = 2
    else:
        h[0] = 1
        h[1:(n + 1) // 2] = 2
    return np.abs(np.fft.ifft(X * h))


def load(run_id):
    d = SIMS / f'2026-10-01_{run_id}'
    with h5py.File(d / f'{run_id}.h5', 'r') as f:
        key = next(iter(f['rxs/rx1'].keys()))
        a = np.asarray(f[f'rxs/rx1/{key}']).ravel()
        dt = float(f.attrs['dt'])
    return a, dt


def parse_in(path):
    txt = path.read_text(encoding='utf-8', errors='replace')
    rx = re.search(r'^#rx:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', txt, re.M)
    rx_x, rx_y, rx_z = float(rx.group(1)), float(rx.group(2)), float(rx.group(3))
    hd = re.search(r'^#hertzian_dipole:\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)', txt, re.M)
    tx_x, tx_y, tx_z = float(hd.group(1)), float(hd.group(2)), float(hd.group(3))
    # cover boxes: '#box: x1 y1 z1 x2 y2 12 cover' ; find column containing midpoint
    mx, my = (rx_x + tx_x) / 2, (rx_y + tx_y) / 2
    best = None
    for m in re.finditer(
            r'^#box:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+cover\s*$',
            txt, re.M):
        x1, y1, z1, x2, y2, z2 = map(float, m.groups())
        if x1 <= mx <= x2 and y1 <= my <= y2 and z2 == 12.0:
            if best is None or z1 < best:
                best = z1
    cover_thk = 12.0 - best if best is not None else None
    return dict(rx_y=rx_y, mid_y=my, cover_thk=cover_thk)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {}
    for gname, (prefix, traces) in GROUPS.items():
        rows = []
        for t in traces:
            rid = f'{prefix}-{t}'
            d = SIMS / f'2026-10-01_{rid}'
            if not (d / f'{rid}.h5').exists():
                continue
            cfg = parse_in(CONFIGS / f'{rid}.in')
            a, dt = load(rid)
            t_ns = np.arange(len(a)) * dt * 1e9
            env = hilbert_env(a)
            e_total = float(np.sum(a ** 2))
            row = dict(trace=t, burned=(gname, t) in BURNED,
                       cover_thk_m=cfg['cover_thk'], mid_y=cfg['mid_y'],
                       finite=bool(np.isfinite(a).all()),
                       tail_ratio=float(np.sum(a[t_ns >= t_ns[-1] * 0.9] ** 2) / e_total))
            for name, (lo, hi) in W.items():
                m = (t_ns >= lo) & (t_ns < hi)
                i = int(np.argmax(env[m]))
                row[f'{name}_t_ns'] = round(float(t_ns[m][i]), 2)
                row[f'{name}_env'] = float(f'{float(env[m][i]):.4e}')
                row[f'{name}_e'] = float(np.sum(a[m] ** 2))
            row['surf_over_direct_e'] = row['surface_e'] / row['direct_e']
            row['iface_over_surf_env'] = row['iface_env'] / row['surface_env']
            row['coda_over_iface_env'] = row['coda_env'] / row['iface_env']
            rows.append(row)
        report[gname] = rows

    (OUT / 'verify_metrics.json').write_text(
        json.dumps(report, indent=1, sort_keys=True), encoding='utf-8')

    for gname, rows in report.items():
        print(f'\n== {gname} ==')
        print(' tr  burned cov_m  surf_t  surf_env   coda_t coda_env   if_t   if_env   if/surf  coda/if  tail')
        for r in rows:
            print(f" {r['trace']} {str(r['burned']):5} {r['cover_thk_m']:.2f} "
                  f"{r['surface_t_ns']:7.2f} {r['surface_env']:.3e} "
                  f"{r['coda_t_ns']:7.2f} {r['coda_env']:.3e} "
                  f"{r['iface_t_ns']:7.2f} {r['iface_env']:.3e} "
                  f"{r['iface_over_surf_env']:.4f} {r['coda_over_iface_env']:.2f} "
                  f"{r['tail_ratio']:.5f}")
        ok = [r for r in rows if not r['burned']]
        thk = np.array([r['cover_thk_m'] for r in ok])
        it = np.array([r['iface_t_ns'] for r in ok])
        if len(ok) > 2 and np.std(thk) > 1e-9:
            c = np.corrcoef(thk, it)[0, 1]
            print(f' cover-thk vs iface_t corr = {c:.3f} '
                  f'(expect positive: thicker cover -> later echo)')
        se = [r['surf_over_direct_e'] for r in ok]
        print(f' surface/direct energy ratio: mean={np.mean(se):.4f} '
              f'min={np.min(se):.4f} max={np.max(se):.4f}')
        ienv = [r['iface_env'] for r in ok]
        senv = [r['surface_env'] for r in ok]
        print(f' iface echo env: mean={np.mean(ienv):.3e}; '
              f'iface/surface env: mean={np.mean(np.array(ienv)/np.array(senv)):.4f}')


if __name__ == '__main__':
    main()
