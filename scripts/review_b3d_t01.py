# First-trace review for B3D5CM-C3mR2-BG-CO13-t01 (single A-scan)
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

RUN = Path(r'E:\automation_djh\automation_repo\artifacts\simulations\2026-10-01_B3D5CM-C3mR2-BG-CO13-t01')
h5_path = RUN / 'B3D5CM-C3mR2-BG-CO13-t01.h5'
out_png = Path(r'E:\automation_djh\fig_b3d_t01_first_trace.png')

with h5py.File(h5_path, 'r') as f:
    print('datasets:', list(f.keys()))
    rx_path = 'rxs/rx1/Ex'
    Ex = f[rx_path][:]
    dt = f.attrs['dt']
    nsteps = Ex.shape[0]

t = np.arange(nsteps) * dt * 1e9  # ns
Ex = np.asarray(Ex).ravel()

finite = np.isfinite(Ex).all()
peak = float(np.max(np.abs(Ex)))
peak_t = float(t[np.argmax(np.abs(Ex))])

# window-tail energy: last 10% of window vs full
n_tail = nsteps // 10
tail_ratio = float(np.sum(Ex[-n_tail:] ** 2) / np.sum(Ex ** 2))

# direct-wave window energy (first 20%) vs interface-echo window
head_ratio = float(np.sum(Ex[: nsteps // 5] ** 2) / np.sum(Ex ** 2))

print(f'nsteps={nsteps} dt={dt*1e12:.2f} ps window={t[-1]:.1f} ns')
print(f'finite={finite} peak={peak:.3e} at {peak_t:.2f} ns')
print(f'tail_energy_ratio(last10%)={tail_ratio:.5f}')
print(f'head_energy_ratio(first20%)={head_ratio:.4f}')

fig, ax = plt.subplots(figsize=(11, 4.2), dpi=130)
ax.plot(t, Ex, lw=0.8)
ax.axvspan(t[-1] * 0.9, t[-1], color='orange', alpha=0.2, label='tail 10%')
ax.set_xlabel('t (ns)')
ax.set_ylabel('Ex (V/m)')
ax.set_title('B3D5CM-C3mR2-BG-CO13-t01  Ex @ rx (6, 5.15, 27)')
ax.legend(loc='upper right')
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(out_png)
print('saved:', out_png)

report = dict(nsteps=int(nsteps), dt_ps=float(dt * 1e12), window_ns=float(t[-1]),
              finite=bool(finite), peak=float(peak), peak_t_ns=float(peak_t),
              tail_energy_ratio=tail_ratio, head_energy_ratio=head_ratio)
(RUN / 'first_trace_review.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print('review json written')
