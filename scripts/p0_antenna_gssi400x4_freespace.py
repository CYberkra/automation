"""P0 link check: scaled GSSI-400 x4 antenna in free space, GPU double, no subgrid.

Coarse-grid approximation is DECLARED: the x4-scaled antenna has 8 mm-class
thin details (case walls, PCB) quantized onto the 4 cm main grid without a
subgrid (V4 subgrids are CPU-only). This run validates (a) build+solve,
(b) the official SFCW actual-source chain on a resistive voltage source,
(c) antenna ringing timescale. Not a physics acceptance run.
"""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import gprMax
from antenna_gssi400_x4 import scaled_antenna, coarse_adapt, CASE_SIZE
from hs_capsule_identity import sha256

OUT = ROOT / 'artifacts/research_checks/2026-10-05_antenna_gssi400x4_p0'
DL = 0.04
DOMAIN = (2.4, 2.4, 1.6)
ANT = (1.2, 1.2, 0.4)
TIME_WINDOW = 100e-9


def build():
    scene = gprMax.Scene()
    scene.add(gprMax.Title(name='gssi400x4_p0_freespace'))
    scene.add(gprMax.Domain(p1=DOMAIN))
    scene.add(gprMax.Discretisation(p1=(DL, DL, DL)))
    scene.add(gprMax.TimeWindow(time=TIME_WINDOW))
    objs, log = coarse_adapt(scaled_antenna(*ANT), DL)
    (OUT / 'p0_coarse_adapt_log.json').write_text(json.dumps(log, indent=1, default=str) + '\n',
                                                 encoding='utf-8')
    for obj in objs:
        scene.add(obj)
    return scene


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    fn = OUT / 'p0_freespace'
    scene = build()
    gprMax.run(scenes=[scene], geometry_only=False, outputfile=str(fn) + '.h5',
               gpu=[0], gpu_precision='double', subgrid=False)
    h5 = fn.with_suffix('.h5')
    import h5py
    info = {'h5_sha256': sha256(h5)}
    with h5py.File(h5) as h:
        info['attrs'] = {k: (v.tolist() if hasattr(v, 'tolist') else v)
                         for k, v in h.attrs.items() if k in ('gprMax', 'dt', 'Iterations', 'nx_ny_nz', 'dx_dy_dz')}
        info['groups'] = sorted(h.keys())
        info['srcs'] = sorted(h['srcs'].keys()) if 'srcs' in h else None
        if 'srcs' in h:
            for s in h['srcs']:
                info[f'srcs/{s}'] = {'subgroups': sorted(h[f'srcs/{s}'].keys()),
                                     'attrs': {k: v for k, v in h[f'srcs/{s}'].attrs.items()}}
                if 'excitation' in h[f'srcs/{s}']:
                    e = h[f'srcs/{s}/excitation']
                    info[f'srcs/{s}']['excitation_attrs'] = dict(e.attrs)
                    info[f'srcs/{s}']['excitation_samples'] = int(e['samples'].shape[0])
        if 'rxs/rx1/Ey' in h:
            y = h['rxs/rx1/Ey'][:]
            info['rx_dtype'] = str(y.dtype)
            info['rx_finite'] = bool(np.isfinite(y).all())
            info['rx_len'] = int(y.shape[0])
            t = np.arange(y.shape[0]) * float(h.attrs['dt'])
            env = np.abs(y)
            pk = int(np.argmax(env))
            info['rx_peak_time_ns'] = float(t[pk] * 1e9)
            info['rx_peak_abs'] = float(env[pk])
            after = env[env < env[pk] * 0.1]
            info['rx_time_below_10pct_ns'] = float(t[len(env) - len(after)] * 1e9) if len(after) else None
    (OUT / 'p0_h5_structure.json').write_text(json.dumps(info, ensure_ascii=False, indent=1,
                                                         default=str) + '\n', encoding='utf-8')
    print(json.dumps(info, ensure_ascii=False, default=str)[:1500])

    # SFCW chain check on the voltage-source output
    from sfcw_official_loader_v0_2 import FREQ
    from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response
    s = load_source(h5)
    rx = load_receiver(h5, receiver_path='/rxs/rx1', component='Ey')
    take = FREQ[[0, 250, 500]]
    res = direct_frequency_response(s, rx, take, tail_taper_fraction=0.2)
    print('SFCW direct response at 20/105/170 MHz:', res.response)
    print('source_valid:', res.source_valid)
    # antenna port characterisation in our band
    import h5py as h5m
    with h5m.File(h5) as h:
        f = h['ports/port1/frequency'][:]
        s11 = h['ports/port1/S11'][:]
        m = (f >= 20e6) & (f <= 170e6)
        s11_db = 20 * np.log10(np.maximum(np.abs(s11[m]), 1e-30))
        band = {'f_mhz': f[m] / 1e6, 's11_db': s11_db}
        info_port = {'s11_db_min': float(s11_db.min()), 's11_db_max': float(s11_db.max()),
                     's11_db_median': float(np.median(s11_db)),
                     'valid_frequency_range_Hz': h['ports/port1'].attrs['ValidFrequencyRange'].tolist(),
                     'mesh_frequency_limit_Hz': float(h['ports/port1'].attrs['MeshFrequencyLimit'])}
        print('port S11 in 20-170MHz:', json.dumps(info_port))
        (OUT / 'p0_port_s11_band.json').write_text(json.dumps(
            {'port': info_port, 'f_mhz': band['f_mhz'].tolist(),
             's11_db': band['s11_db'].tolist()}, indent=1) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
