"""Deep-domain planar background versus independent air + layered line reference.

Analytic air subtraction leaves direct-wave discretization error in the residual;
this is not a same-grid air-paired reflection or an isolated target validation.
"""
import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import c, mu_0
from scipy.special import hankel2
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response
from analyze_vertical_refinement import metrics

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
reference = Path('artifacts/research_checks/2026-09-25_yz_layers_reference/spectra.npz')
with np.load(reference) as d:
    f, layer = d['frequency_Hz'], d['covered']
air = -2*np.pi*f*mu_0/4*hankel2(0, 2*np.pi*f/c*1.3)
total = air+layer
saved = dict(frequency_Hz=f, analytic_air=air, analytic_layer=layer, analytic_total=total)
inputs = [reference, Path(__file__)]
rows = {}
for mode, shape in [('FINE', [1,2560,4000]), ('ZFINE', [1,2560,8000]), ('YFINE', [1,5120,8000])]:
    path = Path(f'artifacts/research_checks/2026-09-25_DEP_BG_{mode}/DEP_BG_{mode}.h5')
    inputs.append(path)
    src = load_source(path)
    rx = load_receiver(path, receiver_path='name:measurement', component='Ex')
    with h5py.File(path) as h:
        assert np.array_equal(h.attrs['nx_ny_nz'], shape)
        assert np.allclose(h[rx.path].parent.attrs['Position'][1:], [16.65,45], rtol=0, atol=1e-12)
        assert h[rx.path].dtype == np.float64
    assert src.dt == rx.dt and src.time_offset == src.dt/2 and rx.time_offset == 0
    assert np.all(np.isfinite(rx.samples))
    n = min(len(rx.samples), int(np.floor(1200e-9/rx.dt))+1)
    row = {}
    for taper in (200,400):
        result = direct_frequency_response(src, replace(rx, samples=rx.samples[:n]), f,
            tail_taper_fraction=(round(taper*1e-9/rx.dt)-.25)/n)
        assert np.all(result.source_valid) and np.all(np.isfinite(result.response))
        h = result.response
        saved[f'{mode}_{taper}'] = h
        row[str(taper)] = dict(total=metrics(h,total), analytic_air_subtracted=metrics(h-air,layer),
            total_relative_L2=float(np.linalg.norm(h-total)/np.linalg.norm(total)),
            error_relative_to_layer_L2=float(np.linalg.norm(h-total)/np.linalg.norm(layer)))
    row['tail_relative_to_layer_L2'] = float(np.linalg.norm(saved[mode+'_200']-saved[mode+'_400'])/np.linalg.norm(layer))
    rows[mode] = row
result = dict(rows=rows, solver_invoked=False, target_validated=False, physical_acceptance_threshold=None,
    caveat='Analytic-air-subtracted diagnostic includes numerical direct-wave error; total errors can be hidden by strong direct wave. Does not validate 20m target differences.',
    inputs={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs})
(a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
np.savez_compressed(a.output/'arrays.npz',**saved)
print(json.dumps(rows))
