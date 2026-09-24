"""Constructed role-aware diagnostics; not a flight line or physical label set."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import c, mu_0
from scipy.special import hankel2
from research_operator_contract import catalogue, apply_configuration, ConfigUnavailable
from research_evaluation_contract import waveform_metrics, energy_metrics

p = argparse.ArgumentParser()
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
template = Path('artifacts/research_checks/2026-09-25_three_grid_mechanisms/arrays.npz')
reference = Path('artifacts/research_checks/2026-09-25_yz_layers_reference/spectra.npz')
with np.load(template) as d:
    saved = {k: d[k] for k in ('time_s', 'common', 'localized', 'coarse_target', 'fine_target', 'zfine_target')}
with np.load(reference) as d:
    f, layer = d['frequency_Hz'], d['covered']
t = np.arange(8192)/(8192*300e3)
keep = t <= 1200e-9
assert np.array_equal(t[keep], saved['time_s'])
window = np.hanning(len(f))
def wave(s):
    return (2*np.real(np.fft.ifft(s*window, n=8192)*8192/window.sum()*np.exp(2j*np.pi*f[0]*t)))[keep]
shape = (sum(keep), 32)
direct = -2*np.pi*f*mu_0/4*hankel2(0, 2*np.pi*f/c*1.3)
A = np.broadcast_to(wave(direct)[:, None], shape).copy()
L = np.broadcast_to(wave(layer)[:, None], shape).copy()
saved.update(direct=A, layer=L)
for role, low, high in [('interference', 0, 50), ('layer', 90, 240), ('target', 440, 600)]:
    saved[role+'_mask'] = np.broadcast_to(((t[keep]>=low*1e-9)&(t[keep]<=high*1e-9))[:, None], shape).copy()
def contrast(y, truth, role):
    result = waveform_metrics(y, truth, saved[role+'_mask'], reference_kind='paired_contrast', state='isolated', scope='contrast')
    assert result['available'] and not result['absolute_preservation_eligible']
    return result['metrics']
rows = []
for config in [x for x in catalogue() if x['end_gain']==1]:
    cid = config['id']
    for grid in ('coarse', 'fine', 'zfine'):
        for profile in ('common', 'localized'):
            row = dict(config_id=cid, grid=grid, profile=profile)
            try:
                S = saved[grid+'_target'][:, None]*saved[profile]
                yA = apply_configuration(A, cid)['output']
                yAL = apply_configuration(A+L, cid)['output']
                yALS = apply_configuration(A+L+S, cid)['output']
                row.update(available=True, interference=energy_metrics(yA, A, saved['interference_mask']), layer=contrast(yAL-yA, L, 'layer'), target=contrast(yALS-yAL, S, 'target'))
                saved[cid+'_direct_output'] = yA
                saved[cid+'_layer_background_output'] = yAL
                saved[grid+'_'+profile+'_'+cid+'_output'] = yALS
            except ConfigUnavailable as e:
                row.update(available=False, reason=str(e))
            rows.append(row)
for row in rows:
    if row['config_id']=='B0_G1_BG':
        assert row['layer']['nrmse']<1e-9 and row['target']['nrmse']<1e-9
        assert abs(row['interference']['energy_ratio']-1)<1e-12
    if row['config_id']=='B3_G1_BG':
        assert row['interference']['energy_ratio']<1e-20 and row['layer']['nrmse']>.999999
        if row['profile']=='common': assert row['target']['nrmse']>.999999
result = dict(rows=rows, solver_invoked=False, training_eligible=False, unique_label=None,
    physical_target_reference_state='numerically_unresolved', constructed_only=True,
    checks=dict(identity_preserves_both_roles=True, full_mean_suppresses_direct_and_deletes_layer=True),
    inputs={str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in (template, reference, Path(__file__), Path('docs/research/2026-09-25_joint_roles.md'))})
(a.output/'results.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
np.savez_compressed(a.output/'arrays.npz', **saved)
print(json.dumps(dict(rows=len(rows), available=sum(r['available'] for r in rows), checks=result['checks'])))
