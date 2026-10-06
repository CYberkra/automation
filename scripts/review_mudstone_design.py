"""Audit mudstone continuum spectra with gprMax V4; no solver or field ingestion."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import gprMax
from gprMax import config
from gprMax.materials import DispersiveMaterial

from review_material_design_v0_1 import epsilon

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / 'configs/research/mudstone_design_v0_1.json'
ACQUISITION = ROOT / 'configs/research/line9_acquisition_design_v0_1.json'
C = 299792458.0


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out):
    if out.exists():
        raise ValueError('Refuse to overwrite a material review output')
    if gprMax.__version__ != '4.0.0':
        raise ValueError('This audit requires gprMax4.0.0 material API')
    d = json.loads(DESIGN.read_text(encoding='utf-8'))
    companion_path = ROOT / d['companion_material_design']
    companion = json.loads(companion_path.read_text(encoding='utf-8'))
    acq = json.loads(ACQUISITION.read_text(encoding='utf-8'))
    b = d['band_hz']
    f = b['min'] + np.arange(b['count']) * b['step']
    if f[-1] != b['max'] or b != companion['band_hz']:
        raise ValueError('Band or companion design mismatch')
    materials = {**companion['materials'], **d['materials']}
    for name, variant in d['diagnostic_variants_not_solver_batches'].items():
        p = copy.deepcopy(materials[variant['base']]) if 'base' in variant else {}
        p.update(variant)
        materials[name] = p
    results, checks = {}, {}
    for i, (name, p) in enumerate(materials.items()):
        e = epsilon(f, p)
        native = DispersiveMaterial(i+2, name)
        native.type = 'debye'
        native.er, native.se = p['epsilon_infinity'], p['sigma_DC_S_m']
        native.poles = len(p['debye'])
        native.deltaer = [q['delta_epsilon'] for q in p['debye']]
        native.tau = [q['tau_s'] for q in p['debye']]
        official = np.array([native.calculate_er(freq) for freq in f])
        native_error = float(np.linalg.norm(e-official)/np.linalg.norm(official))
        real = np.full(f.shape, p['epsilon_infinity'], dtype=float)
        loss = p['sigma_DC_S_m']/(2*np.pi*f*config.e0)
        for pole in p['debye']:
            x = 2*np.pi*f*pole['tau_s']
            real += pole['delta_epsilon']/(1+x*x)
            loss += pole['delta_epsilon']*x/(1+x*x)
        if native_error > 1e-14 or not np.allclose(e, real-1j*loss, rtol=1e-14, atol=1e-14):
            raise ValueError(f'{name}: constitutive API/formula mismatch')
        if not np.all(np.isfinite(e)) or not np.all(loss>0) or np.any(np.diff(real)>1e-12):
            raise ValueError(f'{name}: finite/passive/monotone check failed')
        if any(q['delta_epsilon']<=0 or q['tau_s']<=0 for q in p['debye']):
            raise ValueError(f'{name}: invalid passive Debye parameter')
        # Independent stable real-valued expression for complex sqrt absorption.
        alpha_independent = 2*np.pi*f/C * loss / np.sqrt(2*(np.hypot(real,loss)+real))
        alpha = -(2*np.pi*f/C*np.sqrt(e)).imag
        if not np.allclose(alpha, alpha_independent, rtol=1e-13):
            raise ValueError(f'{name}: independent absorption check failed')
        rows = []
        for freq in [20,40,95,170]:
            ep = complex(epsilon(freq*1e6, p))
            a = float(-(2*np.pi*freq*1e6/C*np.sqrt(ep)).imag)
            rows.append(dict(frequency_MHz=freq, epsilon_real=ep.real, epsilon_loss=-ep.imag,
                             sigma_effective_S_m=2*np.pi*freq*1e6*config.e0*(-ep.imag),
                             alpha_Np_m=a, amplitude_1_e_m=1/a,
                             absorption_roundtrip_2_8m_layer_dB=-20/np.log(10)*5.6*a,
                             absorption_roundtrip_1_7m_layer_dB=-20/np.log(10)*3.4*a))
        results[name] = rows
        checks[name] = dict(native_relative_L2=native_error,
                           independent_alpha_max_absolute_error=float(np.max(abs(alpha-alpha_independent))))
    m = materials['mudstone']
    ep = epsilon(np.array([20,95,170])*1e6, m)
    if abs(ep[1].real-12)>1e-12 or abs(ep[0].real-ep[2].real-3)>1e-12:
        raise ValueError('Mudstone centre/span constraints failed')
    if acq['start_profile_x_m'] != 220 or acq['end_profile_x_m'] != 25 or acq['length_m'] != 195:
        raise ValueError('User acquisition direction/extent mismatch')
    profile_x = np.array([220,196.705832087,25])
    distance_s = 220-profile_x
    if not np.all(np.diff(distance_s)>0) or not np.allclose(220-distance_s,profile_x):
        raise ValueError('Profile/acquisition coordinate round-trip failed')
    rf = np.array([20,40,95,170])*1e6
    reflection = {}
    for a,b in [('cover','mudstone'),('mudstone','rock')]:
        na, nb = np.sqrt(epsilon(rf,materials[a])),np.sqrt(epsilon(rf,materials[b]))
        r = (na-nb)/(na+nb)
        reflection[f'{a}_to_{b}'] = dict(abs_R=abs(r).tolist(), phase_deg=np.angle(r,deg=True).tolist(),
            definition='Electric-field Fresnel reflection at local normal incidence; not Bscan SNR or a power transmission ratio')
    out.mkdir(parents=True)
    commands = []
    for name in ['cover','mudstone','rock']:
        p = materials[name]
        commands.append(f"#material: {p['epsilon_infinity']:.17g} {p['sigma_DC_S_m']:.17g} 1 0 {name}")
        poles = ' '.join(f"{q['delta_epsilon']:.17g} {q['tau_s']:.17g}" for q in p['debye'])
        if poles:
            commands.append(f"#add_dispersion_debye: {len(p['debye'])} {poles} {name}")
    (out/'materials.txt').write_text('\n'.join(commands)+'\n',encoding='utf-8')
    summary = dict(status='PASS_ANALYTIC_RESEARCH_DESIGN_NOT_SOLVER_OR_SITE_VALIDATION',
        calls_solver=False,calls_training=False,gprMax_version=gprMax.__version__,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [DESIGN,ACQUISITION,companion_path]},
        script_sha256=sha(Path(__file__)),
        shared_epsilon_script_sha256=sha(ROOT/'scripts/review_material_design_v0_1.py'),
        local_materials_py_sha256=sha(Path(gprMax.__file__).parent/'materials.py'),
        checks=checks,table=results,local_normal_reflection=reflection,
        acquisition=dict(length_m=195,profile_x_m=profile_x.tolist(),distance_s_m=distance_s.tolist()),
        thickness_path_scope='2.8m and1.7m are single borehole layer illustrations, not constant layer thickness along the survey; absorption only, no geometrical spreading or interface transmission',
        limitations=d['limitations'])
    (out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':summary['status'],'checks':checks,'mudstone':results['mudstone'],
                      'local_normal_reflection':reflection},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    main(parser.parse_args().out)
