"""Analytic material design audit against installed V4; never calls a solver."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import gprMax
from gprMax import config
from gprMax.materials import DispersiveMaterial

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / 'configs/research/material_design_cover_rock_v0_1.json'
C = 299792458.0


def epsilon(f, p):
    w = 2*np.pi*np.asarray(f)
    e = p['epsilon_infinity'] - 1j*p['sigma_DC_S_m']/(w*config.e0)
    for pole in p['debye']:
        e = e + pole['delta_epsilon']/(1+1j*w*pole['tau_s'])
    return e


def main(out):
    if out.exists():
        raise ValueError('Refuse to overwrite a material review output')
    if gprMax.__version__ != '4.0.0':
        raise ValueError('The audited local V4.0.0 material API is required')
    d = json.loads(DESIGN.read_text(encoding='utf-8'))
    f = np.arange(501)*.3e6+20e6
    all_p = dict(d['materials'], **d['literature_controls_not_adopted'])
    results = {}
    checks = {}
    for n, (name, p) in enumerate(all_p.items()):
        native = DispersiveMaterial(n+2, name)
        native.type = 'debye'; native.er = p['epsilon_infinity']; native.se = p['sigma_DC_S_m']
        native.poles = len(p['debye'])
        native.deltaer = [q['delta_epsilon'] for q in p['debye']]
        native.tau = [q['tau_s'] for q in p['debye']]
        e = epsilon(f, p)
        official = np.array([native.calculate_er(v) for v in f])
        error = np.linalg.norm(e-official)/np.linalg.norm(official)
        # Explicit independent real/loss expressions (no complex division).
        real = np.full(f.shape, p['epsilon_infinity'], dtype=float)
        loss = p['sigma_DC_S_m']/(2*np.pi*f*config.e0)
        for q in p['debye']:
            x = 2*np.pi*f*q['tau_s']
            real += q['delta_epsilon']/(1+x*x)
            loss += q['delta_epsilon']*x/(1+x*x)
        if error > 1e-14 or not np.allclose(e, real-1j*loss, rtol=1e-14):
            raise ValueError('Constitutive frequency-response check failed')
        if not (np.all(np.isfinite(e)) and np.all(-e.imag > 0) and np.all(np.diff(e.real) < 0)):
            raise ValueError('Expected finite, lossy, monotone band response')
        values = epsilon(np.array([20, 40, 95, 170])*1e6, p)
        rows = []
        for freq, v in zip([20,40,95,170], values):
            alpha = -(2*np.pi*freq*1e6/C*np.sqrt(v)).imag
            rows.append(dict(frequency_MHz=freq, epsilon_real=float(v.real), epsilon_loss=float(-v.imag),
                             alpha_Np_m=float(alpha), amplitude_1_e_m=float(1/alpha),
                             absorption_37m_dB=float(-20/np.log(10)*37*alpha)))
        results[name] = rows
        checks[name] = dict(native_calculate_er_relative_L2=float(error))
    r = d['materials']['rock']
    rock = epsilon(np.array([20,95,170])*1e6, r)
    if abs(rock[1].real-9)>1e-12 or abs(rock[0].real-rock[2].real-1.5)>1e-12:
        raise ValueError('User center/span constraints failed')
    old = (ROOT/d['controlled_comparison']['baseline_input']).read_text(encoding='utf-8')
    if '#material: 11 0.001 1 0 cover' not in old or '#add_dispersion_debye: 1 0.5 6.4567e-09 cover' not in old:
        raise ValueError('Baseline cover identity differs')
    ep_cover = epsilon(np.array([20,40,95,170])*1e6, d['materials']['cover'])
    ep_rock = epsilon(np.array([20,40,95,170])*1e6, r)
    reflection = (np.sqrt(ep_cover)-np.sqrt(ep_rock))/(np.sqrt(ep_cover)+np.sqrt(ep_rock))
    out.mkdir(parents=True)
    commands = []
    for name, p in d['materials'].items():
        commands.append(f"#material: {p['epsilon_infinity']:.17g} {p['sigma_DC_S_m']:.17g} 1 0 {name}")
        poles = ' '.join(f"{q['delta_epsilon']:.17g} {q['tau_s']:.17g}" for q in p['debye'])
        commands.append(f"#add_dispersion_debye: {len(p['debye'])} {poles} {name}")
    (out/'materials.txt').write_text('\n'.join(commands)+'\n', encoding='utf-8')
    result = dict(status='PASS_FIXED_MATERIAL_DESIGN_NOT_SOLVER_OR_SITE_VALIDATION',
                  calls_solver=False, calls_training=False, gprMax_version=gprMax.__version__,
                  parameters_sha256=hashlib.sha256(DESIGN.read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  local_materials_py_sha256=hashlib.sha256((Path(gprMax.__file__).parent/'materials.py').read_bytes()).hexdigest(),
                  checks=checks, center_span_absolute_tolerance=1e-12, table=results,
                  local_normal_cover_rock_abs_R=abs(reflection).tolist(),
                  local_normal_cover_rock_R_phase_deg=np.angle(reflection,deg=True).tolist(),
                  scope='Continuum analytic frequency response only; no FDTD constitutive-update convergence certification')
    (out/'summary.json').write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    main(parser.parse_args().out)
