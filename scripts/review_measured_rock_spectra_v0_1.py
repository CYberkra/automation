"""Reconstruct published CPCM fits, not raw measurements or a gprMax run."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

EPS0 = 8.8541878128e-12
C = 299792458.0
ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / 'docs/research/2026-10-06_measured_dielectric_sources.json'
BASE = ROOT / 'artifacts/research_checks/2026-10-05_slope_bscan_2d_fine2m_r1/slope_rough_s00/profile.in'


def cpcm(frequency, p):
    """Bore et al. (2025), Eq.23; exp(+jwt), epsilon=real-j*loss."""
    w = 2 * np.pi * np.asarray(frequency, dtype=float)
    return (p['epsilon_infinity']
            + p['delta_P1'] / (1 + 1j * w * p['tau_P1_s'])
            + p['delta_P2'] / (1 + (1j * w * p['tau_P2_s']) ** p['c_P2'])
            - 1j * p['sigma0_S_m'] / (w * EPS0)
            * (1 + p['m_P3'] / (1 - p['m_P3'])
               * (1 - 1 / (1 + (1j * w * p['tau_P3_s']) ** p['c_P3']))))


def baseline(frequency, er, sigma, delta=0., tau=1.):
    w = 2 * np.pi * np.asarray(frequency, dtype=float)
    return er + delta / (1 + 1j * w * tau) - 1j * sigma / (w * EPS0)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out):
    if out.exists():
        raise ValueError('Refusing to overwrite an existing research output')
    source = json.loads(LEDGER.read_text(encoding='utf-8'))
    text = BASE.read_text(encoding='utf-8')
    for line in ['#material: 9 0.001 1 0 rock',
                 '#material: 11 0.001 1 0 cover',
                 '#add_dispersion_debye: 1 0.5 6.4567e-09 cover']:
        if line not in text.splitlines():
            raise ValueError('Current baseline material identity changed')
    states = source['published_rock_parameters_SI']
    f = np.arange(501, dtype=float) * 0.3e6 + 20e6
    eps = {'current_cover': baseline(f, 11, .001, .5, 6.4567e-9),
           'current_rock': baseline(f, 9, .001)}
    eps.update({p['id']: cpcm(f, p) for p in states})
    # Closed-form special cases check the signs and the two representations.
    p = dict(states[0], delta_P1=0., delta_P2=0., m_P3=0.)
    if not np.allclose(cpcm(f, p), baseline(f, p['epsilon_infinity'], p['sigma0_S_m']), rtol=1e-14):
        raise ValueError('DC conductivity special case failed')
    p = dict(states[0], delta_P1=0., delta_P2=0., c_P3=1.)
    w = 2*np.pi*f
    expected = (p['epsilon_infinity'] + p['sigma0_S_m']*p['m_P3']*p['tau_P3_s']
                / ((1-p['m_P3'])*EPS0*(1+1j*w*p['tau_P3_s']))
                - 1j*p['sigma0_S_m']/(w*EPS0))
    if not np.allclose(cpcm(f, p), expected, rtol=1e-14):
        raise ValueError('Cole-Cole conductivity to Debye special case failed')
    rows = []
    for name, e in eps.items():
        if not (np.all(np.isfinite(e)) and np.all(e.real > 0) and np.all(e.imag < 0)):
            raise ValueError('Nonfinite or nonpassive spectrum')
        reflection = ((np.sqrt(eps['current_cover']) - np.sqrt(e))
                      / (np.sqrt(eps['current_cover']) + np.sqrt(e)))
        alpha = -(2*np.pi*f/C*np.sqrt(e)).imag
        for i in range(f.size):
            rows.append([name, f[i], e[i].real, -e[i].imag, alpha[i],
                         abs(reflection[i]), np.angle(reflection[i], deg=True)])
    out.mkdir(parents=True)
    with (out/'spectra.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['model', 'frequency_Hz', 'epsilon_real', 'epsilon_loss',
                         'alpha_Np_m', 'normal_cover_interface_abs_R', 'normal_cover_interface_phase_deg'])
        writer.writerows(rows)
    selected = [r for r in rows if any(abs(r[1]-x*1e6) < .001 for x in (20, 40.1, 95, 170))]
    summary = dict(status='PUBLISHED_FIT_RECONSTRUCTION_NOT_RAW_DATA_FIT',
                   source_git_base='8e8a387', calls_solver=False, calls_training=False,
                   field_fitted=False, new_Debye_parameters_fitted=False,
                   input_sha256=sha(BASE), ledger_sha256=sha(LEDGER), script_sha256=sha(Path(__file__)),
                   grid=dict(start_MHz=20, stop_MHz=170, step_MHz=.3, count=501),
                   checks=['baseline_input_identity', 'DC_limit_special_case',
                           'conductivity_Debye_equivalence_special_case', 'finite_passive_band'],
                   representative_rows=selected,
                   notes=['40.1 MHz is the nearest actual sweep bin to 40 MHz.',
                          'Reflection is local normal-incidence planar cover-to-rock electric-field R.',
                          'SS-A01 theta3 is diagnostic only: table tau_P3=13.3ms exceeds stated 1ms bound.',
                          'Published full-band RMS is not the band error or raw measurement uncertainty.'])
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    for color_id, (name, e) in enumerate(eps.items()):
        color = plt.get_cmap('tab10')(color_id)
        style = '--' if name.startswith('current') else '-'
        label = name + (' (diagnostic only)' if name.endswith('theta3') else '')
        axes[0].plot(f/1e6, e.real, style, color=color, label=label)
        axes[1].plot(f/1e6, -e.imag, style, color=color)
        if name != 'current_cover':
            r = (np.sqrt(eps['current_cover'])-np.sqrt(e))/(np.sqrt(eps['current_cover'])+np.sqrt(e))
            axes[2].plot(f/1e6, abs(r), style, color=color)
    axes[0].set_ylabel('Relative permittivity: real part')
    axes[1].set_ylabel('Relative permittivity: effective loss')
    axes[2].set_ylabel('|R|, normal cover-to-rock interface')
    for ax in axes:
        ax.set_xlabel('Frequency (MHz)')
        ax.grid(alpha=.25)
    axes[0].legend(fontsize=8)
    fig.suptitle('Published sandstone fitted curves; no raw-data refit or FDTD validation')
    fig.savefig(out/'published_rock_spectra.png', dpi=160)
    plt.close(fig)
    print(json.dumps(dict(status=summary['status'], output=str(out), samples=len(rows))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    main(parser.parse_args().out)
