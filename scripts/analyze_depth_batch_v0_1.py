"""Depth-factor analysis: full forward vs 95MHz-style analytic attenuation extrapolation.

Isolates the interface response (flat case minus shared full-cover background;
target/background definition declared), then compares each depth against a
spectrally resolved extrapolation from the 3 m case (per-frequency two-way
attenuation from the archived material plus cylindrical spreading) and against
the old monochromatic 95 MHz extrapolation style. Reports arrival-time
residuals, amplitude ratios, complex errors and the window-tail check.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from analyze_uav_local_free_space import FREQ, transfer
from green_halfspace_hed_v0_1 import C_SI, epsr_cover
from hs_capsule_identity import sha256

DEPTHS = {'depth3': 3., 'depth10': 10., 'depth18_5': 18.5}
CASES = ('fullcover', 'depth3', 'depth10', 'depth18_5')
AIR_NS = 16./C_SI*1e9  # 8 m down + 8 m up


def attenuation_np_m(freqs):
    """One-way amplitude attenuation constant of the archived cover (Np/m)."""
    e = epsr_cover(freqs)
    k = 2*np.pi*freqs/C_SI*np.sqrt(e)
    return -k.imag if (k.imag < 0).all() else k.imag*0  # decay: Im(k)<0 with e^{+iwt}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('fresh analysis output required')
    a.out.mkdir(parents=True)
    c = json.loads((a.study/'execution_contract.json').read_text('utf-8'))
    audit = json.loads((a.study/'completed_verification.json').read_text('utf-8'))
    if audit['status'] != 'PASS' or audit['contract_sha256'] != sha256(a.study/'execution_contract.json'):
        raise ValueError('completed raw identity required')
    identities = {r['id']: r['raw_sha256'] for r in audit['groups']}
    waves, spectra = {}, {}
    for g in c['groups']:
        raw = Path(g['input']).with_suffix('.h5')
        if sha256(raw) != identities[g['id']]:
            raise ValueError('raw changed')
        value, source, receiver = transfer(raw, 1, 'Ey', g['spacing_m'])
        take = np.array([0, 83, 167, 250, 333, 417, 500])
        sy = source.samples.copy(); ry = receiver.samples.copy()
        n = round(20e-9/receiver.dt)
        ry[-n:] *= .5*(1+np.cos(np.linspace(0, np.pi, n)))
        sft = np.exp(-2j*np.pi*FREQ[take, None]*(source.time_offset+source.dt*np.arange(len(sy))))@sy
        rft = np.exp(-2j*np.pi*FREQ[take, None]*(receiver.time_offset+receiver.dt*np.arange(len(ry))))@ry
        independent = rft/sft/g['spacing_m']
        dft_error = float(np.linalg.norm(value[take]-independent)/np.linalg.norm(independent))
        if dft_error > 1e-9:
            raise ValueError('independent direct DFT mismatch')
        waves[g['id']] = (receiver.dt*np.arange(len(receiver.samples))*1e9, receiver.samples)
        spectra[g['id']] = value
    alpha = attenuation_np_m(FREQ)
    n_eff = np.sqrt(epsr_cover(FREQ).real)
    iso = {case: spectra[case]-spectra['fullcover'] for case in DEPTHS}
    rows = []
    for case, depth in DEPTHS.items():
        # spectrally resolved extrapolation from the 3 m case: two-way extra loss + cylindrical spreading
        extra = depth-3.
        r1 = 8./np.sqrt(epsr_cover(np.array([95e6])).real[0]) + 2*3.
        spread = np.sqrt((8./n_eff+2*3.)/(8./n_eff+2*depth))  # 2D line source ~ 1/sqrt(r)
        extrap = iso['depth3']*spread*np.exp(-2*alpha*extra)
        got = iso[case]
        err = float(np.linalg.norm(got-extrap)/np.linalg.norm(got)) if np.linalg.norm(got) else float('nan')
        # old monochromatic 95 MHz style: single-frequency attenuation of the band-integrated response
        a95 = attenuation_np_m(np.array([95e6]))[0]
        mono = float(np.exp(-2*a95*extra))
        full_band = float(np.linalg.norm(got)/np.linalg.norm(iso['depth3']))
        # arrival-time residual: isolated time-domain peak vs ray prediction at 95 MHz
        t, w = waves[case]
        _, wb = waves['fullcover']
        n95 = float(np.sqrt(epsr_cover(np.array([95e6])).real[0]))
        pred = AIR_NS + 2*depth*n95/C_SI*1e9
        diff = w-wb
        env = np.abs(diff)
        lo = max(0, int(np.searchsorted(t, pred-40)))
        hi = min(len(t), int(np.searchsorted(t, pred+40)))
        peak_t = float(t[lo+int(np.argmax(env[lo:hi]))])
        tail = float(np.sqrt(np.mean(w[int(len(w)*0.9):]**2))/abs(w).max())
        rows.append({'case': case, 'depth_m': depth,
                     'isolated_over_fullcover_band_L2': float(np.linalg.norm(got)/np.linalg.norm(spectra['fullcover'])),
                     'isolated_over_depth3_band': full_band,
                     'mono95_extrapolation_factor': mono,
                     'mono95_vs_fullforward_dB': float(20*np.log10(mono/full_band)) if full_band else None,
                     'spectral_extrapolation_relative_L2_error': err,
                     'spectral_extrapolation_max_phase_deg': float(abs(np.angle(got/extrap)).max()*180/np.pi),
                     'predicted_arrival_ns': float(pred), 'isolated_peak_ns': peak_t,
                     'arrival_residual_ns': peak_t-float(pred),
                     'window_tail_rms_over_peak': tail})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    labels = {'fullcover': '全覆盖背景', 'depth3': '平界面3m', 'depth10': '平界面10m', 'depth18_5': '平界面18.5m'}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for case in CASES:
        t, w = waves[case]
        axes[0, 0].plot(t, w, lw=.7, label=labels[case])
        axes[1, 0].plot(FREQ/1e6, 20*np.log10(abs(spectra[case])/abs(spectra[case]).max()+1e-30), lw=.8, label=labels[case])
    axes[0, 0].set(title='总场 A-scan（单站，2D线源族，非完整B-scan）', xlabel='时间(ns)', ylabel='Ey(V/m)', xlim=(0, 650))
    axes[1, 0].set(title='源归一化复频响（各自峰值归一 dB）', xlabel='频率(MHz)', ylabel='dB')
    for case in DEPTHS:
        axes[0, 1].plot(FREQ/1e6, abs(iso[case]), lw=.8, label=labels[case])
    axes[0, 1].set(title='隔离界面响应 |E界面/(Iℓ)|（平界面−全覆盖背景，声明目标对比）', xlabel='频率(MHz)')
    depths = [DEPTHS[k] for k in DEPTHS]
    band = [r['isolated_over_depth3_band'] for r in rows]
    mono = [r['mono95_extrapolation_factor'] for r in rows]
    axes[1, 1].plot(depths, 20*np.log10(band), 'o-', label='完整正演（频带L2幅比）')
    axes[1, 1].plot(depths, 20*np.log10(mono), 's--', label='旧95MHz单频衰减外推')
    axes[1, 1].set(title='界面响应随深度衰减：正演 vs 单频外推（相对3m基准）', xlabel='界面深度(m)', ylabel='dB')
    for ax in axes.flat:
        ax.grid(alpha=.2); ax.legend(fontsize=7)
    fig.suptitle('深度因素批：3/10/18.5m平界面+共享全覆盖背景（同域同材料同网格，仅深度变）')
    fig.savefig(a.out/'depth_factor_comparison.png', dpi=140)
    plt.close(fig)
    summary = {'status': 'PASS', 'contract_sha256': sha256(a.study/'execution_contract.json'),
               'analysis_code_sha256': sha256(__file__), 'rows': rows,
               'isolation_definition': 'interface response = flat-interface case minus shared full-cover background (declared target/background pair, not a universal clean reference)',
               'scope': '2D line-source mechanism screening; depths are mechanism points, not Line9 truth; not 3D absolute amplitude or field validation.'}
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'rows': rows}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
