"""Independent Sommerfeld half-space comparison of actual-clock native V4 outputs."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from analyze_uav_local_free_space import FREQ, dipole_field, transfer
from green_halfspace_hed_v0_1 import epsr_cover, reflected_ex, selftest as reference_selftest
from hs_capsule_identity import sha256

RHO = 1.3      # inline baseline (m)
ZSUM = 16.0    # tx + rx height above surface (m)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('fresh analysis output required')
    a.out.mkdir(parents=True)
    reference_selftest(a.out/'reference_selftest')
    ref = json.loads((a.out/'reference_selftest'/'selftest.json').read_text('utf-8'))
    if ref['status'] != 'PASS':
        raise ValueError('independent reference self-test failed; no comparison authorized')
    c = json.loads((a.study/'execution_contract.json').read_text('utf-8'))
    audit = json.loads((a.study/'completed_verification.json').read_text('utf-8'))
    if audit['status'] != 'PASS' or audit['contract_sha256'] != sha256(a.study/'execution_contract.json'):
        raise ValueError('completed raw identity required')
    identities = {r['id']: r['raw_sha256'] for r in audit['groups']}
    direct = dipole_field(np.array([RHO, 0., 0.]), 'x')[:, 0]
    reflected = reflected_ex(FREQ, RHO, ZSUM, epsr_cover(FREQ))
    exact = {'freespace': direct, 'halfspace': direct + reflected}
    metrics, products, native = [], {}, {}
    for g in c['groups']:
        raw = Path(g['input']).with_suffix('.h5')
        if sha256(raw) != identities[g['id']]:
            raise ValueError('raw changed')
        value, source, receiver = transfer(raw, 1, 'Ex', g['spacing_m'])
        want = exact[g['id']]
        phase = np.angle(value/want)*180/np.pi
        error = float(np.linalg.norm(value-want)/np.linalg.norm(want))
        # Independent seven-frequency direct sum, with actual staggered clocks.
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
        metrics.append({'case': g['id'], 'raw_sha256': sha256(raw),
                        'complex_relative_L2_to_reference': error,
                        'maximum_phase_error_deg': float(abs(phase).max()),
                        'magnitude_error_dB_minmax': [float((20*np.log10(abs(value/want))).min()),
                                                      float((20*np.log10(abs(value/want))).max())],
                        'independent_DFT_relative_L2': dft_error,
                        'passes_declared_diagnostic': bool(error <= .05 and abs(phase).max() <= 5.)})
        products[g['id']] = value
        native[g['id']] = (receiver.dt*np.arange(len(receiver.samples))*1e9, receiver.samples)
    # FDTD isolated reflection (ground minus free-space control) vs Sommerfeld reflection.
    fdtd_reflection = products['halfspace'] - products['freespace']
    refl_error = float(np.linalg.norm(fdtd_reflection-reflected)/np.linalg.norm(reflected))
    refl_phase = float(abs(np.angle(fdtd_reflection/reflected)).max()*180/np.pi)
    cross = {'fdtd_isolated_reflection_vs_sommerfeld_relative_L2': refl_error,
             'fdtd_isolated_reflection_vs_sommerfeld_max_phase_deg': refl_phase,
             'reflection_over_direct_median_dB': float(np.median(20*np.log10(abs(reflected/direct))))}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    labels = {'halfspace': '半空间(覆盖层)', 'freespace': '自由空间对照'}
    for gid in ('halfspace', 'freespace'):
        time, wave = native[gid]
        axes[0, 0].plot(time, wave, label=labels[gid]+' FDTD')
        axes[0, 1].plot(FREQ/1e6, abs(products[gid]), label=labels[gid]+' FDTD')
        axes[0, 1].plot(FREQ/1e6, abs(exact[gid]), '--', label=labels[gid]+' 解析参考')
        axes[1, 0].plot(FREQ/1e6, np.angle(products[gid]/exact[gid])*180/np.pi, label=labels[gid])
    axes[1, 1].plot(FREQ/1e6, abs(fdtd_reflection), label='FDTD差分(半空间−自由空间)')
    axes[1, 1].plot(FREQ/1e6, abs(reflected), '--', label='Sommerfeld反射积分')
    axes[0, 0].set(title='单道总场 A-scan；100MHz Ricker，Iℓ峰值1A·m', xlabel='时间(ns)', ylabel='Ex(V/m)')
    axes[0, 1].set(title='复频响模值：Ex/(Iℓ)，单位(V/m)/(A·m)', xlabel='频率(MHz)')
    axes[1, 0].set(title='相对解析参考的相位误差', xlabel='频率(MHz)', ylabel='度')
    axes[1, 1].set(title='地层反射场：FDTD差分 vs 独立谱积分', xlabel='频率(MHz)')
    for ax in axes.flat:
        ax.grid(alpha=.2); ax.legend(fontsize=7)
    fig.suptitle('第二轮：8m航高 x 偶极在半空间上的点源验证；单站单道，非完整测线')
    fig.savefig(a.out/'halfspace_comparison.png', dpi=140)
    plt.close(fig)
    summary = {'status': 'PASS' if all(m['passes_declared_diagnostic'] for m in metrics) else 'DIAGNOSTIC_FAILED',
               'contract_sha256': sha256(a.study/'execution_contract.json'),
               'analysis_code_sha256': sha256(__file__),
               'reference_code_sha256': ref['code_sha256'], 'reference_selftest': ref['status'],
               'metrics': metrics, 'reflection_isolation': cross,
               'scope': 'Half-space continuum agreement within declared finite-grid diagnostic limits; not antenna, aperture or field validation.'}
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
