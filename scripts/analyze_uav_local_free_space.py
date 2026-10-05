"""Independent continuous dipole comparison of actual-clock native V4 outputs."""
import argparse
import json
from pathlib import Path
import numpy as np
import h5py
from hs_capsule_identity import sha256
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response

C = 299792458.
MU = 4*np.pi*1e-7
EPS = 1/(MU*C*C)
FREQ = np.linspace(20e6, 170e6, 501)


def dipole_field(displacement, polarisation):
    r = np.linalg.norm(displacement)
    n = np.asarray(displacement)/r
    p = np.eye(3)['xyz'.index(polarisation)]
    k = 2*np.pi*FREQ/C
    transverse = p-n*np.dot(n, p)
    quasistatic = 3*n*np.dot(n, p)-p
    return (np.exp(-1j*k*r)/(4*np.pi*EPS*2j*np.pi*FREQ))[:, None]*(
        k[:, None]**2/r*transverse+(1/r**3+1j*k[:, None]/r**2)*quasistatic)


def transfer(raw, index, component, dl):
    s = load_source(raw)
    r = load_receiver(raw, receiver_path=f'/rxs/rx{index}', component=component)
    fraction = (round(20e-9/r.dt)-.25)/len(r.samples)
    result = direct_frequency_response(s, r, FREQ, tail_taper_fraction=fraction)
    if not result.source_valid.all():
        raise ValueError('invalid source spectrum')
    # Hertzian current element length changes with mesh: report E/(I*dl).
    return result.response/dl, s, r


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('fresh analysis output required')
    c = json.loads((a.study/'execution_contract.json').read_text('utf-8'))
    audit = json.loads((a.study/'completed_verification.json').read_text('utf-8'))
    if audit['status'] != 'PASS' or audit['contract_sha256'] != sha256(a.study/'execution_contract.json'):
        raise ValueError('completed raw identity required')
    identities = {r['id']: r['raw_sha256'] for r in audit['groups']}
    metrics, products, native = [], {}, {}
    for g in c['groups']:
        raw = Path(g['input']).with_suffix('.h5')
        if sha256(raw) != identities[g['id']]:
            raise ValueError('raw changed')
        dl = g['spacing_m']; pol = g['polarisation']
        for i, position in enumerate(g['receivers_m'], 1):
            value, source, receiver = transfer(raw, i, 'E'+pol, dl)
            exact = dipole_field(np.array(position)-np.array(g['tx_m']), pol)[:, 'xyz'.index(pol)]
            phase = np.angle(value/exact)*180/np.pi
            error = float(np.linalg.norm(value-exact)/np.linalg.norm(exact))
            # Independent seven-frequency direct sum, with actual staggered clocks.
            take = np.array([0, 83, 167, 250, 333, 417, 500])
            sy = source.samples.copy(); ry = receiver.samples.copy()
            n = round(20e-9/receiver.dt)
            ry[-n:] *= .5*(1+np.cos(np.linspace(0, np.pi, n)))
            sft = np.exp(-2j*np.pi*FREQ[take, None]*(source.time_offset+source.dt*np.arange(len(sy))))@sy
            rft = np.exp(-2j*np.pi*FREQ[take, None]*(receiver.time_offset+receiver.dt*np.arange(len(ry))))@ry
            independent = rft/sft/dl
            dft_error = float(np.linalg.norm(value[take]-independent)/np.linalg.norm(independent))
            if dft_error > 1e-9:
                raise ValueError('independent direct DFT mismatch')
            row = {'case': g['id'], 'receiver': i, 'raw_sha256': sha256(raw),
                   'complex_relative_L2_to_continuum': error,
                   'maximum_phase_error_deg': float(abs(phase).max()),
                   'magnitude_error_dB_minmax': [float((20*np.log10(abs(value/exact))).min()),
                                                float((20*np.log10(abs(value/exact))).max())],
                   'independent_DFT_relative_L2': dft_error,
                   'passes_declared_diagnostic': bool(error <= .05 and abs(phase).max() <= 5.)}
            metrics.append(row)
            products[g['id'], i] = (value, exact)
            native[g['id'], i] = (receiver.dt*np.arange(len(receiver.samples))*1e9, receiver.samples)
    symmetry = {}
    for dl in sorted({g['spacing_m'] for g in c['groups']}):
        x, y = f'd{dl:g}_x_forward', f'd{dl:g}_y_forward'
        if (x, 1) in products and (y, 2) in products:
            symmetry[f'd{dl:g}_axial'] = float(np.linalg.norm(products[x, 1][0]-products[y, 2][0])/np.linalg.norm(products[x, 1][0]))
            symmetry[f'd{dl:g}_broadside'] = float(np.linalg.norm(products[x, 2][0]-products[y, 1][0])/np.linalg.norm(products[x, 2][0]))
    for stem, i in [('reverse_x', 1), ('reverse_y', 2)]:
        fwd, rev = ('d0.05_y_forward', i), ('d0.05_y_'+stem, 1)
        if fwd in products and rev in products:
            symmetry[stem+'_reciprocity'] = float(np.linalg.norm(products[fwd][0]-products[rev][0])/np.linalg.norm(products[fwd][0]))
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    forward = next(g for g in c['groups'] if g['polarisation'] == 'y' and g['direction'] == 'forward')
    for i, label in [(1, '侧向：基线垂直偶极'), (2, '轴向：基线平行偶极')]:
        value, exact = products[forward['id'], i]
        axes[0, 0].plot(FREQ/1e6, abs(value), label=label+' FDTD')
        axes[0, 0].plot(FREQ/1e6, abs(exact), '--', label=label+' 精确解')
        axes[0, 1].plot(FREQ/1e6, np.angle(value/exact)*180/np.pi, label=label)
        time, wave = native[forward['id'], i]
        axes[1, 0].plot(time, wave, label=label)
    matrix = np.stack([native[forward['id'], i][1] for i in range(1, 5)], axis=1)
    axes[1, 1].imshow(matrix, aspect='auto', cmap='RdBu_r', vmin=-abs(matrix).max(), vmax=abs(matrix).max(),
                      extent=[.5, 4.5, time[-1], 0])
    axes[1, 1].set(xticks=[1, 2, 3, 4], xticklabels=['+x', '+y', '−x', '−y'],
                   title='四接收点道集；单站多方向，非完整测线', ylabel='接收时间(ns)')
    axes[0, 0].set(title='复频响模值：E/(Iℓ)，单位(V/m)/(A·m)', xlabel='频率(MHz)')
    axes[0, 1].set(title='相对精确解的相位误差', xlabel='频率(MHz)', ylabel='度')
    axes[1, 0].set(title='原始总场；100MHz Ricker，Iℓ峰值1A·m', xlabel='时间(ns)', ylabel='Ey(V/m)', xlim=(0, 70))
    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        ax.grid(alpha=.2); ax.legend(fontsize=7)
    fig.suptitle('第一轮：自由空间三维偶极自检；收发1.3m，未拟合幅度或相位')
    fig.savefig(a.out/'free_space_comparison.png', dpi=140); plt.close(fig)
    summary = {'status': 'PASS' if all(m['passes_declared_diagnostic'] for m in metrics) else 'DIAGNOSTIC_FAILED',
               'contract_sha256': sha256(a.study/'execution_contract.json'), 'analysis_code_sha256': sha256(__file__),
               'metrics': metrics, 'symmetry_reciprocity_relative_L2': symmetry,
               'scope': 'Continuum agreement within declared finite-grid diagnostic limits; not antenna, ground or field validation.'}
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
