"""Layout/polarization/scenario analysis of the 12-trace 3D batch (design doc section B).

Total-field first (signed A-scans, declared event windows), then source-normalized
spectra. Layout effects are reported against this batch's own numerical floor
(pre-first-arrival cross-scene differences and the 1e-12 DFT identity), not
against a preset dB threshold.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from analyze_uav_local_free_space import FREQ, transfer
from hs_capsule_identity import sha256

C = 299792458.
# Declared windows (ns), from independent path predictions: direct 1.3 m air = 4.34 ns;
# surface wave 16 m air = 53.4 ns; interface ~3 m cover (n~4.3) two-way added ~84 ns.
WINDOWS = {'pre_arrival': (0., 4.0), 'direct': (3., 9.), 'surface': (46., 62.), 'interface': (120., 160.)}
PREDICTED_NS = {'direct': 1.3/C*1e9, 'surface': 16./C*1e9, 'interface': (16.+2*2.95*4.303)/C*1e9}
SCENARIOS = ('fullcover', 'flat', 'rough')
SC_LABEL = {'fullcover': '全覆盖层', 'flat': '平界面', 'rough': '0.8m起伏界面'}
LAY_LABEL = {'inline': '沿线布局', 'cross': '横线布局'}


def window_stats(time, wave):
    out = {}
    for name, (a, b) in WINDOWS.items():
        m = (time >= a) & (time < b)
        seg = wave[m]
        out[name] = {'peak': float(abs(seg).max()) if m.any() else 0.,
                     'rms': float(np.sqrt(np.mean(seg**2))) if m.any() else 0.}
    return out


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
    waves, spectra, metrics = {}, {}, []
    for g in c['groups']:
        raw = Path(g['input']).with_suffix('.h5')
        if sha256(raw) != identities[g['id']]:
            raise ValueError('raw changed')
        pol = g['polarisation']
        value, source, receiver = transfer(raw, 1, 'E'+pol, g['spacing_m'])
        # Independent seven-frequency direct sum with actual staggered clocks.
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
        time = receiver.dt*np.arange(len(receiver.samples))*1e9
        waves[g['id']] = (time, receiver.samples)
        spectra[g['id']] = value
        metrics.append({'case': g['id'], 'scenario': g['scenario'], 'layout': g['layout'],
                        'polarisation': pol, 'raw_sha256': sha256(raw),
                        'windows': window_stats(time, receiver.samples),
                        'copolar_spectrum_peak': float(abs(value).max()),
                        'independent_DFT_relative_L2': dft_error})
    # This-batch numerical floor: pre-first-arrival cross-scene difference per combo,
    # relative to the direct-window peak; plus pairwise same-combo cross-scene L2.
    floor_rows, layout_rows, pol_rows = [], [], []
    combos = sorted({(m['layout'], m['polarisation']) for m in metrics})
    for layout, pol in combos:
        ids = {s: f'{s}_{layout}_{pol}' for s in SCENARIOS}
        ref_direct = max(metrics_by_id(metrics, i)['windows']['direct']['peak'] for i in ids.values())
        for s1 in SCENARIOS:
            for s2 in SCENARIOS:
                if s1 >= s2:
                    continue
                t1, w1 = waves[ids[s1]]; t2, w2 = waves[ids[s2]]
                m = t1 < WINDOWS['pre_arrival'][1]
                floor_rows.append({'combo': f'{layout}_{pol}', 'pair': f'{s1}-{s2}',
                                   'pre_arrival_maxdiff_over_direct_peak':
                                   float(abs(w1[m]-w2[m]).max()/ref_direct)})
                full = float(np.linalg.norm(w1-w2)/np.linalg.norm(w1))
                post = t1 >= WINDOWS['surface'][0]
                floor_rows[-1]['cross_scene_total_relative_L2_post46ns'] = float(
                    np.linalg.norm((w1-w2)[post])/np.linalg.norm(w1[post]))
        for s in SCENARIOS:
            other = 'cross' if layout == 'inline' else 'inline'
            i1, i2 = f'{s}_{layout}_{pol}', f'{s}_{other}_{pol}'
            v1, v2 = spectra[i1], spectra[i2]
            t1, w1 = waves[i1]; t2, w2 = waves[i2]
            layout_rows.append({'scenario': s, 'polarisation': pol, 'pair': f'{layout}-{other}',
                                'copolar_spectrum_relative_L2': float(np.linalg.norm(v1-v2)/np.linalg.norm(v1)),
                                'ascan_relative_L2': float(np.linalg.norm(w1-w2)/np.linalg.norm(w1)),
                                'ascan_relative_L2_post46ns': float(
                                    np.linalg.norm((w1-w2)[t1 >= 46])/np.linalg.norm(w1[t1 >= 46]))})
    for s in SCENARIOS:
        for layout in ('inline', 'cross'):
            i1, i2 = f'{s}_{layout}_x', f'{s}_{layout}_y'
            t1, w1 = waves[i1]; t2, w2 = waves[i2]
            d1 = metrics_by_id(metrics, i1)['windows']['direct']['peak']
            d2 = metrics_by_id(metrics, i2)['windows']['direct']['peak']
            pol_rows.append({'scenario': s, 'layout': layout,
                             'direct_peak_x_over_y': float(d1/d2) if d2 else None,
                             'surface_peak_ratio_x_over_y': float(
                                 metrics_by_id(metrics, i1)['windows']['surface']['peak']/d2)})
    max_floor = max(r['pre_arrival_maxdiff_over_direct_peak'] for r in floor_rows)
    for row in layout_rows:
        row['exceeds_batch_floor'] = bool(row['ascan_relative_L2_post46ns'] > 10*max_floor)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    combo_titles = [('inline', 'x'), ('inline', 'y'), ('cross', 'x'), ('cross', 'y')]
    fig, axes = plt.subplots(3, 4, figsize=(16, 9), layout='constrained', sharex=True)
    allwaves = np.concatenate([w for _, w in waves.values()])
    ylim = abs(allwaves).max()*1.05
    for r, s in enumerate(SCENARIOS):
        for col, (layout, pol) in enumerate(combo_titles):
            time, wave = waves[f'{s}_{layout}_{pol}']
            ax = axes[r, col]
            ax.plot(time, wave, lw=.7)
            ax.set_ylim(-ylim, ylim)
            ax.set_xlim(0, 200)
            for name, (a0, b0) in WINDOWS.items():
                if name != 'pre_arrival':
                    ax.axvspan(a0, b0, alpha=.08, color={'direct': 'g', 'surface': 'orange', 'interface': 'r'}[name])
            ax.grid(alpha=.2)
            if r == 0:
                ax.set_title(f'{LAY_LABEL[layout]} + {pol}极化', fontsize=10)
            if col == 0:
                ax.set_ylabel(SC_LABEL[s]+'\nEx/Ey(V/m)', fontsize=9)
    fig.suptitle('12道三维布局×极化×场景批：单站总场 A-scan（共享纵轴；绿=直耦窗、橙=地表窗、红=界面窗预测位置；单站非完整B-scan）')
    fig.savefig(a.out/'layout_pol_ascans.png', dpi=130)
    plt.close(fig)
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.6), layout='constrained')
    for r, s in enumerate(SCENARIOS):
        ax = axes[r]
        for layout, pol in combo_titles:
            v = spectra[f'{s}_{layout}_{pol}']
            ax.plot(FREQ/1e6, abs(v), lw=.8, label=f'{LAY_LABEL[layout]}-{pol}')
        ax.set(title=SC_LABEL[s]+'：共极化复频响模值', xlabel='频率(MHz)')
        ax.grid(alpha=.2); ax.legend(fontsize=7)
    labels = [f"{r['scenario']}-{r['polarisation']}" for r in layout_rows if r['pair'] == 'inline-cross']
    vals = [r['ascan_relative_L2_post46ns'] for r in layout_rows if r['pair'] == 'inline-cross']
    axes[3].bar(range(len(vals)), vals)
    axes[3].axhline(max_floor*10, color='r', ls='--', label='10×本批到达前差分底')
    axes[3].axhline(max_floor, color='r', ls=':', label='本批到达前差分底')
    axes[3].set_xticks(range(len(vals)), labels, rotation=45, fontsize=7)
    axes[3].set(title='沿线vs横线：A-scan相对L2(46ns后)', yscale='log')
    axes[3].grid(alpha=.2); axes[3].legend(fontsize=7)
    fig.suptitle('布局效应对照：同场景同极化仅改收发基线方向（数值底来自本批到达前跨场景差分）')
    fig.savefig(a.out/'layout_pol_spectra.png', dpi=130)
    plt.close(fig)
    summary = {'status': 'PASS',
               'contract_sha256': sha256(a.study/'execution_contract.json'),
               'analysis_code_sha256': sha256(__file__),
               'declared_windows_ns': WINDOWS, 'predicted_arrivals_ns': PREDICTED_NS,
               'metrics': metrics, 'batch_floor': {'max_pre_arrival_cross_scene_over_direct': max_floor,
                                                   'detail': floor_rows},
               'layout_effect': layout_rows, 'polarization_observation': pol_rows,
               'scope': 'Single centre station, ideal dipole; layout factors from 3D internal comparison only; not antenna, aperture or field validation.'}
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'max_floor': max_floor,
                      'layout_effect': layout_rows}, ensure_ascii=False, indent=2))


def metrics_by_id(metrics, case):
    return next(m for m in metrics if m['case'] == case)


if __name__ == '__main__':
    main()
