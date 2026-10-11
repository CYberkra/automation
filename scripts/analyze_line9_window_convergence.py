"""Window convergence crop-vs-full analysis for the 2026-10-12 4800 ns batch.

For the three worst-tail stations (X170.0/170.5/171.0) of the overnight line:
- q_full: official direct SFCW response from the new 4800 ns solves.
- q_crop: the same 4800 ns records truncated to the 2400 ns sample count.
- q_night: response from the original 2400 ns solves (hash-checked against the
  overnight delivery analysis).
Metrics isolate (a) pure truncation effect (q_crop vs q_full) and
(b) truncation equivalence to a native short run (q_night vs q_crop).
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np

import gprMax
from gprMax.toolboxes.SFCW import processing as sf
from hs_capsule_identity import sha256 as sha
from line9_v401_version_controls import save

FREQ = 20e6+np.arange(501)*300000.
CROP_SAMPLES = 40703
STATIONS = {170.0: '17000', 170.5: '17050', 171.0: '17100'}


def relative(x, y):
    return float(np.linalg.norm(x-y)/np.linalg.norm(y))


def response_2400(path, g):
    s = sf.load_source(path)
    r = sf.load_receiver(path, '/rxs/rx1', 'Ez')
    q_full = sf.direct_frequency_response(s, r, FREQ, tail_taper_fraction=0)
    rc = replace(r, samples=r.samples[:CROP_SAMPLES])
    sc = replace(s, samples=s.samples[:CROP_SAMPLES])
    q_crop = sf.direct_frequency_response(sc, rc, FREQ, tail_taper_fraction=0)
    assert q_full.source_valid.all() and q_crop.source_valid.all()
    return (replace(q_full, response=q_full.response/.025),
            replace(q_crop, response=q_crop.response/.025))


def band_errors(delta, full):
    bands = {'20_60': (FREQ >= 20e6) & (FREQ < 60e6),
             '60_120': (FREQ >= 60e6) & (FREQ < 120e6),
             '120_170': (FREQ >= 120e6) & (FREQ <= 170e6)}
    return {name: float(np.linalg.norm(delta[m])/np.linalg.norm(full[m])) for name, m in bands.items()}


def main(a):
    assert not a.out.exists() and gprMax.__version__ == '4.0.1'
    a.out.mkdir(parents=True)
    c = json.loads((a.window_exec/'execution_contract.json').read_text('utf-8'))
    v = json.loads((a.window_exec/'completed_verification.json').read_text('utf-8'))
    assert v['completed']
    night = json.loads((a.night_analysis/'analysis.json').read_text('utf-8'))
    night_sha = {g['id']: g['native_sha256'] for g in night['native_audits']}
    window_sha = {g['id']: g['native_sha256'] for g in v['groups']}
    rows = []
    spectra = {}
    for x, tag in STATIONS.items():
        wid, nid = f'window_x{tag}_H1', f'line_x{tag}_H1'
        wpath = a.window_exec/'prepared/cases'/wid/'profile.h5'
        npath = a.night_raw/f'{nid}.h5'
        assert sha(wpath) == window_sha[wid] and sha(npath) == night_sha[nid]
        with h5py.File(wpath) as h:
            assert h.attrs['Iterations'] == 81404
            tail_db = None
        q_full, q_crop = response_2400(wpath, None)
        # tail metric of the 4800 ns record (official definition from the toolbox)
        tail_db = sf.direct_frequency_response(sf.load_source(wpath), sf.load_receiver(wpath, '/rxs/rx1', 'Ez'),
                                               FREQ, tail_taper_fraction=0).receiver_tail_relative_db
        q_night, _ = response_2400(npath, None)
        d_cf = q_crop.response-q_full.response
        rows.append(dict(profile_x_m=x,
                         window_raw_sha256=window_sha[wid], night_raw_sha256=night_sha[nid],
                         window_tail_final_5pct_db=tail_db,
                         crop_vs_full_relative_L2=float(np.linalg.norm(d_cf)/np.linalg.norm(q_full.response)),
                         crop_vs_full_bands=band_errors(d_cf, q_full.response),
                         crop_vs_full_max_tone_relative=float(np.max(np.abs(d_cf)/np.abs(q_full.response))),
                         night_vs_crop_relative_L2=relative(q_night.response, q_crop.response)))
        spectra[x] = (q_full, q_crop, q_night)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    fig, axs = plt.subplots(2, 1, figsize=(11, 9), layout='constrained')
    for x, (q_full, q_crop, q_night) in spectra.items():
        rel = np.abs(q_crop.response-q_full.response)/np.abs(q_full.response)
        axs[0].semilogy(FREQ/1e6, rel, label=f'X{x:g} 截断vs全长')
    axs[0].set(xlabel='频率 / MHz', ylabel='单音调相对差异',
               title='4800 ns 全长 vs 截 2400 ns：501 音调复谱相对差异（官方 direct 链）')
    axs[0].legend(); axs[0].grid(alpha=.3)
    for x, (q_full, q_crop, q_night) in spectra.items():
        t = q_full.receiver.times*1e9
        env = np.abs(q_full.receiver.samples)
        env = env/env.max()
        axs[1].semilogy(t, env, label=f'X{x:g} 接收道包络(归一)')
    axs[1].axvspan(2400, 4800, color='tab:red', alpha=.08, label='2400 ns 以外的延长段')
    axs[1].axhline(10**(-60/20), color='k', ls=':', label='官方 −60 dB 建议')
    axs[1].set(xlabel='时间 / ns', ylabel='相对幅度', ylim=[1e-6, 2],
               title='4800 ns 原始接收道尾部衰减（Ez，fp64）')
    axs[1].legend(fontsize=8); axs[1].grid(alpha=.3)
    fig.suptitle('时窗收敛验证：夜间批尾部最差 3 道（X170.0/170.5/171.0），4800 ns 延长窗对照')
    fig.savefig(a.out/'window_convergence.png', dpi=140)
    plt.close(fig)

    result = dict(status='COMPLETE_WINDOW_CONVERGENCE_3STATIONS',
                  script_sha256=sha(__file__), official_processing_sha256=sha(sf.__file__),
                  window_contract_sha256=sha(a.window_exec/'execution_contract.json'),
                  crop_samples=CROP_SAMPLES, full_samples=81404, stations=rows,
                  method='official_direct crop-vs-full; q_night hash-checked against overnight delivery analysis',
                  limits='3 worst-tail stations only; not a line-wide convergence certificate.')
    save(a.out/'analysis.json', result)
    print(json.dumps(result, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['window_exec', 'night_raw', 'night_analysis', 'out']:
        p.add_argument('--'+key.replace('_', '-'), type=Path, required=True)
    main(p.parse_args())
