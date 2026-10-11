"""Stage A antenna batch analysis: 3-tone patterns, port impedance, ground effect.

Inputs: retrieved execution dir of the Stage A batch (A1/A1b free space,
A2 flat cover ground pair). All quantities from the official NTFF/port
outputs inside the native profile.h5 files; no new post-processing physics.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256 as sha
from line9_v401_version_controls import save

TONES = [20e6, 95e6, 170e6]
FF = 'ntff/stage_a_surface/frequency/stage_a_pattern/far_field/full_sphere'


def load_case(exec_dir, cid):
    p = exec_dir/'prepared/cases'/cid/'profile.h5'
    out = {'sha256': sha(p)}
    with h5py.File(p) as h:
        out['freqs'] = h['ntff/stage_a_surface/frequency/stage_a_pattern/frequencies'][:] if 'ntff' in h else None
        if 'ntff' in h:
            g = h[FF]
            out.update(theta=g['theta'][:], phi=g['phi'][:],
                       gain=g['fields/gain_dbi'][:], directivity=g['fields/directivity_dbi'][:],
                       maxdir=g['maximum_directivity_dbi'][:])
        port = h['ports/feed']
        out['port'] = {k: port[k][:] for k in
                       ['frequency', 'Zin', 'S11', 'Vtotal', 'time', 'valid_Zin', 'valid_S11']}
    return out


def at_tones(freqs, values):
    idx = [int(np.argmin(np.abs(freqs-f))) for f in TONES]
    return idx, values[idx]


def cut(gain, theta, phi, phi_deg):
    m = np.isclose(phi, phi_deg)
    order = np.argsort(theta[m])
    return theta[m][order], gain[:, m][:, order]


def main(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    cases = {cid: load_case(a.exec_dir, cid) for cid in ['A1_free_y', 'A1b_free_x', 'A2_ground_pair']}
    summary = {}
    for cid, c in cases.items():
        pf, pz = c['port']['frequency'], c['port']['Zin']
        ps11 = c['port']['S11']
        idx, z = at_tones(pf, pz)
        _, s11 = at_tones(pf, np.abs(ps11))
        row = dict(raw_sha256=c['sha256'],
                   Zin_at_tones={f'{f/1e6:g}MHz': [float(zi.real), float(zi.imag)] for f, zi in zip(TONES, z)},
                   S11_abs_at_tones={f'{f/1e6:g}MHz': float(s) for f, s in zip(TONES, s11)})
        if c.get('gain') is not None:
            imax = c['gain'].argmax(axis=1)
            row['max_gain_dbi'] = {f'{f/1e6:g}MHz': float(c['gain'][k].max()) for k, f in enumerate(TONES)}
            row['max_gain_direction_deg'] = {f'{f/1e6:g}MHz': [float(c['theta'][i]), float(c['phi'][i])]
                                             for f, i in zip(TONES, imax)}
            row['nadir_gain_dbi'] = {f'{f/1e6:g}MHz': float(c['gain'][k][np.isclose(c['theta'], 180)].max())
                                     for k, f in enumerate(TONES)}
            row['zenith_gain_dbi'] = {f'{f/1e6:g}MHz': float(c['gain'][k][np.isclose(c['theta'], 0)].max())
                                      for k, f in enumerate(TONES)}
        summary[cid] = row

    # A2: passive Rx gap record
    with h5py.File(a.exec_dir/'prepared/cases/A2_ground_pair/profile.h5') as h:
        gap = h['rxs/rx1/Ey'][:]  # rxs/rx1 = first #rx in file = passive gap probe
        dt = float(h.attrs['dt'])
    gap_spec = np.fft.rfft(gap)
    gap_freqs = np.fft.rfftfreq(len(gap), dt)
    gi = [int(np.argmin(np.abs(gap_freqs-f))) for f in TONES]
    summary['A2_ground_pair']['rxgap_tone_amplitude'] = {f'{f/1e6:g}MHz': float(abs(gap_spec[i]))
                                                         for f, i in zip(TONES, gi)}
    summary['A2_ground_pair']['rxgap_peak_Vpm'] = float(np.abs(gap).max())

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    fig, axs = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
    for k, f in enumerate(TONES):
        ax = axs[0, k]
        for cid, style, lbl in [('A1_free_y', '-', 'A1 臂沿y(正交基线)'), ('A1b_free_x', '--', 'A1b 臂沿x')]:
            c = cases[cid]
            t, g = cut(c['gain'], c['theta'], c['phi'], 90 if cid == 'A1_free_y' else 0)
            ax.plot(t, g[k], style, label=lbl+' E面')
            t, g = cut(c['gain'], c['theta'], c['phi'], 0 if cid == 'A1_free_y' else 90)
            ax.plot(t, g[k], ':' if cid == 'A1_free_y' else '-.', label=lbl+' H面')
        ax.axvline(180, color='gray', lw=.5)
        ax.set(title=f'{f/1e6:g} MHz 远场增益（官方 NTFF）', xlabel='theta / °（180=正下方）',
               ylabel='gain / dBi', ylim=[-30, 5])
        ax.grid(alpha=.3); ax.legend(fontsize=7)
        ax.invert_xaxis()
    for cid, style, lbl in [('A1_free_y', '-', 'A1 自由空间'), ('A2_ground_pair', '--', 'A2 平地覆盖层上5m')]:
        c = cases[cid]
        pf = c['port']['frequency']/1e6
        m = pf <= 400
        axs[1, 0].plot(pf[m], np.abs(c['port']['Zin'])[m], style, label=lbl)
        axs[1, 1].plot(pf[m], 20*np.log10(np.abs(c['port']['S11'])[m]+1e-12), style, label=lbl)
    t_ns = cases['A1_free_y']['port']['time']*1e9
    axs[1, 2].plot(t_ns, cases['A1_free_y']['port']['Vtotal'], '-', label='A1 馈口总电压', lw=.8)
    axs[1, 2].plot(t_ns, cases['A2_ground_pair']['port']['Vtotal'], '--', label='A2 馈口总电压', lw=.8)
    ax2 = axs[1, 2].twinx()
    ax2.plot(np.arange(len(gap))*dt*1e9, gap, color='tab:red', lw=.6, label='A2 被动Rx间隙Ey')
    ax2.set_ylabel('Rx 间隙 Ey / (V/m)', color='tab:red')
    axs[1, 2].set_ylabel('馈口电压 / V')
    axs[1, 0].set(title='输入阻抗 |Zin|', xlabel='MHz', ylabel='Ω'); axs[1, 0].legend(fontsize=8); axs[1, 0].grid(alpha=.3)
    axs[1, 1].set(title='|S11|（50Ω参考）', xlabel='MHz', ylabel='dB'); axs[1, 1].legend(fontsize=8); axs[1, 1].grid(alpha=.3)
    axs[1, 2].set(title='时域振铃与双天线耦合', xlabel='ns')
    axs[1, 2].grid(alpha=.3)
    fig.suptitle('阶段A 天线表征：等长中馈领结近似（全长3.35m，60°张角，无加载；未聚焦原天线=保守下限）\n'
                 '官方 NTFF 全场方向图 + 馈口端口量；网格 2.5cm，FP64，impulse 源')
    fig.savefig(a.out/'stage_a_patterns_ports.png', dpi=140)
    plt.close(fig)

    result = dict(status='COMPLETE_STAGE_A_CHARACTERIZATION', script_sha256=sha(__file__),
                  cases=summary, tones_hz=TONES,
                  method='official NTFF far-field arrays + voltage-source port outputs from native h5',
                  limits='Unfocused centre-fed bowtie approximation of the folded device antenna; '
                         'patterns are free-space (A1/A1b); A2 flat cover only. Near-field illumination '
                         'of the ground is characterized by the probe lines, not NTFF.')
    save(a.out/'analysis.json', result)
    print(json.dumps(result, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exec-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args())
