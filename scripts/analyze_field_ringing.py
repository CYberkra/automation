"""Field ringing characterization across all Yingshan measured lines (36 dBm).

Per line: q (90-125 MHz bandpass envelope, t0+45-70 / direct peak),
envelope decay profile, free damped-sine fit (F0, tau) per subsampled trace,
fixed 107 MHz/73 ns template R2, ringing-window spectrum, adjacent-trace
coherence, lateral texture. Aggregated comparison against the house V6 model.

Time axis convention: provisional dt = T/(N-1) (per prior audit caveat).
Output: artifacts_check/ringing_field/ringing_field_analysis.json
        fig_ringing_field_overview.png
"""
from pathlib import Path
import json

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import butter, hilbert, sosfiltfilt

DATA = Path(r'E:\automation_djh\real_data_yingshan\营山测线数据')
OUTDIR = Path(r'E:\automation_djh\artifacts_check\ringing_field')
FIG = Path(r'E:\automation_djh\fig_ringing_field_overview.png')
OUTDIR.mkdir(parents=True, exist_ok=True)

LINES = ['Line3', 'Line6', 'Line7', 'Line9', 'LineL1', 'LineX1']
WINS = [(25, 45), (45, 70), (70, 100), (100, 140), (140, 200), (200, 300)]
V6_F0, V6_TAU = 107.0, 73.0  # MHz, ns
FIT_LO, FIT_HI = 25, 220  # fit window rel t0


def load_csv(p):
    with p.open('r', encoding='utf-8-sig') as f:
        hdr = [f.readline() for _ in range(4)]
    vals = [float(s.split('=')[1].split(',')[0]) for s in hdr]
    ns, tw, nr, dx = int(vals[0]), vals[1], int(vals[2]), vals[3]
    raw = np.loadtxt(p, delimiter=',', skiprows=4, usecols=(3,), dtype=float)
    x = raw.reshape(nr, ns)
    dt = tw / (ns - 1)
    t = np.arange(ns) * dt
    return x, t, dt, dict(ns=ns, tw=tw, nr=nr, dx=dx)


def damped_sine(tt, amp, fmhz, tau, phi):
    return amp * np.exp(-tt / tau) * np.sin(2 * np.pi * fmhz * tt + phi)


def analyze_line(name):
    x, t, dt, hdr = load_csv(DATA / f'{name}origin(36).csv')
    i0 = np.argmax(np.abs(x[:, t <= 60]), axis=1)
    pk = np.abs(x[np.arange(len(x)), i0])
    valid = pk > 0
    x, i0, pk = x[valid], i0[valid], pk[valid]
    nrm = x / pk[:, None]
    sos = butter(4, [0.090, 0.125], btype='band', fs=1 / dt, output='sos')
    env = np.abs(hilbert(sosfiltfilt(sos, nrm, axis=1), axis=1))
    rel = t[None, :] - t[i0][:, None]

    prof = np.column_stack([
        np.nanmedian(np.where((rel >= lo) & (rel <= hi), env, np.nan), axis=1)
        for lo, hi in WINS])
    q = prof[:, 1]

    # fit window (rel t0)
    offs = np.arange(int(np.ceil(FIT_LO / dt)), int(np.floor(FIT_HI / dt)) + 1)
    seg = nrm[np.arange(len(x))[:, None], i0[:, None] + offs]
    tt = offs * dt
    # fixed template R2 (107 MHz, 73 ns)
    Bfix = np.column_stack([
        np.exp(-tt / V6_TAU) * np.sin(2 * np.pi * V6_F0 / 1000 * tt),
        np.exp(-tt / V6_TAU) * np.cos(2 * np.pi * V6_F0 / 1000 * tt)])
    cf = seg @ np.linalg.pinv(Bfix).T
    pred = cf @ Bfix.T
    den = np.sum((seg - seg.mean(axis=1, keepdims=True)) ** 2, axis=1)
    r2_fix = 1 - np.sum((seg - pred) ** 2, axis=1) / den
    # free fit on subsampled traces
    step = max(1, len(x) // 150)
    f0s, taus, r2_free, ok = [], [], [], 0
    for r in range(0, len(x), step):
        s = seg[r] - seg[r].mean()
        if np.sqrt((s ** 2).mean()) < 1e-4:
            continue
        try:
            popt, _ = curve_fit(
                damped_sine, tt, s,
                p0=[s.std(), 107.0, 73.0, 0.0],
                bounds=([0, 60, 5, -np.pi], [5 * s.std(), 160, 400, np.pi]),
                maxfev=1200)
            resid = s - damped_sine(tt, *popt)
            r2_free.append(1 - np.sum(resid ** 2) / np.sum(s ** 2))
            f0s.append(popt[1]); taus.append(popt[2]); ok += 1
        except Exception:
            continue
    # spectrum of ringing window
    ft = np.fft.rfft((seg - seg.mean(axis=1, keepdims=True)) * np.hanning(len(offs)),
                     n=2048, axis=1)
    freqs = np.fft.rfftfreq(2048, dt) * 1000
    power = np.abs(ft) ** 2
    band = (freqs >= 20) & (freqs <= 170)
    dom = freqs[band][np.argmax(power[:, band], axis=1)]
    frac_narrow = power[:, (freqs >= 90) & (freqs <= 125)].sum(axis=1) / power[:, band].sum(axis=1)
    # adjacent-trace correlation 30-100 ns
    s2 = x[:, (t >= 30) & (t <= 100)]
    s2 = s2 - s2.mean(axis=1, keepdims=True)
    dn = np.linalg.norm(s2[:-1], axis=1) * np.linalg.norm(s2[1:], axis=1)
    cc = np.sum(s2[:-1] * s2[1:], axis=1)[dn > 0] / dn[dn > 0]
    # lateral texture
    e = np.asarray([np.abs(hilbert(row[(z >= 40) & (z <= 100)])).mean()
                    for row, z in zip(nrm, rel)])
    lat = float(np.sqrt(np.mean(np.diff(e / e.mean()) ** 2)))
    # q along line (33-trace block medians)
    blocks = [float(np.median(q[j:j + 33])) for j in range(0, len(q) - 32, 33)]
    np.save(OUTDIR / f'{name}_q.npy', q)

    # median-trace system ring-down fit: geology cancels in the median,
    # the trace-invariant system ringing survives
    med_trace = np.median(nrm, axis=0)
    i0m = int(np.median(i0))
    offm = np.arange(int(np.ceil(22 / dt)), int(np.floor(300 / dt)) + 1)
    segm = med_trace[i0m + offm]
    tm_ = offm * dt
    sm = segm - segm.mean()
    try:
        popt, _ = curve_fit(damped_sine, tm_, sm,
                            p0=[sm.std(), 107.0, 73.0, 0.0],
                            bounds=([0, 60, 5, -np.pi],
                                    [5 * sm.std(), 160, 400, np.pi]),
                            maxfev=4000)
        resid = sm - damped_sine(tm_, *popt)
        med_fit = dict(f0_mhz=round(float(popt[1]), 1),
                       tau_ns=round(float(popt[2]), 1),
                       r2=round(float(1 - np.sum(resid ** 2) / np.sum(sm ** 2)), 3))
    except Exception:
        med_fit = None
    # envelope decay slope tau from profile medians (45-70 -> 200-300)
    pm = np.nanmedian(prof, axis=0)
    tau_slope = float((200 + 300) / 2 - (45 + 70) / 2) / float(
        np.log(pm[1] / pm[5])) if pm[1] > 0 and pm[5] > 0 else None

    return dict(
        name=name, n_traces=len(x), dt_ns=dt, header=hdr,
        t0_ns=np.percentile(t[i0], [10, 50, 90]).round(2).tolist(),
        q_dist={f'p{p}': float(np.percentile(q, p)) for p in (10, 25, 50, 75, 90)},
        q_blocks=blocks,
        profile_medians=pm.round(5).tolist(),
        tau_env_slope_ns=round(tau_slope, 1) if tau_slope else None,
        median_trace_fit=med_fit,
        r2_fixed_107_73={f'p{p}': float(np.percentile(r2_fix, p)) for p in (10, 50, 90)},
        free_fit=dict(n_ok=ok, n_attempted=len(range(0, len(x), step)),
                      f0_mhz=np.percentile(f0s, [25, 50, 75]).round(1).tolist(),
                      tau_ns=np.percentile(taus, [25, 50, 75]).round(1).tolist(),
                      r2=np.percentile(r2_free, [25, 50, 75]).round(3).tolist()),
        dom_freq_mhz=np.percentile(dom, [25, 50, 75]).round(1).tolist(),
        frac_90_125=np.percentile(frac_narrow, [25, 50, 75]).round(3).tolist(),
        adj_corr=np.percentile(cc, [25, 50, 75]).round(3).tolist(),
        lateral_diff_rms=round(lat, 4),
    )


def v6_synth_profile(dt):
    tt = np.arange(0, 320, dt)
    g = np.exp(-tt / V6_TAU) * np.sin(2 * np.pi * V6_F0 / 1000 * tt)
    sos = butter(4, [0.090, 0.125], btype='band', fs=1 / dt, output='sos')
    env = np.abs(hilbert(sosfiltfilt(sos, g)))
    return tt, env / env.max()


def main():
    recs = [analyze_line(n) for n in LINES]
    for r in recs:
        ff = r['free_fit']
        print(f"{r['name']}: n={r['n_traces']} q50={r['q_dist']['p50']:.3f} "
              f"R2fix50={r['r2_fixed_107_73']['p50']:.3f} "
              f"f0={ff['f0_mhz'][1]}MHz IQR[{ff['f0_mhz'][0]},{ff['f0_mhz'][2]}] "
              f"tau={ff['tau_ns'][1]}ns IQR[{ff['tau_ns'][0]},{ff['tau_ns'][2]}] "
              f"cc50={r['adj_corr'][1]} lat={r['lateral_diff_rms']}")
    (OUTDIR / 'ringing_field_analysis.json').write_text(
        json.dumps({'caveat': 'time axis provisional dt=T/(N-1); '
                              'ringing window mixes geology returns',
                    'v6_model': {'f0_mhz': V6_F0, 'tau_ns': V6_TAU},
                    'lines': recs}, ensure_ascii=False, indent=2),
        encoding='utf-8')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = dict(zip(LINES, ['#6f8093', '#c37426', '#2666ac', '#1a7a3c',
                              '#8c2d6e', '#b09116']))

    fig, axes = plt.subplots(2, 3, figsize=(19, 9.5), dpi=130)
    # (a) q distributions
    ax = axes[0, 0]
    qdata = [np.load(OUTDIR / f'{n}_q.npy') for n in LINES]
    bp = ax.boxplot(qdata, tick_labels=LINES, showfliers=False, notch=True)
    for patch, n in zip(bp['boxes'], LINES):
        patch.set(color=colors[n], lw=1.6)
    l9 = next(r for r in recs if r['name'] == 'Line9')
    ax.axhspan(l9['q_dist']['p25'], l9['q_dist']['p75'], color='#2666ac',
               alpha=.12, label='V6 target band (Line9 p25-p75)')
    ax.set_ylabel('q = 90-125 MHz env median (t0+45~70) / direct peak')
    ax.set_title('(a) ringing level q per line')
    ax.grid(axis='y', alpha=.25); ax.legend(fontsize=8)

    # (b) envelope decay profiles + V6 synth
    ax = axes[0, 1]
    centers = [(lo + hi) / 2 for lo, hi in WINS]
    for n in LINES:
        r = next(x for x in recs if x['name'] == n)
        ax.plot(centers, r['profile_medians'], 'o-', color=colors[n],
                label=n, lw=1.3, ms=4)
    tt, env6 = v6_synth_profile(recs[0]['dt_ns'])
    pm6 = [np.median(env6[(tt >= lo) & (tt <= hi)]) for lo, hi in WINS]
    ax.plot(centers, pm6, 'k--', lw=1.8, label='V6 synthetic (107 MHz, 73 ns)')
    ax.set_yscale('log')
    ax.set_xlabel('window center rel. t0 (ns)')
    ax.set_ylabel('median bandpass env / direct peak')
    ax.set_title('(b) ring-down shape')
    ax.grid(alpha=.25, which='both'); ax.legend(fontsize=8)

    # (c) median-trace system fit: F0 vs tau
    ax = axes[0, 2]
    for n in LINES:
        r = next(x for x in recs if x['name'] == n)
        mf = r['median_trace_fit']
        if mf:
            ax.plot(mf['f0_mhz'], mf['tau_ns'], 'o', color=colors[n], ms=10,
                    label=f"{n} (R2={mf['r2']})")
        if r.get('tau_env_slope_ns'):
            ax.annotate(f"{n} slope-tau={r['tau_env_slope_ns']:.0f}ns",
                        xy=(mf['f0_mhz'], mf['tau_ns']),
                        xytext=(mf['f0_mhz'] + 1.2, mf['tau_ns']),
                        fontsize=7, color=colors[n])
    ax.plot(V6_F0, V6_TAU, 'k*', ms=20, label='V6 model')
    ax.set_xlabel('ring-down F0 (MHz, median-trace fit)')
    ax.set_ylabel('tau (ns)')
    ax.set_title('(c) system ring-down (median trace, t0+22~300 ns)')
    ax.grid(alpha=.25); ax.legend(fontsize=7.5)

    # (d) template R2 fixed vs free
    ax = axes[1, 0]
    rf = [next(r for r in recs if r['name'] == n)['r2_fixed_107_73']['p50']
          for n in LINES]
    rw = [next(r for r in recs if r['name'] == n)['free_fit']['r2'][1]
          for n in LINES]
    xs = np.arange(6)
    ax.bar(xs - .18, rf, width=.36, label='fixed 107 MHz / 73 ns')
    ax.bar(xs + .18, rw, width=.36, label='free (F0, tau)')
    ax.set_xticks(xs, LINES); ax.set_ylim(0, 1)
    ax.set_ylabel('R2 (t0+25~220 ns)')
    ax.set_title('(d) how adequate is the V6 template?')
    ax.grid(axis='y', alpha=.25); ax.legend(fontsize=8)

    # (e) coherence/texture
    ax = axes[1, 1]
    cc = [next(r for r in recs if r['name'] == n)['adj_corr'][1] for n in LINES]
    lat = [next(r for r in recs if r['name'] == n)['lateral_diff_rms']
           for n in LINES]
    ax.plot(LINES, cc, 'o-', label='adjacent-trace corr (30-100 ns)')
    ax.set_ylim(0, 1)
    ax2 = ax.twinx()
    ax2.plot(LINES, lat, 's--', color='#c0392b',
             label='lateral texture diff_rms')
    ax.set_ylabel('adjacent corr'); ax2.set_ylabel('diff_rms', color='#c0392b')
    ax.set_title('(e) system coherence vs texture')
    ax.grid(alpha=.25)
    l1, lb1 = ax.get_legend_handles_labels()
    l2, lb2 = ax2.get_legend_handles_labels()
    ax.legend(l1 + l2, lb1 + lb2, fontsize=8, loc='upper right')

    # (f) q along line
    ax = axes[1, 2]
    for n in LINES:
        r = next(x for x in recs if x['name'] == n)
        b = r['q_blocks']
        ax.plot(np.linspace(0, 1, len(b)), b, '-', color=colors[n], lw=1.2,
                label=f'{n} ({r["n_traces"]} tr)')
    ax.set_xlabel('relative position along line')
    ax.set_ylabel('q block median (33 traces)')
    ax.set_title('(f) ringing stability along line')
    ax.grid(alpha=.25); ax.legend(fontsize=7)

    fig.suptitle('Yingshan field ringing overview: 6 lines / 36 dBm vs house V6'
                 ' model (107 MHz, 73 ns)  |  dt provisional T/(N-1)',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(FIG, dpi=140)
    print('saved:', FIG)
    print('saved:', OUTDIR / 'ringing_field_analysis.json')


if __name__ == '__main__':
    main()
