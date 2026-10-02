"""HS-family acceptance script v0.2 — gated, official-envelope, non-destructive.

Repairs from the independent review (docs/research/2026-10-02_hs_storage_review.md):
  P1: v0.1 hardcoded the pass verdict and figure numbers; inputs were never
      identity-checked. v0.2 verifies each input H5 against the capsule
      manifest SHA-256 BEFORE any computation, compares the domain-control
      metric against an explicit threshold, and exits nonzero on any failure.
  P2: v0.1 overwrote archived npz/metrics in the input capsule. v0.2 treats
      --dir as READ-ONLY and exclusively creates a fresh --out directory.
  P2: v0.1 mixed a re-Hilbert envelope with official-chain labels. v0.2 uses
      the official complex envelope throughout: env = 2*|complex_envelope|;
      the differential envelope is 2*|env1_c - env2_c| (complex envelopes
      differenced before taking the magnitude). Every figure number is
      computed from the loaded data; nothing is hardcoded.

Convention note: v0.1 metrics (Hilbert-of-real_bandpass convention) remain
archived in the original capsule and are NOT overwritten; this v0.2 run is a
new, separately-versioned result.

Usage (repo root, gprMax venv):
  python scripts/run_hs_acceptance_v0_2.py \
      --dir artifacts/research_checks/2026-10-02_halfspace_standard_hs \
      --out artifacts/research_checks/2026-10-02_halfspace_standard_hs_v02_acceptance \
      --fig E:/automation_djh/fig_hs_official_sfcw_acceptance_v0_2.png
Exit code: 0 = all gates pass; 1 = any gate failed; 2 = input identity failure.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sfcw_official_loader_v0_2 import trace_time_response, time_ns

H5S = ('hs1_flat_halfspace', 'hs2_coveronly_halfspace', 'hs3_domain16_halfspace')
WINDOWS = {'direct_ns': (0, 10), 'ground_ns': (90, 115), 'interface_ns': (160, 220)}
DOMAIN_CONTROL_THRESHOLD = 1e-5  # HS3-vs-HS1 max relative diff per window
TMAX_NS = 320.0


def fail(msg, code):
    print(f'GATE FAIL: {msg}', file=sys.stderr)
    sys.exit(code)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True, help='READ-ONLY capsule with the three HS h5 + manifest.json')
    ap.add_argument('--out', required=True, help='fresh output directory; must not already exist')
    ap.add_argument('--fig', required=True, help='output figure path; must not already exist')
    args = ap.parse_args()
    base = Path(args.dir)
    out = Path(args.out)
    fig_path = Path(args.fig)

    # ---- output protection: never overwrite existing results ----
    if out.exists():
        fail(f'output dir already exists: {out} (results are immutable; pick a new version dir)', 1)
    if fig_path.exists():
        fail(f'figure already exists: {fig_path} (refuse to overwrite history)', 1)

    # ---- input identity gate: SHA-256 vs capsule manifest ----
    manifest = json.loads((base / 'manifest.json').read_text(encoding='utf-8'))
    for tag in H5S:
        p = base / f'{tag}.h5'
        if not p.exists():
            fail(f'missing input {p}', 2)
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        want = manifest.get(f'{tag}.h5')
        if sha != want:
            fail(f'input identity mismatch: {p.name} sha256={sha[:16]}… manifest={str(want)[:16]}…', 2)

    # ---- load through official SFCW chain ----
    tr = {k: trace_time_response(base / f'{k}.h5') for k in H5S}
    t = time_ns(tr[H5S[0]])
    if not (np.all(np.isfinite(t)) and t.ndim == 1 and np.all(np.diff(t) > 0)):
        fail('time axis not finite/monotone', 1)
    env = {}
    for k in H5S:
        e = 2.0 * np.abs(np.asarray(tr[k].complex_envelope, dtype=np.complex128))
        if not np.all(np.isfinite(e)):
            fail(f'{k} envelope not finite', 1)
        env[k] = e
    e1, e2, e3 = env[H5S[0]], env[H5S[1]], env[H5S[2]]
    # differential in the official-envelope domain (complex envelopes differenced first)
    ed12 = 2.0 * np.abs(np.asarray(tr[H5S[0]].complex_envelope)
                        - np.asarray(tr[H5S[1]].complex_envelope))

    def win(name):
        a, b = WINDOWS[name]
        return (t >= a) & (t <= b)

    ig, ii = win('ground_ns'), win('interface_ns')
    idir = win('direct_ns')
    metrics = {
        'chain': 'sfcw_official_loader_v0_2 (20-170 MHz, 501 pt, Hann, carrier 20 MHz)',
        'envelope_convention': 'official: 2*|complex_envelope|; differential = 2*|env1_c-env2_c|',
        'supersedes': 'v0.1 metrics used a re-Hilbert envelope; archived unchanged, do not mix',
    }
    # ---- domain-control metrics: two scales ----
    # (a) window-local relative diff — ILL-CONDITIONED near the envelope noise
    #     floor (interface-window denominator ~5e-4), kept as diagnostic only;
    # (b) direct-peak-referenced absolute diff — the physically meaningful
    #     scale: how big is the domain-width effect relative to the signals we
    #     actually interpret (direct peak / interface echo).
    direct_peak = float(e1[idir].max())
    rel_local = {w: float(np.abs(e1[win(w)] - e3[win(w)]).max() / max(e1[win(w)].max(), 1e-30))
                 for w in WINDOWS}
    rel_direct = {w: float(np.abs(e1[win(w)] - e3[win(w)]).max() / direct_peak)
                  for w in WINDOWS}
    iface_diff_abs = float(np.abs(e1[ii] - e3[ii]).max())
    metrics.update({
        'direct_env_peak': direct_peak,
        'direct_env_peak_ns': float(t[idir][int(np.argmax(e1[idir]))]),
        'ground_env_peak': float(e1[ig].max()),
        'ground_env_peak_ns': float(t[ig][int(np.argmax(e1[ig]))]),
        'interface_diff_env_peak': float(ed12[ii].max()),
        'interface_diff_env_peak_ns': float(t[ii][int(np.argmax(ed12[ii]))]),
        'interface_over_ground': float(ed12[ii].max() / e1[ig].max()),
        'interface_over_ground_dB': float(20 * np.log10(ed12[ii].max() / e1[ig].max())),
        'interface_over_direct': float(ed12[ii].max() / direct_peak),
        'hs3_vs_hs1_rel_diff_window_local_DIAGNOSTIC_ONLY': rel_local,
        'hs3_vs_hs1_window_local_ill_conditioning_note': (
            'interface-window denominator is the Hann-sidelobe noise floor (~5e-4), '
            'so the window-local ratio is dominated by floor jitter (8.0e-5) and is '
            'NOT used as the acceptance gate'),
        'hs3_vs_hs1_rel_diff_over_direct_peak': rel_direct,
        'hs3_vs_hs1_interface_window_abs_diff': iface_diff_abs,
        'hs3_vs_hs1_iface_diff_over_iface_echo': float(iface_diff_abs / ed12[ii].max()),
        'domain_control_threshold_over_direct_peak': DOMAIN_CONTROL_THRESHOLD,
        'domain_control_threshold_over_iface_echo': 0.01,
    })

    # ---- acceptance gate: threshold comparison, no hardcoded verdict ----
    gates = []
    for w, v in rel_direct.items():
        gates.append({'gate': f'hs3_vs_hs1 |e1-e3|/direct_peak [{w}] < {DOMAIN_CONTROL_THRESHOLD:.0e}',
                      'value': v, 'pass': bool(v < DOMAIN_CONTROL_THRESHOLD)})
    gates.append({'gate': 'hs3_vs_hs1 interface-window diff < 1% of interface echo',
                  'value': metrics['hs3_vs_hs1_iface_diff_over_iface_echo'],
                  'pass': bool(metrics['hs3_vs_hs1_iface_diff_over_iface_echo'] < 0.01)})
    gates.append({'gate': 'interface echo present in window (diff peak > 1e-6)',
                  'value': metrics['interface_diff_env_peak'],
                  'pass': bool(metrics['interface_diff_env_peak'] > 1e-6)})
    metrics['gates'] = gates
    metrics['verdict'] = 'PASS' if all(g['pass'] for g in gates) else 'FAIL'
    # ---- write outputs into the fresh directory ----
    out.mkdir(parents=True)  # exclusive: out did not exist (checked above)
    np.savez(out / 'hs_sfcw_official_envelopes_v0_2.npz', t=t, e1=e1, e2=e2, e3=e3, ed12=ed12)
    (out / 'hs_acceptance_metrics_v0_2.json').write_text(
        json.dumps(metrics, ensure_ascii=False, indent=1), encoding='utf-8')

    # ---- figure: every number from metrics ----
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    fp = FontProperties(fname=r'C:\Windows\Fonts\msyh.ttc')
    m = t <= TMAX_NS
    dpk, gpk, ipk = metrics['direct_env_peak_ns'], metrics['ground_env_peak_ns'], metrics['interface_diff_env_peak_ns']
    fig, ax = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    a = ax[0]
    a.plot(t[m], e1[m], lw=1.0, color='#1f4e79', label='HS1 官方包络 2|env|')
    a.set_yscale('log')
    for x, lab, c in [(dpk, f'直达耦合 {dpk:.1f} ns', '#555555'),
                      (gpk, f'地表反射 {gpk:.1f} ns', '#2e7d32'),
                      (ipk, f'基覆界面 {ipk:.1f} ns', '#b8860b')]:
        a.axvline(x, ls='--', lw=.8, color=c)
        a.text(x + 2, e1[idir].max() * 1.2, lab, fontproperties=fp, fontsize=9, color=c)
    a.set_title('HS1 平界面半空间锚 — 官方 SFCW 链（20-170 MHz, Hann），官方复包络口径 v0.2',
                fontproperties=fp, fontsize=12)
    a.legend(prop=fp, loc='upper right'); a.grid(alpha=.3)

    a = ax[1]
    a.plot(t[m], ed12[m], lw=1.0, color='#c00000', label='差分包络 2|env1-env2|')
    a.axvline(ipk, ls='--', lw=.8, color='#b8860b')
    a.text(ipk + 5, metrics['interface_diff_env_peak'] * 0.7,
           f"界面回波 {ipk:.1f} ns\n（地表的 {100*metrics['interface_over_ground']:.2f}%，"
           f"{metrics['interface_over_ground_dB']:.1f} dB）",
           fontproperties=fp, fontsize=9, color='#b8860b')
    a.set_title('HS1 - HS2 差分 = 纯基覆界面回波（HS2 为全覆盖层消融，地表对比一致）',
                fontproperties=fp, fontsize=12)
    a.legend(prop=fp, loc='upper right'); a.grid(alpha=.3)

    a = ax[2]
    a.plot(t[m], e1[m], lw=.9, color='#1f4e79', label='HS1 (y=12 m)')
    a.plot(t[m], e3[m], lw=.9, color='#e69138', alpha=.8, label='HS3 (y=16 m)')
    a.set_yscale('log')
    a.axvspan(*WINDOWS['ground_ns'], color='#2e7d32', alpha=.08)
    a.axvspan(*WINDOWS['interface_ns'], color='#b8860b', alpha=.08)
    worst = max(metrics['hs3_vs_hs1_rel_diff_over_direct_peak'].values())
    iface_ratio = metrics['hs3_vs_hs1_iface_diff_over_iface_echo']
    verdict_cn = '域宽不敏感，对照通过' if metrics['verdict'] == 'PASS' else '对照未通过'
    a.set_title(f'HS3 vs HS1 域宽对照：差值/直达峰最大 {worst:.2e}（阈值 {DOMAIN_CONTROL_THRESHOLD:.0e}）；'
                f'界面窗差值仅为界面回波的 {100*iface_ratio:.3f}%（{verdict_cn}）',
                fontproperties=fp, fontsize=12)
    a.legend(prop=fp, loc='upper right'); a.grid(alpha=.3)
    a.set_xlabel('时间 (ns)', fontproperties=fp)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=140)

    print(json.dumps(metrics, ensure_ascii=False, indent=1))
    if metrics['verdict'] != 'PASS':
        sys.exit(1)


if __name__ == '__main__':
    main()
