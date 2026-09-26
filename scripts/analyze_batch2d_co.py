"""Post-run CPU audit + diagnostics of the batch2d_v1_co common-offset batch.

For each of the 11 trace positions (t01..t11, Rx y = 12.65 + 1.0*k, Tx = Rx - 1.30):
  * layout: single rx group, Position/GridPosition match the frozen trace geometry, Ex float64 20352
  * anchor t05 bit-for-bit vs archived single-trace mother h5
  * pairing diagnostics on diff = TGT - BG (same trace position):
      - in_band_energy_ratio_DIFF_over_BG on the 501-pt analysis grid 20-170 MHz
        (diagnostic energy ratio; NOT detectability/SNR/threshold)
      - envelope peak time within 250-450 ns
      - windowed energy ratio on 300-380 ns (cavity arrival) and 400-1200 ns (tail)
  * translation-invariance probes: pairwise bit-equality within BG set, within TGT set,
    and within diff set (flat-layered medium is y-invariant except the target box y14-18)

Deterministic output (no timestamps, sort_keys); two runs must be byte-identical.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
N_TRACES = 11
N_SAMPLES = 20352
DT = 5.896635841874211e-11
ANCHOR_K = 5  # 1-based; t05 == mother geometry
RX_Y = [12.65 + k for k in range(N_TRACES)]
BG = 'B2D-C3m-BG-CO11-t{:02d}'
TG = 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02-CO11-t{:02d}'
ARCHIVE = 'artifacts/research_checks/2026-09-26_{}/{}.h5'
FREQS = np.linspace(20e6, 170e6, 501)


def load(rid):
    with h5py.File(ROOT / 'artifacts/research_checks' / f'2026-09-26_{rid}' / f'{rid}.h5', 'r') as f:
        g = f['rxs/rx1']
        pos = np.asarray(g.attrs.get('Position'), dtype=float)
        return g['Ex'][:].astype(np.float64), pos, f['srcs/src1'].attrs['Position']


def layout_ok(pos, src_pos, k):
    issues = []
    if abs(pos[1] - RX_Y[k]) > 1e-9 or abs(pos[2] - 45.0) > 1e-12:
        issues.append(f'rx Position y={pos[1]} expected {RX_Y[k]}')
    if abs(src_pos[1] - (RX_Y[k] - 1.30)) > 1e-9:
        issues.append(f'src y={src_pos[1]} expected {RX_Y[k] - 1.30}')
    return issues


def band_energy(x):
    spec = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(N_SAMPLES, DT)
    idx = np.argmin(np.abs(freqs[:, None] - FREQS[None, :]), axis=0)
    return float(np.sum(np.abs(spec[idx]) ** 2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path,
                   default=ROOT / 'artifacts/research_checks/2026-09-26_batch2d_v1_co_analysis/results.json')
    a = p.parse_args()

    bg, tg, diffs = [], [], []
    rows = []
    for k in range(N_TRACES):
        rid_b, rid_t = BG.format(k + 1), TG.format(k + 1)
        b, pos_b, src_b = load(rid_b)
        t, pos_t, src_t = load(rid_t)
        issues = layout_ok(pos_b, src_b, k) + layout_ok(pos_t, src_t, k)
        d = t - b
        bg.append(b); tg.append(t); diffs.append(d)
        t_ns = np.arange(N_SAMPLES) * DT * 1e9
        win_arr = (t_ns >= 250) & (t_ns < 450)
        env_peak = float(t_ns[win_arr][np.argmax(np.abs(d)[win_arr])])
        w_cav = (t_ns >= 300) & (t_ns < 380)
        w_tail = (t_ns >= 400) & (t_ns < 1200)
        eb = band_energy(b)
        rows.append({
            'trace': f't{k + 1:02d}', 'rx_y_m': RX_Y[k],
            'over_target_y14_18': 14.0 <= RX_Y[k] <= 18.0,
            'in_band_ratio_DIFF_over_BG': band_energy(d) / eb,
            'window_energy_ratio_300_380ns': float(np.sum(d[w_cav] ** 2)) / float(np.sum(b[w_cav] ** 2)),
            'window_energy_ratio_400_1200ns': float(np.sum(d[w_tail] ** 2)) / float(np.sum(b[w_tail] ** 2)),
            'envelope_peak_ns_250_450': env_peak,
            'layout_issues': issues,
        })

    # anchor bit-for-bit vs archived single-trace mothers
    anchor = {}
    for rid, mother in ((BG.format(ANCHOR_K), 'B2D-C3m-BG'),
                        (TG.format(ANCHOR_K), 'B2D-C3m-D10m-W4m-T0.5m-E20-S0.02')):
        with h5py.File(ROOT / ARCHIVE.format(mother, mother), 'r') as f:
            ref = f['rxs/rx1/Ex'][:]
        cur = bg[ANCHOR_K - 1] if 'BG' in rid else tg[ANCHOR_K - 1]
        anchor[rid] = bool(np.array_equal(cur, ref))

    # translation-invariance probes (bit-equality of same-role traces across positions)
    def bit_eq_matrix(set_):
        m = np.zeros((N_TRACES, N_TRACES), dtype=bool)
        for i in range(N_TRACES):
            for j in range(N_TRACES):
                m[i, j] = np.array_equal(set_[i], set_[j])
        return m
    m_bg = bit_eq_matrix(bg)
    m_tg = bit_eq_matrix(tg)
    m_df = bit_eq_matrix(diffs)

    over = [r['in_band_ratio_DIFF_over_BG'] for r in rows if r['over_target_y14_18']]
    outside = [r['in_band_ratio_DIFF_over_BG'] for r in rows if not r['over_target_y14_18']]
    result = {
        'batch': 'batch2d_v1_co',
        'n_traces': N_TRACES,
        'anchor_t05_bit_for_bit': anchor,
        'all_layout_ok': not any(r['layout_issues'] for r in rows),
        'translation_invariance': {
            'bg_traces_all_bit_equal': bool(m_bg.all()),
            'tgt_traces_all_bit_equal': bool(m_tg.all()),
            'diff_traces_all_bit_equal': bool(m_df.all()),
            'bg_bit_equal_fraction': float(m_bg.sum()) / N_TRACES ** 2,
            'tgt_bit_equal_fraction': float(m_tg.sum()) / N_TRACES ** 2,
            'diff_bit_equal_fraction': float(m_df.sum()) / N_TRACES ** 2,
            'note': 'flat-layered medium is y-invariant except target box y14-18; bit-equality '
                    'across positions probes whether observed variation is arithmetic-level or physical',
        },
        'diagnostic_summary': {
            'in_band_ratio_over_target_mean': float(np.mean(over)),
            'in_band_ratio_outside_mean': float(np.mean(outside)),
            'envelope_peak_ns_range': [min(r['envelope_peak_ns_250_450'] for r in rows),
                                       max(r['envelope_peak_ns_250_450'] for r in rows)],
            'metric_disclaimer': 'in_band_energy_ratio is a diagnostic energy ratio on the '
                                 '501-point analysis grid; not detectability, not SNR, not a '
                                 'physical threshold',
        },
        'traces': rows,
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'all_layout_ok': result['all_layout_ok'],
                      'anchor': anchor,
                      'bg_all_equal': result['translation_invariance']['bg_traces_all_bit_equal'],
                      'tgt_all_equal': result['translation_invariance']['tgt_traces_all_bit_equal'],
                      'diff_all_equal': result['translation_invariance']['diff_traces_all_bit_equal'],
                      'over_mean': result['diagnostic_summary']['in_band_ratio_over_target_mean'],
                      'outside_mean': result['diagnostic_summary']['in_band_ratio_outside_mean'],
                      'output': str(a.output)}, sort_keys=True))


if __name__ == '__main__':
    main()
