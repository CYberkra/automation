"""Read-only numerical check of a cached controller ICZT module.

No instrument imports/commands, field fitting, solver, or training. The fixture
is an ideal frequency-flat delayed response, not an antenna calibration.
Private-source derivatives are restricted to ignored local_checks.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt

ROOT = Path(__file__).resolve().parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--controller-dir', type=Path, required=True)
    ap.add_argument('--out-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.out_dir.resolve()
    assert (ROOT / 'artifacts/local_checks').resolve() in out.parents
    out.mkdir(parents=True, exist_ok=True)
    repo = args.controller_dir.resolve()
    source = repo / 'src/lib/iczt_processor.py'
    spec = importlib.util.spec_from_file_location('controller_iczt_read_only', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    f = np.linspace(20e6, 170e6, 501)
    t = np.linspace(0, 700e-9, 501)
    delay = 14e-9
    h = np.exp(-2j * np.pi * f * delay)
    sos = butter(4, [.09, .125], btype='band', fs=1/1.4, output='sos')
    cases = []
    for window in ['rect', 'kaiser', 'kaiser8', 'hann', 'blackman']:
        x = module.iczt(h, f[0], f[-1], t, window=window)
        direct = module.iczt(h, f[0], f[-1], t, window=window, method='direct')
        rel = float(np.linalg.norm(x-direct)/np.linalg.norm(direct))
        assert rel < 1e-9
        # Independently reconstruct the sum at start, early peak, middle, end.
        win = module._get_window(window, len(h))
        for k in [0, 10, 250, 500]:
            scalar = sum(complex(a)*complex(np.exp(2j*np.pi*fr*t[k])) for a, fr in zip(h*win, f))
            assert abs(scalar-direct[k])/np.max(np.abs(direct)) < 1e-12
        real = x.real
        early = t*1e9 <= 60
        peak_idx = np.argmax(np.abs(real[early]))
        peak = abs(real[peak_idx])
        env = np.abs(hilbert(sosfiltfilt(sos, real/peak)))
        tr = (t-t[peak_idx])*1e9
        q = float(np.median(env[(tr >= 45) & (tr <= 70)]))
        rms_ratio = float(np.sqrt(np.mean(real[(tr >= 40) & (tr <= 140)]**2))/peak)
        cases.append({'window': window, 'bluestein_vs_direct_relative_l2': rel,
                      'q_90_125_mhz': q, 'fullband_tail_rms_over_early_peak': rms_ratio,
                      'complex_magnitude_min': float(np.abs(x).min()),
                      'real_negative_samples': int(np.sum(real < 0)),
                      'magnitude_negative_samples': int(np.sum(np.abs(x) < 0))})
    example = repo / 'src/从P9371B导出保存的一道A-Scan参考数据.csv'
    rows = list(csv.reader(example.open(encoding='utf-8-sig')))
    idx = next(i for i, row in enumerate(rows) if row and row[0]=='Time(s)')
    numeric = np.array([[float(row[0]),float(row[1])] for row in rows[idx+1:] if len(row)==2])
    assert numeric.shape == (501,2)
    assert np.allclose(numeric[:,0], np.linspace(0,900e-9,501), rtol=0, atol=1e-20)
    hashes = {name:sha(repo/name) for name in ['src/lib/iczt_processor.py','src/lib/workers.py',
              'src/lib/vna_controller.py','src/lib/trace_writers.py',str(example.relative_to(repo))]}
    result = {'controller_commit': subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
              'source_hashes': hashes, 'fixture': 'H(f)=exp(-j2pi f 14ns); 20-170 MHz; 501 tones; 0-700 ns, 501 time samples',
              'scope': 'mechanism check only, not historical Yingshan settings or antenna recovery',
              'df_hz':float(f[1]-f[0]), 'dt_ns':float((t[1]-t[0])*1e9),
              'cases':cases,
              'repository_reference_trace': {'date_from_header':'2026-01-10',
                   'channel_from_header': rows[idx][1], 'samples':int(len(numeric)),
                   'dt_ns':float((numeric[1,0]-numeric[0,0])*1e9),
                   'end_ns':float(numeric[-1,0]*1e9),
                   'negative_samples':int(np.sum(numeric[:,1]<0)),
                   'qualification':'unknown scene/settings; not calibration or Yingshan trace'},
              'solver_executed':False,'training_executed':False,'hardware_connected':False}
    target = out/'transform_checks.json'
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({'cases':cases, 'reference_trace':result['repository_reference_trace'],
                      'result_sha256':sha(target)},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
