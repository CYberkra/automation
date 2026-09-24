"""Analytic hypothetical propagation budgets; no FDTD or measured material data.

Run: python scripts/study_layer_loss_budget.py --output <new JSON path>
Only Python standard library is required. Output is never overwritten.
"""
import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path

EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
C0 = 299792458.0


def propagation(f_hz, er, sigma):
    # E~exp(j*w*t-gamma*z); gamma^2=j*w*mu*(sigma+j*w*eps).
    w = 2 * math.pi * f_hz
    return cmath.sqrt(1j * w * MU0 * (sigma + 1j * w * EPS0 * er))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit('Refusing to overwrite existing evidence')
    rows = []
    for thickness in (1.0, 3.0, 5.0):
        for sigma in (0.001, 0.01, 0.05):
            for f_mhz in (20, 80, 170):
                cover = propagation(f_mhz * 1e6, 16, sigma)
                rock = propagation(f_mhz * 1e6, 9, 0.001)
                cover_db = 20 / math.log(10) * 2 * thickness * cover.real
                rock_db = 20 / math.log(10) * 2 * (20-thickness) * rock.real
                rows.append(dict(cover_m=thickness, cover_er=16,
                                 cover_sigma_S_m=sigma, rock_er=9,
                                 rock_sigma_S_m=0.001, f_MHz=f_mhz,
                                 cover_two_way_absorption_dB=cover_db,
                                 rock_two_way_absorption_dB=rock_db,
                                 total_two_way_absorption_dB=cover_db+rock_db))
    checks = {
        'lossless_alpha_zero': propagation(80e6, 9, 0).real == 0,
        'lossless_beta': math.isclose(propagation(80e6, 9, 0).imag,
                                     2*math.pi*80e6*3/C0, rel_tol=1e-8),
        'low_loss_limit': math.isclose(propagation(80e6, 9, 1e-7).real,
                                      1e-7/2*math.sqrt(MU0/(EPS0*9)), rel_tol=1e-8),
        'positive_loss_all_scenarios': all(r['total_two_way_absorption_dB'] > 0 for r in rows),
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    result = {
        'schema': 'analytic-layer-budget/1', 'date': '2026-09-24',
        'status': 'hypothetical sensitivity only; not fitted or measured material properties',
        'solver_executed': False, 'field_data_used': False,
        'model': 'normal ray, nonmagnetic homogeneous layers, constant er/sigma; absorption only',
        'excluded': ['interface coefficients', 'geometric spreading', 'antenna coupling',
                     'scattering', 'dispersion relaxation spectrum', 'noise and receiver dynamic range'],
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'checks': checks, 'rows': rows,
        'lossless_20m_two_way_ns': {str(er): 40*math.sqrt(er)/C0*1e9 for er in (4, 9, 16)},
        'ideal_resolution_scale_m': {str(er): C0/(2*150e6*math.sqrt(er)) for er in (4, 9, 16)},
        'ricker80_infinite_continuous_spectrum_relative_peak_dB': {
            str(f): 20*math.log10((f/80)**2*math.exp(1-(f/80)**2)) for f in (20, 80, 170)},
        'resolution_note': 'v/(2B) scale only; actual window, effective bandwidth and loss alter resolution',
        'ricker_note': 'ideal infinite continuous waveform only; finite sampled V4 source needs separate verification',
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False)
        stream.write('\n')
    print(json.dumps({'rows': len(rows), 'checks': checks, 'output': str(args.output)}))


if __name__ == '__main__':
    main()
