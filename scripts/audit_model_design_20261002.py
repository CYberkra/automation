"""Read-only design audit of archived controls; never runs FDTD or changes contracts.

Run with NumPy and h5py: python scripts/audit_model_design_20261002.py
The installed official reconstruction functions are executed verbatim via AST
extraction, avoiding this checkout's nonportable V4 interpreter. This is CPU
post-processing and synthetic representation verification, not a solver rerun.
"""
from __future__ import annotations

import ast
import argparse
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CONTROLS = ROOT / 'artifacts/research_checks/2026-10-02_benchmark3d_r2_acceptance'
OUT = ROOT / 'artifacts/research_checks/2026-10-02_model_design_audit.json'
OFFICIAL = ROOT / 'artifacts/local_checks/gprmax_v4_gpu_env/Lib/site-packages/gprMax/toolboxes/SFCW/processing.py'
F = np.linspace(20e6, 170e6, 501)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def official_reconstructor():
    source = OFFICIAL.read_text(encoding='utf-8')
    tree = ast.parse(source)
    names = {'spectral_window', 'reconstruct_time_response'}
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert len(nodes) == 2
    module = ast.Module(body=[ast.ImportFrom(module='__future__', names=[ast.alias(name='annotations')], level=0)] + nodes,
                        type_ignores=[])
    ns = {'np': np, 'TimeResponse': SimpleNamespace}
    exec(compile(ast.fix_missing_locations(module), str(OFFICIAL), 'exec'), ns)
    return ns['reconstruct_time_response']


def envelope(x):
    n = len(x)
    h = np.zeros(n)
    h[0] = 1
    h[1:(n + 1) // 2] = 2
    if n % 2 == 0:
        h[n // 2] = 1
    return np.abs(np.fft.ifft(np.fft.fft(x) * h))


def peaks(x, t):
    return {name: float(np.max(np.abs(x[(t >= lo) & (t < hi)])))
            for name, lo, hi in [('direct', 0, 50), ('surface', 85, 140),
                                 ('interface', 165, 215), ('late', 320, 370)]}


def main():
    global OFFICIAL
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--official-processing', type=Path, default=OFFICIAL,
                        help='Path to the installed V4 SFCW processing.py; read-only')
    args = parser.parse_args()
    OFFICIAL = args.official_processing
    reconstruct = official_reconstructor()
    doc = {'schema': 'model-design-audit/1', 'date': '2026-10-02',
           'review_base_commit': 'eb9d1fd', 'solver_runs': 0, 'training_runs': 0,
           'test_families_read': [], 'official_processing_sha256': sha(OFFICIAL),
           'limitations': ['only archived t07 controls available locally; full 13-trace delivered batch absent',
                           'new 501-tone products are audit diagnostics, not historical instrument calibration',
                           'B minus C changes the entire cover, including surface reflection; not a clean interface reference',
                           'B minus E changes material slab width AND PML position; not a PML-only control',
                           'raw and reconstructed amplitudes have different normalization; compare ratios within one representation']}
    paths = ['scripts/study_t3_damage_ladder.py', 'scripts/run_reward_protocol_b2_pilot_v0_1.py',
             'scripts/freeze_s1s3_reference_window_v0_1.py', 'scripts/freeze_benchmark3d_r2_co.py',
             'scripts/augment_sim2real_v0_6.py', 'configs/research/task_definition_contract_v0.1.json',
             'artifacts/research_checks/2026-10-02_benchmark3d_r2_acceptance/plot_sfcw_bscan.py']
    doc['reviewed_source_sha256'] = {p: sha(ROOT / p) for p in paths}
    synthetic = np.zeros(501, dtype=complex)
    synthetic[200] = 1.0  # actual requested tone: 80 MHz
    tr = reconstruct(SimpleNamespace(frequency=F, response=synthetic), window='rectangular', zero_pad_factor=8)
    legacy = np.abs(tr.complex_envelope) * np.cos(2 * np.pi * 95e6 * tr.time + np.angle(tr.complex_envelope))
    half_official = tr.real_bandpass / 2
    expected = np.cos(2 * np.pi * 80e6 * tr.time) / 501
    doc['signed_reconstruction'] = {
        'synthetic_tone_MHz': 80, 'legacy_tone_MHz': 155,
        'extra_frequency_shift_MHz': 75, 'legacy_band_MHz': [95, 245],
        'official_half_amplitude_relative_max_error': float(np.max(np.abs(half_official - expected)) / np.max(np.abs(expected))),
        'legacy_vs_official_half_relative_L2': float(np.linalg.norm(legacy - half_official) / np.linalg.norm(half_official)),
        'affected_scripts': ['study_t3_damage_ladder.py', 'run_reward_protocol_b2_pilot_v0_1.py',
                             'freeze_s1s3_reference_window_v0_1.py'],
        'note': 'amplitude factor two can be a declared convention; the 75 MHz frequency displacement cannot'}
    arrays, spectra, inputs = {}, {}, {}
    for key, name in [('A', 'asis'), ('B', 'flat9p00'), ('C2', 'rockonly'), ('E', 'flat_bigy')]:
        file = CONTROLS / f'ctl3d{key}_t07_{name}.h5'
        with h5py.File(file) as h:
            y = np.asarray(h['rxs/rx1/Ex'])
            dt = float(h.attrs['dt'])
            src = h['srcs/src1/excitation']
            x = np.asarray(src['samples'])
            so = float(src.attrs['TimeSampleOffset'])
            ro = float(h['rxs/rx1/Ex'].attrs['TimeSampleOffset'])
            assert y.dtype == np.dtype('float64') and np.isfinite(y).all()
            assert np.count_nonzero(x) == 1 and x[0] == 1
        t = np.arange(len(y)) * dt
        e = np.exp(-2j * np.pi * F[:, None] * t[None, :])
        # dt cancels between source and receiver. The actual half-step source
        # time origin is retained, matching engineering-convention source division.
        H = (e @ y) * np.exp(2j * np.pi * F * (so - ro))
        probes = [0, 1, 77, 200, 250, 333, 499, 500]
        independent = []
        for j in probes:
            angles = 2 * np.pi * F[j] * (t + ro - so)
            independent.append(complex(math.fsum((y * np.cos(angles)).tolist()),
                                       -math.fsum((y * np.sin(angles)).tolist())))
        dft_error = float(np.max(np.abs(H[probes] - independent)) / np.max(np.abs(H)))
        assert dft_error < 1e-9
        count = round(60e-9 / dt)
        taper = np.ones(len(y))
        taper[-count:] = .5 * (1 + np.cos(np.linspace(0, np.pi, count)))
        HT = (e @ (y * taper)) * np.exp(2j * np.pi * F * (so - ro))
        r = reconstruct(SimpleNamespace(frequency=F, response=H), window='hann', zero_pad_factor=8)
        rt = reconstruct(SimpleNamespace(frequency=F, response=HT), window='hann', zero_pad_factor=8)
        inverse_indices = [0, 120, 208, 230, 410]
        inverse = np.array([np.sum(H * r.weights * np.exp(2j * np.pi * (F - F[0]) * r.time[j])) / len(F)
                            for j in inverse_indices])
        inverse_error = float(np.max(np.abs(inverse - r.complex_envelope[inverse_indices])) / np.max(np.abs(r.complex_envelope)))
        assert inverse_error < 1e-10
        freq = np.fft.rfftfreq(len(y), dt)
        W = np.zeros_like(freq)
        mb = (freq >= F[0]) & (freq <= F[-1])
        W[mb] = np.hanning(int(mb.sum()))
        band = np.fft.irfft(np.fft.rfft(y) * W, n=len(y))
        arrays[key] = (y, band, r.complex_envelope, rt.complex_envelope)
        spectra[key] = H
        inputs[key] = {'path': str(file.relative_to(ROOT)).replace('\\', '/'), 'sha256': sha(file),
                       'input_sha256': sha(file.with_suffix('.in')), 'dtype': str(y.dtype),
                       'samples': len(y), 'dt_s': dt, 'source_time_offset_s': so,
                       'receiver_time_offset_s': ro, 'fft_spacing_MHz': 1 / (len(y) * dt) / 1e6,
                       'fft_band_bin_count': int(mb.sum()),
                       'independent_8point_DFT_relative_max_error': dft_error,
                       'independent_5point_inverse_relative_max_error': inverse_error,
                       'tail_last_10pct_energy_ratio': float(np.sum(y[int(.9 * len(y)):] ** 2) / np.sum(y ** 2)),
                       'tail_taper_60ns_spectral_relative_L2': float(np.linalg.norm(H - HT) / np.linalg.norm(H))}
        inputs[key]['tail_taper_60ns_hann_weighted_spectral_relative_L2'] = float(np.linalg.norm((H - HT) * r.weights) / np.linalg.norm(H * r.weights))
    doc['controls'] = inputs
    doc['representations'] = {}
    for ix, name in enumerate(['raw_impulse', 'fft_bandpass_approximation', 'exact_501_tones_hann', 'exact_501_tones_hann_60ns_taper']):
        a = {k: (envelope(v[ix]) if ix < 2 else np.abs(v[ix])) for k, v in arrays.items()}
        times = t * 1e9 if ix < 2 else r.time * 1e9
        values = {k: peaks(v, times) for k, v in a.items()}
        for label, k1, k2 in [('A_minus_B', 'A', 'B'), ('B_minus_C2', 'B', 'C2'), ('B_minus_E', 'B', 'E')]:
            diff = arrays[k1][ix] - arrays[k2][ix]
            values[label] = peaks(envelope(diff) if ix < 2 else diff, times)
        surface = values['B']['surface']
        doc['representations'][name] = {'peaks': values,
            'A_minus_B_interface_over_B_surface': values['A_minus_B']['interface'] / surface,
            'B_minus_C2_interface_over_B_surface': values['B_minus_C2']['interface'] / surface,
            'B_minus_E_interface_over_B_surface': values['B_minus_E']['interface'] / surface,
            'B_minus_E_over_A_minus_B_in_interface_window': values['B_minus_E']['interface'] / values['A_minus_B']['interface']}
    doc['materials'] = []
    c, eps0 = 299792458.0, 8.8541878128e-12
    for f in [20e6, 30e6, 40e6, 50e6, 95e6, 107e6, 170e6]:
        w = 2 * np.pi * f
        er = 18.017 + 7.878 / (1 + 1j * w * 6.4567e-9) - 1j * .003 / (w * eps0)
        k = w / c * np.sqrt(er)
        alpha = -k.imag
        doc['materials'].append({'frequency_MHz': f / 1e6, 'epsilon_real': float(er.real),
            'sigma_effective_S_m': float(-er.imag * w * eps0), 'attenuation_Np_m': float(alpha),
            'attenuation_dB_m': float(20 / np.log(10) * alpha), 'two_way_3m_amplitude': float(np.exp(-6 * alpha)),
            'phase_wavelength_m': float(2 * np.pi / k.real), 'cells_per_wavelength_5cm': float(2 * np.pi / k.real / .05)})
    from freeze_benchmark3d_r2_co import make_field, build_3d, build_2d, DX, PML2D_5CM
    eta = make_field()
    field_sha = sha(ROOT / 'configs/research/benchmark3d_r2_co/interface_field.npz')
    text3 = build_3d(eta, 6, field_sha)[1].decode()
    text2 = build_2d(eta, 6, DX, PML2D_5CM, 'B2D5CM', field_sha)[1].decode()
    def height_profile(text):
        z = []
        for line in text.splitlines():
            if line.startswith('#box:') and line.endswith(' cover'):
                s = line.split()
                if s[1] == 'inf' or (float(s[1]) <= 6 < float(s[4])):
                    z.append(float(s[3]))
        return np.array(z)
    z3, z2 = height_profile(text3), height_profile(text2)
    assert z3.shape == z2.shape
    delta = z3 - z2
    doc['geometry'] = {'material_extent': 'rock [1,11]x[1,11]x[1,12]; cover tiled inside same xy bounds; air elsewhere',
        'PML_inner_faces_m': [1, 11, 1, 11, 1, 32],
        'control_E_changes': 'y material extent [1,11] -> [1,15] AND domain 12 -> 16; Tx/Rx unchanged',
        '2D_3D_profile_mismatch_bins': int(np.count_nonzero(np.abs(delta) > 1e-10)),
        '2D_3D_profile_max_abs_m': float(np.max(np.abs(delta))),
        'profile_note': '3D averages 5x5 cells per tile; 2D averages only one x row, so dimension is not the sole changed factor',
        'kernel_sigma_m': 1.25, 'ideal_stationary_correlation_1e_length_m': 2.5,
        'kernel_note': 'white-noise Gaussian smoothing gives C(r)=exp(-r^2/(4 sigma^2)); finite reflected realization still needs empirical ACF audit',
        'single_realization_seed': 20261002, 'scan_length_m': 3.0,
        'scan_length_over_declared_CL': 1.2, 'interface_horizontal_tile_m': .25}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, ensure_ascii=False, sort_keys=True, indent=2) + '\n', encoding='utf-8', newline='\n')
    print('audit_sha256', sha(OUT))
    print('signed_reconstruction', doc['signed_reconstruction'])
    for name, row in doc['representations'].items():
        print(name, {k: v for k, v in row.items() if k != 'peaks'})
    print('geometry', doc['geometry'])


if __name__ == '__main__':
    main()
