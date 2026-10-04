"""Per-frequency air-cone diagnostic of archived 8 m Ey snapshots (CPU only).

Finite space/time windows give spectral leakage. |FFT(Ey)|^2 is a field
statistic, not directional energy flux or a reflected-power fraction.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hs4_analysis_metrics import air_propagation_mask
from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
CAP = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c'
C_AIR = 299792458.0
ROWS = {'cover_z11.975': 26, 'air_z16.025': 53, 'air_z19.925': 79}
TWIN = (90.0, 170.0)
FBAND = (60e6, 140e6)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--capsule', type=Path, default=CAP)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new analysis directory required; historical results must not be overwritten')
    cp = args.capsule/'execution_contract.json'
    c = json.loads(cp.read_text('utf-8'))
    v = json.loads((args.capsule/'completed_verification.json').read_text('utf-8'))
    if v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(cp):
        raise ValueError('completed capsule with matching contract required')
    groups = {g['id']: g for g in v['groups']}
    pair = [groups[k] for k in ('mid8_rough', 'mid8_halfspace')]
    iterations = np.asarray(c['snapshot_iterations'])
    for group in pair:
        files = group['snapshots']
        if [s['iteration'] for s in files] != iterations.tolist():
            raise ValueError('archived snapshot timeline differs from contract')
        missing = [s['file'] for s in files
                   if not (args.capsule/group['id']/s['file']).is_file()]
        if missing:
            raise FileNotFoundError(f"{group['id']}: {len(missing)} original frames missing; "
                                    'run on the machine retaining the frames; GIF is insufficient')
        if sha256(args.capsule/group['id']/'profile.h5') != group['raw_sha256']:
            raise ValueError('native receiver identity differs')
    t_ns = iterations * c['dt_s'] * 1e9
    selected = np.flatnonzero((t_ns >= TWIN[0]) & (t_ns <= TWIN[1]))
    if len(selected) < 3 or not np.all(np.diff(iterations[selected]) == np.diff(iterations[selected])[0]):
        raise ValueError('uniform selected time sampling required')
    gathers = {name: [] for name in ROWS}
    for index in selected:
        fields = []
        for group in pair:
            entry = group['snapshots'][index]
            path = args.capsule/group['id']/entry['file']
            if sha256(path) != entry['sha256']:
                raise ValueError(f'snapshot identity differs: {path}')
            with h5py.File(path) as h:
                field = h['Ey'][:]
                if (int(h.attrs['iteration']) != int(iterations[index])
                        or float(h.attrs['time']) != iterations[index] * c['dt_s']
                        or list(field.shape) != c['snapshot_shape_xyz']
                        or field.dtype != np.float64 or not np.isfinite(field).all()):
                    raise ValueError('snapshot metadata/shape/dtype/finite check failed')
                fields.append(field[:, 0, :])
        diff = fields[0] - fields[1]
        for name, iz in ROWS.items():
            gathers[name].append(diff[:, iz])
    spacing = np.asarray(c['snapshot_spacing_m'])
    origin = np.asarray(c['snapshot_extent_m'][:3])
    if not np.array_equal(spacing, [.15, .05, .15]) or not np.array_equal(origin, [13.5, 0., 8.]):
        raise ValueError('this diagnostic requires the reviewed ROI/row layout')
    ts = t_ns[selected] * 1e-9
    nx = c['snapshot_shape_xyz'][0]
    freqs = np.fft.fftshift(np.fft.fftfreq(len(ts), d=ts[1]-ts[0]))
    kxs = np.fft.fftshift(np.fft.fftfreq(nx, d=spacing[0])) * 2*np.pi
    fsel = (freqs >= FBAND[0]) & (freqs <= FBAND[1])
    cone = air_propagation_mask(freqs[fsel], kxs)
    taper = np.hanning(len(ts))[:, None] * np.hanning(nx)[None, :]
    spectra, sums, outside, arrays = {}, {}, {}, {}
    for name in ROWS:
        g = np.asarray(gathers[name])
        band = np.abs(np.fft.fftshift(np.fft.fft2(g*taper))[fsel])**2
        spectra[name] = np.log10(np.maximum(band.sum(axis=0), 1e-300))
        sums[name] = float(band.sum())
        if sums[name] == 0:
            raise ValueError('zero spectral-square denominator')
        outside[name] = float(band[~cone].sum() / sums[name])
        arrays[name + '_gather'] = g
        arrays[name + '_field_spectral_square'] = band
    ref = 'cover_z11.975'
    result = {
        'status': 'COMPLETED_FIELD_SPECTRUM_DIAGNOSTIC',
        'code_sha256': sha256(__file__), 'contract_sha256': sha256(cp),
        'metrics_code_sha256': sha256(ROOT/'scripts/hs4_analysis_metrics.py'),
        'selected_pairs_hash_verified': len(selected),
        'selected_frequency_hz': freqs[fsel].tolist(),
        'air_kx_limit_per_frequency_rad_per_m': (2*np.pi*freqs[fsel]/C_AIR).tolist(),
        'field_spectral_square_sum': sums,
        'field_spectral_square_relative_to_cover_dB':
            {k: float(10*np.log10(value/sums[ref])) for k, value in sums.items()},
        'field_spectral_square_fraction_outside_air_cone': outside,
        'window_ns': list(TWIN), 'fband_Hz': list(FBAND),
        'rows_z_m': {k: float(origin[2]+spacing[2]*(iz+.5)) for k, iz in ROWS.items()},
        'direction_separated': False, 'energy_flux_computed': False,
        'physical_attribution_certified': False,
        'limitations': ['Finite time/space Hann windows cause leakage and censor the field.',
                        'Same fixed time window samples different wave packets at different heights.',
                        'Ey spectral-square fractions are not incident/reflected power percentages.',
                        'No E/H directional decomposition or independent numerical-error budget.']}
    args.out.mkdir(parents=True)
    np.savez_compressed(args.out/'spectral_arrays.npz', time_ns=t_ns[selected],
                        frequency_hz=freqs[fsel], kx_rad_per_m=kxs, air_cone=cone, **arrays)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for name, p in spectra.items():
        ax.plot(kxs, p-p.max(), label=name)
    for f in (FBAND[0], 95e6, FBAND[1]):
        limit = 2*np.pi*f/C_AIR
        ax.axvline(limit, ls=':', label=f'air limit at {f/1e6:g} MHz')
        ax.axvline(-limit, ls=':')
    ax.set(xlabel='kx (rad/m)', ylabel='log10 band-integrated |Ey|² (relative to row max)',
           title='Finite-window field-spectrum diagnostic; fractions use per-frequency air limits')
    ax.grid(alpha=.3); ax.legend(); fig.tight_layout()
    fig.savefig(args.out/'surface_refraction_filter.png', dpi=110); plt.close(fig)
    (args.out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
