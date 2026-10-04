"""Matched-carrier and conductivity comparisons with identical edge policy."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hs_capsule_identity import sha256
from hs4_v4_material_controls import CENTRE
from analyze_hs4_v4_factor_controls import compare
from analyze_hs4_height_wavefield import response
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--capsule', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new analysis directory required')
    cp = args.capsule/'execution_contract.json'
    c = json.loads(cp.read_text('utf-8'))
    v = json.loads((args.capsule/'completed_verification.json').read_text('utf-8'))
    events = [json.loads(s) for s in (args.capsule/'execution.jsonl').read_text('utf-8').splitlines()]
    if (v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(cp)
            or events[-1]['status'] != 'COMPLETED' or events[-1]['traces'] != 4
            or events[-1]['verification_sha256'] != sha256(args.capsule/'completed_verification.json')):
        raise ValueError('complete verified material capsule required')
    verify_official_runtime()
    refs = {}
    for role in ('rough', 'halfspace'):
        path = CENTRE/f'centre_{role}'/'profile.h5'
        row = json.loads((path.parent/'audit.json').read_text('utf-8'))[0]
        if sha256(path) != row['sha256']:
            raise ValueError('reference identity differs')
        refs[role] = response(path)
    native = {}
    for row in v['groups']:
        p = args.capsule/row['id']/'profile.h5'
        if sha256(p) != row['raw_sha256']:
            raise ValueError('material output identity differs')
        r = response(p)
        if (not np.array_equal(r.frequency, refs['rough'].frequency)
                or r.source.quantity != refs['rough'].source.quantity
                or r.source.spatial_scale != refs['rough'].source.spatial_scale):
            raise ValueError('frequency/source semantics differ')
        native[row['id']] = r
    model = refs['rough']
    spectra = {'debye_n': refs['rough'].response-refs['halfspace'].response}
    for tag in ('matched95', 'low_sigma'):
        spectra[tag] = native[f'{tag}_centre_rough'].response-native[f'{tag}_centre_halfspace'].response
    products = {k: reconstruct_time_response(replace(model, response=array), window='hann', zero_pad_factor=8) for k, array in spectra.items()}
    time = products['debye_n'].time*1e9
    comparisons = {'debye_vs_matched95': compare(products['debye_n'], products['matched95'], time),
                   'matched95_vs_low_sigma': compare(products['matched95'], products['low_sigma'], time),
                   'debye_vs_low_sigma': compare(products['debye_n'], products['low_sigma'], time)}
    visibility = {}
    for tag, r in (('debye_n', refs['rough']), ('matched95', native['matched95_centre_rough']), ('low_sigma', native['low_sigma_centre_rough'])):
        total = reconstruct_time_response(r, window='hann', zero_pad_factor=8)
        mask = (time >= 160)&(time <= 220)
        ratio = float(abs(products[tag].complex_envelope[mask]).max()/abs(total.complex_envelope).max())
        visibility[tag] = {'interface_event_over_full_total_peak': ratio,
                           'interface_event_over_full_total_peak_dB': float(20*np.log10(ratio))}
    frequencies = model.frequency
    omega = 2*np.pi*frequencies
    original = 18.017+7.878/(1+1j*omega*6.4567e-9)+.003/(1j*omega*c['epsilon0_F_per_m'])
    matched = c['constant_eps_real']+c['constant_sigma_equivalent95_S_per_m']/(1j*omega*c['epsilon0_F_per_m'])
    low = c['constant_eps_real']+.003/(1j*omega*c['epsilon0_F_per_m'])
    result = {'status': 'COMPLETED_MATCHED_EDGE_POLICY_MATERIAL_DIAGNOSTIC',
              'code_sha256': sha256(__file__), 'contract_sha256': sha256(cp), 'completed_cases': 4,
              'constant_eps_real': c['constant_eps_real'], 'constant_sigma_equivalent95_S_per_m': c['constant_sigma_equivalent95_S_per_m'],
              'matched95_complex_permittivity_error': v['carrier_complex_permittivity_matching_error'],
              'fixed_window_comparisons': comparisons, 'visibility_diagnostic': visibility,
              'source_all501_valid': True, 'physical_attribution_certified': False,
              'limitations': ['Matching complex epsilon at one carrier is not full-band equivalence.',
                              'Changing conductivity alters propagation and reflection, not recoverable time gain alone.',
                              'One ideal line-source 2D centre station; no 3D/field/model-selection claim.']}
    args.out.mkdir(parents=True)
    arrays = {'time_ns': time, 'frequency_hz': frequencies, 'debye_epsilon_complex': original,
              'matched95_epsilon_complex': matched, 'low_sigma_epsilon_complex': low}
    for tag, p in products.items():
        arrays[tag+'_complex_envelope'] = p.complex_envelope
        arrays[tag+'_signed_waveform'] = p.real_bandpass
        arrays[tag+'_frequency_response'] = spectra[tag]
    for key, r in native.items():
        arrays[key+'_native_frequency_response'] = r.response
    np.savez_compressed(args.out/'comparison_arrays.npz', **arrays)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout='constrained')
    mask = (time >= 140)&(time <= 240)
    for tag, p in products.items():
        axes[0, 0].plot(time[mask], p.real_bandpass[mask], label=tag)
        axes[0, 1].plot(time[mask], 2*abs(p.complex_envelope[mask]), label=tag)
    axes[0, 0].set(xlabel='actual receiver time (ns)', ylabel='signed response (transfer units)')
    axes[0, 1].set(xlabel='actual receiver time (ns)', ylabel='2|complex envelope| (transfer units)')
    for name, epsilon in (('debye_n', original), ('matched95', matched), ('low_sigma', low)):
        axes[1, 0].plot(frequencies/1e6, epsilon.real, label=name)
        axes[1, 1].plot(frequencies/1e6, -epsilon.imag, label=name)
    axes[1, 0].set(xlabel='frequency (MHz)', ylabel='relative permittivity real part')
    axes[1, 1].set(xlabel='frequency (MHz)', ylabel='relative permittivity loss part')
    for ax in axes.flat:
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[1, 0].axvline(95, ls=':', color='grey'); axes[1, 1].axvline(95, ls=':', color='grey')
    fig.suptitle('Local V4/double, 15 m centre; same spatial grid and non-averaged cover interfaces')
    fig.savefig(args.out/'material_factor_ascans.png', dpi=130); plt.close(fig)
    (args.out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
